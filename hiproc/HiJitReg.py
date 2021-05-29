#!/usr/bin/env python
"""HiJitReg registers color CCDs to corresponding red CCDs by using the
ISIS tool hijitreg to perform a deconvolution of jittered image data.

This program corrects for spacecraft jitter and prepares
images for coregistration. Using the ISIS program hijitreg it creates
a reference grid of control points (extension .control.pvl) from
the RED product, and then attempts to locate those points within
the BG and IR products.

The configuration file for HiJitReg (HiJitReg.conf) describes the control
point grid density (currently it's a 4-column, 200-line grid) and
correlation tolerance. It also specifies the sizes of the search
window and pattern window. The pattern window "walks" through the
search window in order to locate the local maximum. The calculated
translation is then recorded in the control net file for slithering.

Some adaptability is built-in to HiJitReg. It will add columns if
a channel is missing, and will triple the grid point density if
less than 25% of the points register on the first pass. It will
also increase the size of the search window if more than two points
have the edge of their pattern box close to or beyond the edge of
the search box ("edgy" points).

It also uses a smoothing algorithm to ignore points that are
out-of-bounds or are poorly registered. A JitterPlot is created
showing the results.

HiJigReg also works on one or both of the "color" sets:

- RED4 - BG12 - IR10
- RED5 - BG13 - IR11

HiJitReg must be run after HiColorInit, and it creates a .hislither.pvl
file which is then submitted to the HiSlither_Pipeline.

Data Flow
---------
Input Products:

- RED4 and 5 ``balance.cub`` files which are the result of HiccdStitch.
- BG and IR ``precolor.cub`` files which are the result of HiColorInit.

Output Products:

- creates regdef.pvl and flat.tab and control.pvl files for each
    BG and IR cube provided.

"""

# Copyright 2004-2020, Arizona Board of Regents on behalf of the Lunar and
# Planetary Laboratory at the University of Arizona.
#   - Orignal Perl program.
#
# Copyright 2020, Ross A. Beyer (rbeyer@seti.org)
#   - Elements of this Python program are are based on the original Perl
#     but the logic here is rewritten from scratch to emulate functionality.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# This program is based on HiColor version 5.4.2 (2020/02/14),
# This program is based on these Perl programs:
# - HiJitReg.pm ($Revision: 1.58 $ $Date: 2020/04/28 16:56:16 $)
# - JitStats ($Revision: 1.11 $ $Date: 2020/04/28 16:56:16 $)
# - JitStats.pm ($Revision: 1.17 $ $Date: 2020/02/14 22:46:49 $)
# by Guy McArthur as an employee of the University of Arizona.

import argparse
import collections
import csv
import itertools
import logging
import math
import os
import pkg_resources
import re
import statistics
from datetime import datetime
from pathlib import Path

import pvl

import kalasiris as isis
import hiproc.util as util
import hiproc.HiColorInit as hicolor

logger = logging.getLogger(__name__)


def arg_parser():
    parser = argparse.ArgumentParser(
        description=__doc__,
        parents=[util.parent_parser()],
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-c",
        "--conf",
        required=False,
        type=argparse.FileType('r'),
        default=pkg_resources.resource_stream(
            __name__,
            'data/HiJitReg.conf'
        ),
        help="Path to the HiJitReg config file.  Defaults to "
             "HiJitReg.conf distributed with the library."
    )
    parser.add_argument(
        "cubes",
        metavar="balance.cub and balance.precolor.cub files",
        nargs="+",
        help="Either one or both sets of RED .balance.cub  and IR/BG "
             ".balance.precolor.cub files. However, that's tedious to type,"
             "so you could just type in *.balance*cub here, and the "
             "program will sort out what it needs."
    )
    return parser


def main():
    args = arg_parser().parse_args()

    util.set_logger(args.verbose, args.logfile, args.log)

    with util.main_exceptions(args.verbose):
        successful_ccds = HiJitReg(
            args.cubes,
            pvl.load(args.conf),
            keep=args.keep
        )

    print("Successful CCDs are:")
    for c in successful_ccds:
        print("\t{}".format(str(c)))
    return


class JitterCube(hicolor.HiColorCube, collections.abc.MutableMapping):
    """A class for collecting and analyzing jitter statistics."""

    def __init__(
        self,
        arg,
        config=pvl.load(
            pkg_resources.resource_stream(
                __name__,
                'data/HiJitReg.conf'
            ),
        ),
        matchccd=None,
    ):
        if isinstance(arg, hicolor.HiColorCube):
            super().__init__(arg.path)
        else:
            super().__init__(arg)

        # self = copy.deepcopy(hi_color_cube)
        self.dictionary = dict()
        # HiColorCube isn't a dictionary, but we'll give access to its
        #   members as if it were:
        self.dictionary["bin"] = self.bin
        self.dictionary["tdi"] = self.tdi
        self.dictionary["lines"] = self.lines
        self.dictionary["samps"] = self.samps

        self.dictionary["CanSlither"] = False

        self.IgnoredPoints = set()

        # Just assuming that all of these will be in the self.dictionary
        # self.RegisterCount = None
        # self.AvgSampleOffset = None
        # self.AvgLineOffset = None
        # self.STDSampleOffset = None
        # self.STDLineOffset = None
        # self.SuspectCount = None
        # self.MatchedCount = None
        # self.RegisterCount = None
        # self.SearchSamples = None
        # self.SearchLines = None
        # self.EdgyCount = None
        # self.MatchedLineCount = None
        # self.Tolerance = None
        # self.Columns = None
        # self.Rows = None
        # self.canSlither = None
        # self.PatternSamples = None
        # self.PatternLines = None
        # self.SearchSamples = None
        # self.SearchLines = None

        if isinstance(config, (Path, str)):
            self.conf = pvl.load(str(config))
        elif isinstance(config, collections.abc.Mapping):
            self.conf = config
        else:
            raise TypeError(
                f"The value for *config* was neither an os.PathLike nor a "
                f"Mapping object, it was a {type(config)}"
            )

        if matchccd is None:
            self.matchccd = hicolor.CCD_Corresponence[self.get_ccd()]
        else:
            self.matchccd = matchccd

        sm = self.conf["Smoothing"]
        self.dictionary["ExcludeLimit"] = sm["Exclude_Limit"]
        self.dictionary["BadnessLimit"] = sm["Badness_Limit"]
        self.dictionary["BoxcarLength"] = sm["Boxcar_Length"]

        self.cnet_path = self.get_cnet_path(self)
        self.regdef_path = self.get_regdef_path(self)
        self.flattab_path = self.get_flattab_path(self)

    def __getitem__(self, key):
        return self.dictionary[key]

    def __setitem__(self, key, value):
        self.dictionary[key] = value

    def __delitem__(self, key):
        del self.dictionary[key]

    def __iter__(self):
        return iter(self.dictionary)

    def __len__(self):
        return len(self.dictionary)

    @staticmethod
    def get_pair_name(cube, matchccd=None):
        if matchccd is None:
            if hasattr(cube, "matchccd"):
                matchccd = cube.matchccd
            else:
                matchccd = hicolor.CCD_Corresponence[cube.get_ccd()]
        pair_name = "{}_{}-{}".format(
            str(cube.get_obsid()), matchccd, cube.get_ccd()
        )
        return pair_name

    @staticmethod
    def _get_path(cube, suffix):
        pair = JitterCube.get_pair_name(cube)
        return cube.path.parent / (pair + suffix)

    @staticmethod
    def get_cnet_path(cube):
        return JitterCube._get_path(cube, ".control.pvl")

    @staticmethod
    def get_regdef_path(cube):
        return JitterCube._get_path(cube, ".regdef.pvl")

    @staticmethod
    def get_flattab_path(cube):
        return JitterCube._get_path(cube, ".flat.tab")

    def reset(self):
        self.IgnoredPoints.clear()
        self.parseRegDefs(self.regdef_path)
        self.parseFlatTab(self.flattab_path)
        self.parseCNetPVL(self.cnet_path)

    def parseRegDefs(self, path=None):
        """Parse the register definition file to obtain the search and pattern
        sizes."""
        if path is None:
            path = self.regdef_path

        p = pvl.load(str(path))
        self["PatternSamples"] = p["AutoRegistration"]["PatternChip"][
            "Samples"
        ]
        self["PatternLines"] = p["AutoRegistration"]["PatternChip"]["Lines"]
        self["SearchSamples"] = p["AutoRegistration"]["SearchChip"]["Samples"]
        self["SearchLines"] = p["AutoRegistration"]["SearchChip"]["Lines"]
        return

    def parseFlatTab(self, path=None):
        """Parses the flat file to obtain jitter registration result
        statistics."""
        if path is None:
            path = self.flattab_path

        with open(path, "r") as f:
            flat = f.read()

            match = re.search(r"#\s+Line Spacing:\s+(\S+)", flat)
            self["LineSpacing"] = float(match.group(1))

            match = re.search(r"#\s+Columns, Rows:\s+(\d+)\s+(\d+)", flat)
            self["Columns"] = int(match.group(1))
            self["Rows"] = int(match.group(2))

            match = re.search(r"#\s+Corr. Tolerance:\s+(\S+)", flat)
            self["Tolerance"] = float(match.group(1))

            match = re.search(r"#\s+Total Registers:\s+(\d+) of (\S+)", flat)
            self["MatchedCount"] = int(match.group(1))
            self["RegisterCount"] = int(match.group(2))

            match = re.search(r"#\s+Number Suspect:\s+(\S+)", flat)
            self["SuspectCount"] = int(match.group(1))

            match = re.search(
                r"#\s+Average Sample Offset:\s+(\S+)\s+StdDev:\s+(\S+)", flat
            )
            self["AvgSampleOffset"] = float(match.group(1))
            self["STDSampleOffset"] = float(match.group(2))

            match = re.search(
                r"#\s+Average Line Offset:\s+(\S+)\s+StdDev:\s+(\S+)", flat
            )
            self["AvgLineOffset"] = float(match.group(1))
            self["STDLineOffset"] = float(match.group(2))

            dialect = csv.Dialect
            dialect.delimiter = " "
            dialect.skipinitialspace = True
            dialect.quoting = csv.QUOTE_NONE
            dialect.lineterminator = "\n"

            reader = csv.DictReader(
                itertools.filterfalse(
                    lambda x: x.startswith("#") or x.isspace() or len(x) == 0,
                    flat.splitlines(),
                ),
                dialect=dialect,
            )

            if "EdgyCount" not in self:
                self["EdgyCount"] = 0

            lineCount = 0
            for row in reader:
                # how many pixels in x is the edge of the pattern box from
                # the reg point
                deltaSamp = (
                    abs(float(row["RegSamp"]) - int(row["MatchSamp"]))
                    + self["PatternSamples"] / 2
                )
                # how many pixels in y is the edge of the pattern box from
                # the reg point
                deltaLine = (
                    abs(float(row["RegLine"]) - int(row["MatchLine"]))
                    + self["PatternLines"] / 2
                )

                # if the edge of the pattern box is more than two pixels away
                # from the search box, increment the count of marginal
                # control points
                if (deltaSamp > (self["SearchSamples"] / 2 - 2)) or (
                    deltaLine > (self["SearchLines"] / 2 - 2)
                ):
                    self["EdgyCount"] += 1

                    logger.info(
                        "Marginal register {} lines, {} samples to "
                        "edge.".format(
                            deltaLine - self["SearchLines"] / 2,
                            deltaSamp - self["SearchSamples"] / 2,
                        )
                    )
                lineCount += 1
            self["MatchedLineCount"] = lineCount
        return

    def parseCNetPVL(self, path=None):
        """Parses the control net output from hijitreg, performs smoothing,
        and sets the array of ignorable points based on smoothing and
        badness ("goodness of fit").
        """

        if path is None:
            path = self.cnet_path

        p = pvl.load(str(path))

        count = [0] * self["Columns"]

        self.control_measures = self._get_control_measures(p)

        lineCount = 0
        for i, cm in enumerate(self.control_measures):
            offset = int(i - self["BoxcarLength"] / 2)
            length = int(self["BoxcarLength"])

            if offset < 0:
                offset = 0
                length = int(self["BoxcarLength"] / 2 + i)

            if self["BoxcarLength"] > (len(self.control_measures) - i):
                # offset not changed
                length = int(
                    self["BoxcarLength"] / 2 + (len(self.control_measures) - i)
                )

            boxcar = map(
                lambda x: x["ErrorMagnitude"],
                self.control_measures[offset : offset + length],
            )

            median = statistics.median(boxcar)
            delta = abs(cm["ErrorMagnitude"] - median)

            if (
                cm["GoodnessOfFit"] > self["BadnessLimit"]
                or delta > self["ExcludeLimit"]
            ):
                self["MatchedCount"] -= 1
                self.IgnoredPoints.add(cm["PointId"])
                logger.info(
                    "Ignorable point {} with ".format(cm["PointId"])
                    + "badness {} and ".format(cm["GoodnessOfFit"])
                    + f"smoothing delta {delta}"
                )
            else:
                if "Row" in cm:
                    lineCount += 1

                if "Column" in cm:
                    count[cm["Column"]] += 1

            if len(tuple(filter(lambda x: x > 3, count))) >= 3:
                self["CanSlither"] = True
        self["MatchedLineCount"] = lineCount
        return

    @staticmethod
    def _get_control_measures(pvl) -> list:
        control_measures = list()

        # Original Perl issue: there were two "conditions" for
        # extracting information, one, labeled "<3.4" was to find
        # a ControlMeasure with a Reference = False key.  The other
        # labeled ">=3.4" was a ControlMeasure with MeasureType =
        # Candidate.  However, this condition really just ended
        # the line-by-line parsing, because the "Candidate"
        # ControlMeasure was the second one in the ControlPoint.
        # The proper logic is to get information from the
        # ControlMeasure that meets the conditions as implemented
        # below.

        for cp in pvl["ControlNetwork"].getlist("ControlPoint"):
            if "PointId" not in cp or "ControlMeasure" not in cp:
                continue
            for cm in cp.getlist("ControlMeasure"):
                if (
                    "MeasureType" not in cm
                    # or 'Reference' not in cm
                    or "GoodnessOfFit" not in cm
                    or "LineResidual" not in cm
                    or "SampleResidual" not in cm
                ):
                    continue

                if cm["MeasureType"] == "RegisteredPixel":
                    cm["ErrorMagnitude"] = math.hypot(
                        cm["SampleResidual"].value, cm["LineResidual"].value
                    )
                    # Tack on a few extra values here, and then append
                    cm["PointId"] = cp["PointId"]
                    match = re.search(
                        r"Row\s+(\d+)\s+Column\s+(\d+)", cp["PointId"]
                    )
                    if match:
                        cm["Row"] = int(match.group(1))
                        cm["Column"] = int(match.group(2))

                    control_measures.append(cm)
        return control_measures

    def filterCNetPVL(self, path=None):
        """Filters the CNET file and adds Ignored point information."""

        if len(self.IgnoredPoints) == 0:
            return

        if path is None:
            path = self.cnet_path

        p = pvl.load(str(path))

        cn = pvl.PVLModule()

        badness = 0
        for (k, v) in p["ControlNetwork"].items():
            if k == "ControlPoint":
                if (
                    v["PointId"] in self.IgnoredPoints
                    and "Ignore" not in v.keys()
                ):
                    v.append("Ignore", True)
                    badness += 1
                    logger.info("Ignoring point {}".format(v["PointId"]))
            cn.append(k, v)

        logger.info(f"{badness} point(s) ignored.")

        new_pvl = pvl.PVLModule(ControlNetwork=cn)

        with open(path, "w") as stream:
            pvl.dump(new_pvl, stream, encoder=pvl.encoder.ISISEncoder())


def HiJitReg(cube_paths: list, conf: dict, keep=False) -> list:
    cubes = list(map(hicolor.HiColorCube, cube_paths))
    red4, red5, ir10, ir11, bg12, bg13 = hicolor.separate_ccds(cubes)

    ccds = list()
    for c in red4, red5, ir10, ir11, bg12, bg13:
        if c is not None:
            ccds.append(str(c))

    logger.info(f"HiJitReg start: {', '.join(map(str, ccds))}")

    successful_ccds = list()
    if red4 is not None:
        for c in [ir10, bg12]:
            if c is not None:
                if jitter_iter(red4, c, conf, keep=keep):
                    logger.info(f"Iterations completed for {c}")
                    successful_ccds.append(c)
    if red5 is not None:
        for c in [ir11, bg13]:
            if c is not None:
                if jitter_iter(red5, c, conf, keep=keep):
                    logger.info(f"Iterations completed for {c}")
                    successful_ccds.append(c)

    # Not going to check to make sure that at most one pair fails.

    if bg12 not in successful_ccds and bg13 not in successful_ccds:
        raise RuntimeError("Registration failed for both BG halves.")

    logger.info(f"HiJitReg done: {', '.join(map(str, successful_ccds))}")
    return successful_ccds


def jitter_iter(
    red: hicolor.HiColorCube,
    color: hicolor.HiColorCube,
    conf: dict,
    keep=False,
) -> bool:
    """Iterates through hijitreg for the color cube."""

    temp_token = datetime.now().strftime("HiJitReg-%y%m%d%H%M%S")

    bin_ratio = color.bin / red.bin

    jit_param = dict()
    conf_ar = conf["AutoRegistration"]
    jit_param["GROUP"] = "HiJitReg"
    jit_param["COLS"] = conf_ar["ControlNet"]["Control_Cols"]
    jit_param["ROWS"] = conf_ar["ControlNet"]["Control_Lines"]
    jit_param["TOLERANCE"] = conf_ar["Algorithm"]["Tolerance"]
    jit_param["PATTERN_SAMPLES"] = conf_ar["PatternChip"]["Samples"]
    jit_param["PATTERN_LINES"] = conf_ar["PatternChip"]["Lines"]
    jit_param["SEARCH_SAMPLES"] = conf_ar["SearchChip"]["Samples"]
    jit_param["SEARCH_LINES"] = conf_ar["SearchChip"]["Lines"]
    jit_param["SEARCHLONGER_SAMPLES"] = conf_ar["SearchLongerChip"]["Samples"]
    jit_param["SEARCHLONGER_LINES"] = conf_ar["SearchLongerChip"]["Lines"]

    if bin_ratio > 3:
        jit_param["TOLERANCE"] -= conf_ar["Algorithm"]["INCREMENT"]

    channels = isis.getkey_k(color.path, "Instrument", "StitchedProductIds")

    coverage = 1.0

    if len(channels) < 2:
        coverage /= 2
        jit_param["COLS"] += jit_param["COLS"] / 2

    # A two-step process with completely different outcomes at each step,
    # so we can't really make a loop.
    step = 1

    logger.info(f"Attempting hijitreg iteration #{step} for {color}")

    color_jitter = JitterCube(color, conf)

    run_HiJitReg(red.path, color_jitter, jit_param, temp_token, keep=keep)

    ret = Analyze_Flat(color_jitter, step, coverage)

    if ret == -1:
        # edgy or suspect points only
        if jit_param["SEARCH_LINES"] == jit_param["SEARCHLONGER_LINES"]:
            return True
        else:
            # use larger search box for all subsequent iterations
            # (other CCDs too)
            jit_param["SEARCH_SAMPLES"] = jit_param["SEARCHLONGER_SAMPLES"]
            jit_param["SEARCH_LINES"] = jit_param["SEARCHLONGER_LINES"]
    elif ret == 0:
        # not enough points found
        # increase grid density
        jit_param["ROWS"] = jit_param["ROWS"] * 2
        if len(channels) >= 2:
            jit_param["COLS"] += 2
            coverage /= 2
    else:
        return True

    step += 1
    logger.info(f"Attempting hijitreg iteration #{step} for {color}")

    # second pass
    run_HiJitReg(red.path, color_jitter, jit_param, temp_token, keep=keep)

    # analyze output again
    ret = Analyze_Flat(color_jitter, step, coverage)

    if ret == 0:
        logger.info(f"Jitter registration failed for {color}")
        return False
    elif ret < 0:
        logger.info("!!! Validation Required !!!")
        return True
    else:
        return True


def run_HiJitReg(
    red_path: os.PathLike,
    color: JitterCube,
    params: dict,
    temptoken: str,
    keep=False,
):
    """Examine output of control net and/or flat file to automatically remove
    out-of-bound points."""

    file_status = "OVERWRITE"
    if color.regdef_path.exists():
        file_status = pvl.load(str(color.regdef_path))["AutoRegistration"][
            params["GROUP"]
        ]["File_Status"]

    if file_status == "KEEP":
        logger.info("Using existing regdef file due to KEEP file status.")
    else:
        logger.info(f"Writing new regdef file {color.regdef_path}")
        logger.info(params)
        write_regdef(color.regdef_path, params)

    tmp_control = color.cnet_path.with_suffix(".net")
    isis.hijitreg(
        red_path,
        match=color.path,
        regdef=color.regdef_path,
        rows=params["ROWS"],
        columns=params["COLS"],
        flat=color.flattab_path,
        cnet=tmp_control,
    )
    isis.cnetbin2pvl(tmp_control, to=color.cnet_path)
    if not keep:
        tmp_control.unlink()
    return


def write_regdef(out_path: os.PathLike, parameters: dict):
    """Writes PVL file that will be given to HiJitReg."""
    out_p = Path(out_path)

    pvl_text = """Object = AutoRegistration

  Version = 2

  Group = {GROUP}
     File_Status  = "OVERWRITE"
     Control_Cols = {COLS}
     Control_Rows = {ROWS}
  End_Group

  Group = Algorithm
    Name      = MaximumCorrelation
    Tolerance = {TOLERANCE}
  End_Group

  Group = PatternChip
    Samples = {PATTERN_SAMPLES}
    Lines   = {PATTERN_LINES}
  End_Group

  Group = SearchChip
    Samples = {SEARCH_SAMPLES}
    Lines   = {SEARCH_LINES}
  End_Group

End_Object

"""

    out_p.write_text(pvl_text.format(**parameters))
    return


def Analyze_Flat(
    cube: JitterCube, step: int, fraction: float, hijitreg=True
) -> int:
    cube.reset()
    cube.filterCNetPVL()

    logger.info(
        "Matched Registers     = {} of {}".format(
            cube["MatchedCount"], cube["RegisterCount"]
        )
    )
    logger.info("Average Sample Offset = {}".format(cube["AvgSampleOffset"]))
    logger.info("Average Line Offset   = {}".format(cube["AvgLineOffset"]))
    logger.info("Edgy Count            = {}".format(cube["EdgyCount"]))
    logger.info("Suspect Points        = {}".format(cube["SuspectCount"]))

    if cube["AvgSampleOffset"] is None or cube["AvgLineOffset"] is None:
        logger.warning("No points met the correlation tolerance.")
        return 0

    if hijitreg and cube["CanSlither"] is False:
        logger.warning(
            "Too few correlated lines found for cubic slither fit."
        )
        return 0

    good_fraction = (cube["MatchedCount"] - cube["SuspectCount"]) / cube[
        "RegisterCount"
    ]

    if good_fraction < 0.5 * fraction and step <= 1:
        logger.info(
            f"Too few correlated points ({good_fraction}) "
            "found at this tolerance."
        )
        return 0

    elif hijitreg and good_fraction < 0.25 * fraction and step <= 2:
        logger.info(
            f"Too few correlated points ({good_fraction}) "
            "found at this tolerance."
        )
        return -1

    if hijitreg and cube["EdgyCount"] > 2 and good_fraction > 0.8 * fraction:
        logger.info("More than two edgy points with search box size.")
        return -1

    if cube["SuspectCount"] > 3:
        logger.info("More than three suspect points with search box size.")
        return -1

    return 1
