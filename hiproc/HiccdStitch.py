#!/usr/bin/env python
"""Stitch together the same HiRISE color-filter CCDs from an observation to
create a single mosaicked image file.

This program:

- Optionally performs cubenorm processing on each CCD using a training area
    that has been selected by the user. Cubenorm processing
    is required whenever radiometric problems, such as vertical
    striping, remain after the radiometric calibration step in the
    HiCal Pipeline. The cubenorm processing is described in the
    ISIS program cubenorm.
- Balance the CCD products, created by HiStitch, to radiometrically match.
    A multiplicative constant is applied to each CCD to force the
    overlapping areas of adjacent CCDs to be identical. The resulting
    balanced CCD products, have the naming convention ``*.balance.cub``.
    These products are used by subsequent pipeline processing for
    creating color products, RDR products, and DTM products.
- Join (stitch together) the CCD products to form an image of the entire
    observation.

HiccdStitch is the fourth step in the HiRISE processing chain, after
HiStitch. If an observation has RED, IR and BG CCDs, then this program would
need to be run once for each set.


Data Flow
---------
Input Products:

- ``.cub`` files of the same color (RED, IR, or BG) which are the result of
    HiStitch.
- ``.json`` files that start with the CCD ID of each source .cub file.

Output Products:

- A stitched ``.cub`` file for that color.
- A ``.balance.cub`` file for each input cube file.
- A ``.json`` file with summary information about the stitched image.

"""

# Copyright 2006-2020, Arizona Board of Regents on behalf of the Lunar and
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
# This program is based on HiccdStitch version 2.4.8 (2020/11/03),
# and on the Perl HiccdStitch program ($Revision: 1.45 $
#                                      $Date: 2020/11/03 23:44:15 $)
# by Eric Eliason as an employee of the University of Arizona.

import argparse
import csv
import itertools
import json
import logging
import operator
import os
import pkg_resources
import pvl
import statistics
import subprocess
import sys
from pathlib import Path

import kalasiris as isis
import hiproc.hirise as hirise
import hiproc.util as util

logger = logging.getLogger(__name__)


class HiccdStitchCube(hirise.CCDID):
    """A class for HiRISE CCD IDs with additional capabilities for
    HiccdStitch."""

    def __init__(self, pathlike, cubenormstep=False):

        self.path = Path(pathlike)
        self.nextpath = self.path
        super().__init__(hirise.get_CCDID_fromfile(self.path))
        self.ns = int(isis.getkey_k(self.path, "Dimensions", "Samples"))
        self.nl = int(isis.getkey_k(self.path, "Dimensions", "Lines"))
        self.bin = int(isis.getkey_k(self.path, "Instrument", "Summing"))
        self.cubenormstep = cubenormstep
        self.cubenorm_stddev = None
        self.sl_cubenorm = None
        self.nl_cubenorm = None
        self.ss_balance_left = None
        self.ns_balance_left = None
        self.ss_balance_right = None
        self.ns_balance_right = None
        self.sl_balance = None
        self.nl_balance = None
        self.smag = None
        self.lmag = None
        self.ls_path = None
        self.rs_path = None
        self.lm_path = None
        self.rm_path = None
        self.rstats = None
        self.lstats = None
        self.correction = None
        self.hical_status = None
        self.snr_list = list()

    def __repr__(self):
        return f"{self.__class__.__name__}('{self.path}')"

    def set_cubenorm_lines(self, skiptop, skipbot, sline, eline, minbin):
        # For the cubenorm process, determine starting line and number of lines
        start_line = skiptop / self.bin + 1
        end_line = self.nl - skipbot / self.bin
        if end_line < start_line:
            start_line = 1
            end_line = self.nl

        if sline is not None and sline > 0:
            start_line = int(sline * (minbin / self.bin))
            if start_line > self.nl:
                raise ValueError(
                    f"The effective starting line ({sline}) of "
                    "the training area is outside the range of "
                    f"the image ({self.path})"
                )

        if eline is not None and eline > 0:
            end_line = int(eline * (minbin / self.bin))
            if end_line > self.nl:
                raise ValueError(
                    f"The effective ending line ({eline}) of "
                    "the training area is outside the range of "
                    f"the image ({self.path})"
                )

        if end_line < start_line:
            raise RuntimeError(
                "Something wrong with the calculated training "
                "area: end_line < start_line for cubenorm "
                "trainging area.  The calculated start, end "
                f"lines are: {start_line}, {end_line}."
            )

        self.sl_cubenorm = int(start_line)
        self.nl_cubenorm = int(end_line - start_line + 1)
        return (self.sl_cubenorm, self.nl_cubenorm)

    def set_balance(
        self, skiptop: int, skipbot: int, area: dict, minbinlines: int
    ):
        # For balance cube process, determine the crop area for left
        # and right crop areas and the scaling needed for each crop are.
        (skip, samps) = area[self.bin]

        self.ss_balance_left = int(skip + 1)
        self.ns_balance_left = int(samps)
        self.ss_balance_right = int(self.ns + 1 - skip - samps)
        self.ns_balance_right = int(samps)

        sl_bal = skiptop / self.bin + 1
        nl_bal = minbinlines / self.bin - skipbot / self.bin - sl_bal + 1
        if nl_bal < 1:
            sl_bal = 1
            nl_bal = minbinlines / self.bin
        self.sl_balance = int(sl_bal)
        self.nl_balance = int(nl_bal)

        return

    def set_ls_path(self, p: Path):
        self.ls_path = p
        self.lm_path = p

    def set_rs_path(self, p: Path):
        self.rs_path = p
        self.rm_path = p

    def gather_from_db(self, dbs: list = None):
        """There is some data that we need to pull from the Channel DB files,
        which this method finds and extracts."""
        bad_flag = False
        if dbs is None:
            dbs = list()
            # Scan the parent directory of cube for any .json DB files that
            # match the CCD name
            for p in self.path.parent.glob(str(self) + "*.json"):
                with open(p, "r") as f:
                    dbs.append(json.load(f))

        for db in dbs:
            if "hical_status" in db:
                if "BadCal" == db["hical_status"]:
                    bad_flag = True
                else:
                    self.hical_status = db["hical_status"]

            if "IMAGE_SIGNAL_TO_NOISE_RATIO" in db:
                self.snr_list.append(float(db["IMAGE_SIGNAL_TO_NOISE_RATIO"]))

        if bad_flag:
            self.hical_status = "BadCal"
        return


def arg_parser():
    parser = argparse.ArgumentParser(
        description=__doc__,
        parents=[util.parent_parser()],
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-o", "--output",
        required=False,
        default=".HiccdStitch.cub",
        help="The name of the output .cub file to write.  Optionally, if "
             "it starts with a '.' it is considered a suffix"
             "and will be added to the Observation ID of the "
             "input files. Default: %(default)s",
    )
    parser.add_argument(
        "-c",
        "--conf",
        required=False,
        type=argparse.FileType('r'),
        default=pkg_resources.resource_stream(
            __name__,
            'data/HiccdStitch.conf',
        ),
        help="Path to the HiccdStitch config file.  Defaults to "
             "HiccdStitch.conf distributed with the library."
    )
    parser.add_argument(
        "--db",
        required=False,
        default=".HiCat.json",
        help="The .json file to output.  Optionally, if it "
        "starts with a '.' it is considered an extension "
        "and will be swapped with the output file's extension "
        "to determine the .json filename to use. Default: %(default)s",
    )
    parser.add_argument(
        "--sline",
        required=False,
        default=None,
        type=int,
        help="If given, will be used as the starting line to crop the image "
             "at in order to create a training area to use for cubenorm "
             "processing."
    )
    parser.add_argument(
        "--eline",
        required=False,
        default=None,
        type=int,
        help="If given, will be used as the ending line for cubenorm cropping, "
             "see --sline for more information."
    )
    parser.add_argument(
        "--cubenorm",
        required=False,
        nargs="+",
        type=bool,
        help="To engage cubenorm processing a list of true or false values "
             "(could be 0 or 1) that must match the number of input cubes"
             "must be given to indicate which of the input cubes should "
             "have cubenorm processing applied.  If you only had four input"
             "cubes, then ``--cubenorm 0 0 1 1`` would not run cubenorm "
             "processing on the first two cubes, but would run it on the last"
             "two, etc.  The default is not to run this processing on any."
    )
    parser.add_argument(
        "cubes",
        metavar="cub-file",
        nargs="+",
        help="Cubes to assemble, which are presumably the output of HiStitch. "
             "They must all be from the same detectors, so either all RED, "
             "all IR, or all BG cubes."
    )
    return parser


def main():
    # The Original Perl took a .pvl file as input which mostly just had the
    # filenames of the ccd files to stitch together.  We'll just take those
    # on the command line and into args.cubes.

    args = arg_parser().parse_args()

    util.set_logger(args.verbose, args.logfile, args.log)

    # outcub_path = set_outcube(args.output, pid0)

    if args.cubenorm is not None:
        if len(args.cubenorm) != len(args.cubes):
            logger.critical(
                f"The number of cubes ({len(args.cubes)}) and "
                "the number of cubenorm flags given "
                f"({len(args.cubenorm)}) did not match. Exiting."
            )
            sys.exit()
    else:
        args.cubenorm = [False] * len(args.cubes)

    # Perl: GetConfigurationParameters()
    conf = pvl.load(args.conf)
    conf_check(conf)

    # We may not need to read anything from the DB, only write to it?
    # Nope, we need to read in each of the CHANNEL!(!) HiCat files to get
    # their hical_status value.  Hmm.  Rather than devise a
    # command line strategy for manually loading them, I think we'll just
    # use the args.cubes filenames to find them.

    # Perl: GetPVLParameters()
    # Gets some options and items from the pvl file that the.
    #    Need to recreate?
    # Looks like there's a True/False for each input about Cubenorming or not
    # Then a HiccdStitch/Start_line, HiccdStitch/End_Line,
    # and a HiccdStitch/Reduce_Factor
    # If those aren't present then then default to 0, 0, -9999.  Probably
    # should all be 'None's.
    # Upon inspection of the HiStitch_Next_Pipe, there is no logic to set these
    # values, they are simply always FALSE for each CCD, 0, 0, and -9999.
    # So maybe this was put in for manual use, but not a 'normal' part of the
    # pipeline?
    # Adding as command line arguments

    if len(args.cubes) == 1:
        # Assume this is a filepath to a file of cube names
        cubes = list(
            map(
                HiccdStitchCube,
                Path(args.cubes[0]).read_text().splitlines(),
                args.cubenorm,
            )
        )
    else:
        cubes = list(map(HiccdStitchCube, args.cubes, args.cubenorm))

    for c in cubes:
        c.gather_from_db()

    with util.main_exceptions(args.verbose):
        (db, outpath) = HiccdStitch(
            cubes,
            args.output,
            conf,
            args.sline,
            args.eline,
            keep=args.keep,
        )

    db_path = util.path_w_suffix(args.db, outpath)

    with open(db_path, "w") as f:
        json.dump(db, f, indent=0, sort_keys=True)

    return


def HiccdStitch(
    cubes: list,
    out_path: os.PathLike,
    conf: dict,
    sline=None,
    eline=None,
    keep=False
) -> tuple:
    logger.info(f"HiccdStitch start: {', '.join(map(str, cubes))}")

    # Perl: GetConfigurationParameters()
    # conf = pvl.load(str(conf_path))
    conf_check(conf)

    out_p = set_outpath(out_path, cubes)

    # Perl: GetImageDims()
    cubes = GetImageDims(cubes, conf, sline, eline)

    # This string will get placed in the filename for all of our
    # temporary files. It will (hopefully) prevent collisions with
    # existing files and also allow for easy clean-up if keep=True
    # temp_token = datetime.now().strftime('HiccdStitch-%y%m%d%H%M%S')

    to_delete = isis.PathSet()

    # Perl: MakeList()
    # This makes a file, but doesn't really do anything with it, because
    # a different file is actually made later to give to hiccdstitch

    for i, c in enumerate(cubes):
        if c.cubenormstep:
            cubes[i] = CubeNormStep(c, conf["HiccdStitch"], keep=keep)

    if conf["HiccdStitch"]["HiccdStitch_Balance"]:
        cubes = BalanceStep(cubes, conf["HiccdStitch"], keep=keep)

    for c in cubes:
        SpecialProcessingFlags(c)

    cubes.sort()

    logger.info(
        "The Original Perl looked for a custom file for "
        "hiccdstitch's shiftdef parameter, but the default ISIS "
        "file seems better, so this isn't implemented."
    )

    with isis.fromlist.temp([str(c.nextpath) for c in cubes]) as f:
        isis.hiccdstitch(
            fromlist=f,
            to=out_p,
            interp=conf["HiccdStitch"]["HiccdStitch_Interpolation"],
        )

    SNR_Check(cubes, conf["HiccdStitch"]["HiccdStitch_SNR_Threshold"])

    logger.info("HiccdStitch done.")
    if not keep:
        to_delete.unlink()

    # Afterwards inserts into CCD_Processing_Statistics table.
    db = {"OBSERVATION_ID": str(cubes[0].get_obsid())}
    for c in cubes:
        ccd_db = {
            "CCDID": str(c),
            "RADIOMETRIC_MATCHING_CORRECTION": c.correction,
            "LEFT_OVERLAP_AVERAGE": c.lstats,
            "RIGHT_OVERLAP_AVERAGE": c.rstats,
            "CUBENORM_COLUMN_CORRECTION_STANDARD_DEVIATION": c.cubenorm_stddev,
        }
        db[c.get_ccd()] = ccd_db

    logger.info(f"HiccdStitch done: {out_p}")
    return db, out_p


def conf_check(conf: dict) -> dict:
    """Various checks on parameters in the configuration."""

    util.conf_check_strings(
        "HiccdStitch_Clean",
        ("DELETE", "KEEP"),
        conf["HiccdStitch"]["HiccdStitch_Clean"],
    )

    util.conf_check_bool(
        "HiccdStitch_Version_Enable",
        conf["HiccdStitch"]["HiccdStitch_Version_Enable"],
    )

    util.conf_check_bool(
        "HiccdStitch_Balance", conf["HiccdStitch"]["HiccdStitch_Balance"]
    )

    util.conf_check_bounds(
        "HiccdStitch_Normalization_Minimum",
        (-16384.0, 16364.0),
        conf["HiccdStitch"]["HiccdStitch_Normalization_Minimum"],
    )

    util.conf_check_bounds(
        "HiccdStitch_Normalization_Maximum",
        (-16384.0, 16364.0),
        conf["HiccdStitch"]["HiccdStitch_Normalization_Maximum"],
    )

    util.conf_check_strings(
        "HiccdStitch_Cubenorm_Method",
        ("DIVIDE", "SUBTRACT"),
        conf["HiccdStitch"]["HiccdStitch_Cubenorm_Method"],
    )

    util.conf_check_strings(
        "HiccdStitch_Interpolation",
        ("NEAREST", "BILINEAR", "CUBIC"),
        conf["HiccdStitch"]["HiccdStitch_Interpolation"],
    )

    util.conf_check_strings(
        "HiccdStitch_Balance_Correction",
        ("MULTIPLY", "ADD"),
        conf["HiccdStitch"]["HiccdStitch_Balance_Correction"],
    )

    util.conf_check_strings(
        "HiccdStitch_Balance_Method",
        ("AVERAGE", "MEDIAN"),
        conf["HiccdStitch"]["HiccdStitch_Balance_Method"],
    )

    util.conf_check_count(
        "HiccdStitch_Bin01_Area",
        2,
        "value",
        conf["HiccdStitch"]["HiccdStitch_Bin01_Area"],
    )

    util.conf_check_count(
        "HiccdStitch_Bin02_Area",
        2,
        "value",
        conf["HiccdStitch"]["HiccdStitch_Bin02_Area"],
    )

    util.conf_check_count(
        "HiccdStitch_Bin03_Area",
        2,
        "value",
        conf["HiccdStitch"]["HiccdStitch_Bin03_Area"],
    )

    util.conf_check_count(
        "HiccdStitch_Bin04_Area",
        2,
        "value",
        conf["HiccdStitch"]["HiccdStitch_Bin04_Area"],
    )

    util.conf_check_count(
        "HiccdStitch_Bin08_Area",
        2,
        "value",
        conf["HiccdStitch"]["HiccdStitch_Bin08_Area"],
    )

    util.conf_check_count(
        "HiccdStitch_Bin16_Area",
        2,
        "value",
        conf["HiccdStitch"]["HiccdStitch_Bin16_Area"],
    )


def set_outpath(out: os.PathLike, cubes) -> Path:
    # Check that they all have the same Obs ID:
    if len(set(map(lambda c: c.get_obsid(), cubes))) != 1:
        raise ValueError(
            "These cube files don't all have the same "
            f"Observation ID: {cubes}"
        )

    if str(out).startswith("."):
        return cubes[0].path.parent / Path(
            str(cubes[0].get_obsid()) + "_" + cubes[0].ccdname + out
        )
    else:
        return Path(out)


def GetImageDims(cubes: list, in_conf: dict, sline, eline) -> list:
    """Gathers information for each image."""
    # The functionality is different in this function than the original
    # Perl because I've placed more logic into the HiccdStitchCube class.
    conf = in_conf["HiccdStitch"]

    # cubes should be a list of HiccdStitchCube objects.
    minbin = min(map(lambda c: c.bin, cubes))
    minbinlines = min(map(lambda c: c.bin * c.nl, cubes))

    for c in cubes:
        c.set_cubenorm_lines(
            conf["HiccdStitch_Bin1_Skip_Top_Lines"],
            conf["HiccdStitch_Bin1_Skip_Bot_Lines"],
            sline,
            eline,
            minbin,
        )

        area = make_area_dict(conf)
        c.set_balance(
            conf["HiccdStitch_Bin1_Skip_Top_Lines"],
            conf["HiccdStitch_Bin1_Skip_Bot_Lines"],
            area,
            minbinlines,
        )

        if c.bin == minbin:
            ns_balance_scale = c.ns_balance_left
            nl_balance_scale = c.nl_balance

    for i, c in enumerate(cubes):
        cubes[i].smag = ns_balance_scale / cubes[i].ns_balance_left
        cubes[i].lmag = nl_balance_scale / cubes[i].nl_balance

    return cubes


def make_area_dict(conf: dict) -> dict:
    d = {
        1: conf["HiccdStitch_Bin01_Area"],
        2: conf["HiccdStitch_Bin02_Area"],
        3: conf["HiccdStitch_Bin03_Area"],
        4: conf["HiccdStitch_Bin04_Area"],
        8: conf["HiccdStitch_Bin08_Area"],
        16: conf["HiccdStitch_Bin16_Area"],
    }
    return d


def CubeNormStep(cube, hconf: dict, keep=False) -> HiccdStitchCube:
    to_del = isis.PathSet()

    # crop removing the top and bottom portion of the image
    crop_p = to_del.add(cube.nextpath.with_suffix(".crop.cub"))
    isis.crop(
        cube.nextpath,
        to=crop_p,
        line=cube.sl_cubenorm,
        nlines=cube.nl_cubenorm,
    )

    # run cubenorm to get statistics of the cropped area
    stats_p = to_del.add(cube.nextpath.with_suffix(".cubenorm.tab"))
    isis.cubenorm(
        crop_p,
        stats=stats_p,
        format_="TABLE",
        direction="COLUMN",
        normalizer="AVERAGE",
        MODE=hconf["HiccdStitch_Cubenorm_Method"],
        PRESERVE=True,
    )

    stats_filtered_p = to_del.add(cube.nextpath.with_suffix(".cubenorm2.tab"))

    # Original Perl: write cubenorm_stdev to the DB (again?) but we'll just
    cube.cubenorm_stddev = AnalyzeStats(stats_p, stats_filtered_p)

    # Original Perl: if HiccdStitch: cubenorm_stdev is written to the DB here
    # Original Perl: if HiccdStitchC: write cubenorm_stdev to a PVL file?

    # run cubenorm again, this time make the correction to the file
    next_path = cube.nextpath.with_suffix(".cubenorm.cub")
    to_s = "{}+SignedWord+{}:{}".format(
        next_path,
        hconf["HiccdStitch_Normalization_Minimum"],
        hconf["HiccdStitch_Normalization_Maximum"],
    )
    isis.cubenorm(
        cube.nextpath,
        to=to_s,
        fromstats=stats_filtered_p,
        statsource="TABLE",
        direction="COLUMN",
        normalizer="AVERAGE",
        MODE=hconf["HiccdStitch_Cubenorm_Method"],
        PRESERVE=True,
    )

    cube.nextpath = next_path

    if not keep:
        to_del.unlink()

    return cube


def AnalyzeStats(
    cubenorm_out_p: os.PathLike, filtered_p: os.PathLike
) -> float:
    # Does more than just 'analyze,' also performs some correction.
    valid_points = list()
    averages = list()
    medians = list()
    other_cols = list()
    # header = list()
    with open(cubenorm_out_p) as csvfile:
        reader = csv.DictReader(csvfile, dialect=isis.cubenormfile.Dialect)
        # header = reader.fieldnames
        for row in reader:
            valid_points.append(int(row.pop("ValidPoints")))
            averages.append(float(row.pop("Average")))
            medians.append(float(row.pop("Median")))
            other_cols.append(row)

    # Calculate the average of the column averages
    average = sum(map(operator.mul, valid_points, averages)) / sum(
        valid_points
    )

    # Treat the left and right halves independently
    half = int((len(valid_points) + 1) / 2)

    maxvp_1sthalf = max(valid_points[:half])
    for i, vp in enumerate(valid_points[:half]):
        if vp < maxvp_1sthalf:
            averages[i] = average
            medians[i] = average
    maxvp_2ndhalf = max(valid_points[half:])
    for i, vp in enumerate(valid_points[half:]):
        if vp < maxvp_2ndhalf:
            averages[i] = average
            medians[i] = average

    # Write the results of the filtered cubenorm data
    with open(filtered_p, "w") as csvfile:
        writer = isis.cubenormfile.DictWriter(csvfile)
        writer.writeheader()
        for (d, vp, av, md) in zip(
            other_cols, valid_points, averages, medians
        ):
            d["ValidPoints"] = str(vp)
            d["Average"] = "{:f}".format(av)
            d["Median"] = "{:f}".format(md)
            writer.writerow(d)

    # compute standard deviation of averages
    return statistics.variance(averages)


def BalanceStep(cubes, conf, keep=False) -> list:
    to_del = isis.PathSet()

    # Sort the cubes so that they are in CCD order
    cubes.sort()

    cubes, to_delete = crop_and_scale(cubes)
    to_del.update(to_delete)

    # Original Perl: Generate CCD number array for each CCD file, not needed
    # now, since we can just query the object.

    # Original Perl: Used $0 (the program name) instead of the index 0
    # here, but I've arranged things with the set_ls_path() and set_rs_path()
    # to just make these a full copy so there's no
    # need to mess with a conditional assignment and also streamlines
    # the following logic, since they're identical except when you
    # explicitly change them.

    # The third step is to mask the left and right overlap areas. We
    # want to zap pixels where there is not common coverage.
    for i, c in enumerate(cubes):
        if i + 1 < len(cubes) and int(cubes[i].ccdnumber) + 1 == int(
            cubes[i + 1].ccdnumber
        ):
            cubes[i].rm_path = to_del.add(
                c.nextpath.with_suffix(".right.mask.cub")
            )
            cubes[i + 1].lm_path = to_del.add(
                cubes[i + 1].nextpath.with_suffix(".left.mask.cub")
            )

            for f, m, t in zip(
                [cubes[i].rs_path, cubes[i + 1].ls_path],
                [cubes[i + 1].ls_path, cubes[i].rs_path],
                [cubes[i].rm_path, cubes[i + 1].lm_path],
            ):
                isis.mask(
                    f,
                    mask=m,
                    to=t,
                    preserve="INSIDE",
                    min_=conf["HiccdStitch_Normalization_Minimum"],
                    max_=conf["HiccdStitch_Normalization_Maximum"],
                )

    # The fourth step is to get image statistics for left and right
    # overlap areas of each CCD image.
    cubes = get_stats(cubes)

    # Look for a break in joining CCDs, defined by a break in the CCD number,
    # or the right or left statistics are undefined, due to an all null channel
    #
    # In the original Perl there was a loop to determine if there was a break,
    # but then nothing was done with that information?  Oh, it was used
    # differently: the code past that point develops a series of sequences
    # from $first to $last.  If there are no breaks, then it only runs a
    # single sequence.  If there are breaks, it runs the sequences it finds.
    #
    # Here's the pythonic version:
    cubes.sort()
    for (offset, group) in get_group_i(cubes):
        logger.info("Correction before redistribution.")
        for ccd in group:
            i = ccd + offset
            cubes[i].correction = get_correction(
                cubes[i],
                cubes[i - 1],
                conf["HiccdStitch_Balance_Correction"],
                i,
            )
            logger.info(
                f"CCDID: {cubes[i]}, correction: {cubes[i].correction}"
            )

        normalization = get_normalization(
            cubes, group, offset, conf["HiccdStitch_Control_CCD"]
        )

        logger.info("Correction after redistribution.")
        for ccd in group:
            i = ccd + offset
            cubes[i].correction /= normalization
            logger.info(
                f"CCDID: {cubes[i]}, correction: {cubes[i].correction}, "
                f"left: {cubes[i].lstats}, right: {cubes[i].rstats}"
            )

            # In the original Perl, they wrote out to the DB here, but we'll
            # do it later.  There was also a distinction that if it was
            # HiccdStitchC that the data was written out to a PVL file.  Not
            # sure why.

        # Create the balance cubes
        for ccd in group:
            i = ccd + offset
            balance_path = cubes[i].nextpath.with_suffix(".balance.cub")
            make_balance(cubes[i], conf, balance_path)
            cubes[i].nextpath = balance_path

    if not keep:
        to_del.unlink()

    return cubes


def crop_and_scale(cubes: list) -> list:
    to_del = isis.PathSet()
    for i, c in enumerate(cubes):
        # First step in balancing process is to crop out the left and
        # right overlap areas of each CCD
        lc_path = to_del.add(c.nextpath.with_suffix(".left.crop.cub"))
        rc_path = to_del.add(c.nextpath.with_suffix(".right.crop.cub"))

        for t, s, n in zip(
            [lc_path, rc_path],
            [c.ss_balance_left, c.ss_balance_right],
            [c.ns_balance_left, c.ns_balance_right],
        ):
            isis.crop(
                c.nextpath,
                to=t,
                sample=s,
                nsamples=n,
                line=c.sl_balance,
                nlines=c.nl_balance,
            )

        # Second step is to scale all of the croped files to have the
        # same lines and samples, needed for mask step
        if c.smag == 1 and c.lmag == 1:
            cubes[i].set_ls_path(lc_path)
            cubes[i].set_rs_path(rc_path)
        else:
            cubes[i].set_ls_path(
                to_del.add(c.nextpath.with_suffix(".left.scale.cub"))
            )
            cubes[i].set_rs_path(
                to_del.add(c.nextpath.with_suffix(".right.scale.cub"))
            )
            for lc, ls in zip(
                [lc_path, rc_path], [cubes[i].ls_path, cubes[i].rs_path]
            ):
                isis.enlarge(
                    lc, to=ls, sscale=c.smag, lscale=c.lmag, interp="CUBIC"
                )
    return (cubes, to_del)


def get_group_i(cubes: list) -> list:
    """Given a list of CCDIDs, return a list of lists where each list contains
    continuous CCD number range indexes that have good statistics."""
    cubes.sort()
    good_indexes = list()
    for i, c in enumerate(cubes):
        if c.rstats is not None:
            if i + 1 < len(cubes):
                if cubes[i + 1].lstats is not None:
                    good_indexes.append(int(c.ccdnumber))
            else:
                good_indexes.append(int(c.ccdnumber))

    good_indexes.sort()
    # print(f'good indexes: {good_indexes}')
    pairs = list()
    for offset, group in itertools.groupby(
        good_indexes, lambda x, c=itertools.count(): next(c) - x
    ):
        pairs.append((offset, list(group)))
    return pairs


def get_stats(cubes: list) -> list:
    for i, c in enumerate(cubes):
        cubes[i].lstats = float(
            pvl.loads(isis.stats(c.lm_path).stdout)["Results"]["Average"]
        )
        cubes[i].rstats = float(
            pvl.loads(isis.stats(c.rm_path).stdout)["Results"]["Average"]
        )

    return cubes


def get_correction(this_c, prev_c, Balance_corr, is_not_first=True) -> float:
    if "MULTIPLY" == Balance_corr:
        if is_not_first:
            return prev_c.correction * (prev_c.rstats / this_c.lstats)
        else:
            return 1
    else:
        if is_not_first:
            logger.warning(
                "Original Perl issue: non-MULTIPLY value "
                "for HiccdStitch_Balance_Correction leads "
                "to a correction value that is calculated "
                "strangely."
            )
            # I think the below should be:
            # prev_c.correction + (prev_c.rstats - this_c.lstats)
            return prev_c.correction + (prev_c.rstats - this_c.rstats)
        else:
            return 0


def get_normalization(
    cubes: list, group: list, offset: int, control_ccds: list
) -> float:
    # If a CCD is a control CCD, defined in HiccdStitch_Control_CCD
    # keyword in the confguration file then use it as the normalization,
    # otherwise use the average
    normalization = None
    corrections = list()
    for ccd in group:
        i = ccd + offset
        corrections.append(cubes[i].correction)
        if int(cubes[i].ccdnumber) in control_ccds:
            normalization = cubes[i].correction

    if normalization is None:
        normalization = statistics.mean(corrections)

    return normalization


def make_balance(cube, conf, balance_path):
    a_param = None
    c_param = None
    if conf["HiccdStitch_Balance_Correction"] == "MULTIPLY":
        a_param = cube.correction
        c_param = 0
    elif conf["HiccdStitch_Balance_Correction"] == "ADD":
        a_param = 1
        c_param = cube.correction

    to_s = "{}+SignedWord+{}:{}".format(
        balance_path,
        conf["HiccdStitch_Normalization_Minimum"],
        conf["HiccdStitch_Normalization_Maximum"],
    )
    isis.algebra(
        cube.nextpath, to=to_s, operator="unary", a=a_param, c=c_param
    )
    logger.info(f"Created {balance_path}")


def SpecialProcessingFlags(cube: HiccdStitchCube):
    """Set the special processing flags in the ISIS label."""
    status = "NOMINAL"

    if cube.hical_status == "BadCal":
        status = "BADCAL"

    if cube.cubenormstep:
        status = "CUBENORM"

    try:
        isis.getkey_k(cube.nextpath, "Instrument", "Special_Processing_Flag")
        option = "MODKEY"
    except subprocess.CalledProcessError:
        option = "ADDKEY"

    isis.editlab(
        cube.nextpath,
        option=option,
        grpname="Instrument",
        keyword="Special_Processing_Flag",
        value=status,
    )


def SNR_Check(cubes: list, snr_threshold: float):
    snr_list = list()
    for c in cubes:
        logger.info(f"{c}, binning: {c.bin}, SNRs: {c.snr_list}")
        if c.bin == 1:
            if len(c.snr_list) > 0:
                snr_list.extend(c.snr_list)
            else:
                logger.info(
                    "Without SNR values, this won't be included in "
                    "the check."
                )
        else:
            logger.info(
                "Since binning is more than 1, this CCD won't be "
                " included in the check."
            )

    if len(snr_list) == 0:
        return

    snr_average = statistics.mean(snr_list)
    if snr_average < snr_threshold:
        logger.warning(
            f"The average SNR is {snr_average}, which is greater "
            f"than the threshold: {snr_threshold}"
        )
