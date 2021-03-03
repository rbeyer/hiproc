#!/usr/bin/env python
"""Start from .img files and run them through a complete processing chain."""

# Copyright 2020-2021, Ross A. Beyer (rbeyer@seti.org)
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

# Need to think about endpoints, etc.
# Minimal 'run' from channels to single file is to stop at HiccdStitch
#
# If you knew you only wanted to stop at HiNoProj, could do this next,
# but in order to run HiPrecisionInit to decide HiNoProj or HiJACK, need
# slither files.
#
# For the remaining 'end points' these next three must be done:
#   - HiColorInit
#   - HiJitReg
#   - HiSlither
#
# "Color" end:
# Then run HiColorNorm and HiBeautify
#
# "Red" end:
# HiPrecisionInit, and then HiNoProj or HiJACK.
#
# So, stops are: HiccdStitch, "Color", and "red"


import argparse
import itertools
import logging
from pathlib import Path

import hiproc.hirise as hirise
import hiproc.util as util
import hiproc.EDR_Stats as EDR_Stats
import hiproc.HiCal as HiCal
import hiproc.HiStitch as HiStitch
import hiproc.HiccdStitch as HiccdStitch
import hiproc.HiColorInit as HiColorInit
import hiproc.HiJitReg as HiJitReg
import hiproc.HiSlither as HiSlither
import hiproc.HiColorNorm as HiColorNorm
import hiproc.HiBeautify as HiBeautify
import hiproc.HiPrecisionInit as HiPrecisionInit
import hiproc.HiNoProj as HiNoProj
import hiproc.HiJACK as HiJACK

logger = logging.getLogger(__name__)


class ChannelCube(hirise.ChannelID):
    def __init__(self, pathlike, db=None):
        self.path = Path(pathlike)
        self.nextpath = self.path
        super().__init__(hirise.get_ChannelID_fromfile(self.path))
        self.db = db


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, parents=[util.parent_parser()]
    )
    parser.add_argument(
        "-c", "--color", action="store_true", help="Perform Color processing."
    )
    parser.add_argument(
        "-p",
        "--precision",
        action="store_true",
        help="Perform HiPrecision processing to result "
        "in NOPROJ.cub that might be processed through "
        "HiJACK.",
    )
    parser.add_argument(
        "-j",
        "--jack",
        action="store_true",
        help="Do not test if HiJACK is needed, but force " "HiJACK to run.",
    )
    parser.add_argument(
        "img",
        metavar="some.img-file",
        nargs="+",
        help="More than one can be listed here.",
    )
    parser.add_argument(
        "--conf_dir",
        required=False,
        type=Path,
        default=Path(__file__).resolve().parent.parent / "data",
    )
    parser.add_argument(
        "--db",
        required=False,
        default=".HiCat.json",
        help="The .json file to use.  Optionally, if it "
        "starts with a '.' it is considered an extension "
        "and will be swapped with the input file's "
        "extension to find the .json file to use.",
    )

    args = parser.parse_args()

    util.set_logger(args.verbose, args.logfile, args.log)

    oid = hirise.get_ObsID_fromfile(args.img[0])
    parent = Path(args.img[0]).parent
    if args.img == 1:
        imgs = parent.glob(f"{oid}*.img") + parent.glob(f"{oid}*.IMG")
    else:
        imgs = args.img

    db_list = None
    try:
        get_cubes(f"{oid}*balance.cub", parent)
    except FileNotFoundError:
        chancubes = edr2stitch(imgs, args.conf_dir, keep=args.keep)
        db_list = [x.db for x in chancubes]

    if args.color:
        color(oid, parent, args.conf_dir, parent, db_list, keep=args.keep)

    if args.precision or args.jack:
        jack = False
        if args.jack:
            jack = True
        else:
            slithers = list(parent.glob(f"*.slither.txt"))

            # HiPrecisionInit to determine if you need to HiNoProj or HiJACK
            #   takes *slither.txt
            HiJACK_flags = HiPrecisionInit.start(
                slithers, args.conf_dir / "HiPrecisionInit.conf"
            )

            if any(HiJACK_flags):
                jack = True

        precision(oid, args.conf_dir, parent, hijack=jack, keep=args.keep)


def get_cubes(glob: str, parent: Path) -> list:
    cubes = list(parent.glob(glob))
    if len(cubes) == 0:
        raise FileNotFoundError(
            "Did not find any files that match " f"{glob} in {parent}."
        )
    return cubes


def edr2stitch(images, conf_dir, keep=False):
    chids = list()
    for i in images:
        out_edr = util.path_w_suffix(".EDR_Stats.cub", i)

        # EDR_Stats
        db = EDR_Stats.EDR_Stats(
            i, out_edr, conf_dir / "EDR_Stats_gains_config.pvl", keep=keep
        )

        # HiCal
        out_hical = util.path_w_suffix(".HiCal.cub", out_edr)

        db = HiCal.start(
            out_edr,
            out_hical,
            db,
            HiCal.conf_setup(conf_dir / "HiCal.conf", "NoiseFilter.conf"),
            conf_dir / "HiCal.conf",
            None,
            None,
            keep=keep,
        )

        chids.append(ChannelCube(out_hical, db))

    # HiStitch
    # get Channel pairs
    cids = list()
    for chid1, chid2 in get_CCDpairs(chids):
        (db, o_path) = HiStitch.start(
            chid1.nextpath,
            chid2.nextpath,
            chid1.db,
            chid2.db,
            ".HiStitch.cub",
            conf_dir / "HiStitch.conf",
            keep=keep,
        )
        cid = HiccdStitch.HiccdStitchCube(o_path)
        cid.gather_from_db(db)
        cids.append(cid)

    # HiccdStitch, makes balance cubes
    # need to separate by color:
    color_groups = get_color_groups(cids)
    for color_group in color_groups.values():
        db, out_stitch = HiccdStitch.start(
            color_group,
            conf_dir / "HiccdStitch.conf",
            ".HiccdStitch.cub",
            sline=None,
            eline=None,
            keep=keep,
        )
    # HiColorInit
    #   takes *balance.cub
    #   creates *[IR|BG]*.balance.precolor.cub
    #   Can then run JitPlot on these *.balance.precolor.cub
    HiColorInit.start([c.nextpath for c in cids], ".precolor.cub", keep=keep)

    # HiJitReg
    #   takes tmp/*balance.cub tmp/*balance.precolor.cub
    #   creates *regdef.pvl and *flat.tab files
    for_jitreg = list()
    for color, balcubes in color_groups.items():
        if color == "RED":
            for c in balcubes:
                for_jitreg.append(c.nextpath)
        else:
            for c in balcubes:
                for_jitreg.append(c.nextpath.with_suffix(".precolor.cub"))

    HiJitReg.start(for_jitreg, conf_dir / "HiJitReg.conf", keep=keep)

    # HiSlither
    #   takes same as HiJitReg (and assumes its products are available.
    #   creates *slither.txt, *slither.cub, and *COLOR[4|5].cub
    #   Can then run SliterStats on the *slither.txt
    HiSlither.start(for_jitreg)

    return chids


def color(
    obsid, path: Path, conf_dir: Path, parent: Path, db_list: list, keep=False
):
    colors = get_cubes(f"{obsid}_COLOR*.cub", parent)

    # HiColorNorm - only for color
    #   takes *COLOR[4|5].cub
    #   creates *UNFILTERED_COLOR[4|5].cub and *COLOR[4|5].HiColorNorm.cub
    HiColorNorm.start(
        colors, "_COLOR.cub", conf_dir / "HiColorNorm.conf", db_list, keep=keep
    )

    # HiBeautify - only for color
    #   takes tmp/*.HiColorNorm.cub
    #   creates *IRB.cub and *RGB.cub
    for x in colors:
        x.with_suffix(".HiColorNorm.cub")
    HiBeautify.start(
        [x.with_suffix(".HiColorNorm.cub") for x in colors],
        conf_dir / "HiBeautify.conf",
    )
    return


def precision(obsid, conf_dir: Path, parent: Path, hijack=False, keep=False):

    if hijack:
        bal_cubs = get_cubes(f"{obsid}*balance.cub", parent)
        # HiJACK - starts with balance cubes, needs all the colors
        #   takes all color balance.cubs
        #   creates a whole slow of files in HiJACK/
        HiJACK.start(bal_cubs, conf_dir, outdir=(parent / "HiJACK"), keep=keep)
    else:
        red_bal_cubs = get_cubes(f"{obsid}_RED*balance.cub", parent)
        # HiNoProj - alternate to HiccdStitch, starts with balance cubes
        #   takes only REDS: *RED*balance.cub
        #   creates PSP_010502_2090_RED.NOPROJ.cub
        HiNoProj.start(red_bal_cubs, conf_dir / "HiNoProj.conf", keep=keep)

    return


def get_CCDpairs(chids: list) -> list:
    """Return a list of two-tuples which contain pairs of ChannelCubes
    that belong to the same CCD from *chids*.
    """
    chids.sort()
    pairs = list()
    for k, g in itertools.groupby(chids, lambda x: x.ccdnumber):
        pairs.append(list(g))

    return pairs


def get_color_groups(cids: list) -> list:
    """Return a list of lists where each list contains
    ChannelCubes that are from the same COLOR of detector
    from *cids*.
    """
    cids.sort()
    color_groups = dict()
    for k, g in itertools.groupby(cids, lambda x: x.ccdname):
        color_groups[k] = list(g)

    return color_groups
