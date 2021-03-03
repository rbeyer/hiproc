#!/usr/bin/env python
"""Resamples images from individual CCDs to ideal camera geometry and
mosaicks them into a single cube."""

# Copyright 2008-2020, Arizona Board of Regents on behalf of the Lunar and
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
# This program is based on HiPrecision version 3.2.3 (2020/06/23),
# and on the Perl HiNoProj program ($Revision: 1.10 $
#                                    $Date: 2020/02/15 00:36:19 $)
# by Audrie Fennema as an employee of the University of Arizona.

import argparse
import itertools
import logging
import os
import re
import shutil
import subprocess
from datetime import datetime
from pathlib import Path

import pvl

import kalasiris as isis
import hiproc.util as util
import hiproc.HiColorInit as hci
import hiproc.HiColorNorm as hcn

logger = logging.getLogger(__name__)


class Cube(hci.HiColorCube):
    def __init__(self, pathlike):
        super().__init__(pathlike)
        self.next_path = None
        self.line_offset = 0
        self.samp_offset = 0


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, parents=[util.parent_parser()]
    )
    parser.add_argument(
        "-o", "--output", required=False, default="_RED.NOPROJ.cub"
    )
    parser.add_argument(
        "-c",
        "--conf",
        required=False,
        default=Path(__file__).resolve().parent.parent
        / "data"
        / "HiNoProj.conf",
    )
    parser.add_argument("-b", "--base_ccd_number", required=False, default=5)
    parser.add_argument("cubes", metavar="balance.cub-files", nargs="+")

    args = parser.parse_args()

    util.set_logger(args.verbose, args.logfile, args.log)

    try:
        start(
            args.cubes,
            args.conf,
            args.output,
            args.base_ccd_number,
            keep=args.keep,
        )
    except ValueError as err:
        print(err)
    except subprocess.CalledProcessError as err:
        print("Had an ISIS error:")
        print(err.cmd)
        print(err.stdout)
        print(err.stderr)
        raise err
    return


def start(
    cube_paths: os.PathLike,
    conf_path: os.PathLike,
    output="_RED.NOPROJ.cub",
    base_ccd_number=5,
    keep=False,
):

    cubes = list(map(Cube, cube_paths))
    cubes.sort()

    if not all(c.ccdname == "RED" for c in cubes):
        raise ValueError("Not all of the input files are RED CCD files.")

    sequences = list()
    for k, g in itertools.groupby(
        (int(c.ccdnumber) for c in cubes),
        lambda x, c=itertools.count(): next(c) - x,
    ):
        sequences.append(list(g))

    if len(sequences) != 1:
        raise ValueError(
            "The given cubes are not a single run of sequential "
            "HiRISE CCDs, instead there are "
            f"{len(sequences)} groups with these "
            f"CCD numbers: {sequences}."
        )

    base_ccd = list(
        filter(lambda x: x.ccdnumber == str(base_ccd_number), cubes)
    )
    if len(base_ccd) != 1:
        raise ValueError(
            f"The base ccd, number {base_ccd_number}, "
            "is not one of the given cubes."
        )

    conf = pvl.load(str(conf_path))
    conf_check(conf["HiNoProj"])

    outcub_path = hcn.set_outpath(output, cubes)

    HiNoProj(cubes, base_ccd[0], outcub_path, conf["HiNoProj"], keep=keep)


def conf_check(conf: dict) -> None:
    """Various checks on parameters in the configuration."""

    t = "Shape"
    util.conf_check_strings(t, ("ELLIPSOID", "SYSTEM", "USER"), conf[t])

    return


def is_polar(cubes, pole_tolerance: float, temp_token: str) -> bool:

    abs_lats = list()
    for c in cubes:
        temp_p = c.path.with_suffix(f".{temp_token}.spice.cub")
        shutil.copyfile(c.path, temp_p)

        isis.spiceinit(temp_p)

        cpvl = pvl.loads(isis.camrange(temp_p).stdout)

        abs_lats.append(
            abs(float(cpvl["UniversalGroundRange"]["MinimumLatitude"]))
        )
        abs_lats.append(
            abs(float(cpvl["UniversalGroundRange"]["MaximumLatitude"]))
        )

        temp_p.unlink()

    if max(abs_lats) > float(pole_tolerance):
        return True
    else:
        return False


def handmos_side(cubes, base_cube, out_p: os.PathLike, left=True):
    """Runs handmos to add cubes which are to one side or the other
       of the *base_cube*."""
    # handmos left side
    side = 1
    priority = "beneath"
    if left:
        side = -1
        priority = "ontop"
    ssm = 1
    slm = 1
    for c in cubes[(cubes.index(base_cube) + side) :: side]:
        slm -= side * round(c.line_offset)
        ssm -= side * round(c.samp_offset)
        isis.handmos(
            c.next_path,
            mosaic=out_p,
            priority=priority,
            outline=slm,
            outsample=ssm,
            outband=1,
        )
    return


def copy_and_spice(
    inpath: os.PathLike, outpath: os.PathLike, conf: dict, polar=False
):
    shutil.copyfile(inpath, outpath)

    spiceinit_args = {
        "from": outpath,
        "shape": conf["Shape"],
        "cksmithed": True,
        "spksmithed": True,
        "spkrecon": True,
        "spkpredicted": False,
    }

    if conf["Shape"] == "USER":
        if polar:
            spiceinit_args["model"] = conf["Polar_Shape_Model_Path"]
        else:
            spiceinit_args["model"] = conf["Shape_Model_Path"]

    isis.spiceinit(**spiceinit_args)

    return


def get_offsets(
    cube: os.PathLike, match: os.PathLike, flat: os.PathLike
) -> tuple:
    isis.hijitreg(cube, match=match, flat=flat)

    with open(flat, "r") as f:
        flat_text = f.read()

        match = re.search(r"#\s+Average Line Offset:\s+(\S+)", flat_text)
        avg_line_offset = float(match.group(1))

        match = re.search(r"#\s+Average Sample Offset:\s+(\S+)", flat_text)
        avg_samp_offset = float(match.group(1))

    return (avg_line_offset, avg_samp_offset)


def add_offsets(
    cubes: list, base_ccdnumber: int, temp_token: str, keep=False
) -> tuple:
    flats = isis.PathSet()
    for i, c in enumerate(cubes[:-1]):
        pair = "{}-{}".format(cubes[i + 1].get_ccd(), cubes[i].get_ccd())

        flat_p = flats.add(
            c.path.with_suffix(f".{temp_token}.{pair}.flat.tab")
        )

        (avg_line_offset, avg_samp_offset) = get_offsets(
            cubes[i].next_path, cubes[i + 1].next_path, flat_p
        )

        j = i
        if int(c.ccdnumber) >= base_ccdnumber:
            j = i + 1

        cubes[j].line_offset = avg_line_offset
        cubes[j].samp_offset = avg_samp_offset

    if not keep:
        flats.unlink()
        return (cubes, None)
    else:
        return (cubes, flats)


def fix_labels(
    cubes: list, path: os.PathLike, matched_cube: str, prodid: str
) -> None:
    isis.editlab(
        path,
        option="modkey",
        grpname="Archive",
        keyword="ProductId",
        value=prodid,
    )
    isis.editlab(
        path,
        option="modkey",
        grpname="Instrument",
        keyword="MatchedCube",
        value=str(matched_cube),
    )

    # Fix ck kernel in InstrumentPointing in RED label
    # This doesn't seem to be needed, maybe it was HiROC-specific.

    #  Add SourceProductIds to Archive group in label
    logger.info(
        "Original Perl just assumes that both channels are included "
        "in the balance cube."
    )
    source_ids = list()
    for c in cubes:
        source_ids.append(
            isis.getkey_k(c.path, "Instrument", "StitchedProductIds")
        )

    isis.editlab(
        path,
        option="ADDKEY",
        grpname="Archive",
        keyword="SourceProductId",
        value="({})".format(", ".join(source_ids)),
    )
    return


def HiNoProj(
    cubes: list, base_cube, outcub_path: os.PathLike, conf: dict, keep=False
):

    temp_token = datetime.now().strftime("HiNoProj-%y%m%d%H%M%S")
    to_del = isis.PathSet()
    out_p = Path(outcub_path)

    polar = False
    if conf["Shape"] == "USER":
        polar = is_polar(cubes, conf["Pole_Tolerance"], temp_token)

    for c in cubes:

        temp_p = to_del.add(c.path.with_suffix(f".{temp_token}.spiced.cub"))
        copy_and_spice(c.path, temp_p, conf, polar)

        isis.spicefit(temp_p)

        c.next_path = to_del.add(
            c.path.with_suffix(f".{temp_token}.noproj.cub")
        )
        c.path = temp_p

    for c in cubes:
        isis.noproj(
            c.path, match=base_cube.path, to=c.next_path, source="frommatch"
        )

    # Run hijitreg on adjacent noproj'ed ccds to get average line/sample offset
    (cubes, _) = add_offsets(
        cubes, int(base_cube.ccdnumber), temp_token, keep=keep
    )

    # Mosaic noproj'ed ccds using average line/sample offset
    shutil.copyfile(base_cube.next_path, out_p)
    logger.info(
        "Original Perl hard codes this file copy from RED5, even if "
        "another cube is selected as the base_ccd."
    )

    handmos_side(cubes, base_cube, out_p, left=True)
    handmos_side(cubes, base_cube, out_p, left=False)

    isis.editlab(
        out_p,
        option="addkey",
        grpname="Instrument",
        keyword="ImageJitterCorrected",
        value=0,
    )
    fix_labels(
        cubes,
        out_p,
        base_cube,
        "{}_{}".format(str(cubes[0].get_obsid()), cubes[0].ccdname),
    )

    if not keep:
        to_del.unlink()

    return
