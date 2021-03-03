#!/usr/bin/env python
"""HiSlither uses the output of HiJitReg to perform the registration of
   color CCDs to red CCDs and output a color cube."""

# Copyright 2007-2020, Arizona Board of Regents on behalf of the Lunar and
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
# This program is based on HiColor version 5.4.3 (2020/04/28),
# and on the Perl HiSlither program ($Revision: 1.43 $
#                                    $Date: 2020/04/28 16:56:16 $)
# by Guy McArthur as an employee of the University of Arizona.

import argparse
import logging
import math
from datetime import datetime

import kalasiris as isis
import hiproc.util as util
import hiproc.HiColorInit as hicolor
import hiproc.HiJitReg as hjr

logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, parents=[util.parent_parser()]
    )
    parser.add_argument(
        "cubes",
        metavar="RED balance.cub and color balance.precolor.cub files",
        nargs="+",
    )

    args = parser.parse_args()

    util.set_logger(args.verbose, args.logfile, args.log)

    start(args.cubes, keep=args.keep)
    return


def start(cube_paths: list, keep=False):
    cubes = list(map(hicolor.HiColorCube, cube_paths))
    (red4, red5, ir10, ir11, bg12, bg13) = hicolor.separate_ccds(cubes)

    if cube_check(red4, ir10, bg12):
        HiSlither(red4, ir10, bg12, keep=keep)

    if cube_check(red5, ir11, bg13):
        HiSlither(red5, ir11, bg13, keep=keep)

    return


def cube_check(
    red: hicolor.HiColorCube, ir: hicolor.HiColorCube, bg: hicolor.HiColorCube
) -> bool:
    if red is None and ir is None and bg is None:
        return False

    if red is None or bg is None:
        have_ccds = [filter(lambda x: x is not None, [red, ir, bg])]
        if len(have_ccds) == 1:
            have_str = str(have_ccds[0])
        else:
            have_str = "{} and {}".format(*have_ccds)
        logger.info(
            "For a slither set, we need at least a RED and BG CCD."
            + f"We only have: {have_str}"
        )
        return False

    return True


def HiSlither(
    red: hicolor.HiColorCube,
    ir: hicolor.HiColorCube,
    bg: hicolor.HiColorCube,
    keep=False,
):
    if red is None:
        raise TypeError("Expected a RED HiColorCube, but got None.")
    if bg is None:
        raise TypeError("Expected a BG HiColorCube, but got None.")

    logger.info("Beginning HiSlither for {}".format(str(red)))

    temp_token = datetime.now().strftime("HiSlither-%y%m%d%H%M%S")

    if ir is None:
        # If the IR band is missing, it is filled with a "dummy" band
        # consisting entirely of null pixels. This band is created by
        # masking the red band with itself. The dummy band allows for
        # RGB (synthetic B) to be created by HiBeautify even if the
        # IR band is missing.
        logger.info("Missing IR, creating a null band.")
        ir_slither_p = make_dummy_IR(red, bg)
        ir = hicolor.HiColorCube(ir_slither_p)
    else:
        ir_slither_p = run_slither(ir)

    bg_slither_p = run_slither(bg)

    logger.info("Slither finished, merging stacks.")
    logger.info(
        "Creating {}-{}-{} mosaic".format(
            ir.get_ccd(), red.get_ccd(), bg.get_ccd()
        )
    )

    stacked_path = red.path.parent / "{}_COLOR{}.cub".format(
        red.get_obsid(), red.ccdnumber
    )
    tmp_stack = stacked_path.with_suffix(f".{temp_token}.cub")
    isis.hicubeit(ir=ir_slither_p, red=red.path, bg=bg_slither_p, to=tmp_stack)

    trim_args = {"from": tmp_stack, "to": stacked_path}
    # trim 5 px in bin 1, 3 in bin 2, 2 in bin 4
    trim = math.ceil(5 / red.bin)
    if red.ccdnumber == "4":
        trim_args["right"] = trim
    else:
        trim_args["left"] = trim
    isis.trim(**trim_args)

    if not keep:
        tmp_stack.unlink()

    return stacked_path


def get_slither_path(cube):
    return cube.path.parent / "{}.slither.cub".format(str(cube))


def run_slither(cube):
    cnet_path = hjr.JitterCube.get_cnet_path(cube)
    s = get_slither_path(cube)
    t = s.with_suffix(".txt")
    isis.slither(
        cube.path,
        to=s,
        control=cnet_path,
        results=t,
        dir="REVERSE",
        spline="CUBIC",
        interp="BILIN",
    )
    return s


def make_dummy_IR(red, bg):
    bg_slither_path = get_slither_path(bg)
    ir_name = bg_slither_path.name.replace(
        bg.get_ccd(), "IR" + str(int(bg.ccdnumber) - 2)
    )
    ir_path = bg_slither_path.parent / ir_name

    if ir_path.exists():
        raise FileExistsError(
            "{} exists ".format(str(ir_path))
            + "and we don't want to overwrite it."
        )

    ir_ccd = "IR" + str(int(red.ccdnumber) + 6)

    isis.mask(red, mask=red, to=ir_path, preserve="outside")
    isis.editlab(
        ir_path,
        options="modkey",
        grpname="Instrument",
        keyword="CpmmNumber",
        value=int(red.ccdnumber) + 2,
    )
    isis.editlab(
        ir_path,
        options="modkey",
        grpname="Instrument",
        keyword="CcdId",
        value=ir_ccd,
    )
    isis.editlab(
        ir_path,
        options="modkey",
        grpname="Archive",
        keyword="ProductID",
        value="{}_{}".format(red.get_obsid(), ir_ccd),
    )
    isis.editlab(
        ir_path,
        options="modkey",
        grpname="BandBin",
        keyword="Name",
        value="NearInfrared",
    )
    isis.editlab(
        ir_path,
        options="modkey",
        grpname="BandBin",
        keyword="Center",
        value=900,
        units="NANOMETERS",
    )
    isis.editlab(
        ir_path,
        options="modkey",
        grpname="BandBin",
        keyword="Width",
        value=200,
        units="NANOMETERS",
    )
    return ir_path
