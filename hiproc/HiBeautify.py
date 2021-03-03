#!/usr/bin/env python
"""Creates special color products ("pretty pictures")."""

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
# This program is based on HiColor version 5.4.3 2020/04/28
# and on the Perl HiBeautify program ($Revision: 1.92 $
#                                     $Date: 2020/04/28 16:56:16 $)
# by Guy McArthur as an employee of the University of Arizona.

import argparse
import logging
import os
from datetime import datetime
from pathlib import Path

import pvl

import kalasiris as isis
import hiproc.util as util
import hiproc.HiColorNorm as hcn

logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        parents=[util.parent_parser()],
        conflict_handler="resolve",
    )
    parser.add_argument(
        "-o_irb", "--output_irb", required=False, default="_IRB.cub"
    )
    parser.add_argument(
        "-o_rgb", "--output_rgb", required=False, default="_RGB.cub"
    )
    parser.add_argument(
        "-c",
        "--conf",
        required=False,
        default=Path(__file__).resolve().parent.parent
        / "data"
        / "HiBeautify.conf",
    )
    # parser.add_argument('-f', '--frost', action='store_true',
    #                     help='Use the frost/ice color stretch, and disable '
    #                     'auto-detection of frost/ice.')
    # parser.add_argument('--nofrost', action='store_true',
    #        help='Do not use the frost/ice color stretch, and disable '
    #                     'auto-detection of frost/ice.')
    parser.add_argument(
        "cubes",
        metavar=(
            "the COLOR4.HiColorNorm.cub and " "COLOR5.HiColorNorm.cub files"
        ),
        nargs=2,
    )

    args = parser.parse_args()

    util.set_logger(args.verbose, args.logfile, args.log)

    start(
        args.cubes, args.conf, args.output_irb, args.output_rgb, keep=args.keep
    )


def start(
    cube_paths: list,
    conf_path: os.PathLike,
    out_irb="_IRB.cub",
    out_rgb="_RGB.cub",
    keep=False,
):

    # GetConfigurationParameters()
    conf = pvl.load(str(conf_path))

    cubes = list(map(hcn.ColorCube, cube_paths))

    cubes.sort()

    outirb_path = hcn.set_outpath(out_irb, cubes)
    outrgb_path = hcn.set_outpath(out_rgb, cubes)

    # HiBeautify(cubes, outcub_paths, conf, frost=args.frost, keep=args.keep)
    HiBeautify(cubes, (outirb_path, outrgb_path), conf, keep=keep)


def HiBeautify(cubes: list, outcub_paths: list, conf: dict, keep=False):
    temp_token = datetime.now().strftime("HiBeautify-%y%m%d%H%M%S")
    irb_out_p = Path(outcub_paths[0])
    rgb_out_p = Path(outcub_paths[1])
    to_del = isis.PathSet()

    # Create an IRB mosaic from the HiColorNorm halves.

    # If a half is missing, we create a mosaic with the proper width and place
    # the half in it at the proper location.
    # Logic ofr image_midpoint and total_width come from original HiColorInit
    # irbmerged_p = to_del.add(out_p.with_suffix(f'.{temp_token}_IRB.cub'))
    image_midpoint = int((2000 / cubes[0].red_bin) + 1)
    outsample = image_midpoint
    if cubes[0].ccdnumber == "4":
        outsample = 1
    total_width = int(2 * cubes[0].samps - (48 / cubes[0].red_bin))
    import subprocess

    try:
        isis.handmos(
            cubes[0].path,
            mosaic=irb_out_p,
            outline=1,
            outsample=outsample,
            outband=1,
            create="Y",
            nlines=cubes[0].lines,
            nsamp=total_width,
            nbands=3,
        )
    except subprocess.CalledProcessError as err:
        print(err.stdout)
        print(err.stderr)

    if len(cubes) == 1:
        logger.info("Warning, missing one half!")
    else:
        logger.info("Using both halves")
        isis.handmos(
            cubes[1].path,
            mosaic=irb_out_p,
            outline=1,
            outsample=image_midpoint,
            outband=1,
        )

    # Nothing is actually done to the pixels here regarding the FrostStats, so
    # I'm just going to skip them here.
    #
    # # Determine if Frost/ICE may be present using FrostStats module.
    # frost = None
    # if args.frost:
    #     frost = True
    #     logging.info('Frost override: disabling auto-detection and using '
    #                  'the frost/ice color stretch')
    # if args.nofrost:
    #     frost = False
    #     logging.info('Frost override: disable auto-detection and not using '
    #                  'the frost/ice color stretch')
    # if frost is None:
    #     pass
    #     # get frost stats

    # Subtract the unaltered RED band from the high pass filtered BG for
    # synthetic blue.
    logger.info("Creating synthetic B, subtracting RED from BG")
    rgbsynthb_p = to_del.add(irb_out_p.with_suffix(f".{temp_token}_B.cub"))
    isis.algebra(
        f"{irb_out_p}+3",
        from2=f"{irb_out_p}+2",
        op="subtract",
        to=rgbsynthb_p,
        A=conf["Beautify"]["Synthetic_A_Coefficient"],
        B=conf["Beautify"]["Synthetic_B_Coefficient"],
    )

    # HiBeautify gathers and writes a bunch of statistics to PVL that is
    # important to the HiRISE GDS, but not relevant to just producing pixels
    # so I'm omitting it.
    #
    # # Determine the min and max DN values of each band (RED, IR, BG, B) we're
    # # working with.
    # (upper, lower) = conf['Beautify']['Stretch_Exclude_Lines']
    # if upper == 0 and lower == 0:
    #     synthbcrp_p = rgbsynthb_p
    #     irbmrgcrp_p = irbmerged_p
    # else:
    #     synthbcrp_p = to_del.add(
    #           out_p.with_suffix(f'.{temp_token}_Bx.cub'))
    #     irbmrgcrp_p = to_del.add(
    #           out_p.with_suffix(f'.{temp_token}_IRBx.cub'))
    #
    #     for (f, t) in ((rgbsynthb_p, synthbcrp_p),
    #                    (irbmerged_p, irbmrgcrp_p)):
    #         logging.info(isis.crop(f, to=t, propspice=False,
    #                                line=(1 + upper),
    #                                nlines=(
    #                                   cubes[0].lines - lower + upper)).args)
    #
    # stats = dict()
    # stats['B'] = Get_MinMax(synthbcrp_p,
    #                         conf['Beautify']['Stretch_Reduction_Factor'],
    #                         temp_token, keep=keep)
    #
    # for band in cubes[0].band.keys():
    #     stats[band] = Get_MinMax('{}+{}'.format(str(irbmrgcrp_p),
    #                                             cubes[0].band[band]),
    #                              conf['Beautify']['Stretch_Reduction_Factor'],
    #                              temp_token, keep=keep)

    # Create an RGB cube using the RED from the IRB mosaic,
    # the BG from the IRB mosaic and the synthetic B that we just made.
    isis.cubeit_k(
        [f"{irb_out_p}+2", f"{irb_out_p}+3", rgbsynthb_p], to=rgb_out_p
    )

    if not keep:
        to_del.unlink()

    return


# Turns out that HiBeautify doesn't really need this functionality.  This was
# just development code and hasn't been tested.
#
# def Get_MinMax(cube, scale=1,
#                temp_token=datetime.now().strftime('%y%m%d%H%M%S'),
#                keep=False) -> tuple:
#     '''Reduces an Isis cub by a specified scaling factor.  Runs hist to
#        return the pixel values from the histogram at the ordinal
#        positions passed (or min & max if unspecified).
#
#        Returns (min, max, avg).
#     '''
#     # Taken from HiArch, not HiBeautify, but seems to be called here first.
#
#     in_p = Path(cube)
#     to_del = isis.PathSet()
#
#     if scale > 1:
#         reduce_p = to_del.add(in_p.with_suffix(f'{temp_token}.reduce.cub'))
#         logging.info(isis.reduce(in_p, to=reduce_p, sscale=scale,
#                                  lscale=scale, validper=1,
#                                  vper_replace='nearest', mode='SCALE').args)
#
#         in_p = reduce_p
#
#     elif scale < 1:
#         raise ValueError(
#             f'The value for scale must be >= 1, but was {scale}.')
#
#     stats = isis.stats_k(in_p)
#
#     if not keep:
#         to_del.unlink()
#
#     if int(stats['NullPixels']) >= int(stats['TotalPixels']):
#         return (-1, -1, 0)
#
#     return(int(stats['Minimum']), int(stats['Maximum']),
#            int(stats['Average']))
