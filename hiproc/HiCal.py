#!/usr/bin/env python
"""Generate radiometrically corrected HiRISE Channel products.

These are the functionalities that the EDR_Stats pipeline performs that are
reproduced here:

- Perform a radiometric calibration correction and conversion to "I/F" units
  using the ISIS ``hical`` program.
- Perform furrow correction. The image columns at the channel join
  are checked for furrows. Pixels in the furrow region (at the channel
  join) whose DN values have gone above a threshold are set to the
  null pixel value as these pixels can not be calibrated. If a
  RED-filter image experiences furrowing, a rare event usually caused
  by an improperly commanded observation, then the furrowed pixels
  will be permanently set to null pixels for all of the HiRISE standard
  and extras products. For BG and IR-filter images that make up the
  color products, the furrowed pixels will be interpolated in the
  HiColorNorm pipeline step.
- If an observation is determined to have furrows then an entry
  with the key "zapped" is set to "true" in the output .json file.
- Due to HiRISE instrument instability problems, the HiCal pipeline
  performs a noise reduction procedure to reduce the number of bad
  pixels in an image observation. The noise correction is applied
  when the standard deviation of the dark pixel or mask regions
  exceed a threshold, or if the number of LIS pixels exceeds a
  threshold.  A newer version of that processing can be engaged
  if --bitflipwidth is set to a non-zero value.
- Sometimes operating at lower temperatures can result in LIS pixels
  in the calibration areas, but not in the image area that can result
  in an erroneous calibration.  Some new processing to deal with that
  can be engaged if --lisfix is set to a value less than 1.
- A high-pass filter "cubenorm" step is applied to the calibrated
  image. Due to camera instabilities, residual vertical striping
  often exists in the imaging that is corrected by this empirical
  method. The average standard deviation of this change is calculated
  and stored in the .json file.
- The ISIS program "hidestripe" is applied to suppress horizontal
  stripping seen whenever an observation is acquired using mixed-binning
  commanding. The standard deviation of this change is calculated
  and stored in the .json file.


HiCal is the second step in the HiRISE processing chain, after EDR_Stats and
results in data that is read for the next step: HiStitch.  If you are
running HiCal individually without all EDRs from the Observation in the
same directory, please be aware of the --bin2 and --bin4 options.


Data Flow
---------
Input Products:

- A ``.cub`` file that has been converted from an EDR.
- A ``.json`` file that contains summary information about the .cub file.

Output Products:

- A calibrated ``.cub`` file for each input ``.cub`` file processed.
- Modifies each provided ``.json`` file.

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
# This program is based on HiCal version 4.1.1 (2021/05/03),
# and on the Perl HiCal program ($Revision: 1.58 $
#                                $Date: 2021/05/03 22:33:55 $)
# and on the Perl NoiseFilter program ($Revision: 1.19
#                                      $Date: 2021/04/22 22:38:39 $)
# by Eric Eliason and Richard Leis as employees of the University of Arizona.

import argparse
import collections
import concurrent.futures
import copy
import csv
import json
import logging
import math
import os
import pkg_resources
import re
import statistics
import sys
import textwrap
import warnings
from datetime import datetime
from pathlib import Path

import numpy as np

import kalasiris as isis
import kalasiris.version as isisversion
import pvl.new as pvl

import hiproc.bitflips as bf
import hiproc.lisfix as lisfix
import hiproc.hirise as hirise
import hiproc.util as util

logger = logging.getLogger(__name__)

tempfile_re = re.compile(r"HiCal-\d{12}\.")


def arg_parser():
    parser = argparse.ArgumentParser(
        description=__doc__,
        parents=[util.parent_parser()],
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("-o", "--output", required=False, default=".HiCal.cub")
    parser.add_argument(
        "-c",
        "--conf",
        required=False,
        type=Path,
        default=Path(
            pkg_resources.resource_filename(
                __name__,
                'data/HiCal.conf'
            )
        ),
        help="Configuration file to use, Defaults to "
             "HiCal.conf distributed with the library."
    )
    parser.add_argument(
        "--db",
        required=False,
        default=".HiCat.json",
        help="The .json file to use.  Optionally, if it "
        "starts with a '.' it is considered an extension "
        "and will be swapped with the input file's "
        "extension to find the .json file to use. "
        "Default: %(default)s",
    )
    parser.add_argument(
        "-b",
        "--bitflipwidth",
        required=False,
        type=int,
        default=None,
        help="The number of medstd widths for bit-flip cleaning. Leaving it "
             "unset will result in the setting being determined by the value "
             "of LIS_Mask in the configuration file.  Explicitly setting it to "
             "zero implies no bit-flip cleaning and will use the old "
             "ANALYZECUBENORMSTATS behavior. Explicitly setting it to 5 is "
             "identical to the BITFLIPS setting.",
    )
    parser.add_argument(
        "-t",
        "--lis_tolerance",
        type=float,
        default=None,
        help="If the fraction of LIS pixels in each of the Reverse-Clock and "
             "Buffer areas are higher than this value, lis-fix processing "
             "will occur for that area. Leaving it unset will result in "
             "the value being determined by the value of LIS_Fix in the "
             "configuration fie.  Explicitly setting it to 0.4 is the same "
             "as LIS_Fix = TRUE, and explicitly setting it to 1.0 is the "
             "same as LIS_Fix = FALSE."
    )
    # parser.add_argument('--hgfconf', required=False, default='HiGainFx.conf')
    parser.add_argument(
        "--nfconf",
        required=False,
        type=Path,
        default=Path(
            pkg_resources.resource_filename(
                __name__,
                'data/NoiseFilter.conf'
            )
        ),
        help="Defaults to NoiseFilter.conf distributed with the library."
    )
    parser.add_argument(
        "--bin2",
        required=False,
        action="store_true",
        default=None,
        help="By default HiCal will look at all of the HiRISE .cub files in "
             "the directory with the input file(s) in it, and determine "
             "whether any of the EDRs from this observation was bin2 and will "
             "internally take care of setting this switch, but if you're just "
             "processing some RED channels and not all of them, then although "
             "IR10 may be binned, and you didn't convert it to a .cub file or "
             "even download it, then it isn't there for this program to find, "
             "and explicitly setting this switch is needed."
    )
    parser.add_argument(
        "--nobin2",
        required=False,
        action="store_false",
        default=None,
        dest="bin2",
        help="Forcibly sets the --bin2 switch to off, even if it would "
             "otherwise be automatically determined to be on."
    )
    parser.add_argument(
        "--bin4",
        required=False,
        action="store_true",
        default=None,
        help="Same as the --bin2 switch, but for setting whether any EDRs "
             "from this observation are bin4."
    )
    parser.add_argument(
        "--nobin4",
        required=False,
        action="store_false",
        default=None,
        dest="bin4",
        help="Forcibly sets the --bin4 switch to off, even if it would "
             "otherwise be automatically determined to be on."
    )
    parser.add_argument(
        "-n", "--info",
        action="store_true",
        help="If given, the program will not run, but will instead report "
             "information about the input files and various things that "
             "the program would have done, if it was run regularly."
    )
    parser.add_argument(
        "--max_workers",
        default=None,
        type=int,
        help="If more than one image is provided to work on, this program "
             "will engage multiprocessing to parallelize the work.  This "
             "multiprocessing will default to the number of processors on the "
             "machine.  If you want to throttle this to use less resources on "
             "your machin, indicate the number of processors you want to use."
    )
    parser.add_argument(
        "cube",
        metavar="cube_file(s)",
        nargs="+",
        help="More than one can be listed here.",
    )
    return parser


def main():
    args = arg_parser().parse_args()

    util.set_logger(args.verbose, args.logfile, args.log)

    if len(args.cube) > 1 and (
        not args.output.startswith(".") or not args.db.startswith(".")
    ):
        logger.critical(
            "With more than one input cube file, the --output "
            "and --db must start with a period, and one of them "
            f"does not: {args.output} {args.db}"
        )
        sys.exit()

    try:
        conf = conf_setup(pvl.load(args.conf), pvl.load(args.nfconf))
    except (TypeError, NotADirectoryError, FileNotFoundError) as err:
        logger.critical(err)
        sys.exit()

    if args.bitflipwidth is None:
        d = {"ANALYZECUBENORMSTATS": 0, "BITFLIPS": 5}
        args.bitflipwidth = d[conf["HiCal"]["LIS_Mask"]]

    if args.lis_tolerance is None:
        d = {True: 0.4, False: 1.0}
        args.lis_tolerance = d[conf["HiCal"]["LIS_Fix"]]

    with util.main_exceptions(args.verbose, logger=logger):
        if args.info:
            for c in args.cube:
                db_path = util.pid_path_w_suffix(args.db, c)
                with open(db_path, "r") as f:
                    db = json.load(f)

                if args.info:
                    print(f"HiCal Information for {c}:")
                    print(textwrap.indent(
                        "\n".join(
                            norun_info(c, db, conf, args.bin2, args.bin4)
                        ),
                        "  "
                    ))
            return

        if len(args.cube) == 1:
            # No need to fire up the multiprocessing for just one image.
            db_path = util.pid_path_w_suffix(args.db, args.cube[0])
            with open(db_path, "r") as f:
                db = json.load(f)

            out_cube = util.path_w_suffix(args.output, args.cube[0])
            db = HiCal(
                args.cube[0],
                out_cube,
                db,
                conf,
                args.conf,
                args.bin2,
                args.bin4,
                args.bitflipwidth,
                args.lis_tolerance,
                keep=args.keep,
            )
            with open(db_path, "w") as f:
                json.dump(db, f, indent=0, sort_keys=True)

            logger.info(f"Wrote {db_path}")
        else:
            with concurrent.futures.ProcessPoolExecutor(
                max_workers=args.max_workers
            ) as executor:
                future_dbs = dict()
                for c in args.cube:
                    db_path = util.pid_path_w_suffix(args.db, c)
                    with open(db_path, "r") as f:
                        db = json.load(f)

                    out_cube = util.path_w_suffix(args.output, c)

                    f = executor.submit(
                        HiCal,
                        c,
                        out_cube,
                        db,
                        conf,
                        args.conf,
                        args.bin2,
                        args.bin4,
                        args.bitflipwidth,
                        args.lis_tolerance,
                        keep=args.keep,
                    )
                    future_dbs[f] = db_path

                for future in concurrent.futures.as_completed(future_dbs):
                    db_path = future_dbs[future]
                    with open(db_path, "w") as f:
                        json.dump(future.result(), f, indent=0, sort_keys=True)

                    logger.info(f"Wrote {db_path}")
    return


def conf_setup(conf: dict, nfconf: dict) -> dict:
    # Get Configuration Parameters
    conf_check(conf)

    # # If the sub-conf arguments aren't 'findable', look for them in
    # # the main conf directory.
    # # hgf_path = util.get_path(Path(args.hgfconf), Path(args.conf).parent)
    # nf_path = util.get_path(
    #     Path(nfconf_path),
    #     (
    #         Path(conf_path).parent,
    #         Path(__file__).resolve().parent.parent / "data",
    #     ),
    # )

    # Merge the configuration files together into a single dict
    # conf['HiGainFx'] = pvl.load(str(hgf_path))['HiGainFx']
    # conf["NoiseFilter"] = pvl.load(str(nf_path))["NoiseFilter"]
    conf["NoiseFilter"] = nfconf["NoiseFilter"]

    return conf


def norun_info(cube: os.PathLike, db: dict, conf: dict, bin2: bool, bin4: bool):

    lines = list()
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("ignore")
        in_cube, cid, ccdchan, lis_per = setup(cube, db, conf)
        if len(w) > 0:
            for wa in w:
                lines.append(f"Would have warned: {wa.message}")

    lines.append(
        f'Image LIS: {db["LOW_SATURATED_PIXELS"]} pixels, '
        f'{lis_per:.2f}%'
    )

    # Non-imaging areas:
    label = pvl.load(cube)
    specialpix = getattr(
        isis.specialpixels, label["IsisCube"]["Core"]["Pixels"]["Type"]
    )
    t_name = "HiRISE Calibration Image"
    hci_dict = isis.cube.get_table(in_cube, t_name)
    cal_vals = np.array(hci_dict["Calibration"])
    cal_image = np.ma.masked_outside(
        cal_vals, specialpix.Min, specialpix.Max
    )
    revclk_liscount = np.ma.count_masked(cal_image[:20, :])
    revclk_lisfrac = revclk_liscount / cal_image[:20, :].size
    lines.append(
        f"Reverse-Clock LIS: {revclk_liscount} pixels, "
        f"{revclk_lisfrac:.2%}"
    )

    hca_dict = isis.cube.get_table(cube, "HiRISE Calibration Ancillary")
    calbuf_vals = np.array(hca_dict["BufferPixels"])
    calbuf = np.ma.masked_outside(
        calbuf_vals, specialpix.Min, specialpix.Max
    )
    calbuf_liscount = np.ma.count_masked(calbuf[:20, :])
    calbuf_lisfrac = calbuf_liscount / calbuf[:20, :].size
    lines.append(
        f"Reverse-Clock Buffer LIS: {calbuf_liscount} pixels, "
        f"{calbuf_lisfrac:.2%}"
    )

    ha_dict = isis.cube.get_table(cube, "HiRISE Ancillary")
    buffer_vals = np.array(ha_dict["BufferPixels"])
    buffer = np.ma.masked_outside(
        buffer_vals, specialpix.Min, specialpix.Max
    )
    buffer_liscount = np.ma.count_masked(buffer)
    buffer_lisfrac = buffer_liscount / buffer.size
    lines.append(
        f"Buffer LIS: {buffer_liscount} pixels, "
        f"{buffer_lisfrac:.2%}"
    )

    b = 1, 2, 4, 8, 16
    flags = set_flags(
        conf["HiCal"], db, ccdchan, b.index(int(db["BINNING"]))
    )
    lines.append(str(flags))

    destripe_filter = check_destripe(in_cube, int(db["BINNING"]), bin2, bin4)
    lines.append(f"Would the destripe filter be run: {destripe_filter}")

    return lines


def setup(cube: os.PathLike, db: dict, conf: dict):
    in_cube = Path(cube)
    cid = hirise.get_ChannelID_fromfile(in_cube)
    if str(cid) != db["PRODUCT_ID"]:
        raise ValueError(
            "The Product ID in the file ({}) does not match "
            "the one in the database ({}).".format(str(cid), db["PRODUCT_ID"])
        )
    ccdchan = (cid.get_ccd(), cid.channel)
    lis_per = (
        float(db["LOW_SATURATED_PIXELS"])
        / (int(db["IMAGE_LINES"]) * int(db["LINE_SAMPLES"]))
        * 100.0
    )
    # if (
    #     ccdchan == ("IR10", "1")
    #     and lis_per > conf["HiCal"]["HiCal_Bypass_IR10_1"]
    # ):
    #     raise UserWarning("Bypassing IR10_1.")
    if lis_per > conf["HiCal"]["HiCal_Bypass"]:
        raise UserWarning(
            f"The image area has a LIS percent ({lis_per}) greater than "
            f"{conf['HiCal']['HiCal_Bypass']}, so this image will not be "
            f"processed through HiCal."
        )

    return in_cube, cid, ccdchan, lis_per


def HiCal(
    in_cube: os.PathLike,
    out_cube: os.PathLike,
    db: dict,
    conf: dict,
    conf_path: os.PathLike,
    bin2: bool,
    bin4: bool,
    bitflipwidth=0,
    lis_tolerance=1.0,
    keep=False,
) -> dict:
    # The original Perl Setup00() builds data structures that we don't
    # need.
    # The original Perl Setup01() set up data routing and did filename
    # checking that we don't need here.

    # The original Perl Setup03 queried the HiCat.Planned_Observations
    # table to get all of the possible binning values for CCDs in this
    # CCD's observation.  We now accomplish the same thing with
    # check_destripe() farther down.

    # The original Perl ProcessingSwitches() set various flags.  That
    # functionality is now broken up and spread out to the check_destripe()
    # function, the set_flags() function, and the if statement below that
    # checks HiCal_Bypass_IR10_1
    in_cube, cid, ccdchan, lis_per = setup(in_cube, db, conf)
    logger.info(f"{cid}: HiCal start: {in_cube}")

    destripe = check_destripe(in_cube, int(db["BINNING"]), bin2, bin4)

    # The original Perl SetHiCalVersion() placed the HiCal version into the
    # db, but since we're not interested in persistence, we can ignore it.

    # Allows for indexing in lists ordered by bin value.
    b = 1, 2, 4, 8, 16

    in_c = Path(in_cube)

    # Keep from having to write out ['HiCal'] all the time.
    hconf = conf["HiCal"]

    # This string will get placed in the filename for all of our
    # temporary files. It will (hopefully) prevent collisions with
    # existing files and also allow for easy clean-up if keep=True
    temp_token = datetime.now().strftime("HiCal-%y%m%d%H%M%S")

    logger.info(f"{cid}: destripe: {destripe}")

    flags = set_flags(
        hconf, db, ccdchan, b.index(int(db["BINNING"])), bitflipwidth > 0
    )
    logger.info(f"{cid}: {flags}")

    # Start processing cube files
    to_delete = isis.PathSet()
    furrows_found = False
    if int(db["BINNING"]) > 1 and float(db["IMAGE_MEAN"]) > 7000.0:
        furrow_cube = to_delete.add(
            in_c.with_suffix(f".{temp_token}.ffix.cub")
        )
        furrows_found = furrow_nulling(
            in_c, furrow_cube, int(db["BINNING"]), ccdchan
        )
        next_cube = furrow_cube
    else:
        next_cube = to_delete.add(in_c.with_suffix(f".{temp_token}.cub"))
        logger.info(f"{cid}: Symlink {in_c} to {next_cube}")
        next_cube.symlink_to(in_c.resolve())
    logger.info(f"{cid}: Furrow Fix application: {furrows_found}")

    # Run hical
    lis_per = (
        int(db["LOW_SATURATED_PIXELS"])
        / (int(db["IMAGE_LINES"]) * int(db["LINE_SAMPLES"]))
        * 100.0
    )
    hical_file = to_delete.add(next_cube.with_suffix(".hical.cub"))
    hical_status = run_hical(
        next_cube,
        hical_file,
        conf,
        conf_path,
        lis_per,
        float(db["IMAGE_BUFFER_MEAN"]),
        int(db["BINNING"]),
        flags.noise_filter,
        bitflipwidth=bitflipwidth,
        lis_tolerance=lis_tolerance,
        keep=keep,
    )
    next_cube = hical_file

    if furrows_found:
        lpfz_file = to_delete.add(next_cube.with_suffix(".lpfz.cub"))
        isis.lowpass(
            next_cube,
            to=lpfz_file,
            lines=3,
            samples=3,
            minopt="COUNT",
            minimum=5,
            filter_="OUTSIDE",
        )
        next_cube = lpfz_file

    # # Perform gain-drift correction
    # if(db['BINNING'] != '8'):  # There is no gain fix for bin8 imaging
    #     higain_file = to_delete.add(next_cube.with_suffix('.fx.cub'))
    #     hgf_path = util.get_path(Path('HiGainFx.conf'),
    #                              Path(conf_path).parent)
    #     conf['HiGainFx'] = pvl.load(str(hgf_path))['HiGainFx']
    #     HiGainFx(next_cube, higain_file,
    #              conf['HiGainFx']['HiGainFx_Coefficient_Path'],
    #              conf['HiGainFx']['HiGainFx_Version'])
    #     next_cube = higain_file

    # Perform the high-pass filter cubenorm steps

    (sl, nl) = set_lines(
        int(hconf["HiCal_Bin1_Skip_Top_Lines"]),
        int(hconf["HiCal_Bin1_Skip_Bot_Lines"]),
        int(db["BINNING"]),
        int(db["IMAGE_LINES"]),
    )
    crop_file = to_delete.add(next_cube.with_suffix(".crop.cub"))
    isis.crop(next_cube, to=crop_file, line=sl, nlines=nl)

    # In the original Perl, this file was kept for the next step, HiStitch.
    # However, since FurrowCheck() is now performed here, it does not need
    # to be kept.
    stats_file = to_delete.add(out_cube.with_suffix(".cubenorm.tab"))
    stats_fix_file = to_delete.add(next_cube.with_suffix(".cubenorm_fix.tab"))
    isis.cubenorm(crop_file, stats=stats_file, format_="TABLE")

    (std_final, zapped) = Cubenorm_Filter(
        stats_file,
        stats_fix_file,
        boxfilter=5,
        pause=True,
        divide=flags.divide,
        chan=int(ccdchan[1]),
    )

    # Now perform the cubnorm_plus correction
    div_or_sub = {True: "DIVIDE", False: "SUBTRACT"}
    cubenorm_args = {
        "direction": "COLUMN",
        "statsource": "TABLE",
        "normalize": "AVERAGE",
        "preserve": "FALSE",
        "mode": div_or_sub[flags.divide],
    }
    cubenorm_file = to_delete.add(next_cube.with_suffix(".cn.cub"))
    to_s = "{}+SignedWord+{}:{}".format(
        cubenorm_file,
        hconf["HiCal_Normalization_Minimum"],
        hconf["HiCal_Normalization_Maximum"],
    )
    isis.cubenorm(
        next_cube, fromstats=stats_fix_file, to=to_s, **cubenorm_args
    )
    next_cube = cubenorm_file

    # NoiseFilter() - external Noise_Filter
    if flags.noise_filter:
        noisefilter_file = to_delete.add(next_cube.with_suffix(".nf.cub"))
        NoiseFilter(
            next_cube,
            output=noisefilter_file,
            conf=conf["NoiseFilter"],
            minimum=hconf["HiCal_Normalization_Minimum"],
            maximum=hconf["HiCal_Normalization_Maximum"],
            zapc=flags.zapcols,
            keep=keep,
        )
        next_cube = noisefilter_file

    # Hidestripe() - isis.[hidestripe,hipass,lowpass,algebra]
    diff_std_dev = None
    if destripe:
        hidestripe_file = to_delete.add(next_cube.with_suffix(".hd.cub"))
        diff_std_dev = Hidestripe(
            next_cube,
            hidestripe_file,
            int(db["BINNING"]),
            hconf["HiCal_Normalization_Minimum"],
            hconf["HiCal_Normalization_Maximum"],
            hconf["HiCal_Hidestripe_Correction"],
            int(db["LINE_SAMPLES"]),
            keep=keep,
        )
        next_cube = hidestripe_file

    # Create final output file
    to_delete.remove(next_cube)
    logger.info(f"{cid}: Rename {next_cube} to {out_cube}.")
    next_cube.rename(out_cube)

    if not keep:
        to_delete.unlink()

    db["HIGH_PASS_FILTER_CORRECTION_STANDARD_DEVIATION"] = std_final
    db["DESTRIPED_DIFFERENCE_STANDARD_DEVIATION"] = diff_std_dev
    db["zapped"] = zapped
    db["hical_status"] = hical_status
    # The zapped flag is not in the original code, I'm putting it in
    # to the DB to use in a later step.

    logger.info(f"{cid}: HiCal done: {out_cube}")
    return db


def FurrowCheck(vpnts: list, channel: int) -> bool:
    # This function was brought forward from HiStitch, because
    # it is more appropriate to perform the check here, and then
    # store the result for later.
    zap = False
    i = -1 * channel

    if vpnts[i] != max(vpnts):
        zap = True

    return zap


def conf_check(conf: dict) -> None:
    """Various checks on parameters in the configuration."""

    util.conf_check_strings(
        "HiCal_Clean_EDR_Stats",
        ("DELETE", "KEEP"),
        conf["HiCal"]["HiCal_Clean_EDR_Stats"],
    )

    util.conf_check_strings(
        "HiCal_HPF_Cubenorm",
        ("DIVIDE", "SUBTRACT"),
        conf["HiCal"]["HiCal_HPF_Cubenorm"],
    )

    util.conf_check_count(
        "HiCal_Noise_Processing",
        28,
        "CCD/channel",
        conf["HiCal"]["HiCal_Noise_Processing"],
    )

    util.conf_check_bounds(
        "HiCal_Bypass",
        (0.0, 100.0),
        conf["HiCal"]["HiCal_Bypass"],
    )

    # util.conf_check_bounds(
    #     "HiCal_Bypass_IR10_1",
    #     (0.0, 100.0),
    #     conf["HiCal"]["HiCal_Bypass_IR10_1"],
    # )

    util.conf_check_count(
        "HiCal_Noise_Bin_DarkPixel_STD",
        5,
        "possible bin mode",
        conf["HiCal"]["HiCal_Noise_Bin_DarkPixel_STD"],
    )

    util.conf_check_count(
        "HiCal_Noise_Bin_Mask_STD",
        5,
        "possible bin mode",
        conf["HiCal"]["HiCal_Noise_Bin_Mask_STD"],
    )

    util.conf_check_bounds(
        "HiCal_Normalization_Minimum",
        (-16384.0, 16364.0),
        conf["HiCal"]["HiCal_Normalization_Minimum"],
    )

    util.conf_check_bounds(
        "HiCal_Normalization_Maximum",
        (-16384.0, 16384.0),
        conf["HiCal"]["HiCal_Normalization_Maximum"],
    )

    util.conf_check_strings(
        "HiCal_Hidestripe_Correction",
        ("ADD", "MULTIPLE"),
        conf["HiCal"]["HiCal_Hidestripe_Correction"],
    )

    util.conf_check_bounds(
        "HiCal_Minimum_Percent",
        (0.0001, 2.0),
        conf["HiCal"]["HiCal_Minimum_Percent"],
    )

    util.conf_check_bounds(
        "HiCal_Maximum_Percent",
        (98.0, 99.9999),
        conf["HiCal"]["HiCal_Maximum_Percent"],
    )

    # There are checks for HiCal_Jpeg_Quality, HiCal_Thumb_Samples, and
    # HiCal_Browse_Samples in the original Perl that aren't relevant here,
    # ignoring.  There are also empty string checks, but if there's nothing
    # in a parameter, it will error when it gets used in the code, so no need
    # to check here.
    return


def check_destripe(
    cube: os.PathLike, mybinning: int, bin2=None, bin4=None
) -> bool:
    """Determines whether to run the destripe filtering based on whether *any*
    of the CCDs in an observation were bin 2 or bin 4.

    Without a database, we must either rely on the incoming bin2 and bin4
    arguments, or scour the directory that the input cube file is in, and
    read the cube labels present to try and find the information.
    """
    destripe_filter = False

    if bin2 is None and mybinning == 2:
        bin2 = True

    if bin4 is None and mybinning == 4:
        bin4 = True

    warn_message = None
    if (mybinning == 1 or mybinning == 2) and (bin2 is None or bin4 is None):
        binnings = get_bins_fromfiles(cube)
        powered_count = len(
            list(
                filter(
                    lambda x: x == "On",
                    isis.getkey_k(cube, "Instrument", "PoweredCpmmFlag").split(
                        ","
                    ),
                )
            )
        )
        if bin2 is None and "2" in binnings.values():
            bin2 = True

        if bin4 is None and "4" in binnings.values():
            bin4 = True

        if len(binnings) < powered_count:
            warn_message = (
                f"Expecting to find {powered_count} CCD files, but found"
                f"{len(binnings)}, as a consequence, cannot determine if the "
                "missing CCDs might have been bin 2 or bin 4.  "
                f"hidestripe will not be applied to this image ({cube})."
            )
            warn_message.format(powered_count, len(binnings))

    bin2 = bool(bin2)
    bin4 = bool(bin4)

    if (1 == mybinning and (bin2 or bin4)) or (2 == mybinning and bin4):
        destripe_filter = True
    elif warn_message is not None:
        logger.warning(warn_message)

    return destripe_filter


def get_bins_fromfiles(cube: os.PathLike) -> dict:
    """Extract summing values from all cubes in the input cube's directory
    that belong to the same Observation.
    """
    bins = dict()
    oid = hirise.get_ObsID_fromfile(cube)
    for path in Path(cube).parent.glob("*.cub"):
        if tempfile_re.search(str(path)) is not None:
            # Want to exclude files like "HiCal-%y%m%d%H%M%S"
            # which during multiprocessing could come and go between
            # the call to glob() above and the call to getkey_k()
            # below.
            continue
        try:
            p = hirise.get_ChannelID_fromfile(path)
            if (oid.phase, oid.orbit_number, oid.target) == (
                p.phase,
                p.orbit_number,
                p.target,
            ):
                ccd = p.get_ccd()
                if ccd in bins:
                    continue
                else:
                    bins[ccd] = isis.getkey_k(path, "Instrument", "Summing")
        except ValueError:
            # Couldn't extract a Product ID from a cube, ignore.
            continue
        except AttributeError:
            # Product ID doesn't have a CCD, ignore.
            continue
    return bins


def set_flags(
    conf: dict, db: dict, ccdchan: tuple, bindex: int, newlogic=True
) -> collections.namedtuple:
    """Set various processing flags based on various configuration
    parameters."""
    HiCalFlags = collections.namedtuple(
        "HiCalFlags", ["noise_filter", "zapcols", "divide"]
    )

    noise_filter = False
    if (process_this(ccdchan, conf["HiCal_Noise_Processing"])) and (
        (
            float(db["IMAGE_DARK_STANDARD_DEVIATION"])
            >= float(conf["HiCal_Noise_Bin_DarkPixel_STD"][bindex])
        )
        or (
            float(db["CAL_MASK_STANDARD_DEVIATION"])
            >= float(conf["HiCal_Noise_Bin_Mask_STD"][bindex])
        )
        or (
            int(db["LOW_SATURATED_PIXELS"])
            >= int(conf["HiCal_Noise_LIS_Count"])
        )
    ):
        noise_filter = True

    zapcols = False
    if int(db["LOW_SATURATED_PIXELS"]) >= int(conf["HiCal_Noise_LIS_Count"]):
        zapcols = True
        if newlogic:
            noise_filter = True

    if "DIVIDE" == conf["HiCal_HPF_Cubenorm"]:
        divide = True
    elif "SUBTRACT" == conf["HiCal_HPF_Cubenorm"]:
        divide = False
    else:
        raise KeyError(
            "The HiCal_HPF_Cubenorm keyword can be DIVIDE or "
            "SUBTRACT, but was "
            "{}".format(conf["HiCal_HPF_Cubenorm"])
        )

    return HiCalFlags(noise_filter, zapcols, divide)


def set_lines(
    skip_top: int, skip_bottom: int, binning: int, image_lines: int
) -> collections.namedtuple:
    """Determine the right values for the start lines and number of lines."""
    Lines = collections.namedtuple("Lines", ["start", "number"])

    sl = skip_top / binning + 1
    nl = image_lines - (skip_top + skip_bottom) / binning
    if nl < 2000 / binning:
        sl = 1
        nl = image_lines

    # The original Perl code doesn't force sl and nl back into ints, but we do.
    return Lines(int(sl), int(nl))


def run_hical(
    in_cube: os.PathLike,
    hical_cub: os.PathLike,
    conf: dict,
    conf_path: os.PathLike,
    lis_per: float,
    image_buffer_mean: float,
    binning: int,
    noiseclean: bool,
    bitflipwidth=0,
    lis_tolerance=1.0,
    keep=False,
) -> str:

    to_d = isis.PathSet()
    next_cube = Path(in_cube)
    status = "Standard"

    to_s = "{}+SignedWord+{}:{}".format(
        hical_cub,
        conf["HiCal"]["HiCal_Normalization_Minimum"],
        conf["HiCal"]["HiCal_Normalization_Maximum"],
    )
    hical_args = {"to": to_s, "units": "IOF"}
    if conf["HiCal"]["HiCal_ISIS_Conf"] != "DEFAULT":
        dirs = (
            Path(conf_path).parent,
            Path(
                pkg_resources.resource_filename(
                    __name__,
                    'data/'
                )
            )
        )
        if lis_per < 5 and image_buffer_mean > 0:
            hical_args["conf"] = util.get_path(
                conf["HiCal"]["HiCal_ISIS_Conf"], dirs
            )
        else:
            hical_args["conf"] = util.get_path(
                conf["HiCal"]["HiCal_ISIS_Conf_Noise"], dirs
            )
            status = "BadCal"

    if noiseclean:
        mask_cube = to_d.add(next_cube.with_suffix(".mask.cub"))
        mask(
            next_cube,
            mask_cube,
            conf["NoiseFilter"]["NoiseFilter_Raw_Min"],
            conf["NoiseFilter"]["NoiseFilter_Raw_Max"],
            binning,
            width=bitflipwidth,
            keep=keep,
        )
        next_cube = mask_cube

    lisfix_cube = next_cube.with_suffix(".lisfix.cub")
    if (
        lisfix.fix(next_cube, lisfix_cube, lis_tolerance)
    ):
        to_d.add(lisfix_cube)
        next_cube = lisfix_cube

    isis.hical(next_cube, **hical_args)

    if not keep:
        to_d.unlink()
        cid = hirise.get_ChannelID_fromfile(hical_cub)
        Path(hical_cub).with_name(str(cid)).with_suffix(".hical.log").unlink()

    return status


def process_this(ccdchan: tuple, flag_list: list) -> int:
    if len(flag_list) != 28:
        raise IndexError("The list must have 28 elements")
    ccd_number = int(hirise.get_ccdnumber(ccdchan[0]))
    i = (2 * ccd_number) + int(ccdchan[1])
    return int(flag_list[i])


def furrow_nulling(
    cube: os.PathLike,
    out_cube: os.PathLike,
    binning: int,
    ccdchan: tuple,
    keep=False,
) -> bool:
    furrows_found = False

    out_path = Path(out_cube)

    to_del = isis.PathSet()
    fcrop_file = to_del.add(out_path.with_suffix(".fcrop.cub"))

    isis.mask(cube, mask=cube, to=out_path, minimum=1000000, maximum=1000000)

    furrow_values = furrow_setup(ccdchan[0], binning)
    chan_samp = chan_samp_setup(int(ccdchan[1]), binning)

    # Crop out the portion of the image that will not be furrow checked.
    # ^- that's what the original said, but this is really cropping out
    #    the part that will be furrow corrected.
    isis.crop(cube, to=fcrop_file, samp=chan_samp.ssamp, nsamp=chan_samp.nsamp)

    isis.handmos(
        fcrop_file,
        mosaic=out_path,
        insamp=1,
        outsamp=chan_samp.ssamp,
        create="no",
    )

    # for each column subject to furrowing, crop out each column,
    # run mask to null pixels above the furrow threshold, and
    # mosaic into the output image.
    for (s, furrow_v) in zip(chan_samp.samp, furrow_values):
        fscrop_file = to_del.add(out_path.with_suffix(f".crop{s}.cub"))
        isis.crop(cube, to=fscrop_file, samp=s, nsamp=1)

        fsmask_file = to_del.add(out_path.with_suffix(f".mask{s}.cub"))
        isis.mask(
            fscrop_file,
            mask=fscrop_file,
            to=fsmask_file,
            minimum=0,
            maximum=furrow_v,
        )

        isis.handmos(fsmask_file, mosaic=out_path, insamp=1, outsamp=s)

    else:
        # finally, check to see if any furrow pixels were zapped.
        # The original Perl code checked on the fist iteration, this
        # mechanism checks on the last, not sure if there was anything
        # magical about the first, or maybe we should do it each time?
        fscrop_stats = pvl.loads(isis.stats(fscrop_file).stdout)
        fsmask_stats = pvl.loads(isis.stats(fsmask_file).stdout)

        if int(fsmask_stats["Results"]["NullPixels"]) > int(
            fscrop_stats["Results"]["NullPixels"]
        ):
            furrows_found = True
            # Perform simple 3x3 LPFZ filter to clean up the edges
            # along the furrow
            trim_file = out_path.with_suffix(".trim.cub")
            isis.trimfilter(
                out_path, to=trim_file, lines=3, samples=3, minopt="COUNT"
            )
            trim_file.rename(out_path)

    if not keep:
        to_del.unlink()

    return furrows_found


def mask(
    in_cube: os.PathLike,
    out_cube: os.PathLike,
    noisefilter_min: float,
    noisefilter_max: float,
    binning: int,
    width=5,
    cid=None,
    keep=False,
) -> None:
    """mask out unwanted pixels"""
    logger.debug(mask.__doc__)
    to_del = isis.PathSet()
    out_path = Path(out_cube)
    temp_cube = to_del.add(out_path.with_suffix(".mask_temp.cub"))

    if cid is None:
        cid = hirise.get_ChannelID_fromfile(in_cube)

    isis.mask(
        in_cube,
        mask=in_cube,
        to=temp_cube,
        minimum=noisefilter_min,
        maximum=noisefilter_max,
        preserve="INSIDE",
        spixels="NONE",
    )
    if width == 0:
        cubenorm_stats_file = to_del.add(temp_cube.with_suffix(".cn.stats"))
        isis.cubenorm(temp_cube, stats=cubenorm_stats_file)
        (mindn, maxdn) = analyze_cubenorm_stats(
            cubenorm_stats_file, binning, cid=cid
        )
        isis.mask(
            temp_cube,
            mask=temp_cube,
            to=out_path,
            minimum=mindn,
            maximum=maxdn,
            preserve="INSIDE",
            spixels="NONE",
        )
    else:
        # The analyze_cubenorm_stats() function is meant to make sure we
        # don't blow away valid data, with the philosophy that it is better to
        # let in a little bad data in order to keep the good.  However, in
        # images with a high percentage of LIS%, there is so much bad that it
        # swamps the good, so *this* approach attempts to find the median
        # standard deviation that might be an approximation to the 'real'
        # standard devation of the data without the bit-flip noise.
        logger.info(
            f"{cid}: More severe handling of images with LIS pixels engaged."
        )
        bf.clean_cube(temp_cube, out_path, width=width, keep=keep)

    if not keep:
        to_del.unlink()

    return


def analyze_cubenorm_stats(
    statsfile: os.PathLike, binning: int, cid=None
) -> tuple:
    with open(statsfile) as csvfile:
        valid_points = list()
        std_devs = list()
        mins = list()
        maxs = list()
        reader = csv.DictReader(csvfile, dialect=isis.cubenormfile.Dialect)
        for row in reader:
            valid_points.append(int(row["ValidPoints"]))
            std_devs.append(float(row["StdDev"]))
            mins.append(int(row["Minimum"]))
            maxs.append(int(row["Maximum"]))

    if cid is None:
        cid = hirise.get_ChannelID_fromfile(statsfile)

    maxvp = max(valid_points)
    logger.info(f"{cid}: Maximum count of valid pixels: {maxvp}")

    # Original note:
    # # Get the median standard deviation value for all columns that have
    # # the maximum valid pixel count
    # #
    # # 2016-12-02 Note: this may not pick the median standard
    # # deviation value but the value from an index 0.95 times the number
    # # of entries in the sorted cubenorm statistics file.
    #
    # That seems to be exactly what it does.  I think the term 'median'
    # in the original (apparently pre-2016) comment is wrong.
    # The variable is called 'facstd' and when it is used below, it refers
    # to this being 'medstd + tol' which indicates that it really is
    # meant to be a factor above the median, which it is.
    std_w_maxvp = list()
    for (vp, std) in zip(valid_points, std_devs):
        if vp == maxvp:
            std_w_maxvp.append(std)

    std_w_maxvp.sort()
    facstd = std_w_maxvp[int((len(std_w_maxvp) - 1) * 0.95)]
    logger.info(
        f"{cid}: 95th percentile standard deviation of all "
        + "columns ({}) that have the ".format(len(std_w_maxvp))
        + "maximum valid pixel count: {}".format(facstd)
    )

    # Original note:
    # # find the minimum of minimums and the maximum of maximums for any
    # # column whose std is less than or equal to $medstd + $tol;
    min_w_maxvp = list()
    max_w_maxvp = list()
    for (vp, std, mi, ma) in zip(valid_points, std_devs, mins, maxs):
        if vp >= (maxvp * 0.9) and std < facstd:
            min_w_maxvp.append(mi)
            max_w_maxvp.append(ma)

    min_w_maxvp.sort()
    max_w_maxvp.sort()

    mindn = min_w_maxvp[int((len(min_w_maxvp) - 1) * 0.05)]
    maxdn = max_w_maxvp[int((len(max_w_maxvp) - 1) * 0.95)]

    if 1 == binning:
        mindn *= 0.7
        maxdn *= 1.3
    elif 2 == binning:
        mindn *= 0.6
        maxdn *= 1.4
    else:
        mindn *= 0.5
        maxdn *= 1.5

    return mindn, maxdn


def HiGainFx(
    cube: os.PathLike,
    outcube: os.PathLike,
    coef_path: os.PathLike,
    version: str,
    keep=False,
) -> None:
    """Perform a Gain-Drift correction on a HiRISE Channel Image."""
    # Eric Eliason (2019 Sept):
    # The HiCal pipeline is apparently double correcting for the
    # Gain Line Drift.   After looking at the source code of the HiCal
    # pipeline, hical program, HiGainFx and the various log files
    # sent to me, I too have come to the same conclusion.  My memory
    # from 10 years ago is starting to come back to me.  Here’s how
    # I remember things:

    # The HiGainFX perl script inserted into the HiCal pipeline was
    # always meant to be a stop-gap measure until a more thorough and
    # robust correction could be developed in the ISIS hical program.
    # Alan D. developed the algorithms for both the stop-gap and
    # full-correction models.  To further complicate things, there
    # were various other intermediary renditions that could be found
    # in the ISIS hicalbeta program but I don’t think hicalbeta was
    # ever used in the HiCal pipeline. I assume hicalbeta eventually
    # became hical.

    warnings.warn(
        "HiGainFx should not be used when Gain Drift correction is "
        "being used in ISIS hical, as HiGainFx is redundant and not "
        "as complete.",
        DeprecationWarning,
    )

    logger.info(HiGainFx.__doc__)
    binning = isis.getkey_k(cube, "Instrument", "Summing")
    ccd = isis.getkey_k(cube, "Instrument", "CcdId")
    chan = isis.getkey_k(cube, "Instrument", "ChannelNumber")

    coef_dir = Path(coef_path)

    if not coef_dir.exists() or not coef_dir.is_dir():
        coef_dir = Path(
            pkg_resources.resource_filename(
                __name__,
                'data/'
            )
        )
        logger.warning(
            "The HiGainFx coefficient directory {} could not be "
            "found, using {} instead.".format(coef_path, coef_dir)
        )

    coef_p = (
        Path(coef_dir)
        / f"HiRISE_Gain_Drift_Correction_Bin{binning}.{version}.csv"
    )

    with open(coef_p) as csvfile:
        reader = csv.DictReader(csvfile, skipinitialspace=True)
        for row in reader:
            if hirise.get_ccdchannel(row["CCD CH"]) == (ccd, chan):
                max_line = row["Max line"]
                a_coef = list()
                for r in (row["R(0)"], row["R(1)"], row["R(2)"]):
                    if r.startswith("-"):
                        # add another '-' for fx syntax
                        a_coef.append(f"-{r}")
                    else:
                        a_coef.append(r)

    eqn = r"\((F1/({0}+({1}*line)+({2}*line*line)))*".format(
        *a_coef
    ) + r"(line<{0}) + (F1*(line>={0})))".format(max_line)

    tfile = outcube.with_suffix(".tempfx.cub")
    isis.fx(f1=cube, to=tfile, mode="CUBES", equation=eqn)
    isis.specadd(tfile, to=outcube, match=cube)
    if not keep:
        tfile.unlink()
    return


def Cubenorm_Filter(
    cubenorm_tab: os.PathLike,
    outfile: os.PathLike,
    pause=False,
    boxfilter=5,
    divide=False,
    chan=None,
) -> tuple:
    """Returns a two-tuple with the standard deviation of the
    filtered average and a boolean indicating whether the maximum
    value of the ValidPoints in *cubenorm_tab* are at one end or
    the other (depending on *chan*).
    """
    logger.debug(
        "Perform a highpass filter on the cubenorm table output of the "
        "columnar average and median values."
    )
    if boxfilter < 3:
        raise ValueError(f"boxfilter={boxfilter} is less than 3")

    if chan is None:
        chan = int(hirise.get_ccdchannel(Path(cubenorm_tab).name)[1])

    # Make a list to receive each column
    valid_points = list()
    averages = list()
    medians = list()
    other_cols = list()
    with open(cubenorm_tab) as csvfile:
        reader = csv.DictReader(csvfile, dialect=isis.cubenormfile.Dialect)
        for row in reader:
            valid_points.append(int(row.pop("ValidPoints")))
            averages.append(float(row.pop("Average")))
            medians.append(float(row.pop("Median")))
            other_cols.append(row)

    zapped = FurrowCheck(valid_points, chan)

    avgflt = Cubenorm_Filter_filter(
        averages,
        boxfilter=boxfilter,
        iterations=50,
        chan=chan,
        pause=pause,
        vpoints=valid_points,
        divide=divide,
    )
    medflt = Cubenorm_Filter_filter(
        medians,
        boxfilter=boxfilter,
        iterations=50,
        chan=chan,
        pause=pause,
        vpoints=valid_points,
        divide=divide,
    )

    with open(outfile, "w") as csvfile:
        writer = isis.cubenormfile.DictWriter(csvfile)
        writer.writeheader()
        for (d, vp, av, md) in zip(other_cols, valid_points, avgflt, medflt):
            d["ValidPoints"] = str(vp)
            d["Average"] = "{:f}".format(av)
            d["Median"] = "{:f}".format(md)
            writer.writerow(d)

    # Calculate the standard deviation value for the filtered average:
    return statistics.stdev(avgflt), zapped


def cut_size(chan: int, length: int) -> collections.namedtuple:
    """Determine the appropriate values for the cut sizes."""
    Cut = collections.namedtuple("Cut", ["left", "right"])
    left = 6
    right = 6

    if 0 == chan:
        if 511 == length:
            left = 40
        elif 255 == length:
            left = 50
    elif 1 == chan:
        if 511 == length:
            right = 40
        elif 255 == length:
            right = 50
    else:
        raise ValueError(f"chan={chan} is not 0 or 1")

    return Cut(left, right)


def Cubenorm_Filter_filter_boxfilter(
    inlist: list, boxfilter: int, iterations=50
) -> list:
    x = copy.deepcopy(inlist)
    ori = copy.deepcopy(inlist)
    hwidth = int(boxfilter / 2)
    frac = [0.25, 0.125]
    for step in range(3):
        for it in range(iterations):
            xflt = copy.deepcopy(x)
            for i, _ in enumerate(x):
                start = i - hwidth
                if start < 0:
                    start = 0
                try:
                    xflt[i] = statistics.mean(
                        filter(lambda y: y != 0, x[start : i + hwidth + 1])
                    )
                except statistics.StatisticsError:
                    xflt[i] = 0
            x = copy.deepcopy(xflt)
        # Zap any columns that are different from the average by more then 25%
        if step < 2:
            for i, (orig, new) in enumerate(zip(ori, x)):
                if (
                    orig != 0
                    and new != 0
                    and (abs(orig - new) / new) > frac[step]
                ):
                    # print(f'orig: {orig}, new: {new}, setting to zero')
                    ori[i] = x[i] = 0
                else:
                    x[i] = orig
                # logger.info(
                # f"frac: {frac[step]}, index: {i}, orig: {orig}, new: {new}, "
                # f"stored: {x[i]}")
    return x


def Cubenorm_Filter_filter(
    inlist: list,
    boxfilter: int,
    iterations: int,
    chan: int,
    pause: bool,
    vpoints: list,
    divide: bool,
) -> list:
    """This performs highpass filtering on the passed list."""
    logger.debug(Cubenorm_Filter_filter.__doc__)

    x = copy.deepcopy(inlist)
    cut = cut_size(chan, len(x))

    # zap the left edge
    x[: cut.left] = [0] * cut.left

    # zap the right edge
    x[(len(x) - cut.right) :] = [0] * cut.right

    # zap the pause point pixels
    if pause and 1024 == len(x):
        for samp, width in zip(util.ch_pause[chan], util.ch_width[chan]):
            zap_slice = util.pause_slicer(samp, width)
            x[zap_slice] = [0] * abs(width)

    # boxfilter
    x = Cubenorm_Filter_filter_boxfilter(x, boxfilter, iterations)

    # Perform the highpass difference of divide the original from lowpass
    maxvp = max(vpoints)
    for i, (orig, vp) in enumerate(zip(inlist, vpoints)):
        if orig != 0 and x[i] != 0 and vp == maxvp:
            if divide:
                x[i] = orig / x[i]
            else:
                x[i] = orig - x[i]
        else:
            x[i] = None

    # Need to patch up any of those None values with neighboring values:
    first_i = (-1, 0)
    patch = (1, -1)

    if x[first_i[chan]] is None:
        x[first_i[chan]] = int(divide)

    if 0 == chan:
        e = reversed(list(enumerate(x)))
    else:
        e = list(enumerate(x))

    for (i, x_val) in e:
        if x_val is None:
            x[i] = x[i + patch[chan]]

    return x


def highlow_destripe(
    high_cube: os.PathLike,
    low_cube: os.PathLike,
    out_cube: os.PathLike,
    conf: dict,
    isisnorm="",
    lnull=True,
    lhrs=True,
    lhis=True,
    llrs=True,
    llis=True,
    keep=False,
) -> None:
    # Perform highpass/lowpass filter vertical destripping
    to_delete = isis.PathSet()
    lpf_cub = to_delete.add(out_cube.with_suffix(".lpf.cub"))
    isis.lowpass(
        low_cube,
        to=lpf_cub,
        line=conf["NoiseFilter_LPF_Line"],
        samp=conf["NoiseFilter_LPF_Samp"],
        minopt="PERCENT",
        replace="NULL",
        minimum=conf["NoiseFilter_LPF_Minper"],
        null=lnull,
        hrs=lhrs,
        his=lhis,
        lrs=llrs,
        lis=llis,
    )

    hpf_cub = to_delete.add(out_cube.with_suffix(".hpf.cub"))
    isis.highpass(
        high_cube,
        to=hpf_cub,
        minopt="PERCENT",
        line=conf["NoiseFilter_HPF_Line"],
        samp=conf["NoiseFilter_HPF_Samp"],
        minimum=conf["NoiseFilter_HPF_Minper"],
    )

    isis.algebra(
        from_=lpf_cub,
        from2=hpf_cub,
        to=out_cube.with_suffix(".cub" + isisnorm),
        operator="ADD",
    )
    if not keep:
        to_delete.unlink()

    return


def getHistVal(histogram: isis.Histogram, conf: dict) -> tuple:
    """Return information about the histogram."""
    lisper = 0
    maxval = None
    if int(histogram["Total Pixels"]) - int(histogram["Null Pixels"]) > 0:
        lisper = (
            int(histogram["Lis Pixels"])
            / (int(histogram["Total Pixels"]) - int(histogram["Null Pixels"]))
            * 100
        )
    cumper = conf["NoiseFilter_HighEnd_Percent"]
    if cumper < 99.0:
        cumper = 99.0

    if lisper > conf["NoiseFilter_Hard_Tolmax"]:
        hard_high_end = conf["NoiseFilter_Hard_HighEnd_Percent"]
        if hard_high_end < 99.9:
            hard_high_end = 99.9
        cumper = hard_high_end

    for row in histogram:
        if float(row.CumulativePercent) > cumper:
            if tuple(isisversion.version_info()[:3]) < (4, 3, 0):
                maxval = float(row.DN)
            else:
                maxval = float(row.MaxExclusive)
            break

    if maxval is None:
        raise ValueError(
            "Did not find a CumulativePercent value greater "
            "than {}".format(cumper)
        )

    return lisper, maxval


def NoiseFilter_noisefilter(
    from_cube: os.PathLike,
    to_cube: os.PathLike,
    flattol: float,
    conf: dict,
    maxval: float,
    tolmin: float,
    tolmax: float,
) -> None:
    """Convenience function for the repeated noisefiltering."""
    isis.noisefilter(
        from_cube,
        to=to_cube,
        flattol=flattol,
        low=conf["NoiseFilter_Minimum_Value"],
        high=maxval,
        tolmin=tolmin,
        tolmax=tolmax,
        sample=conf["NoiseFilter_Noise_Samp"],
        line=conf["NoiseFilter_Noise_Line"],
        toldef="STDDEV",
        replace="NULL",
        lisisnoise=True,
        lrsisnoise=True,
    )


def NoiseFilter_cubenorm_edit(
    in_tab: os.PathLike,
    out_tab: os.PathLike,
    out_tab2: os.PathLike,
    chan: int,
    binning: int,
    conf: dict,
    zapc=False,
) -> None:
    """This function returns nothing, but creates a pair of edited cubenorm
    files.
    """
    # Slightly different values from other function, not entirely sure why.
    # Pause point locations are 1-based pixel numbers, so -1 to get list index.
    # Width values are the number of pixels to affect, including the pause
    # point pixel.
    logger.debug(
        "Zaps the relevant pixels in the cubenorm output and creates an "
        "edited cubenorm file."
    )

    # pause point sample locations:
    ch_pause = {
        0: (1, 252, 515, 778),  # Channel 0
        1: (247, 510, 773, 1024),
    }  # Channel 1

    # Number of pixels to cut from pause point, the sign indicates
    # the direction of cut from the pause point.
    ch_width = {0: (3, 6, 6, 6), 1: (-8, -7, -6, -3)}

    vpnts = list()
    other_cols = list()
    with open(in_tab) as csvfile:
        reader = csv.DictReader(csvfile, dialect=isis.cubenormfile.Dialect)
        for row in reader:
            vpnts.append(int(row.pop("ValidPoints")))
            other_cols.append(row)

    max_vpnts = max(vpnts)
    if max_vpnts <= 0:
        max_vpnts = 1

    # Create a 'unity array'
    # Zap any columns with less then NoiseFilter_Zap_Fraction
    norm = [1] * len(vpnts)
    for i, v in enumerate(vpnts):
        if zapc and v / max_vpnts < float(conf["NoiseFilter_Zap_Fraction"]):
            norm[i] = 0

    # Copy it.
    norm2 = copy.deepcopy(norm)

    # Determine if the pause point pixels need to be zapped
    if binning == 1:
        trigger = NoiseFilter_zaptrigger(
            vpnts,
            max_vpnts,
            conf["NoiseFilter_Nonvalid_Fraction"],
            ch_pause[chan],
            ch_width[chan],
        )

        for samp, width in zip(ch_pause[chan], ch_width[chan]):
            zap_slice = util.pause_slicer(samp, width)
            for i in range(zap_slice.start, zap_slice.stop):
                if trigger:
                    norm[i] = 0
                norm2[i] = 0

    NoiseFilter_cubenorm_writer(out_tab, other_cols, vpnts, norm)
    NoiseFilter_cubenorm_writer(out_tab2, other_cols, vpnts, norm2)

    return


def NoiseFilter_zaptrigger(
    vpnts, max_vpnts, nonvalid_frac, ch_pause, ch_width
) -> bool:
    trigger = False
    for samp, width in zip(ch_pause, ch_width):
        zap_slice = util.pause_slicer(samp, width)
        for i in range(zap_slice.start, zap_slice.stop):
            if vpnts[i] / max_vpnts < nonvalid_frac:
                trigger = True
    if trigger:
        logger.info(f"Fraction of non-valid pixels > {nonvalid_frac}")
        logger.info(f"Pause point pixels will be zapped.")
    return trigger


def NoiseFilter_cubenorm_writer(
    file_path: os.PathLike, other_cols: list, vpnts: list, norm: list
):
    with open(file_path, "w") as csvfile:
        writer = isis.cubenormfile.DictWriter(csvfile)
        writer.writeheader()
        for (d, vp, n) in zip(other_cols, vpnts, norm):
            d["ValidPoints"] = vp
            d["Average"] = n
            d["Median"] = n
            d["StdDev"] = n
            d["Minimum"] = n
            d["Maximum"] = n
            writer.writerow(d)


def NoiseFilter(
    in_cube: os.PathLike,
    output: os.PathLike,
    conf: dict,
    minimum=None,
    maximum=None,
    zapc=False,
    keep=False,
) -> None:
    """NoiseFilter: Perform salt/pepper noise removal."""
    logger.debug(NoiseFilter.__doc__)
    binning = int(isis.getkey_k(in_cube, "Instrument", "Summing"))
    (ccd, chan) = hirise.get_ccdchannel(
        isis.getkey_k(in_cube, "Archive", "ProductId")
    )
    isisnorm = ""
    if minimum is not None and maximum is not None:
        isisnorm = "+SignedWord+{}:{}".format(minimum, maximum)

    to_delete = isis.PathSet()

    h = isis.Histogram(in_cube)
    (LisP, MaxVal) = getHistVal(h, conf)

    cn_tab = to_delete.add(output.with_suffix(".cn.tab"))
    isis.cubenorm(in_cube, stats=cn_tab, format_="TABLE", direction="COLUMN")

    cn2_tab = to_delete.add(output.with_suffix(".cn2.tab"))
    cn3_tab = to_delete.add(output.with_suffix(".cn3.tab"))
    NoiseFilter_cubenorm_edit(
        cn_tab, cn2_tab, cn3_tab, int(chan), binning, conf, zapc
    )

    # Zap the bad columns for the highpass filter
    zap_cub = to_delete.add(output.with_suffix(".zap.cub"))
    isis.cubenorm(
        in_cube,
        to=zap_cub,
        fromstats=cn2_tab,
        statsource="TABLE",
        mode="DIVIDE",
        norm="AVE",
        preserve="FALSE",
    )

    # Zap the bad columns for the lowpass filter
    zap2_cub = to_delete.add(output.with_suffix(".zap2.cub"))
    isis.cubenorm(
        in_cube,
        to=zap2_cub,
        fromstats=cn3_tab,
        statsource="TABLE",
        mode="DIVIDE",
        norm="AVE",
        preserve="FALSE",
    )

    # Perform highpass/lowpass filter vertical destripping
    add_cub = to_delete.add(output.with_suffix(".add.cub"))
    highlow_destripe(
        zap_cub, zap2_cub, add_cub, conf, isisnorm, llis=False, keep=keep
    )

    # Perform the 1st noise filter
    tolmin = float(conf["NoiseFilter_Tolmin"])
    tolmax = float(conf["NoiseFilter_Tolmax"])
    if LisP >= float(conf["NoiseFilter_Hard_Filtering"]):
        tolmin = float(conf["NoiseFilter_Hard_Tolmin"])
        tolmax = float(conf["NoiseFilter_Hard_Tolmax"])
    flattol = float(h["Std Deviation"]) * float(conf["NoiseFilter_Flattol"])
    if flattol < 0.00001:
        flattol = 0.00001

    nf1_cub = to_delete.add(output.with_suffix(".nf1.cub"))
    NoiseFilter_noisefilter(
        add_cub,
        nf1_cub.with_suffix(".cub" + isisnorm),
        flattol,
        conf,
        MaxVal,
        tolmin,
        tolmax,
    )

    # Perform the 2nd noise filter
    nf2_cub = to_delete.add(output.with_suffix(".nf2.cub"))
    NoiseFilter_noisefilter(
        nf1_cub,
        nf2_cub.with_suffix(".cub" + isisnorm),
        flattol,
        conf,
        MaxVal,
        tolmin,
        tolmax,
    )

    # Perform the 3rd noise filter
    nf3_cub = to_delete.add(output.with_suffix(".nf3.cub"))
    NoiseFilter_noisefilter(
        nf2_cub,
        nf3_cub.with_suffix(".cub" + isisnorm),
        flattol,
        conf,
        MaxVal,
        tolmin,
        tolmax,
    )

    # Perform another highpass/lowpass filter now that the
    # data are much cleaner
    add2_cub = to_delete.add(output.with_suffix(".add2.cub"))
    highlow_destripe(
        nf3_cub,
        nf3_cub,
        add2_cub,
        conf,
        isisnorm,
        lnull=False,
        lhrs=False,
        lhis=False,
        llrs=False,
        llis=False,
        keep=keep,
    )

    if ccd.startswith("RED"):
        # Perform LPFZ  filters if we have a RED filter image.
        # For IR and BG filter data, assume that the HiColorNorm pipeline
        # step will interpolate using the BG/RED and IR/RED ratio data.
        lowmin = int(
            int(conf["NoiseFilter_LPFZ_Line"])
            * int(conf["NoiseFilter_LPFZ_Samp"])
            / 3
        )
        lpfz_cub = to_delete.add(output.with_suffix(".lpfz.cub"))
        isis.lowpass(
            add2_cub,
            to=lpfz_cub.with_suffix(".cub" + isisnorm),
            sample=3,
            line=3,
            minopt="COUNT",
            minimum=1,
            filter_="OUTSIDE",
            null=True,
            hrs=False,
            his=True,
            lrs=True,
            lis=True,
        )
        isis.lowpass(
            lpfz_cub,
            to=output.with_suffix(".cub" + isisnorm),
            sample=conf["NoiseFilter_LPFZ_Samp"],
            line=conf["NoiseFilter_LPFZ_Line"],
            minopt="COUNT",
            minimum=lowmin,
            filter_="OUTSIDE",
            null=True,
            hrs=False,
            his=True,
            lrs=True,
            lis=True,
        )
    else:
        to_delete.remove(add2_cub)
        add2_cub.rename(output)

    if not keep:
        to_delete.unlink()
    return


def Hidestripe(
    in_cube: os.PathLike,
    out_cube: os.PathLike,
    binning: int,
    minimum: float,
    maximum: float,
    hidcorr: str,
    line_samples: int,
    keep=False,
) -> float:
    # SignedWord+$HiCal_Normalization_Minimum:$HiCal_Normalization_Maximum
    to_s = "+SignedWord+{}:{}".format(minimum, maximum)
    to_del = isis.PathSet()
    if 1 == binning:
        temp_cub = to_del.add(out_cube.with_suffix(".hd.cub"))
        isis.hidestripe(
            in_cube,
            to=str(temp_cub) + to_s,
            parity="EVEN",
            correction=hidcorr,
        )
        isis.hidestripe(
            temp_cub,
            to=str(out_cube) + to_s,
            parity="ODD",
            correction=hidcorr,
        )
    else:
        lpf_cube = to_del.add(out_cube.with_suffix(".lpf.cub"))
        hpf_cube = to_del.add(out_cube.with_suffix(".hpf.cub"))
        boxsamp = (2 * line_samples) - 1

        isis.lowpass(
            in_cube,
            to=lpf_cube,
            samples=boxsamp,
            lines=3,
            null=False,
            hrs=False,
            his=False,
            lrs=False,
            lis=False,
        )
        isis.highpass(
            in_cube, to=hpf_cube, samples=boxsamp, lines=1, propagate=True
        )
        isis.algebra(
            from_=lpf_cube,
            from2=hpf_cube,
            to=str(out_cube) + to_s,
            operator="ADD",
            a="1.0",
            b="1.0",
        )

    # Standard deviation of difference between cube after noise
    # filter and cube after hidestripe
    diff_cube = to_del.add(out_cube.with_suffix(".diff.cub"))
    isis.algebra(
        from_=out_cube, from2=in_cube, to=diff_cube, operator="subtract"
    )
    diff_pvl = pvl.loads(isis.stats(diff_cube).stdout)["Results"]
    try:
        stddev = float(diff_pvl["StandardDeviation"])
    except (KeyError, ValueError) as err:
        if int(diff_pvl["ValidPixels"]) <= 0:
            raise KeyError(
                "There is no StandardDeviation computed from "
                f"{diff_cube} because there were no ValidPixels. "
                f"Probably because {in_cube} has values outside "
                f"the range: {minimum} to {maximum}."
            )
        else:
            raise err

    if not keep:
        to_del.unlink()

    return stddev


def chan_samp_setup(channel: int, binning: int) -> collections.namedtuple:
    """Returns a named tuple which contains a list and two numbers."""

    if not isinstance(channel, int):
        raise TypeError(
            "channel should be an int, but it "
            "is: {} {}".format(channel, type(channel))
        )

    if not isinstance(binning, int):
        raise TypeError(
            "binning should be an int, but it "
            "is: {} {}".format(binning, type(binning))
        )

    samp = collections.defaultdict(dict)
    # samp[chan][binning]
    samp[0][2] = 1, 2, 3, 4, 5, 6, 7, 8, 9, 10
    samp[1][2] = 512, 511, 510, 509, 508, 507, 506, 505, 504, 503
    samp[0][4] = 1, 2, 3, 4, 5, 6
    samp[1][4] = 256, 255, 254, 253, 252, 251

    ssamp = collections.defaultdict(dict)
    # ssamp[chan][binning]
    ssamp[0][2] = 11
    ssamp[1][2] = 1
    ssamp[0][4] = 7
    ssamp[1][4] = 1

    nsamp = collections.defaultdict(dict)
    # nsamp[chan][binning]
    nsamp[0][2] = 502
    nsamp[1][2] = 502
    nsamp[0][4] = 250
    nsamp[1][4] = 250

    ChanSamp = collections.namedtuple("ChanSamp", ["samp", "ssamp", "nsamp"])

    c = ChanSamp(
        samp[channel][binning],
        ssamp[channel][binning],
        nsamp[channel][binning],
    )
    return c


def furrow_setup(ccd: str, binning: int):
    """Returns the right tuple of furrow."""
    # Expect to call the returned dict like this: d[ccd][binning]
    # Assume each channel for each CCD will have the same threshold
    d = collections.defaultdict(dict)
    d["RED0"][2] = (
        8000,
        8100,
        8700,
        9200,
        9600,
        10000,
        12000,
        12000,
        12000,
        12000,
    )
    d["RED0"][4] = (8000, 9000, 9500, 9900, 9900, 10000)
    d["RED1"][2] = (
        7200,
        7200,
        7800,
        8400,
        9000,
        9500,
        12000,
        12000,
        12000,
        12000,
    )
    d["RED1"][4] = (8000, 8100, 9200, 9600, 9800, 10000)
    d["RED2"][2] = (
        7800,
        7800,
        8400,
        9000,
        9600,
        10000,
        12000,
        12000,
        12000,
        12000,
    )
    d["RED2"][4] = (8000, 8700, 9500, 9800, 9900, 10000)
    d["RED3"][2] = (
        7800,
        8100,
        8300,
        9200,
        9600,
        10000,
        12000,
        12000,
        12000,
        12000,
    )
    d["RED3"][4] = (7900, 9200, 9700, 9900, 10000, 10500)
    d["RED4"][2] = (
        7800,
        7800,
        8300,
        9000,
        9500,
        9900,
        12000,
        12000,
        12000,
        12000,
    )
    d["RED4"][4] = (8000, 8700, 9700, 10000, 10300, 10600)
    d["RED5"][2] = (
        7900,
        8200,
        8600,
        9200,
        9600,
        10000,
        12000,
        12000,
        12000,
        12000,
    )
    d["RED5"][4] = (8000, 9300, 9700, 9900, 10200, 10700)
    d["RED6"][2] = (
        7500,
        7500,
        8100,
        8500,
        9200,
        10000,
        12000,
        12000,
        12000,
        12000,
    )
    d["RED6"][4] = (8000, 8400, 9700, 10000, 10500, 10700)
    d["RED7"][2] = (
        7600,
        8300,
        8900,
        9400,
        9900,
        11000,
        12000,
        12000,
        12000,
        12000,
    )
    d["RED7"][4] = (7700, 9600, 10000, 10200, 11000, 12000)
    d["RED8"][2] = (
        7200,
        7200,
        7900,
        8500,
        9000,
        9400,
        12000,
        12000,
        12000,
        12000,
    )
    d["RED8"][4] = d["RED7"][4]
    d["RED9"][2] = (
        7600,
        8300,
        8600,
        9200,
        9600,
        10000,
        12000,
        12000,
        12000,
        12000,
    )
    d["RED9"][4] = (8000, 8800, 9200, 9400, 9800, 10500)
    d["IR10"][2] = d["RED0"][2]
    d["IR10"][4] = (7600, 8300, 9000, 10000, 10500, 12000)
    d["IR11"][2] = d["RED0"][2]
    d["IR11"][4] = d["IR10"][4]
    d["BG12"][2] = d["RED0"][2]
    d["BG12"][4] = d["IR10"][4]
    d["BG13"][2] = d["RED0"][2]
    d["BG13"][4] = d["IR10"][4]

    return d[ccd][binning]


if __name__ == "__main__":
    main()
