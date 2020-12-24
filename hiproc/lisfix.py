#!/usr/bin/env python
"""Converts LIS values in reverse-clock and buffer areas to
derived values.

For a variety of reasons, but highlighted with the August 2020
change of ADC settings, the HiRISE LUT process can clip the DN
histogram resulting in good data in the image area, but mostly LIS
values in the reverse-clock and buffer areas, which are crucial for
calibration.
"""

# Copyright 2020, Ross A. Beyer (rbeyer@seti.org)
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


import argparse
import logging
import os
import shutil
import subprocess
import sys
import traceback
from pathlib import Path

import numpy as np

import pvl
import kalasiris as isis

import hiproc.util as util
from hiproc.bitflips import apply_special_pixels

logger = logging.getLogger(__name__)


def main():
    try:
        parser = argparse.ArgumentParser(
            description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter,
            parents=[util.parent_parser()],
        )
        parser.add_argument(
            "-o", "--output", required=False, default=".lisfix.cub"
        )
        parser.add_argument("file", help="ISIS Cube file.")

        args = parser.parse_args()

        util.set_logger(logger, args.log, args.logfile)

        out_p = util.path_w_suffix(args.output, args.file)

        fix(
            args.file,
            out_p,
        )
        sys.exit(0)
    except subprocess.CalledProcessError as err:
        print("Had an ISIS error:", file=sys.stderr)
        print(" ".join(err.cmd), file=sys.stderr)
        print(err.stdout, file=sys.stderr)
        print(err.stderr, file=sys.stderr)
        sys.exit(1)
    except Exception as err:
        traceback.print_exc(file=sys.stderr)
        print(err, file=sys.stderr)
        sys.exit(1)


def fix(
    in_path: os.PathLike,
    out_path: os.PathLike,
):
    """The ISIS cube file at *out_path* will be the result of running lis-fix
    processing of the file at *in-path*.

    This function anticipates a HiRISE cube that has the following
    table objects in its label: HiRISE Calibration Ancillary,
    HiRISE Calibration Image, HiRISE Ancillary.

    The HiRISE Calibration Ancillary table contains the BufferPixels
    and DarkPixels from either side of the HiRISE Calibration Image.
    The HiRISE Calibration Image table contains the Reverse-Clock,
    Mask, and Ramp Image areas.  The HiRISE Ancillary table contains
    the BufferPixels and DarkPixels from either side of the main
    Image Area.
    """

    in_p = Path(in_path)
    out_p = Path(out_path)

    label = pvl.load(in_p)
    if "IsisCube" not in label:
        raise ValueError(f"The file at {in_p} is not an ISIS Cube.")

    binning = label["IsisCube"]["Instrument"]["Summing"]
    specialpix = getattr(
        isis.specialpixels, label["IsisCube"]["Core"]["Pixels"]["Type"]
    )

    shutil.copy(in_path, out_path)

    # Fix the reverse-clock area:
    t_name = "HiRISE Calibration Image"
    hci_dict = isis.cube.get_table(in_path, t_name)
    cal_vals = np.array(hci_dict["Calibration"])

    cal_image = np.ma.masked_outside(
        cal_vals, specialpix.Min, specialpix.Max
    )

    mask_lines = int(20 / binning)

    # Get the average slope of the dark ramp
    hca_dict = isis.cube.get_table(in_path, "HiRISE Calibration Ancillary")
    calbufdark_vals = np.array(hca_dict["DarkPixels"])
    calbufdark = np.ma.masked_outside(
        calbufdark_vals, specialpix.Min, specialpix.Max
    )
    dark_slopes = np.ma.apply_along_axis(
        get_ramp_slope, 0, calbufdark, mask_lines
    )

    # CORRECT FOR MASKED LINES
    # ASSUME SLOPE IS DEFINED BY DARK COLUMNS
    dark_slope_mean = np.mean(dark_slopes)
    # FOR # J = 12, SZ(1) - 1 # DO
    # zero.line(j) = zero.line(j) - SL_D * (
    #         ystart(j) + 1 + 20. / info.bin * 103. / 89)
    zero_correction = dark_slope_mean * (1 + 20 / binning * 103 / 89)

    revclk_lisfrac = np.ma.count_masked(
        cal_image[:20, :]
    ) / cal_image[:20, :].size
    logger.info(
        f"Fraction of LIS pixels in Reverse-Clock area: {revclk_lisfrac}"
    )

    # Fix the reverse-clock lines.
    fixed_cal = np.ma.apply_along_axis(
        fix_rev_clock, 0, cal_image, mask_lines, zero_correction
    )

    # write the table back out
    hci_dict["Calibration"] = apply_special_pixels(
        fixed_cal, specialpix
    ).data.tolist()
    isis.cube.overwrite_table(out_p, t_name, hci_dict)

    # if any((buffer_area, dark_area)):
    #     t_name = 'HiRISE Ancillary'
    #     HA_dict = isis.cube.get_table(cube, t_name)
    #     buffer_image = np.array(HA_dict['BufferPixels'])
    #     dark_image = np.array(HA_dict['DarkPixels'])
    #
    #     clean_buffer = clean_buffer_table()
    #     clean_dark = clean_dark_table()
    #
    #     # write the table back out into outcube
    #     if not dryrun:
    #         HA_dict['BufferPixels'] = mask_listoflists(
    #             HA_dict['BufferPixels'], buffer_image, specialpix)
    #         HA_dict['DarkPixels'] = mask_listoflists(
    #             HA_dict['DarkPixels'], buffer_image, specialpix)
    #         isis.cube.overwrite_table(outcube, t_name, HA_dict)

    return


def get_ramp_slope(
    col: np.ma.array,
    mask_lines: int,
):
    """Returns the slope from the ramp pixels."""
    fit = np.ma.polyfit(
            np.arange(col[20 + mask_lines:].size),
            col[20 + mask_lines:],
            1
        )
    return fit[0]


def fix_rev_clock(
    cal_image_col: np.ma.array,
    mask_lines: int,
    zero_correction=0,
    lis_fraction=0.2,
):
    """Returns a np.masked_masked array where the LIS pixels in the
    reverse-clock area are replaced with a fit from the ramp area, if
    the percent of LIS pixels is greater than or equal to *lis_fraction*.

    The size of the mask is variable, based on binning, so the number of
    *mask_lines* must be provided.
    """
    if np.ma.count_masked(cal_image_col[:20]) / 20 >= lis_fraction:
        # Attempt to derive value from the ramp via a linear fit
        fit = np.ma.polyfit(
            np.arange(cal_image_col[20 + mask_lines:].size),
            cal_image_col[20 + mask_lines:],
            1
        )

        cal_image_col[:20] = cal_image_col[:20].filled(
            fill_value=fit[1] - zero_correction
        )
    return cal_image_col


if __name__ == "__main__":
    main()
