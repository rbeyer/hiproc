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
#
# Parts of this were inspired by fix_74.pro, find_zero2REV.pro
# by Alan Delamere, Dec 2020, but the implementation
# here was written from scratch.


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
    tolerance=0.4
):
    """The ISIS cube file at *out_path* will be the result of running lis-fix
    processing of the file at *in-path*.

    If the fraction of LIS pixels in the Calibration Buffer Pixel area or
    the first 20 lines of Dark Pixels is greater than *tolerance* a ValueError
    is raised.

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

    #############################
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
    caldark_vals = np.array(hca_dict["DarkPixels"])
    caldark = np.ma.masked_outside(
        caldark_vals, specialpix.Min, specialpix.Max
    )
    # print(caldark.shape)
    dark_slopes = np.ma.apply_along_axis(
        get_ramp_slope, 0, caldark[:,1:], mask_lines
    )
    # print(dark_slopes.shape)
    # print(dark_slopes)

    # CORRECT FOR MASKED LINES
    # ASSUME SLOPE IS DEFINED BY DARK COLUMNS
    dark_slope_mean = np.mean(dark_slopes)
    logger.info(f"Dark Ramp mean slope: {dark_slope_mean}")
    # FOR # J = 12, SZ(1) - 1 # DO
    # zero.line(j) = zero.line(j) - SL_D * (
    #         ystart(j) + 1 + 20. / info.bin * 103. / 89)
    zero_correction = dark_slope_mean * (20 / binning * 103 / 89)

    revclk_lisfrac = np.ma.count_masked(
        cal_image[:20, :]
    ) / cal_image[:20, :].size
    logger.info(
        f"Fraction of LIS pixels in Reverse-Clock area: {revclk_lisfrac}"
    )

    logger.info("Fixing Reverse-Clock pixels.")
    fixed_cal = np.ma.apply_along_axis(
        fix_rev_clock, 0, cal_image, mask_lines, zero_correction
    )

    logger.info("Writing out Reverse-Clock pixels.")
    hci_dict["Calibration"] = apply_special_pixels(
        fixed_cal, specialpix
    ).data.tolist()
    isis.cube.overwrite_table(out_p, t_name, hci_dict)

    ######################
    # Fix the buffer area:
    calbuf_vals = np.array(hca_dict["BufferPixels"])
    calbuf = np.ma.masked_outside(
        calbuf_vals, specialpix.Min, specialpix.Max
    )

    ha_dict = isis.cube.get_table(in_path, "HiRISE Ancillary")
    dark_vals = np.array(ha_dict["DarkPixels"])
    dark = np.ma.masked_outside(
        dark_vals, specialpix.Min, specialpix.Max
    )
    # print(dark.shape)
    buffer_vals = np.array(ha_dict["BufferPixels"])
    buffer = np.ma.masked_outside(
        buffer_vals, specialpix.Min, specialpix.Max
    )
    # print(buffer)
    # print(buffer.dtype)
    # print(buffer.shape)

    # first_im_line = 19+(20+label["IsisCube"]["Instrument"]["Tdi"])/binning

    if np.ma.count_masked(calbuf) / calbuf.size <= tolerance:
        raise ValueError(
            "Less than 40% of the Calibration Buffer Pixels have a "
            "real value, using zero as the median."
        )
    else:
        refr = int(np.ma.median(calbuf))
        logger.info(f"Median of Calibration Buffer Pixels: {refr}")

    if np.ma.count_masked(dark[:20]) / dark[:20].size <= tolerance:
        raise ValueError(
            "Less than 40% of the first 20 lines of Dark Pixels have a "
            "real value, using zero as the median."
        )
    else:
        refd = int(np.ma.median(dark[:20]))
        logger.info(f"Median of first 20 lines of Dark Pixels: {refd}")

    model_buffer = dark[:, 2:14] - refd + refr
    # print(model_buffer)
    # print(model_buffer.dtype)
    # print(model_buffer.shape)

    logger.info("Fixing Buffer pixels.")
    fixed_buf = np.ma.apply_along_axis(
        fix_buffer,
        1,
        np.ma.concatenate((buffer, model_buffer), 1),
        buffer.shape[1]
    )
    # print(fixed_buf)
    # print(fixed_buf.dtype)
    # print(fixed_buf.shape)

    ha_dict["BufferPixels"] = apply_special_pixels(
        fixed_buf, specialpix
    ).data.tolist()
    logger.info("Writing out Buffer pixels.")
    isis.cube.overwrite_table(out_p, "HiRISE Ancillary", ha_dict)

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
    """Returns a np.masked array where the LIS pixels in the
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


def fix_buffer(
    buf_row: np.ma.array,
    stop,
    lis_fraction=0.5,
):
    """Returns a np.masked array the length of *buf_row[:stop]* with those
    same values, but with the masked pixels in the [:stop] part the row
    replaced with pixels from the [stop:] part, if the percent of masked
    pixels in [stop:] is greater than or equal to *lis_fraction*."""
    if np.ma.count_masked(buf_row[:stop]) / stop >= lis_fraction:
        # print(f"row: {buf_row[:stop]} {buf_row[stop:]}")
        # print(buf_row.data[:stop])
        for i in range(stop):
            # print(f"masks: {buf_row.mask[i]}, {buf_row.mask[i + stop]}")
            # print(f"mask type: {type(buf_row.mask[i])} ")
            # print(f"values: {buf_row[i]}, {buf_row[i + stop]}")
            # The mask values are really np.bool_ values, so compare
            # to ints, not Python True, False objects.
            if buf_row.mask[i] == 1 and buf_row.mask[i + stop] == 0:
                buf_row[i] = buf_row[i + stop]
            # print(f"values: {buf_row[i]}, {buf_row[i + stop]}")
        # print(buf_row.data[:stop])
        # print(buf_row[:stop])

    return buf_row[:stop]


if __name__ == "__main__":
    main()
