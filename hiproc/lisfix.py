#!/usr/bin/env python
"""Converts LIS values in reverse-clock and buffer areas to
derived values.

For a variety of reasons, but highlighted with the August 2020
change of ADC settings, the HiRISE LUT process can clip the DN
histogram resulting in good data in the image area, but mostly LIS
values in the reverse-clock and buffer areas, which are crucial for
calibration.  This program attempts to replace "missing" data in those
areas with values derived from the ramp area, the dark area, and models
of expected reverse-clock medians.
"""

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
#
# Parts of this were inspired by fix_74.pro, find_zero2REV.pro (Dec 2020),
# and fix_74_all_new.pro (Feb 2021) by Alan Delamere, but the implementation
# here was written from scratch.


import argparse
import logging
import os
import shutil
import sys
from collections import abc
from pathlib import Path

import numpy as np
import scipy.signal as signal

import pvl
import kalasiris as isis

import hiproc.util as util
from hiproc.bitflips import apply_special_pixels
from hiproc.hirise import get_ChannelID_fromfile

logger = logging.getLogger(__name__)


def arg_parser():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        parents=[util.parent_parser()],
    )
    parser.add_argument(
        "-o", "--output", required=False, default=".lisfix.cub"
    )
    parser.add_argument("file", help="ISIS Cube file.")
    return parser


def main():
    args = arg_parser().parse_args()

    util.set_logger(args.verbose, args.logfile, args.log)

    out_p = util.path_w_suffix(args.output, args.file)

    with util.main_exceptions(args.verbose):
        fixed = fix(
            args.file,
            out_p,
        )

    if not fixed:
        sys.exit(f"{args.file} did not need lisfix cleaning.")

    return


def fix(
    in_path: os.PathLike,
    out_path: os.PathLike,
    tolerance=0.4
) -> bool:
    """The ISIS cube file at *out_path* will be the result of running lis-fix
    processing of the file at *in-path*, if allowed by *tolerance*.  If the
    *tolerance* setting would result in the application of lis-fix, then
    no file will be created at *out_path* and this function will return
    False.  Otherwise there will be a file there, and this function
    will return True.

    If the fraction of LIS pixels in the Reverse-Clock areas is
    greater than *tolerance*, then all LIS pixels in the rev-clock
    will be converted to a real DN value based on modeling the slope
    of the ramp area.

    If the fraction of LIS pixels in the Buffer Pixel area is greater
    than *tolerance*, this algorithm will attempt to fix those LIS pixels,
    as well.

    If the file at *in_path* isn't an ISIS cube, a ValueError will be raised.

    This function anticipates a HiRISE cube that has the following
    table objects in its label: HiRISE Calibration Ancillary,
    HiRISE Calibration Image, HiRISE Ancillary.  If they are not present
    a KeyError will be raised.

    The HiRISE Calibration Ancillary table contains the BufferPixels
    and DarkPixels from either side of the HiRISE Calibration Image.
    The HiRISE Calibration Image table contains the Reverse-Clock,
    Mask, and Ramp Image areas.  The HiRISE Ancillary table contains
    the BufferPixels and DarkPixels from either side of the main
    Image Area.
    """

    in_p = Path(in_path)
    out_p = Path(out_path)
    fixed = False

    label = pvl.load(in_p)
    if "IsisCube" not in label:
        raise ValueError(f"The file at {in_p} is not an ISIS Cube.")

    binning = label["IsisCube"]["Instrument"]["Summing"]
    specialpix = getattr(
        isis.specialpixels, label["IsisCube"]["Core"]["Pixels"]["Type"]
    )

    # Get Rev-Clock data
    t_name = "HiRISE Calibration Image"
    hci_dict = isis.cube.get_table(in_path, t_name)
    cal_vals = np.array(hci_dict["Calibration"])

    cal_image = np.ma.masked_outside(
        cal_vals, specialpix.Min, specialpix.Max
    )

    mask_lines = int(20 / binning)
    ramp_start = 20 + mask_lines

    # Get the average slope of the dark ramp
    hca_dict = isis.cube.get_table(in_path, "HiRISE Calibration Ancillary")
    caldark_vals = np.array(hca_dict["DarkPixels"])
    caldark = np.ma.masked_outside(
        caldark_vals, specialpix.Min, specialpix.Max
    )
    # print(caldark.shape)
    dark_slopes = np.ma.apply_along_axis(
        get_ramp_slope, 0, caldark[ramp_start:, 1:]
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
        f"Fraction of LIS pixels in Reverse-Clock area: {revclk_lisfrac:.2}"
    )
    if revclk_lisfrac > tolerance:
        fixed = True
    else:
        logger.info(
            f"Less than tolerance ({tolerance}), will not fix Reverse-Clock."
        )

    # Get Buffer Data
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

    buffer_lisfrac = np.ma.count_masked(buffer) / buffer.size
    logger.info(
        f"Fraction of LIS pixels in Buffer area: {buffer_lisfrac:.2}"
    )
    if buffer_lisfrac > tolerance:
        fixed = True
    else:
        logger.info(
            f"Less than tolerance ({tolerance}), will not fix Buffer."
        )

    ramp_lisfrac = np.ma.count_masked(
        cal_image[ramp_start:]
    ) / cal_image[ramp_start:].size

    if ramp_lisfrac > tolerance:
        fixed = False
        logger.warning(
            f"The fraction of LIS pixels in the ramp area ({ramp_lisfrac}) "
            f"is greater than the tolerance ({tolerance}), which prevents "
            f"lisfix processing."
        )

    if fixed:
        shutil.copy(in_path, out_path)
    else:
        return False

    # Fix rev-clock
    if revclk_lisfrac > tolerance:
        rev_model = ramp_rev_clock(
            cal_image[ramp_start:],
            label["IsisCube"]["Instrument"]["ChannelNumber"],
            zero_correction
        ).astype(int)
        # # Mask out all of the rev-clock for testing, so we get wholesale
        # # replacement
        # orig_cal_image = np.ma.copy(cal_image)
        # o_mean = np.ma.mean(orig_cal_image[:20, :], axis=0)
        # y_err_upper = np.ma.max(orig_cal_image[:20, :], axis=0) - o_mean
        # y_err_lower = o_mean - np.ma.min(orig_cal_image[:20, :], axis=0)

        # cal_image[:20, :] = np.ma.masked

        # # Some plotting for debugging:
        # import matplotlib.pyplot as plt
        # plt.plot(o_mean, 'k', label='Mean of Rev-Clock columns')
        # # plt.errorbar(
        # #     np.linspace(0, o_mean.size, num=o_mean.size, endpoint=False),
        # #     o_mean,
        # #     yerr=[y_err_lower, y_err_upper],
        # #     fmt='k',
        # #     label='Mean of Rev-Clock columns'
        # # )
        # for i in range(19):
        #     plt.plot(orig_cal_image[i, :], 'k.')
        # plt.plot(orig_cal_image[19, :], 'k.', label='Original Rev-Clock')
        # plt.plot(rev_model, 'r', label="Fixed")
        # plt.legend(loc='best')
        # plt.show()
        # sys.exit()

        # When we "fix" pixels, we must ensure that we are only fixing the LIS
        # pixels.  The cal_image has all special pixels masked out, so
        # we need to create a new structure that only masks the LIS pixels.
        cal_lismasked = np.ma.masked_equal(cal_vals, specialpix.Lis)
        logger.info("Fixing Reverse-Clock LIS pixels.")
        fixed_cal = np.ma.apply_along_axis(
            fix_rev_clock, 0, np.ma.concatenate((
                cal_lismasked,
                rev_model.reshape((1, rev_model.size))
            ))
        )

        logger.info("Writing out Reverse-Clock pixels.")
        hci_dict["Calibration"] = apply_special_pixels(
            fixed_cal, specialpix
        ).data.tolist()
        isis.cube.overwrite_table(out_p, t_name, hci_dict)

    # Fix the buffer area:
    if buffer_lisfrac > tolerance:
        # first_im_line = 19+(20+label["IsisCube"]["Instrument"]["Tdi"])/binning
        sp_frac = np.ma.count_masked(calbuf[:20, :]) / calbuf[:20, :].size
        if sp_frac <= tolerance:
            refr = int(np.ma.median(calbuf[:20, :]))
            logger.info(
                f"Median of first 20 lines of Calibration Buffer Pixels: {refr}"
            )
        else:
            logger.info(
                f"More than {tolerance:.1%} of the Calibration Buffer "
                f"Pixels have special pixel values: {sp_frac:.2%}."
            )
            isp = pvl.loads(
                isis.catoriglab(in_path).stdout
            )["INSTRUMENT_SETTING_PARAMETERS"]
            adc = label["IsisCube"]["Instrument"]["ADCTimingSetting"]
            if adc == -9999:
                adc = isp["MRO:ADC_TIMING_SETTINGS"]

            cid = get_ChannelID_fromfile(in_path)
            refr = int(
                revclock_model(
                    cid.ccdname + cid.ccdnumber + "_" + cid.channel,
                    binning,
                    label[
                        "IsisCube"
                    ]["Instrument"]["FpaPositiveYTemperature"].value,
                    adc
                )
            )
            logger.info(f"Using a model-based substitute ({refr}).")
            # import hiproc.img as img
            # lut = img.LUT_Table(isp["MRO:LOOKUP_CONVERSION_TABLE"])
            # if refr <= lut.table[1]:
            #     logger.error(f"Using a model-based substitute ({refr}).")
            # else:
            #     logger.error(
            #         f"Model ({refr}) was greater than LUT floor "
            #         f"({lut.table[1]}). Setting to LUT floor."
            #     )
            #     refr = lut.table[1]

        if np.ma.count_masked(dark[:20]) / dark[:20].size <= tolerance:
            refd = int(np.ma.median(dark[:20]))
            logger.info(f"Median of first 20 lines of Dark Pixels: {refd}")
        else:
            refd = int(np.ma.median(fixed_cal[:20, :]))
            logger.info(
                f"More than {tolerance:.1%} of the first 20 lines of "
                f"Dark Pixels have a real "
                f"value: {np.ma.count_masked(dark[:20]) / dark[:20].size} "
                f"Using the rev-clock median ({refd})."
            )

        model_buffer = dark[:, 2:14] - refd + refr
        # print(model_buffer)
        # print(model_buffer.dtype)
        # print(model_buffer.shape)

        # # Mask out all of the buffer for testing, so we get wholesale
        # # replacement
        # buffer.mask = True

        logger.info("Fixing Buffer pixels.")
        buffer_lismasked = np.ma.masked_equal(buffer_vals, specialpix.Lis)
        fixed_buf = np.ma.apply_along_axis(
            fix_buffer,
            1,
            np.ma.concatenate((buffer_lismasked, model_buffer), 1),
            buffer_lismasked.shape[1]
        )
        # print(fixed_buf)
        # print(fixed_buf.dtype)
        # print(fixed_buf.shape)

        ha_dict["BufferPixels"] = apply_special_pixels(
            fixed_buf, specialpix
        ).data.tolist()
        logger.info("Writing out Buffer pixels.")
        isis.cube.overwrite_table(out_p, "HiRISE Ancillary", ha_dict)

    return fixed


def fit_ramp(col: np.ma.array):
    """Returns the slope from the ramp pixels."""
    return np.ma.polyfit(np.arange(col.size), col, 1)


def get_ramp_slope(col: np.ma.array):
    return fit_ramp(col)[0]


def get_ramp_intercept(col: np.ma.array):
    return fit_ramp(col)[1]


def fix_rev_clock(
    cal_image_col: np.ma.array,
    lis_fraction=0.2,
):
    """Returns a np.masked array where the masked pixels in the
    *cal_image_col* are replaced with a fill value from the end
    of the column, if the percent of masked pixels is greater than or
    equal to *lis_fraction*.

    It is assumed that cal_image_col array has one more final entry
    (at position [-1]) than it should, and that value is used to
    fill the masked areas, if needed, and then is not included
    in the array on return.
    """
    if np.ma.count_masked(cal_image_col[:20]) / 20 >= lis_fraction:
        cal_image_col[:20] = cal_image_col[:20].filled(
            fill_value=cal_image_col[-1]
        )
    return cal_image_col[:-1]


def ramp_rev_clock(
    ramp_image: np.ma.array,
    chan: int,
    zero_correction=0,
):
    ramp_model = np.ma.apply_along_axis(
        get_ramp_intercept, 0, ramp_image
    )

    return flatten(ramp_model - zero_correction, chan)


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


def revclock_model(
    channelid: str,
    binning: int,
    temperature: float,
    adc: abc.Sequence
) -> float:
    # The slope and intercept data is from a presentation by Ken
    # Herkenhoff, in July 2020, entitled "Analysis of flight reverse-clocked
    # data" in the file Full\ reverse\ analysis.pptx
    # However, I think the ADC 7, 4 values were added after that date, but
    # sometime before Feb 2021.
    adc = tuple(adc)
    if adc == (5, 4):
        if binning != 1:
            raise KeyError(
                "There is no model for ADC 5,4 data with a binning other "
                "than 1."
            )
        pairs = {
            "RED0_0": (1278, 4.08),
            "RED1_0": (1072, 3.32),
            "RED2_0": (1214, 3.45),
            "RED3_0": (1083, 3.70),
            "RED4_0": (1124, 3.52),
            "RED5_0": (1148, 4.90),
            "RED6_0": (992, 4.14),
            "RED7_0": (1117, 1.90),
            "RED8_0": (1269, 4.56),
            "RED9_0": (1144, 3.33),
            "IR10_0": (1121, 4.92),
            "IR11_0": (1171, 3.35),
            "BG12_0": (986, 3.78),
            "BG13_0": (1027, 4.10),
            "RED0_1": (1341, 5.70),
            "RED1_1": (1114, 4.65),
            "RED2_1": (1093, 4.73),
            "RED3_1": (1187, 5.15),
            "RED4_1": (1078, 4.54),
            "RED5_1": (1335, 5.50),
            "RED6_1": (1154, 5.36),
            "RED7_1": (1355, 5.55),
            "RED8_1": (1398, 4.27),
            "RED9_1": (1255, 5.47),
            "IR10_1": (991, 5.81),
            "IR11_1": (1145, 4.20),
            "BG12_1": (1112, 5.67),
            "BG13_1": (974, 5.32),
        }
        intercept, slope = pairs[channelid]

    elif adc == (7, 4):
        pairs = {
            "RED1_0": {
                1: (1010.4, 0.7),
                2: (930.4, 1.9),
                4: (914.4, 2.8),
            },
            "RED1_1": {
                1: (869.2, 3.5),
                2: (845.3, 3.1),
                4: (699.2, 7.1),
            },
            "RED2_0": {
                1: (1051.2, 1.6),
                2: (962.7, 2.6),
                4: (923.6, 4.1),
            },
            "RED2_1": {
                1: (842.5, 3.9),
                2: (734.9, 5.1),
                4: (644.9, 7.7),
            },
            "RED3_0": {
                1: (958.8, 1.8),
                2: (841.5, 3.7),
                4: (788.1, 5.5),
            },
            "RED3_1": {
                1: (916.8, 3.5),
                2: (759.8, 6.1),
                4: (648.4, 9.1),
            },
        }
        intercept, slope = pairs[channelid][binning]

    else:
        raise KeyError(
            f"Do not have a model for {adc}."
        )

    return temperature * slope + intercept


def flatten(arr: np.array, chan: int):
    # The idea here is to flatten the zero_correction array between the pause
    # points
    sos = signal.butter(1, 0.5, output="sos")
    pl = list()
    for samp, width in zip(util.ch_pause[chan], util.ch_width[chan]):
        sl = util.pause_slicer(samp, width)
        pl.append(sl.start)
        pl.append(sl.stop)

    # # Supports debug plotting below, 2 lines
    # fixed = np.copy(arr)
    # filtered_whole = np.full_like(arr, np.mean(arr))
    for pslice in (
        slice(None, pl[0]),
        slice(pl[1], pl[2]),
        slice(pl[3], pl[4]),
        slice(pl[5], None)
    ):
        m = np.mean(arr[pslice])
        filtered = signal.sosfilt(sos, arr[pslice] - m)
        arr[pslice] = (arr[pslice] - m - filtered) + m
        # # Supports debug plotting below, 2 lines
        # fixed[pslice] = (arr[pslice] - m - filtered) + m
        # filtered_whole[pslice] = filtered + m

    # # Some plotting for debugging:
    # import matplotlib.pyplot as plt
    # plt.plot(arr, 'k', label='Original')
    # # plt.plot(filtered_whole, 'b', label='Filtered')
    # plt.plot(fixed, 'r', label="Fixed")
    # plt.legend(loc='best')
    # plt.show()
    # sys.exit()

    return arr


if __name__ == "__main__":
    main()
