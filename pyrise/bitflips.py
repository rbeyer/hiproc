#!/usr/bin/env python
"""Deal with stuck bits in HiRISE pixels.

If there is some process which is flipping one (or more) of the
HiRISE bits in a 14-bit pixel the wrong way, it will inadvertently
add or subtract a certain DN value to a pixel.  Naturally, it goes
from 2^13 (8192) all the way down to 2.

In practice, these 'bit-flipped' pixels show up as spikes in the
image histogram at DN values approximately centered on the image
Median plus or minus the bit-flip value (8192, 4096, etc.).  These
are not narrow spikes and have variable width, so care must be taken
when dealing with them.  The large bit-flip values often show up
as easily-identifiable islands in the DN histogram, but small
bit-flip values (2, 4, 8, etc.) are often within the variability
of the data and it is more difficult to cleanly isolate them.

This library and program deals with two aspects of the bit-flip
problem: identification and mitigation.


Identification
--------------

This can be as simple as deciding that some bit-flip position (maybe
16 or 32) is the cut-off, and simply saying that the DN window of
'good' pixles is Median plus/minus the cut-off.  However, this is
often a blunt instrument, and can still leave problematic outliers.

An improved strategy is to attempt to identify the 'islands' in the
DN histogram that are centered on the bit-flip positions, and attempt
to find the gaps between those islands.  Again, picking a bit-flip
value cut-off, but this time selecting the low-point between two
bit-flip-value-centered islands as the boundaries of the good-DN-window.
This strategy works okay for test data, but real data of the martian
surface often has complicated histograms, functions for this approach
can be found un bitunflip.py.

The best strategy we have developed so far is to take the above one
step further and treat the histogram curve as its own signal,
identify all of the minima in the histogram, along with a sense of
the standard deviation of the 'good' data in the image, and attempt
to smartly pick the bounds of the good-DN-window.


Mitigation
----------
This library and the program can perform one of two general activities
to mitigate the bit-flip pixels once they have been identified.

1. Unflipping: It can attempt to un-flip the bits in the pixel to
   attempt to return it to what its 'original' DN value might have
   been.  This is limited by however you select the bit-flip cut-off.
   Since you can't identify below that cut-off, you also really can't
   unflip below that cut-off.  Functions for this approach can be found
   in bitunflip.py

2. Masking: This mitigation strategy simply sets these pixels to
   an arbitrary DN value or the ISIS NULL value.  This basically just
   exercises the 'identification' part of this library and leaves
   it up to other image processing activities to decide how to deal
   with this data.
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

# This program can also be run from within the Perl HiCal program's
# Mask() function like so:
#
#   $cmd = "bitflips.py -l INFO -o $mask_file $tmask_file"

import argparse
import logging
import json
import math
import os
import shutil
import subprocess
import sys
import traceback
from pathlib import Path

import numpy as np
from osgeo import gdal_array
from scipy import interpolate
from scipy.signal import find_peaks
from scipy.stats import mstats

import pvl
import kalasiris as isis

import pyrise.img as img
import pyrise.util as util


# This is only here for temporary testing
goal_dn = dict()
with open(Path("/Users/rbeyer/projects/HiRISE/bitflip/dn_ideal.json"), 'r') as f:
    dn_ideal = json.load(f)

for k, v in dn_ideal.items():
    p = Path(v['path']).name
    for key, title in (("img", "Image Area"), ("rev", "Reverse-Clock")):
        if key in v:
            goal_dn[f"{p} {title}"] = v[key]


def main():
    try:
        parser = argparse.ArgumentParser(
            description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter,
            parents=[util.parent_parser()],
        )
        parser.add_argument(
            "-o", "--output", required=False, default=".bitflip.cub"
        )
        parser.add_argument(
            "-w",
            "--width",
            required=False,
            default=5,
            type=int,
            help="The number of medstd widths for bit-flip " "cleaning. "
                 "(default: %(default)s)",
        )
        parser.add_argument(
            "--line",
            required=False,
            action="store_true",
            help="Performs statistics along the line "
            "direction instead of column for the image area.",
        )
        parser.add_argument(
            "-r",
            "--replacement",
            required=False,
            type=int,
            help="By default, the program will replace "
            "identified pixels with an appropriate NULL data "
            "value, but if provided this value will be used "
            "instead.",
        )
        parser.add_argument(
            "-p",
            "--plot",
            required=False,
            action="store_true",
            help="Displays interactive plot for each area.",
        )
        parser.add_argument(
            "--saveplot",
            required=False,
            nargs="?",
            default=False,
            const=True,
            help="Saves plot for each area to a file.  If a directory is "
                 "provided it will be used to save the plots to, otherwise "
                 "the directory of the input file will be used."
        )
        parser.add_argument(
            "-n",
            "--dryrun",
            required=False,
            action="store_true",
            help="Does not produce a cleaned output file.",
        )
        parser.add_argument("file", help="ISIS Cube file or PDS IMG to clean.")

        args = parser.parse_args()

        util.set_logging(args.log, args.logfile)

        out_p = util.path_w_suffix(args.output, args.file)

        clean(
            args.file,
            out_p,
            width=args.width,
            axis=(1 if args.line else 0),
            replacement=args.replacement,
            plot=args.plot,
            saveplot=args.saveplot,
            dryrun=args.dryrun,
            keep=args.keep,
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


def clean(
    in_path: os.PathLike,
    out_path: os.PathLike,
    width=5,
    replacement=None,
    axis=0,
    plot=False,
    saveplot=False,
    dryrun=False,
    keep=False,
):
    """The file at *out_path* will be the result of running bit-flip
    cleaning of the file at *in-path*.

    The file at *out_path* can be an ISIS cube file, or a PDS IMG file.

    **WARNING**: at this time, only the image-area and the reverse-clock
    area are undergoing bit-flip cleaning.

    If *replacement* is not specified, pixels that are identified
    as being beyond the allowable DN window (which may be defined
    differently for each of the image areas), will be set to the
    equivalent ISIS NULL value, otherwise they will be set to this
    value.

    If *keep* is True, then all intermediate files will be preserved,
    otherwise, this function will clean up any intermediary files
    it creates.
    """

    in_p = Path(in_path)
    out_p = Path(out_path)

    if saveplot:
        try:
            saveplot = Path(saveplot)
            if not saveplot.is_dir():
                raise NotADirectoryError(
                    f"{saveplot} either doesn't exists or isn't a directory."
                )
        except TypeError:
            saveplot = in_p.parent

    label = pvl.load(str(in_p))
    if "IsisCube" in label:
        clean_cube(
            in_p, out_p, label, width, replacement, axis, plot, saveplot,
            dryrun, keep
        )
    elif "PDS_VERSION_ID" in label:
        clean_img(
            in_p, out_p, label, width, replacement, axis, plot, saveplot,
            dryrun, keep
        )
    else:
        raise ValueError(
            f"The file at {in_p} is not an ISIS Cube or a " "PDS IMG fie."
        )
    return


def clean_cube(
    in_p: Path,
    out_p: Path,
    label=None,
    width=5,
    replacement=None,
    axis=0,
    plot=False,
    saveplot=False,
    dryrun=False,
    keep=False,
):
    """ISIS Cube version of clean().

    Please see clean() for argument details.
    """

    to_del = isis.PathSet()

    # Bit-flip correct the non-image areas.
    tblcln_p = to_del.add(in_p.with_suffix(".tableclean.cub"))
    clean_tables_from_cube(
        in_p, tblcln_p, width=width, plot=plot, saveplot=saveplot, dryrun=dryrun
    )

    # Now clean the image area.
    if label is None:
        label = pvl.load(str(in_p))
    specialpix = getattr(
        isis.specialpixels, label["IsisCube"]["Core"]["Pixels"]["Type"]
    )
    image = np.ma.masked_outside(
        gdal_array.LoadFile(str(in_p)), specialpix.Min, specialpix.Max
    )

    logging.info(f"Bit-flip cleaning Image area.")
    # These four lines are just informational.
    img_mean = np.ma.mean(image)
    img_mode = mstats.mode(image, axis=None)[0][0]
    d = img_mean - img_mode
    logging.info(f"Mean: {img_mean}, Mode: {img_mode}, diff: {d}")

    (s_min, s_max) = find_smart_window_from_ma(
        image,
        width=width,
        axis=axis,
        plot=plot,
        plottitle=(f"{in_p.name} Image Area" if plot or saveplot else None),
        saveplot=(
            saveplot / in_p.with_suffix(".bf-image.pdf").name
        ) if saveplot else False,
    )

    if not dryrun:
        isis.mask(
            tblcln_p,
            mask=tblcln_p,
            to=out_p,
            minimum=s_min,
            maximum=s_max,
            preserve="INSIDE",
            spixels="NONE",
        )

        if replacement is not None:
            null_p = to_del.add(tblcln_p.with_suffix(".null.cub"))
            shutil.copy(out_p, null_p)
            isis.stretch(null_p, to=out_p, null=replacement)

        if not keep:
            to_del.unlink()

    return


def clean_img(
    in_path: Path,
    out_path: Path,
    label=None,
    width=5,
    replacement=None,
    axis=0,
    plot=False,
    saveplot=False,
    dryrun=False,
    keep=False,
):
    """PDS IMG file version of clean().

    This function is currently quite slow and takes almost a minute to
    process a 50,000 line image.  Since the ability to process IMG files
    isn't a primary task, and this is just a proof-of-concept, we can live
    with slow.  It can always be optimized, if needed.

    Please see clean() for argument details.
    """

    if label is None:
        label = pvl.load(str(in_path))
    if "PDS_VERSION_ID" not in label:
        raise ValueError(
            f"The file at {in_path} does not appear to be a PDS IMG file."
        )

    # Bit-flip correct the non-image areas.
    # This is going to diverge from cubes for the Buffer and Dark pixels.
    # In ISIS-land, they are 'tables' but in a IMG these are part of the
    # image array, and will need to be handled below.
    tblcln_p = in_path.with_suffix(".tableclean.img")
    clean_tables_from_img(
        in_path, tblcln_p, label, width, replacement,
        plot=plot, saveplot=saveplot, dryrun=dryrun
    )

    # Now clean the image area.
    lut = img.LUT_Table(
        label["INSTRUMENT_SETTING_PARAMETERS"]["MRO:LOOKUP_CONVERSION_TABLE"]
    )
    specialpix = lut.specialpixels()
    t_name = "IMAGE"

    # Just doing this takes almost a minute for a 50,000 line image!
    img_arr = img.object_asarray(in_path, t_name)

    # Using GDAL is better, but need to unlut.
    #   Hunh ... The results when I do this just aren't right, and I'm
    #   not quite sure what's going wrong.  These print statements seem
    #   to show well-behaved arrays that are properly unlutted.  However,
    #   when I run this with ESP_061686_1725_RED3_1.IMG, my test image,
    #   it only finds a single 'good' column with 34,279 valid pixels
    #   (instead of 21 with 50,000), but it has a huge medstd: 1207
    #   (instead of ~50).  So this might mean that the unlut + mask
    #   operation is somehow masking more values than it should (lower
    #   valid count), but also leaving some high-value values (to produce
    #   the high medstd).  The unlut and masking are literally the same
    #   functions, so that seems odd.  Maybe gdal_array.LoadFile() is
    #   doing unexpected?  Not sure.  If we really need the speed-up,
    #   we can try and work this.
    # unlut = np.vectorize(lut.unlut)
    # from_gdal = gdal_array.LoadFile(str(in_path))
    # print(from_gdal)
    # img_arr = unlut(from_gdal)
    # print(img_arr)

    image = np.ma.masked_outside(img_arr, specialpix.Min, specialpix.Max)

    logging.info(f"Bit-flip cleaning Image area.")
    (s_min, s_max) = find_smart_window_from_ma(
        image,
        width=width,
        axis=axis,
        plot=plot,
        plottitle=(f"{in_path.name} Image Area" if plot or saveplot else None),
        saveplot=(
            saveplot / in_p.with_suffix(".bf-image.pdf").name
        ) if saveplot else False
    )

    if not dryrun:
        shutil.copy(tblcln_p, out_path)
        clean_image = np.ma.masked_outside(image, s_min, s_max)
        if replacement is not None:
            specialpix.Null = replacement
        img_arr = apply_special_pixels(clean_image, specialpix)
        img.overwrite_object(out_path, "IMAGE", img_arr)

        if not keep:
            tblcln_p.unlink()

    return


def clean_tables_from_cube(
    in_path: Path,
    out_path: Path,
    width=5,
    replacement=None,
    rev_area=True,
    mask_area=False,
    ramp_area=False,
    buffer_area=False,
    dark_area=False,
    plot=False,
    saveplot=False,
    dryrun=False,
):
    """The file at *out_path* will be the result of running bit-flip
    cleaning of the specified non-image areas in table objects within
    the ISIS cube file at *in_path*.

    The various *_area* booleans direct whether these areas will
    undergo bit-flip cleaning.

    **WARNING**: at this time, only the reverse-clock area is
    enabled for undergoing bit-flip cleaning.  Specifying True
    for any others will result in a NotImplementedError.

    This means that pixels that are identified as being beyond the
    allowable DN window (which may be defined differently for each
    of the image areas), will be set to the ISIS NULL value.

    If *keep* is True, then all intermediate files will be preserved,
    otherwise, this function will clean up any intermediary files
    it creates.

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
    for notimpl in (mask_area, ramp_area, buffer_area, dark_area):
        if notimpl is True:
            raise NotImplementedError(
                f"Bit-flip cleaning for {notimpl} is not yet implemented."
            )

    label = pvl.load(str(in_path))
    if "IsisCube" not in label:
        raise ValueError(
            f"The file at {in_path} does not appear to be an ISIS Cube."
        )

    binning = label["IsisCube"]["Instrument"]["Summing"]
    specialpix = getattr(
        isis.specialpixels, label["IsisCube"]["Core"]["Pixels"]["Type"]
    )

    if not dryrun:
        shutil.copy(in_path, out_path)

    if any((rev_area, mask_area, ramp_area)):
        t_name = "HiRISE Calibration Image"
        HCI_dict = isis.cube.get_table(in_path, t_name)
        cal_vals = np.array(HCI_dict["Calibration"])

        cal_image = np.ma.masked_outside(
            cal_vals, specialpix.Min, specialpix.Max
        )

        clean_cal = clean_cal_tables(
            cal_image,
            binning,
            width,
            rev_area,
            mask_area,
            ramp_area,
            plot,
            str(in_path.name),
            (
                (
                saveplot / in_path.with_suffix(".bf-revclk.pdf").name
                ) if saveplot else False
            ),
        )
        if not dryrun:
            # write the table back out
            if replacement is not None:
                specialpix.Null = replacement
            HCI_dict["Calibration"] = apply_special_pixels(
                clean_cal, specialpix
            ).data.tolist()
            isis.cube.overwrite_table(out_path, t_name, HCI_dict)

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


def clean_tables_from_img(
    in_path: Path,
    out_path: Path,
    label=None,
    width=5,
    replacement=None,
    rev_area=True,
    mask_area=False,
    ramp_area=False,
    buffer_area=False,
    dark_area=False,
    plot=False,
    saveplot=False,
    dryrun=False,
):
    """The file at *out_path* will be the result of running bit-flip
    cleaning of the specified non-image areas in table objects within
    the PDS IMG file at *in_path*.

    The various *_area* booleans direct whether these areas will
    undergo bit-flip cleaning.

    **WARNING**: at this time, only the reverse-clock area is
    enabled for undergoing bit-flip cleaning.  Specifying True
    for any others will result in a NotImplementedError.

    This means that pixels that are identified as being beyond the
    allowable DN window (which may be defined differently for each
    of the image areas), will be set to the ISIS NULL value.

    If *keep* is True, then all intermediate files will be preserved,
    otherwise, this function will clean up any intermediary files
    it creates.

    This function anticipates a HiRISE EDR PDS IMG file that has
    the following objects in its label: CALIBRATION_IMAGE,
    LINE_PREFIX_TABLE, and LINE_SUFFIX_TABLE.

    The CALIBRATION_IMAGE contains the Reverse-Clock, Mask, and
    Ramp Image areas and their Buffer and Dark pixels.  The HiRISE
    Ancillary table contains the BufferPixels and DarkPixels from
    either side of the main Image Area.

    The LINE_PREFIX_TABLE contains the Buffer pixels for the image area.

    The LINE_SUFFIX_TABLE contains the Dark pixels for the image area.
    """

    for notimpl in (mask_area, ramp_area, buffer_area, dark_area):
        if notimpl is True:
            raise NotImplementedError(
                f"Bit-flip cleaning for {notimpl} is not yet implemented."
            )

    if label is None:
        label = pvl.load(str(in_path))
    if "PDS_VERSION_ID" not in label:
        raise ValueError(
            f"The file at {in_path} does not appear to be a PDS IMG file."
        )

    binning = label["INSTRUMENT_SETTING_PARAMETERS"]["MRO:BINNING"]
    specialpix = img.LUT_Table(
        label["INSTRUMENT_SETTING_PARAMETERS"]["MRO:LOOKUP_CONVERSION_TABLE"]
    ).specialpixels()

    if not dryrun:
        shutil.copy(in_path, out_path)

    if any((rev_area, mask_area, ramp_area)):
        t_name = "CALIBRATION_IMAGE"
        cal_vals = img.object_asarray(in_path, t_name)

        # prefix = label[t_name]['LINE_PREFIX_BYTES']
        # suffix = label[t_name]['LINE_SUFFIX_BYTES']
        # rev_clock_slice = np.s_[:20, prefix:(-1 * suffix)]
        rev_clock_slice = np.s_[:20, :]

        # print(cal_vals.shape)
        # print(cal_vals[:20, 18:-16])
        # print(cal_vals[:20, 18:-16].shape)
        # print(specialpix)
        cal_image = np.ma.masked_outside(
            cal_vals[rev_clock_slice], specialpix.Min, specialpix.Max
        )

        # print("cal_image:")
        # print(cal_image)
        # print(cal_image.shape)

        clean_cal = clean_cal_tables(
            cal_image,
            binning,
            width,
            rev_area,
            mask_area,
            ramp_area,
            plot,
            str(in_path.name),
            (
                (
                    saveplot / in_p.with_suffix(".bf-revclk.pdf").name
                ) if saveplot else False
            ),
        )
        if not dryrun:
            # write the table back out
            if replacement is not None:
                specialpix.Null = replacement
            cal_vals[rev_clock_slice] = apply_special_pixels(
                clean_cal, specialpix
            )
            img.overwrite_object(out_path, t_name, cal_vals)
    return


def clean_cal_tables(
    cal_image,
    binning,
    width=5,
    rev_area=True,
    mask_area=False,
    ramp_area=False,
    plot=False,
    imagestr=None,
    saveplot=False,
):
    # Deal with the HiRISE Calibration Image first (Reverse-clock, Mask,
    # and Ramp

    for notimpl in (mask_area, ramp_area):
        if notimpl is True:
            raise NotImplementedError(
                f"Bit-flip cleaning for {notimpl} is not yet implemented."
            )

    if rev_area:
        logging.info(f"Bit-flip cleaning Reverse-Clock area.")
        # Apply the cleaning along the lines, rather than along
        # the columns, since the statistics are better, and also
        # restrict the medstd_limit down to 200, since these reverse
        # clock pixels should not be that divergent.
        rev_clean = clean_array(
            cal_image[:20, :],
            width=width,
            axis=1,
            medstd_limit=200,
            plot=plot,
            plottitle=(f"{imagestr} Reverse-Clock" if plot or saveplot else None),
            saveplot=saveplot
        )
        cal_image[:20, :] = rev_clean

    # mask_lines = int(20 / binning)
    # if mask_area:
    #     mask_pixels = cal_image[20:mask_lines, :]

    # if ramp_area:
    #     ramp_pixels = cal_image[20 + mask_lines:, :]

    return cal_image


def clean_array(data: np.ma.array, *args, **kwargs):
    """Returns a numpy masked array whose mask is based on applying
    the smart window bounds from find_smart_window_from_ma() applied
    to *data* with the specified parameters.  All of the parameters
    provided to this function are passed along to
    find_smart_window_from_ma(), see its documentation for details.

    If all of the values in *data* are masked, a ValueError is raised.
    """
    if np.all(data.mask):
        logging.info("All of the values in data are masked.")
        return data

    (w_min, w_max) = find_smart_window_from_ma(data, *args, **kwargs)
    return np.ma.masked_outside(data, w_min, w_max)


def fit_array(data: np.ma.array):
    """Returns an array with the masked elements of *data* replaced
    by fitted values.

    This function uses the scipy.interpolate.griddata algorithm to
    interpolate a 2-D function across the unmasked data to replace
    the masked values of *data*.
    """
    # An early version of this function had an option to run B-spline
    # interpolation along the rows, using the code below.  However,
    # due to the 2D nature of the data, this B-spline interpolation
    # ended up looking 'stripey' in the row direction, so it was not
    # developed further.
    if data.ndim != 2:
        raise ValueError(
            "The provided array does not have two dimensions, "
            f"it has {data.ndim}."
        )

    xx, yy = np.meshgrid(
        np.arange(data.shape[0]), np.arange(data.shape[1]), indexing="ij"
    )
    x1 = xx[~data.mask]
    y1 = yy[~data.mask]
    good = data[~data.mask]

    interp = interpolate.griddata(
        (x1, y1), good.ravel(), (xx, yy), method="nearest"
    )
    return interp


def apply_special_pixels(array: np.ma, specialpix) -> np.ma:
    """Return a modified version of *array* where the values of
    array.data that are masked (array.mask is True) are examined,
    and the array.data value is set to one of the values in the
    *specialpix* namedtuple (which is assumed to conform to the
    kalasiris.specialpixels.SpecialPixels namedtuple interface).

    If the value is already the same as the Null, Lrs, Lis, His,
    or Hrs value of *specialpixels* it is left as-is.  If it is
    less than specialpix.Min it is set to specialpix.Lrs.  If it
    is greater than specialpix.Max, it is set to specialpix.Hrs.
    Finally, if it is any other value, it is set to specialpix.Null.
    """

    def sp_pix(val):
        if val in (
            specialpix.Null,
            specialpix.Lrs,
            specialpix.Lis,
            specialpix.His,
            specialpix.Hrs,
        ):
            return val
        elif val < specialpix.Min:
            return specialpix.Lrs
        elif val > specialpix.Max:
            return specialpix.Hrs
        else:
            return specialpix.Null

    sp_apply = np.frompyfunc(sp_pix, 1, 1)
    array.data[array.mask] = sp_apply(array.data[array.mask])
    return array


def find_smart_window_from_ma(
    data: np.ma.array,
    width=5,
    axis=0,
    medstd_limit=300,
    medstd_fallback=64,
    plot=False,
    plottitle=None,
    saveplot=False
):
    """Returns a two-tuple with the result of find_smart_window().

    This function mostly just does set-up based on the provided
    masked array *data* and the provided *width* and *axis* value
    to calculate the inputs for find_smart_window().  The
    value of *medstd_limit* and *medstd_fallback* are passed along
    to min_max_ex(), please see that documentation for more information.
    """
    # For some reason, I appear to have built this median_limit()
    # function, but in practice it doesn't seem to help, as there
    # doesn't seem to be a problem with letting the median be whatever
    # it actually is.
    # median = median_limit(np.ma.median(data), data)
    median = np.ma.median(data)
    logging.info(f"Median: {median}")
    medstd = median_std_from_ma(data, axis=axis)

    unique, unique_counts = np.unique(data.compressed(), return_counts=True)

    mindn, maxdn, ex = min_max_ex(
        median, medstd, width, medstd_limit, medstd_fallback
    )

    return find_smart_window(
        unique,
        unique_counts,
        mindn,
        maxdn,
        median,
        central_exclude_dn=ex,
        plot=plot,
        plottitle=plottitle,
        saveplot=saveplot
    )


def median_limit(median, data: np.ndarray, limit=4000):
    """Return the 'best' median of the provided data.

    If the extracted median is larger than *limit*, the algorithm
    will attempt to find a better representation of the median.  If
    it cannot, it will return the median, even if larger than
    *limit*.
    """
    logging.info(f"Median: {median}")

    if median > limit and np.any([data < limit]):
        median = np.ma.median(data[data < limit])
        logging.info(
            f"The median was too high (> {limit}), "
            f"found a better one: {median}."
        )

    return median


def median_std_from_ma(data: np.ma, axis=0):
    """On the assumption that there are bit-flips in the *data*,
    attempt to find a value that might represent the standard
    deviation of the 'real' data.  The *data* object must be a
    numpy masked array.

    The value of *axis* determines which way the data are handled.
    The default is 0 to scan vertially to accumulate statistics for
    columns.  In this case, only those columns with the most unmasked
    data are evaluated.  For them, the standard deviation is found
    for each column, and the returned value is the median of those
    standard deviations.  If *axis* is 1, then this is applied to
    rows, not columns.
    """
    valid_points = data.count(axis=axis)
    std_devs = np.std(data, axis=axis)
    return median_std(valid_points, std_devs)


def median_std(valid_points: np.ma, std_devs: np.ma):
    """On the assumption that there are bit-flips in the data
    represented by the two arrays *valid_points* and *std_devs*,
    attempt to find a value that might represent the standard
    deviation of the 'real' data.  The *valid_points* numpy masked
    array should be the count of valid points along some axis (for
    the rows or columns).  The *std_devs* numpy masked array
    should be the same size as *valid_points* and should be the
    standard deviation of those rows or columns.
    """
    maxvp = max(valid_points)
    logging.info(f"Maximum count of valid pixels along axis: {maxvp}")

    std_w_maxvp = np.extract(valid_points == maxvp, std_devs)
    medstd = np.ma.median(std_w_maxvp)
    logging.info(f"Number of rows/columns with that count: {len(std_w_maxvp)}")
    logging.info(
        "Median standard deviation of elements along the axis "
        f"that have the maximum valid pixel count: {medstd}."
    )
    return medstd


def min_max_ex(
    central, medstd, width, medstd_limit=300, medstd_fallback=64
) -> tuple:
    """Return a minimum, maximum, and exclusion value based on the
    provided *central* value, and adding and subtracting the result
    of multiplying the *medstd* by the *width*.

    If *medstd* is larger than *medstd_limit*, then *medstd_fallback*
    will be used instead.
    """
    # Sometimes even the medstd can be too high because the 'good' lines
    # still had too many outliers in it.  What is 'too high' and how do
    # we define it rigorously?  I'm not entirely sure, and I wish there
    # was a better way to determine this, or, alternately an even more robust
    # way to find 'medstd' in the first place.
    # I am going to select an abitrary value based on my experience (300)
    # and then also pick a replacement of 64 which is a complete guess.
    # Finally, in this case, since the 'statistics' are completely broken,
    # I'm also going to set the exclusion to be arbitrary, and lower than
    # the enforced medstd.
    if medstd > medstd_limit:
        medstd = medstd_fallback
        ex = 16
        logging.info(
            "The derived medstd was too big, setting the medstd "
            f"to {medstd} and the exclusion value to {ex}."
        )
    else:
        # We want to ignore minima that are too close to the central value.
        # Sometimes the medstd is a good choice, sometimes 16 DN (which is a
        # minimal bit flip level) is better, so use the greatest:
        ex = max(medstd, 16)

    mindn = central - (width * medstd)
    maxdn = central + (width * medstd)

    return mindn, maxdn, ex


# def mask(in_path: Path, out_path: Path, line=False, plot=True, keep=False):
#     """Attempt to mask out pixels beyond the central DNs of the median
#     based on minima in the histogram.
#
#     This is now superceded by clean_cube()
#     """
#
#     from pyrise.HiCal import analyze_cubenorm_stats2
#
#     to_del = isis.PathSet()
#
#     hist_p = to_del.add(in_path.with_suffix('.hist'))
#     hist = histogram(in_path, hist_p)
#
#     median = math.trunc(float(hist['Median']))
#     logging.info(f'Median: {median}')
#
#     # With terrible noise, the computed median is the median of the
#     # noise, not of the data, so if that happens, try and find a
#     # DN peak that is more reasonable:
#     median_limit = 4000
#     if median > median_limit:
#         trunc_hist = list()
#         for row in hist:
#             if int(row.DN) < median_limit:
#                 trunc_hist.append((int(row.DN), int(row.Pixels)))
#         try:
#             median = max(trunc_hist, key=lambda x: x[1])[0]
#             logging.info(
#                 f'The median was too high, found a better one: {median}.')
#         except ValueError:
#             # No values in the trunc_hist, which means that there are no
#             # DN less than maedian_limit, so we just leave things as they
#             # are.
#             pass
#
#     cubenorm_stats_file = to_del.add(in_path.with_suffix('.cn.stats'))
#     if line:
#         util.log(isis.cubenorm(in_path, stats=cubenorm_stats_file,
#                                direction='line').args)
#     else:
#         util.log(isis.cubenorm(in_path, stats=cubenorm_stats_file).args)
#     (mindn, maxdn) = analyze_cubenorm_stats2(cubenorm_stats_file, median,
#                                              hist, width=5, plot=plot)
#
#     # To bypass the 'medstd' calculations in analyze_cubenorm_stats2(),
#     # this mechanism can be used to set the limits directly.
#     # (mindn, maxdn) = find_smart_window(hist,
#     #                                    math.trunc(float(hist['Minimum'])),
#     #                                    math.trunc(float(hist['Maximum'])),
#     #                                    median, plot=True)
#
#     # util.log(isis.mask(in_path, to=out_path, minimum=maskmin,
#     #                    maximum=maskmax).args)
#
#     if not keep:
#         to_del.unlink()
#     return (mindn, maxdn)


def find_select_idx_name(
    central_idx: int, limit_idx: int, close_to_limit: bool
):
    """Returns a two-tuple which contains an int of value 0 or -1 in
    the first position and a string of value 'min' or 'max' in the second.

    The string is determined based on the relative values of
    *central_idx* to *limit_idx*, and the value of the int is meant
    as an index to a list, based on the three input variables.

    Raises *ValueError* if *central_idx* and *limit_idx* are the same.

    This is primarily a helper function to find_minima_index().
    """
    if limit_idx < central_idx:
        idx_name = "min"
        if close_to_limit is False:
            select_idx = -1
        else:
            select_idx = 0
    elif limit_idx > central_idx:
        idx_name = "max"
        if close_to_limit is False:
            select_idx = 0
        else:
            select_idx = -1
    else:
        raise ValueError

    return select_idx, idx_name


def find_minima_index(
    central_idx: int,
    limit_idx: int,
    minima_idxs: np.ndarray,
    pixel_counts: np.ndarray,
    close_to_limit=True,
) -> int:
    """Searches the pixel_counts at the minima_idxs positions between
    central_idx and limit_idx to find the one with the lowest
    pixel_count value. If multiple lowest values, choose the closest
    to the limit_idx position (the default, if *close_to_limit* is true,
    otherwise picks the index closest to *central_idx* if *close_to_limit*
    is False).

    If no minima_idx values are found in the range, the next value
    in minima_idx past the limit_idx will be returned.
    """
    # print(f'central_idx: {central_idx}')
    # print(f'limit_idx: {limit_idx}')
    info_str = "Looking for {} {} range ..."

    try:
        select_idx, idx_name = find_select_idx_name(
            central_idx, limit_idx, close_to_limit
        )
        logging.info(info_str.format(idx_name, "inside"))

        min_i = min(central_idx, limit_idx)
        max_i = max(central_idx, limit_idx)
        inrange_i = minima_idxs[(minima_idxs >= min_i) * (minima_idxs <= max_i)]
        logging.info(str(inrange_i) + " are the indexes inside the range.")

        try:
            value = min(np.take(pixel_counts, inrange_i))
            logging.info(
                f"{value} is the minimum Pixel count amongst those indexes."
            )

            idx = inrange_i[
                np.asarray(pixel_counts[inrange_i] == value).nonzero()
            ][select_idx]
        except ValueError:
            if limit_idx < central_idx:
                logging.info(info_str.format("min", "outside"))
                idx = max(minima_idxs[minima_idxs < limit_idx])
            elif limit_idx > central_idx:
                logging.info(info_str.format("max", "outside"))
                idx = min(minima_idxs[minima_idxs > limit_idx])

        logging.info(f"{idx} is the minimum index.")
    except ValueError:
        logging.info(
            "Could not find a valid minima."
            # "Could not find a valid minima, returning "
            # f"the limit index: {limit_idx}."
        )
        # idx = limit_idx
        idx = None

    return idx


def find_prominence_index(
    minima_idxs, prominences, min_idx, max_idx, counts,
    min_idx_limit=None, max_idx_limit=None
):
    # scaled = list()
    # for c, p in zip(counts, prominences):
    #     scaled.append(math.log10(c + p) - math.log10(c))
    if min_idx_limit is None:
        min_idx_limit = minima_idxs[0]

    if max_idx_limit is None:
        max_idx_limit = minima_idxs[-1]

    scaled = np.log10(counts + prominences) - np.log10(counts)

    p = list()
    for condition in(
        ((minima_idxs > min_idx) | (minima_idxs < min_idx_limit)),
        ((minima_idxs < max_idx) | (minima_idxs > max_idx_limit))
        # (minima_idxs < max_idx, len(prominences) - 1)

    ):
        # masked = np.ma.masked_where(condition, prominences)
        masked = np.ma.masked_where(condition, scaled)
        if masked.count() > 0:
            prom_i = masked.argmax()
            p.append(minima_idxs[prom_i])
        else:
            p.append(None)

    return p[0], p[1]


def find_smart_window(
    dn: np.ndarray,
    counts: np.ndarray,
    mindn: int,
    maxdn: int,
    centraldn: int,
    central_exclude_dn=0,
    plot=False,
    closest=True,
    plottitle=None,
    saveplot=False
) -> tuple:
    """Returns a minimum and maximum DN value from *dn* which are
       based on using the find_minima_index() function with the
       given *mindn*, *maxdn*, and *centraldn* values.  The *dn*
       array must be a sorted list of unique values, and *counts*
       is the number of times each of the values in *dn* occurs.

       If *central_exclude_dn* is given, the returned minimum and
       maximum DN are guaranteed to be at least *central_exclude_dn*
       away from *centraldn*.  This is useful if you don't want returned
       minimum and maximum DN to be too close to the *centraldn*.

       If *plot* is True, this function will display an interactive
       plot describing its work.  The curve represents the hist
       values, the shaded area marks the window between the given
       mindn and maxdn.  The 'x'es mark all the detected minima,
       and the red dots indicate the minimum and maximum DN values
       that this function will return.

       The value of *closest* is passed on to find_minima_index().

       If *saveplot* is not False, it will be assumed to be a filename
       that the plot should be saved as.  If False, the plot will not
       be saved.
    """
    # hist_list = sorted(hist, key=lambda x: int(x.DN))
    # pixel_counts = np.fromiter((int(x.Pixels) for x in hist_list), int)
    # dn = np.fromiter((int(x.DN) for x in hist_list), int)

    # print(dn)

    # print(f'mindn {mindn}')
    # print(f'maxdn {maxdn}')

    mindn_i = (np.abs(dn - mindn)).argmin()
    maxdn_i = (np.abs(dn - maxdn)).argmin()
    central_i = np.where(dn == centraldn)
    # print(f'centraldn: {centraldn}')
    # print(f'central_exclude_dn: {central_exclude_dn}')
    central_min_i = (np.abs(dn - (centraldn - central_exclude_dn))).argmin()
    # tobemin = np.abs(dn - (centraldn + central_exclude_dn))
    # print(tobemin)
    central_max_i = (np.abs(dn - (centraldn + central_exclude_dn))).argmin()

    minima_i, minprops = find_peaks(
        np.negative(counts),
        threshold=(None, None),
        distance=1,
        prominence=(None, None),
        width=(None, None),
        rel_height=1
    )
    # print("minprops: ")
    # print(minprops)
    maxima_i, maxprops = find_peaks(
        counts,
        threshold=(None, None),
        distance=1,
        prominence=(None, None),
        width=(None, None),
        rel_height=1,
    )
    # print("maxprops: ")
    # print(maxprops)

    # print(f'central_i {central_i}')
    # print(f'central_min_i {central_min_i}')
    # print(f'central_max_i {central_max_i}')
    # print(f'mindn_i {mindn_i}')
    # print(f'maxdn_i {maxdn_i}')

    if len(minima_i) == 0:
        logging.info("No minima found.")
        min_i = 0
        max_i = len(dn) - 1
        min_i_ips = min_in = min_prom = min_i
        max_i_ips = max_in = max_prom = max_i
        max_prom_idx = None
    else:
        min_i = find_minima_index(
            central_min_i, mindn_i, minima_i, counts, close_to_limit=closest
        )
        max_i = find_minima_index(
            central_max_i, maxdn_i, minima_i, counts, close_to_limit=closest
        )

        min_prom, max_prom = find_prominence_index(
            minima_i, minprops["prominences"], central_min_i, central_max_i,
            counts[minima_i]
        )
        min_in, max_in = find_prominence_index(
            minima_i, minprops["prominences"], central_min_i, central_max_i,
            counts[minima_i], min_i, max_i
        )

        max_prom_idx = np.argmax(maxprops["prominences"])
        max_prom_left = maxprops["left_ips"][max_prom_idx]
        max_prom_right = maxprops["right_ips"][max_prom_idx]
        # This reaches out to the "next" minima beyond the ips boundaries
        # # print(minima_i)
        # if max_prom_left not in minima_i and np.size(np.where(minima_i < max_prom_left)):
        #     # find the minima to the immediate left of max_prom_left
        #     max_prom_left = minima_i[np.where(minima_i < max_prom_left)[0][-1]]
        #     # print(f"!!! max_prom_left diff {np.where(minima_i < max_prom_left)[0][-1]}")
        # if max_prom_right not in minima_i and np.size(np.where(minima_i > max_prom_right)):
        #     # find the minima to the immediate left of max_prom_left
        #     max_prom_right = minima_i[np.where(minima_i > max_prom_right)[0][0]]
        #     # print(f"!!! max_prom_right diff {minima_i[np.where(minima_i > max_prom_right)[0][0]]}")
        print(f"initial max_prom_left/right: {max_prom_left}, {max_prom_right}")

        # If the image is very noisy (in some sense it is not clear that
        # it is worth fixing if it is this bad, but we will try our best),
        # the mode (which is what max_prom_idx really is giving above) might be
        # garbage.  If it is outside the mindn / maxdn range, reject it,
        # and pick the next best.
        print(f"mode index: {maxima_i[max_prom_idx]}")
        if not mindn_i <= maxima_i[max_prom_idx] <= maxdn_i:
            # scaled_maxima = np.log10(
            #     counts[maxima_i] + maxprops["prominences"]
            # ) - np.log10(counts[maxima_i])
            # sorted_prom_idxs = np.argsort(scaled_maxima)
            sorted_prom_idxs = np.argsort(maxprops["prominences"])
            print(f"sorted_prom_idxs: {sorted_prom_idxs}")
            print(f"maxima_i[sorted_prom_idxs]: {maxima_i[sorted_prom_idxs]}")

            max_prom_idx_i = -2  # -1 we already know isn't good.
            max_prom_idx = sorted_prom_idxs[max_prom_idx_i]
            while not maxprops["left_ips"][max_prom_idx] <= central_i[0][0] <= maxprops["right_ips"][max_prom_idx]:
                print(
                    f"while not {maxprops['left_ips'][max_prom_idx]} <= {central_i[0][0]} <= {maxprops['right_ips'][max_prom_idx]}"
                )
                max_prom_idx_i -= 1
                max_prom_idx = sorted_prom_idxs[max_prom_idx_i]

            # This is the original, so -1 idx is right:
            if maxima_i[sorted_prom_idxs[-1]] > maxdn_i:
                max_prom_right = maxprops["right_ips"][max_prom_idx]
                if max_prom_left > max_prom_right:
                    max_prom_left = maxprops["left_ips"][max_prom_idx]
            else:
                print("less than max")
                max_prom_left = maxprops["left_ips"][max_prom_idx]
                if max_prom_right < max_prom_left:
                    max_prom_right = maxprops["right_ips"][max_prom_idx]

        # Again, if the image is noisy, and central and mode are separate,
        # but not outside the mindn/maxdn bounds, you could get major prominences
        # between them, which shouldn't be allowed, so this adjusts that.
        if not max_prom_left <= central_i[0][0] <= max_prom_right:
            sorted_prom_idxs = np.argsort(maxprops["prominences"])
            max_prom_idx_i = -2  # -1 we already know isn't good.
            max_prom_idx = sorted_prom_idxs[max_prom_idx_i]
            while not maxprops["left_ips"][max_prom_idx] <= central_i[0][0] <= maxprops["right_ips"][max_prom_idx]:
                max_prom_idx_i -= 1
                max_prom_idx = sorted_prom_idxs[max_prom_idx_i]

            # Adjust the prominences:
            if maxima_i[sorted_prom_idxs[-1]] > central_i[0][0]:
                max_prom_left = maxprops["left_ips"][max_prom_idx]
            else:
                print("less than max")
                max_prom_right = maxprops["right_ips"][max_prom_idx]

            # Don't allow searching between the central DN and the mode.
            if central_min_i > min(central_i[0][0], maxima_i[sorted_prom_idxs[-1]]):
                    central_min_i = min(central_i[0][0], maxima_i[sorted_prom_idxs[-1]])

            if central_max_i < max(central_i[0][0], maxima_i[sorted_prom_idxs[-1]]):
                central_max_i = max(central_i[0][0], maxima_i[sorted_prom_idxs[-1]])

        print(f"max_prom_left/right: {max_prom_left}, {max_prom_right}")
        print(
            f"DN edges of prominence envelope: {dn[int(max_prom_left)]}, "
            f"{dn[int(max_prom_right)]}"
        )

        min_i_ips = find_minima_index(
            central_min_i, max_prom_left,
            minima_i, counts, close_to_limit=closest
        )
        max_i_ips = find_minima_index(
            central_max_i, max_prom_right,
            minima_i, counts, close_to_limit=closest
        )

    print(f"min_prom, max_prom: {min_prom}, {max_prom}")
    print(f"min_in, max_in: {min_in}, {max_in}")
    print(f"min_i_ips, max_i_ips: {min_i_ips}, {max_i_ips}")
    min_prom = 0 if min_prom is None else min_prom
    max_prom = len(dn) - 1 if max_prom is None else max_prom
    min_in = 0 if min_in is None else min_in
    max_in = len(dn) - 1 if max_in is None else max_in
    min_i = 0 if min_i is None else min_i
    max_i = len(dn) - 1 if max_i is None else max_i
    min_i_ips = 0 if min_i_ips is None else min_i_ips
    max_i_ips = len(dn) - 1 if max_i_ips is None else max_i_ips

    new_min_i, new_max_i = pick_index(
        min_i, min_i_ips, max_i, max_i_ips,
        maxima_i[max_prom_idx],
        # maxprops['left_ips'][max_prom_idx], maxprops['right_ips'][max_prom_idx],
        minprops["prominences"],
        counts, minima_i, dn, centraldn, mindn, maxdn
    )

    logging.info(f"indexes: {min_i}, {max_i}")
    logging.info(f"DN window: {dn[min_i]}, {dn[max_i]}")

    print(f"new_indexes: {new_min_i}, {new_max_i}")
    print(f"new DN window: {dn[new_min_i]}, {dn[new_max_i]}")

    if plot or saveplot:
        import matplotlib.pyplot as plt

        plt.ioff()
        fig, (ax0, ax1) = plt.subplots(2, 1, constrained_layout=True)

        indices = np.arange(0, len(counts))
        dn_window = np.fromiter(
            map(
                (lambda i: mindn_i <= i <= maxdn_i),
                (x for x in range(len(counts))),
            ),
            dtype=bool,
        )

        if plottitle is not None:
            fig.suptitle(plottitle)
            eval_bounds(
                plottitle, dn[min_i], dn[max_i], dn[min_prom], dn[max_prom],
                dn[min_in], dn[max_in]
            )

        ax0.set_ylabel("Pixel Count")
        ax0.set_xlabel("DN Index")
        ax0.set_yscale("log")
        ax0.fill_between(indices, counts, where=dn_window, color="lightgray")
        ax0.axvline(x=central_i, c="gray", label=f"Central DN: {centraldn}")
        mode = np.argmax(counts)
        ax0.axvline(x=mode, c="lime", ls="--", label=f"Mode: {dn[mode]}")
        ax0.plot(counts)
        ax0.plot(minima_i, counts[minima_i], "x")
        ax0.vlines(
            x=minima_i, ymin=counts[minima_i] + minprops["prominences"],
            ymax=counts[minima_i], color="C1"
        )
        ax0.hlines(
            y=(-1 * minprops["width_heights"]), xmin=minprops["left_ips"],
            xmax=minprops["right_ips"], color="C1"
        )
        # ax0.hlines(
        #     y=range(1, len(minprops["left_bases"])+1),
        #     xmin=minprops["left_bases"],
        #     xmax=minprops["right_bases"],
        #     color="C1",
        # )
        ax0.plot(maxima_i, counts[maxima_i], "^")
        ax0.vlines(
            x=maxima_i, ymin=counts[maxima_i] - maxprops["prominences"],
            ymax=counts[maxima_i], color="C2"
        )
        ax0.hlines(
            y=maxprops["width_heights"], xmin=maxprops["left_ips"],
            xmax=maxprops["right_ips"], color="C2"
        )
        # ax0.hlines(
        #     # y=range(10, (len(maxprops["left_bases"])+1) * 10, 10),
        #     y=counts[maxima_i],
        #     xmin=maxprops["left_bases"],
        #     xmax=maxprops["right_bases"],
        #     color="C2",
        # )

        ax0.plot(min_i, counts[min_i], "o", c="red")
        ax0.plot(max_i, counts[max_i], "o", c="red")
        ax0.plot(min_prom, counts[min_prom], ">", c="purple")
        ax0.plot(max_prom, counts[max_prom], "<", c="purple")
        ax0.plot(min_in, counts[min_in], "|", c="purple")
        ax0.plot(max_in, counts[max_in], "|", c="purple")
        ax0.plot(min_i_ips, counts[min_i_ips], "x", c="teal")
        ax0.plot(max_i_ips, counts[max_i_ips], "x", c="teal")
        ax0.plot(new_min_i, counts[new_min_i], "+", c="cyan")
        ax0.plot(new_max_i, counts[new_max_i], "+", c="cyan")

        ax0.legend()
        ax1.set_ylabel("Pixel Count")
        ax1.set_xlabel("DN")
        ax1.set_yscale("log")
        ax1.set_ybound(lower=0.5)
        if plottitle in goal_dn:
            ax1.axvline(
                x=goal_dn[plottitle][0],
                c="blue",
                label=f"Ideal DN: {goal_dn[plottitle]}"
            )
            ax1.axvline(
                x=goal_dn[plottitle][1],
                c="blue",
            )
        ax1.axvline(
            x=dn[min_i], c="red", ls="--",
            label=f"DN Limits: {dn[min_i]}, {dn[max_i]}"
        )
        ax1.axvline(x=dn[max_i], c="red", ls="--" )
        ax1.axvline(
            x=dn[new_min_i], c="cyan", ls="dotted",
            label=f"New DN Limits: {dn[new_min_i]}, {dn[new_max_i]}"
        )
        ax1.axvline(x=dn[new_max_i], c="cyan", ls="dotted")

        ax1.scatter(dn, counts, marker=".", s=1)
        ax1.legend()

        if saveplot:
            eval_bounds(
                plottitle, dn[min_i], dn[max_i], dn[min_prom], dn[max_prom],
                dn[min_in], dn[max_in], log_path= saveplot.with_name('bounds.log')
            )
            plt.savefig(saveplot)

        if plot:
            plt.show()

    # return dn[min_i], dn[max_i]
    return dn[new_min_i], dn[new_max_i]


def pick_index(
    min_i, min_i_ips, max_i, max_i_ips,
    max_prom_i, # max_prom_left, max_prom_right,
    prominences,
    counts, minima_i, dn, centraldn, mindn, maxdn
):
    scaled = np.log10(
        counts[minima_i] + prominences
    ) - np.log10(counts[minima_i])
    # print(f"scaled: {scaled}")

    span_factor = 2
    count_thresh = min(5000, counts[max_prom_i] / 100)
    scaled_depth_thresh = 0.3
    print(f"min/max i: {min_i}, {max_i}")
    print(f"min/max ips: {min_i_ips}, {max_i_ips}")
    print(f"count_thresh: {count_thresh}")
    print(f"count min/max i: {counts[min_i]} {counts[max_i]}")
    print(f"count min/max ips: {counts[min_i_ips]} {counts[max_i_ips]}")

    # print(f"{max_prom_left}, {max_prom_right}")
    # print(
    #     f"DN edges of prominence envelope: {dn[int(max_prom_left)]}, "
    #     f"{dn[int(max_prom_right)]}"
    # )

    new_i = [None, None]
    span_left = (centraldn - ((centraldn - mindn) * span_factor))
    span_right = (centraldn + ((maxdn - centraldn) * span_factor))
    print(f"DN edges of {span_factor} span envelope: {span_left}, {span_right}")
    for m, i, ips, dnlim, message in (
        (0, min_i, min_i_ips, mindn, "min"),
        (1, max_i, max_i_ips, maxdn, "max")
    ):
        print(f"Picking the index for {message}")
        if i == ips:
            print("i == ips")
            new_i[m] = i
        else:
            left = min(i, ips)
            right = max(i, ips)
            # print(f"minima_i: {minima_i}")
            # print(f"{left} LR {right}")
            if m == 0:
                potential_idxs = minima_i[(minima_i <= right)]
            else:
                # print(minima_i >= left)
                potential_idxs = minima_i[(minima_i >= left)]
            print(f"potential_idxs: {potential_idxs}")

            below_thresh = potential_idxs[counts[potential_idxs] < count_thresh]
            print(f"below_thresh: {below_thresh}")
            if np.size(below_thresh) == 0:
                below_thresh = potential_idxs
                print(f"fixed below_thresh: {below_thresh}")
                if np.size(below_thresh) == 0:
                    # Whoa, no potentials, either, so just fall back to the pair
                    below_thresh = np.array([i, ips])
            print(f"below_thresh: {below_thresh}")

            if m == 0:
                in_span = below_thresh[dn[below_thresh] >= span_left]
            else:
                in_span = below_thresh[dn[below_thresh] <= span_right]
            print(f"in_span: {in_span}")
            if np.size(in_span) == 0:
                in_span = below_thresh
                print(f"fixed in_span: {in_span}")


        #     print("sorting deep enough:")
        #     scaled_idxs = list()
        #     for span_i in in_span:
        #         scaled_idxs.append(np.where(minima_i == span_i)[0].item())
        #     print(scaled[scaled_idxs])
        #     deep_enough = np.argwhere(scaled[scaled_idxs] >= scaled_depth_thresh)
        #     if np.size(deep_enough) > 0:
        #         print(f"deep_enough: {deep_enough}")
        #         print(f"deep index: {in_span[deep_enough]}")
        #         good_idxs = in_span[deep_enough]
        #     else:
        #         scaled_idxs = list()
        #         for below_i in below_thresh:
        #             scaled_idxs.append(np.where(minima_i == below_i)[0].item())
        #         deep_enough = np.argwhere(scaled[scaled_idxs] >= scaled_depth_thresh)
        #         if np.size(deep_enough) > 0:
        #             print(f"fixed deep_enough: {deep_enough}")
        #             print(f"fixed deep_enough index: {below_thresh[deep_enough]}")
        #             good_idxs = below_thresh[deep_enough]
        #         else:
        #             print(f"depth thresh yielded no minima")
        #             good_idxs = in_span


        #     print(good_idxs)
        #     # print(good_idxs.min())
        #     # print(good_idxs.max())

        #     new_i[m] = best_index(scaled, counts, minima_i, good_idxs.min(), good_idxs.max(), m)
            new_i[m] = best_index(scaled, counts, minima_i, in_span[0], in_span[-1], m)
            # new_i[m] = best_index_frac(counts, in_span, m)

        # elif counts[i] > count_thresh:
        #     print("count i above thresh")
        #     new_i[m] = ips
        # elif counts[ips] > count_thresh:
        #     print("count ips above thresh")
        #     new_i[m] = i
        # else:
        #     # if max_prom_left < i < max_prom_right:
        #     #     print(f"{m} _i inside prominence envelope")
        #     #     new_i[m] = best_index(scaled, minima_i, i, ips)
        #     if not span_left < dn[ips] < span_right:
        #         print(f"DN of {m} ips is outside span, revert to red dot")
        #         new_i[m] = i
        #     else:
        #         print("best_index")
        #         # print(f"scaled: {scaled}")
        #         new_i[m] = best_index(scaled, counts, minima_i, i, ips, m)

    # print(f"--in pick_index new are: {new_i}")
    return new_i[0], new_i[1]


def best_index_frac(counts, idxs, maximum=True):
    # Examines all of the minima in idxs, and selects the minima
    # that excludes noisy pixels based on the area under the curve
    # between idxs.

    if np.size(idxs) == 1:
        return idxs[0]

    frac = 0.0001
    pixel_count = 0
    total = np.sum(counts)
    print(f"total: {total}")
    idxs.sort()
    if not maximum:
        idxs = np.flip(idxs)
    last_idx = idxs[0]
    mi = min(idxs[0], idxs[1])
    ma = max(idxs[0], idxs[1])
    last_count = np.sum(counts[mi:ma])
    best_idx = idxs[0]
    for i in idxs[1:]:
        print(f"indexes: {last_idx, i}")
        # print(f"counts: {counts[last_idx:i]}")
        mi = min(last_idx, i)
        ma = max(last_idx, i)
        c = np.sum(counts[mi:ma])
        print(f"counts: {c}")
        print(f"last_count: {last_count}")
        # print(f"count between idxs: {c}")
        # print(f"fraction: {c / total}")
        pixel_count += c
        # print(f"cumulative count: {pixel_count}")
        # print(f"cumulative fraction: {pixel_count / total}")
        count_ratio = c / last_count
        print(f"last_count / total: {last_count / total}")
        print(f"count_ratio: {count_ratio}")
        if last_count / total > frac and count_ratio < 0.1:
            # if (pixel_count / total) > frac:
            # if pixel_count > 5000:
            best_idx = last_idx
            pixel_count = 0
            # break
        print(f"best_idx: {best_idx}")
        last_idx = i
        last_count = c

    print(f"picked {best_idx}")
    return best_idx


def best_index(scaled, counts, minima_i, i1, i2, maximum=True):
    # looks at all minima between i1 and i2, and
    # and returns the one with the deepest scaled minima
    # print(f"--in best_index")

    if i1 == i2:
        return i1

    # print(scaled)
    # print(minima_i)
    # print(i1)
    # print(i2)
    minima_idx_i1 = np.where(minima_i == i1)[0]
    minima_idx_i2 = np.where(minima_i == i2)[0]
    # print(minima_idx_i1)
    # print(minima_idx_i2)
    if minima_idx_i1.size == 0:
        if i1 < i2:
            s_min_i = 0
            s_max_i = minima_idx_i2[0]
        else:
            s_min_i = minima_idx_i2[0]
            s_max_i = len(minima_i)
    elif minima_idx_i2.size == 0:
        if i1 < i2:
            s_min_i = minima_idx_i1[0]
            s_max_i = len(minima_i)
        else:
            s_min_i = 0
            s_max_i = minima_idx_i1[0]
    else:
        s_min_i = min(minima_idx_i1[0], minima_idx_i2[0])
        s_max_i = max(minima_idx_i1[0], minima_idx_i2[0])
    # print(s_min_i)
    # print(s_max_i)
    # print(scaled[s_min_i:s_max_i + 1])
    # print(np.where(scaled[s_min_i:s_max_i + 1] == np.amax(scaled[s_min_i:s_max_i + 1])))
    deepest_idxs = np.where(
        scaled[s_min_i:s_max_i + 1] == np.amax(scaled[s_min_i:s_max_i + 1])
    )[0] + s_min_i
    # print(deepest_idxs)
    if deepest_idxs.size == 1:
        print("Only one deepest idx")
        idx = minima_i[deepest_idxs[0]]
    else:
        print("Multiple deepest minima")
        low_i = minima_i[deepest_idxs[0]]
        high_i = minima_i[deepest_idxs[-1]]

        print(f"{low_i}, {high_i}")
        pixels = np.sum(counts[low_i:high_i])
        frac = pixels / np.sum(counts)
        print(frac)
        if frac > 0.1:
            print("greater than frac, keep?")
            # lots of pixels, pick outermost
            if maximum:
                idx = high_i
            else:
                idx = low_i
        else:
            print("less than frac, noise?")
            # not may pixels, probably noise?
            if maximum:
                idx = low_i
            else:
                idx = high_i

    # best_idx = np.argmax(scaled[s_min_i:s_max_i + 1])
    # print(best_idx)
    # print(scaled[s_min_i])
    # print(scaled[s_max_i])
    # idx = minima_i[s_min_i + best_idx]
    # print(f"--in best_index, idx is: {idx}")

    return idx


def eval_bounds(
    name, min_i, max_i, min_p, max_p, min_pin, max_pin, log_path=None
):
    # This is only here for temporary testing

    goal_limits = goal_dn.get(name, None)
    lines = list()
    if goal_limits is not None:
        lines.append(f"Goal DN: {goal_limits}")
        lines.append("Min Index | Prominence | Prominence in limits")
        lines.append(f"{min_i} {max_i} | {min_p} {max_p} | {min_pin} {max_pin}")
        min_i_g = min_i == goal_limits[0]
        max_i_g = max_i == goal_limits[1]
        min_p_g = min_p == goal_limits[0]
        max_p_g = max_p == goal_limits[1]
        min_pin_g = min_pin == goal_limits[0]
        max_pin_g = max_pin == goal_limits[1]
        lines.append(f"{min_i_g} {max_i_g} | {min_p_g} {max_p_g} | {min_pin_g} {max_pin_g}")

    else:
        lines.append(f"{name} isn't in list of known bounds")

    if log_path is None:
        print("\n".join(lines))
    else:
        lines.insert(0, name)
        with open(log_path, 'a+t') as f:
            f.write("\n".join(lines))
            f.write("\n\n")


if __name__ == "__main__":
    main()
