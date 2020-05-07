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


def main():
    try:
        parser = argparse.ArgumentParser(
            description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter,
            parents=[util.parent_parser()])
        parser.add_argument('-o', '--output',
                            required=False, default='.bitflip.cub')
        parser.add_argument('-w', '--width', required=False, default=5,
                            help="The number of medstd widths for bit-flip "
                            "cleaning.")
        parser.add_argument('--line', required=False, action='store_true',
                            help="Performs statistics along the line "
                            "direction instead of column for the image area.")
        parser.add_argument('-r', '--replacement', required=False, type=int,
                            help="By default, the program will replace "
                            "identified pixels with an appropriate NULL data "
                            "value, but if provided this value will be used "
                            "instead.")
        parser.add_argument('-p', '--plot', required=False,
                            action='store_true',
                            help="Displays plot for each area.")
        parser.add_argument('-n', '--dryrun', required=False,
                            action='store_true',
                            help="Does not produce a cleaned output file.")
        parser.add_argument('file', help='ISIS Cube file or PDS IMG to clean.')

        args = parser.parse_args()

        util.set_logging(args.log, args.logfile)

        out_p = util.path_w_suffix(args.output, args.file)

        clean(args.file, out_p, width=args.width,
              axis=(1 if args.line else 0), replacement=args.replacement,
              plot=args.plot, dryrun=args.dryrun, keep=args.keep)

        sys.exit(0)
    except subprocess.CalledProcessError as err:
        print('Had an ISIS error:', file=sys.stderr)
        print(' '.join(err.cmd), file=sys.stderr)
        print(err.stdout, file=sys.stderr)
        print(err.stderr, file=sys.stderr)
        sys.exit(1)
    except Exception as err:
        traceback.print_exc(file=sys.stderr)
        print(err, file=sys.stderr)
        sys.exit(1)


def clean(in_path: os.PathLike, out_path: os.PathLike, width=5,
          replacement=None, axis=0, plot=False, dryrun=False, keep=False):
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

    label = pvl.load(str(in_p))
    if 'IsisCube' in label:
        clean_cube(in_p, out_p, label, width, replacement, axis, plot,
                   dryrun, keep)
    elif 'PDS_VERSION_ID' in label:
        clean_img(in_p, out_p, label, width, replacement, axis, plot,
                  dryrun, keep)
    else:
        raise ValueError(f"The file at {in_p} is not an ISIS Cube or a "
                         "PDS IMG fie.")
    return


def clean_cube(in_p: Path, out_p: Path, label=None, width=5,
               replacement=None, axis=0, plot=False, dryrun=False, keep=False):
    """ISIS Cube version of clean().

    Please see clean() for argument details.
    """

    to_del = isis.PathSet()

    # Bit-flip correct the non-image areas.
    tblcln_p = to_del.add(in_p.with_suffix('.tableclean.cub'))
    clean_tables_from_cube(in_p, tblcln_p, width=width, plot=plot,
                           dryrun=dryrun)

    # Now clean the image area.
    if label is None:
        label = pvl.load(str(in_p))
    specialpix = getattr(isis.specialpixels,
                         label['IsisCube']['Core']['Pixels']['Type'])
    image = np.ma.masked_outside(gdal_array.LoadFile(str(in_p)),
                                 specialpix.Min, specialpix.Max)

    logging.info(f"Bit-flip cleaning Image area.")
    # These four lines are just informational.
    img_mean = np.ma.mean(image)
    img_mode = mstats.mode(image, axis=None)[0][0]
    d = img_mean - img_mode
    logging.info(f"Mean: {img_mean}, Mode: {img_mode}, diff: {d}")

    (s_min, s_max) = find_smart_window_from_ma(
        image, width=width, axis=axis,
        plot=(f"{in_p.name} Image Area" if plot else False))

    if not dryrun:
        util.log(isis.mask(tblcln_p, mask=tblcln_p, to=out_p,
                 minimum=s_min, maximum=s_max,
                 preserve='INSIDE', spixels='NONE').args)

        if replacement is not None:
            null_p = to_del.add(tblcln_p.with_suffix('.null.cub'))
            shutil.copy(out_p, null_p)
            util.log(isis.stretch(null_p, to=out_p, null=replacement))

        if not keep:
            to_del.unlink()

    return


def clean_img(in_path: Path, out_path: Path, label=None, width=5,
              replacement=None, axis=0, plot=False, dryrun=False, keep=False):
    """PDS IMG file version of clean().

    This function is currently quite slow and takes almost a minute to
    process a 50,000 line image.  Since the ability to process IMG files
    isn't a primary task, and this is just a proof-of-concept, we can live
    with slow.  It can always be optimized, if needed.

    Please see clean() for argument details.
    """

    if label is None:
        label = pvl.load(str(in_path))
    if 'PDS_VERSION_ID' not in label:
        raise ValueError(
            f"The file at {in_p} does not appear to be a PDS IMG file.")

    # Bit-flip correct the non-image areas.
    # This is going to diverge from cubes for the Buffer and Dark pixels.
    # In ISIS-land, they are 'tables' but in a IMG these are part of the
    # image array, and will need to be handled below.
    tblcln_p = in_path.with_suffix('.tableclean.img')
    clean_tables_from_img(in_path, tblcln_p, label, width, replacement,
                          plot=plot, dryrun=dryrun)

    # Now clean the image area.
    lut = img.LUT_Table(
        label['INSTRUMENT_SETTING_PARAMETERS']['MRO:LOOKUP_CONVERSION_TABLE'])
    specialpix = lut.specialpixels()
    t_name = 'IMAGE'

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

    image = np.ma.masked_outside(img_arr,
                                 specialpix.Min, specialpix.Max)

    logging.info(f"Bit-flip cleaning Image area.")
    (s_min, s_max) = find_smart_window_from_ma(
        image, width=width, axis=axis,
        plot=(f"{in_path.name} Image Area" if plot else False))

    if not dryrun:
        shutil.copy(tblcln_p, out_path)
        clean_image = np.ma.masked_outside(image, s_min, s_max)
        if replacement is not None:
            specialpix.Null = replacement
        img_arr = apply_special_pixels(clean_image, specialpix)
        img.overwrite_object(out_path, 'IMAGE', img_arr)

        if not keep:
            tblcln_p.unlink()

    return


def clean_tables_from_cube(in_path: Path, out_path: Path, width=5,
                           replacement=None,
                           rev_area=True, mask_area=False, ramp_area=False,
                           buffer_area=False, dark_area=False,
                           plot=False, dryrun=False):
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
                f"Bit-flip cleaning for {notimpl} is not yet implemented.")

    label = pvl.load(str(in_path))
    if 'IsisCube' not in label:
        raise ValueError(
            f"The file at {in_path} does not appear to be an ISIS Cube.")

    binning = label['IsisCube']['Instrument']['Summing']
    specialpix = getattr(isis.specialpixels,
                         label['IsisCube']['Core']['Pixels']['Type'])

    if not dryrun:
        shutil.copy(in_path, out_path)

    if any((rev_area, mask_area, ramp_area)):
        t_name = 'HiRISE Calibration Image'
        HCI_dict = isis.cube.get_table(in_path, t_name)
        cal_vals = np.array(HCI_dict['Calibration'])

        cal_image = np.ma.masked_outside(cal_vals,
                                         specialpix.Min, specialpix.Max)

        clean_cal = clean_cal_tables(cal_image, binning, width,
                                     rev_area, mask_area, ramp_area,
                                     (str(in_path.name) if plot else False))
        if not dryrun:
            # write the table back out
            if replacement is not None:
                specialpix.Null = replacement
            HCI_dict['Calibration'] = apply_special_pixels(
                clean_cal, specialpix).data.tolist()
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


def clean_tables_from_img(in_path: Path, out_path: Path, label=None, width=5,
                          replacement=None,
                          rev_area=True, mask_area=False, ramp_area=False,
                          buffer_area=False, dark_area=False,
                          plot=False, dryrun=False):
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
                f"Bit-flip cleaning for {notimpl} is not yet implemented.")

    if label is None:
        label = pvl.load(str(in_path))
    if 'PDS_VERSION_ID' not in label:
        raise ValueError(
            f"The file at {in_path} does not appear to be a PDS IMG file.")

    binning = label['INSTRUMENT_SETTING_PARAMETERS']['MRO:BINNING']
    specialpix = img.LUT_Table(
        label['INSTRUMENT_SETTING_PARAMETERS'][
            'MRO:LOOKUP_CONVERSION_TABLE']).specialpixels()

    if not dryrun:
        shutil.copy(in_path, out_path)

    if any((rev_area, mask_area, ramp_area)):
        t_name = 'CALIBRATION_IMAGE'
        cal_vals = img.object_asarray(in_path, t_name)

        # prefix = label[t_name]['LINE_PREFIX_BYTES']
        # suffix = label[t_name]['LINE_SUFFIX_BYTES']
        # rev_clock_slice = np.s_[:20, prefix:(-1 * suffix)]
        rev_clock_slice = np.s_[:20, :]

        # print(cal_vals.shape)
        # print(cal_vals[:20, 18:-16])
        # print(cal_vals[:20, 18:-16].shape)
        # print(specialpix)
        cal_image = np.ma.masked_outside(cal_vals[rev_clock_slice],
                                         specialpix.Min, specialpix.Max)

        # print("cal_image:")
        # print(cal_image)
        # print(cal_image.shape)

        clean_cal = clean_cal_tables(cal_image, binning, width,
                                     rev_area, mask_area, ramp_area,
                                     (str(in_path.name) if plot else False))
        if not dryrun:
            # write the table back out
            if replacement is not None:
                specialpix.Null = replacement
            cal_vals[rev_clock_slice] = apply_special_pixels(clean_cal,
                                                             specialpix)
            img.overwrite_object(out_path, t_name, cal_vals)
    return


def clean_cal_tables(cal_image, binning, width=5,
                     rev_area=True, mask_area=False, ramp_area=False,
                     plot=False):
    # Deal with the HiRISE Calibration Image first (Reverse-clock, Mask,
    # and Ramp

    for notimpl in (mask_area, ramp_area):
        if notimpl is True:
            raise NotImplementedError(
                f"Bit-flip cleaning for {notimpl} is not yet implemented.")

    mask_lines = int(20 / binning)

    if rev_area:
        logging.info(f"Bit-flip cleaning Reverse-Clock area.")
        rev_clean = clean_array(
            cal_image[:20, :], width=width, axis=1,
            plot=(f"{plot} Reverse-Clock" if plot else False))
        cal_image[:20, :] = rev_clean

    # if mask_area:
    #     mask_pixels = cal_image[20:mask_lines, :]

    # if ramp_area:
    #     ramp_pixels = cal_image[20 + mask_lines:, :]

    return cal_image


def clean_array(data: np.ma.array, width=5, axis=0, plot=False):
    """Returns a numpy masked array whose mask is based on applying
    the smart window bounds from find_smart_window_from_ma() applied
    to *data* with the specified *wdith* and *axis*.

    If all of the values in *data* are masked, a ValueError is raised.
    """
    if np.all(data.mask):
        logging.info("All of the values in data are masked.")
        return data

    (w_min, w_max) = find_smart_window_from_ma(data, width=width, axis=axis,
                                               plot=plot)
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
        raise ValueError("The provided array does not have two dimensions, "
                         f"it has {data.ndim}.")

    xx, yy = np.meshgrid(np.arange(data.shape[0]),
                         np.arange(data.shape[1]), indexing='ij')
    x1 = xx[~data.mask]
    y1 = yy[~data.mask]
    good = data[~data.mask]

    interp = interpolate.griddata((x1, y1), good.ravel(), (xx, yy),
                                  method='nearest')
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
        if val in (specialpix.Null, specialpix.Lrs, specialpix.Lis,
                   specialpix.His, specialpix.Hrs):
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


def find_smart_window_from_ma(data: np.ma.array, width=5, axis=0, plot=False):
    """Returns a two-tuple with the result of find_smart_window().

    This function mostly just does set-up based on the provided
    masked array *data* and the provided *width* and *axis* value
    to calculate the inputs for find_smart_window().
    """
    median = median_limit(np.ma.median(data), data)
    medstd = median_std_from_ma(data, axis=axis)

    unique, unique_counts = np.unique(data.compressed(), return_counts=True)

    mindn, maxdn, ex = min_max_ex(median, medstd, width)

    return find_smart_window(unique, unique_counts, mindn, maxdn, median,
                             central_exclude_dn=ex, plot=plot)


def median_limit(median, data: np.ndarray, limit=4000):
    """Return the 'best' median of the provided data.

    If the extracted median is larger than *limit*, the algorithm
    will attempt to find a better representation of the median.  If
    it cannot, it will return the median, even if larger than
    *limit*.
    """
    logging.info(f"Median: {median}")

    if median > limit and np.any([data < limit]):
        median = np.median(data[data < limit])
        logging.info(f"The median was too high (> {limit}), "
                     f"found a better one: {median}.")

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
    logging.info("Median standard deviation of elements along the axis "
                 f"that have the maximum valid pixel count: {medstd}.")
    return medstd


def min_max_ex(central, medstd, width) -> tuple:
    """Return a minimum, maximum, and exclusion value based on the
    provided *central* value, and adding and subtracting the result
    of multiplying the *medstd* by the *width*.

    If the medstd is too large, the resulting values are based on
    conservative guesses.
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
    if medstd > 300:
        medstd = 64
        ex = 16
        logging.info("The derived medstd was too big, setting the medstd "
                     f"to {medstd} and the exclusion value to {ex}.")
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


def find_select_idx_name(central_idx: int, limit_idx: int,
                         close_to_limit: bool):
    """Returns a two-tuple which contains an int of value 0 or -1 in
    the first posiion and a string of value 'min' or 'max' in the second.

    The string is determined based on the relative values of
    *central_idx* to *limit_idx*, and the value of the int is meant
    as an index to a list, based on the three input variables.

    Raises *ValueError* if *central_idx* and *limit_idx* are the same.

    This is primarily a helper function to find_minima_index().
    """
    select_idx = None
    idx_name = None

    if limit_idx < central_idx:
        idx_name = 'min'
        if close_to_limit is False:
            select_idx = -1
        else:
            select_idx = 0
    elif limit_idx > central_idx:
        idx_name = 'max'
        if close_to_limit is False:
            select_idx = 0
        else:
            select_idx = -1
    else:
        raise ValueError

    return select_idx, idx_name


def find_minima_index(central_idx: int, limit_idx: int,
                      minima_idxs: np.ndarray, pixel_counts: np.ndarray,
                      close_to_limit=True) -> int:
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
    info_str = 'Looking for {} {} range ...'

    try:
        select_idx, idx_name = find_select_idx_name(central_idx, limit_idx,
                                                    close_to_limit)
        logging.info(info_str.format(idx_name, 'inside'))

        min_i = min(central_idx, limit_idx)
        max_i = max(central_idx, limit_idx)
        inrange_i = minima_idxs[
            (minima_idxs >= min_i) * (minima_idxs < max_i)]
        logging.info(str(inrange_i) + ' are the indexes inside the range.')

        value = min(np.take(pixel_counts, inrange_i))
        logging.info(f'{value} is the minimum Pixel count amongst those '
                     'indexes.')

        idx = inrange_i[np.asarray(pixel_counts[inrange_i] ==
                                   value).nonzero()][select_idx]
    except ValueError:
        try:
            if limit_idx < central_idx:
                logging.info(info_str.format('min', 'outside'))
                idx = max(minima_idxs[minima_idxs < limit_idx])
            elif limit_idx > central_idx:
                logging.info(info_str.format('max', 'outside'))
                idx = min(minima_idxs[minima_idxs > limit_idx])
            else:
                raise ValueError

            logging.info(f'{idx} is the minimum index.')
        except ValueError:
            logging.info('Could not find a valid minima, returning '
                         f'the limit index: {limit_idx}.')
            idx = limit_idx

    return idx


def find_smart_window(dn: np.ndarray, counts: np.ndarray,
                      mindn: int, maxdn: int,
                      centraldn: int, central_exclude_dn=0,
                      plot=False, closest=True) -> tuple:
    '''Returns a minimum and maximum DN value from *dn* which are
       based on using the find_minima_index() function with the
       given *mindn*, *maxdn*, and *centraldn* values.  The *dn*
       array must be a sorted list of unique values, and *counts*
       is the number of times each of the values in *dn* occurs.

       If *central_exclude_dn* is given, the returned minimum and
       maximum DN are guaranteed to be at least *central_exclude_dn*
       away from *centraldn*.  This is useful if you don't want returned
       minimum and maximum DN to be too close to the *centraldn*.

       If plot is True, this function will display a plot describing
       its work.  The curve represents the hist values, the shaded
       area marks the window between the given mindn and maxdn.  The
       'x'es mark all the detected minima, and the red dots indicate
       the minimum and maximum DN values that this function will
       return.

       The value of *closest* is passed on to find_minima_index().
    '''
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

    minima_i, _ = find_peaks(np.negative(counts))

    # print(f'central_i {central_i}')
    # print(f'central_min_i {central_min_i}')
    # print(f'central_max_i {central_max_i}')
    # print(f'mindn_i {mindn_i}')
    # print(f'maxdn_i {maxdn_i}')

    min_i = find_minima_index(central_min_i, mindn_i, minima_i, counts,
                              close_to_limit=closest)
    max_i = find_minima_index(central_max_i, maxdn_i, minima_i, counts,
                              close_to_limit=closest)
    logging.info(f'indexes: {min_i}, {max_i}')
    logging.info(f'DN window: {dn[min_i]}, {dn[max_i]}')

    if plot:
        import matplotlib.pyplot as plt

        plt.ioff()
        fig, (ax0, ax1) = plt.subplots(2, 1)

        indices = np.arange(0, len(counts))
        dn_window = np.fromiter(map((lambda i: i >= mindn_i and i <= maxdn_i),
                                    (x for x in range(len(counts)))),
                                dtype=bool)

        if isinstance(plot, str):
            fig.suptitle(plot)

        ax0.set_ylabel('Pixel Count')
        ax0.set_xlabel('DN Index')
        ax0.set_yscale('log')
        ax0.fill_between(indices, counts, where=dn_window,
                         color='lightgray')
        ax0.axvline(x=central_i, c='gray')
        ax0.axvline(x=np.argmax(counts), c='lime', ls='--')
        ax0.plot(counts)
        ax0.plot(minima_i, counts[minima_i], "x")
        ax0.plot(min_i, counts[min_i], "o", c='red')
        ax0.plot(max_i, counts[max_i], "o", c='red')

        ax1.set_ylabel('Pixel Count')
        ax1.set_xlabel('DN')
        ax1.set_yscale('log')
        ax1.set_ybound(lower=0.5)
        ax1.scatter(dn, counts, marker='.', s=1)
        ax1.axvline(x=dn[min_i], c='red')
        ax1.axvline(x=dn[max_i], c='red')

        plt.show()

    return (dn[min_i], dn[max_i])


if __name__ == "__main__":
    main()
