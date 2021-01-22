#!/usr/bin/env python
"""
Perform the jitter derivation for a HiRISE observation.

Outputs are the jitter in the x (sample) and y (line) directions in
pixels, realt (ephemeris time at that translation), average error
between the derived jitter function and the original data, linerate
(line time read from flat files), and TDI. These outputs are
written to a text file. Another output is the pixel smear, also
written to a file.
"""

# Copyright (C) 2013-2020 Arizona Board of Regents on behalf of the Lunar and
# Planetary Laboratory at the University of Arizona.
#   - Original MatLab program written by Aaron Boyd and Sarah Mattson for
#     HiROC as part of the process to describe and correct geometric
#     distortions caused by jitter in HiRISE images.
#     Original version written approximately 6/2008.
#
#  Copyright (c) 2012, United States Government as represented by the
#  Administrator of the National Aeronautics and Space Administration. All
#  rights reserved.
#   - C++ version written by Oleg Alexandrov based on the above MatLab
#     program, resolveJitter4HiJACK.m version 1.4
#
# Copyright 2020, Ross A. Beyer (rbeyer@seti.org)
#   - Elements of this Python program are are based on the C++ version but
#     the logic here is rewritten from scratch to emulate functionality.
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
import csv
import itertools
import logging
import math
import os
import statistics
import sys
from collections import abc
from pathlib import Path

import numpy as np
from scipy.signal import medfilt
from scipy.interpolate import PchipInterpolator

import hiproc.util as util
from hiproc.FlatFile import FlatFile


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     parents=[util.parent_parser()])
    parser.add_argument('image_location', type=Path)
    parser.add_argument("image_id", type=str)
    parser.add_argument("line_interval", type=float)
    parser.add_argument("file_path1", type=Path)
    parser.add_argument("which1", type=int, choices=[-1, 1])
    parser.add_argument("file_path2", type=Path)
    parser.add_argument("which2", type=int, choices=[-1, 1])
    parser.add_argument("file_path3", type=Path)
    parser.add_argument("which3", type=int, choices=[-1, 1])

    args = parser.parse_args()

    util.set_logging(args.log)

    if args.line_interval <= 0:
        raise ValueError("The parameter 'line_interval' must be positive.")

    start(
        args.image_location, args.image_id, args.line_interval,
        args.file_path1, args.which1, args.file_path2, args.which2,
        args.file_path3, args.which3
    )
    return


def start(
    image_location: Path, image_id: str, line_interval: float,
    file_path1: Path, which1: int, file_path2: Path, which2: int,
    file_path3: Path, which3: int
):
    fp1 = set_file_path(image_location, file_path1)
    fp2 = set_file_path(image_location, file_path2)
    fp3 = set_file_path(image_location, file_path3)

    window_size = 11
    window_width = 2

    # The first file to be parsed sets et, et_shift, t0, and duration
    (tdi1, nfft1, linerate1, tt1, et, et_shift, t0, duration, dt1, ddt1, t1,
     offx1, offy1,
     xinterp1, yinterp1, x1, y1, overxx1, overyy1) = parse_file(
        fp1, which1, line_interval, window_size, window_width
    )
    (tdi2, nfft2, linerate2, tt2, _, _, _, _, dt2, ddt2, t2, offx2, offy2,
     xinterp2, yinterp2, x2, y2, overxx2, overyy2) = parse_file(
        fp2, which2, line_interval, window_size, window_width,
        et, et_shift, t0, duration
    )
    (tdi3, nfft3, linerate3, tt3, _, _, _, _, dt3, ddt3, t3, offx3, offy3,
     xinterp3, yinterp3, x3, y3, overxx3, overyy3) = parse_file(
        fp3, which3, line_interval, window_size, window_width,
        et, et_shift, t0, duration
    )

    if np.array_equal(tt1, tt2) and np.array_equal(tt2, tt3):
        tt = tt3
    else:
        raise ValueError("The values of tt are not identical.")

    if nfft1 == nfft2 == nfft3:
        nfft = nfft3
    else:
        raise ValueError("The values of nfft are not identical.")

    if tdi1 == tdi2 == tdi3:
        tdi = tdi3
    else:
        raise ValueError("The values of tdi are not identical.")

    if linerate1 == linerate2 == linerate3:
        linerate = linerate3
    else:
        raise ValueError("The values of linerate are not identical.")

    # starting a loop to find correct phase tol
    # int      k           = 0;
    # double   minAvgError = numeric_limits<double>::max();
    # int      minK        = 0;
    # double   error       = errorTol;
    # ArrayXd  minJitterX, minJitterY;
    # MatrixXd overxxx1, overyyy1, overxxx2, overyyy2, overxxx3, overyyy3;

    # Tolerance for the error
    error_tol = 0.0000000001
    repetitions = 50

    # Tolerance coefficient for the phase difference
    tolcoef = 0.01

    k = 0
    error = error_tol
    min_avg_error = sys.float_info.max
    # min_k = 0
    while error >= error_tol and k < repetitions:
        k += 1

        # setting the phase tolerance
        phasetol = k * tolcoef

        # null the frequencies that cause a problem (really zero them out)
        overxxx1, overyyy1 = make_null_problematic_frequencies(
            phasetol, ddt1, overxx1, overyy1
        )
        overxxx2, overyyy2 = make_null_problematic_frequencies(
            phasetol, ddt2, overxx2, overyy2
        )
        overxxx3, overyyy3 = make_null_problematic_frequencies(
            phasetol, ddt3, overxx3, overyy3
        )

        # Overwrite null frequencies.
        overxxx1 = overwrite_null_freq(overxxx1, overxxx2, overxxx3)
        overxxx2 = overwrite_null_freq(overxxx2, overxxx1, overxxx3)
        overxxx3 = overwrite_null_freq(overxxx3, overxxx1, overxxx2)

        overyyy1 = overwrite_null_freq(overyyy1, overyyy2, overyyy3)
        overyyy2 = overwrite_null_freq(overyyy2, overyyy1, overyyy3)
        overyyy3 = overwrite_null_freq(overyyy3, overyyy1, overyyy2)

        # Adding all frequencies together
        stackedx = np.stack((overxxx1, overxxx2, overxxx3))
        stackedy = np.stack((overyyy1, overyyy2, overyyy3))

        overxxx = stackedx.mean(axis=0)
        overyyy = stackedy.mean(axis=0)

        # take the sum of each row
        overx = overxxx.sum(axis=0)
        overy = overyyy.sum(axis=0)

        jitterx = overx - overx[0]
        jittery = overy - overy[0]

        # checking
        jittercheckx1 = np.interp(
            tt + dt1/duration, tt, jitterx, left=0, right=0
        ) - jitterx
        jitterchecky1 = np.interp(
            tt + dt1/duration, tt, jittery, left=0, right=0
        ) - jittery

        jittercheckx2 = np.interp(
            tt + dt2/duration, tt, jitterx, left=0, right=0
        ) - jitterx
        jitterchecky2 = np.interp(
            tt + dt2/duration, tt, jittery, left=0, right=0
        ) - jittery

        jittercheckx3 = np.interp(
            tt + dt3/duration, tt, jitterx, left=0, right=0
        ) - jitterx
        jitterchecky3 = np.interp(
            tt + dt3/duration, tt, jittery, left=0, right=0
        ) - jittery

        error_vec = 1.0/6.0*(
            abs(xinterp1 - (jittercheckx1 + np.real(x1[0])/2.0)) +
            abs(xinterp2 - (jittercheckx2 + np.real(x2[0])/2.0)) +
            abs(xinterp3 - (jittercheckx3 + np.real(x3[0])/2.0)) +
            abs(yinterp1 - (jitterchecky1 + np.real(y1[0])/2.0)) +
            abs(yinterp2 - (jitterchecky2 + np.real(y2[0])/2.0)) +
            abs(yinterp3 - (jitterchecky3 + np.real(y3[0])/2.0))
        )

        error = error_vec.mean()
        logging.info(f"Error 1: {error}")

        if error < min_avg_error:
            min_avg_error = error
            # min_k = k
            min_jitterx = jitterx
            min_jittery = jittery

    # end while

    k = 0
    min_k = 0
    error = error_tol
    min_avg_error = sys.float_info.max

    # The jitter in the x (sample) and y (line) directions in pixels.
    # This is the jitter with minimum error, after scanning through all
    # frequencies omega.
    # ArrayXd Sample, Line;
    # ArrayXd  minJitterCheckX1, minJitterCheckX2, minJitterCheckX3;
    # ArrayXd  minJitterCheckY1, minJitterCheckY2, minJitterCheckY3;

    # starting a loop to find the correct filter size
    while error >= error_tol and k < repetitions:

        k += 1

        omega = k - 1
        c = omega/(2.0 * nfft)

        jitterxx = filter_data(nfft, c, min_jitterx)
        jitteryy = filter_data(nfft, c, min_jittery)

        jitterxx = jitterxx - jitterxx[0]
        jitteryy = jitteryy - jitteryy[0]

        # checking
        jittercheckx1 = np.interp(
            tt + dt1/duration, tt, jitterxx, left=0, right=0
        ) - jitterxx
        jitterchecky1 = np.interp(
            tt + dt1/duration, tt, jitteryy, left=0, right=0
        ) - jitteryy

        jittercheckx2 = np.interp(
            tt + dt2/duration, tt, jitterxx, left=0, right=0
        ) - jitterxx
        jitterchecky2 = np.interp(
            tt + dt2/duration, tt, jitteryy, left=0, right=0
        ) - jitteryy

        jittercheckx3 = np.interp(
            tt + dt3/duration, tt, jitterxx, left=0, right=0
        ) - jitterxx
        jitterchecky3 = np.interp(
            tt + dt3/duration, tt, jitteryy, left=0, right=0
        ) - jitteryy

        error_vec = 1.0/6.0*(
            abs(xinterp1 - (jittercheckx1 + np.real(x1[0])/2.0)) +
            abs(xinterp2 - (jittercheckx2 + np.real(x2[0])/2.0)) +
            abs(xinterp3 - (jittercheckx3 + np.real(x3[0])/2.0)) +
            abs(yinterp1 - (jitterchecky1 + np.real(y1[0])/2.0)) +
            abs(yinterp2 - (jitterchecky2 + np.real(y2[0])/2.0)) +
            abs(yinterp3 - (jitterchecky3 + np.real(y3[0])/2.0))
        )

        error = error_vec.mean()
        logging.info(f"error 2: {error}")

        if error < min_avg_error:
            min_k = k
            min_avg_error = error
            sample = jitterxx
            line = jitteryy
            min_jitter_check_x1 = jittercheckx1
            min_jitter_check_x2 = jittercheckx2
            min_jitter_check_x3 = jittercheckx3
            min_jitter_check_y1 = jitterchecky1
            min_jitter_check_y2 = jitterchecky2
            min_jitter_check_y3 = jitterchecky3
    # end while

    # The outputs
    print(f" Average error is {min_avg_error} at min index {min_k}")
    print(f"linerate is {linerate}")
    print(f"TDI = {tdi}")

    max_smear_sample, max_smear_line, max_smear_mag = pixel_smear(
        sample, line, et, linerate, tdi, image_location, image_id
    )

    # To do: Remove the py suffix from jitter and smear!

    # Make a text file of the jitter data
    jitter_p = image_location / (image_id + "_jitter_py.txt")
    jitter_text = [f"""# Using image {image_id} the jitter was found with an
# Average Error of {min_avg_error}
# Maximum Cross-track pixel smear {max_smear_sample}
# Maximum Down-track pixel smear {max_smear_line}
# Maximum Pixel Smear Magnitude {max_smear_mag}
#
# Sample                 Line                   ET"""]
    for s, l, e in zip(sample, line, et):
        jitter_text.append(f"""{s}     {l}     {e}""")

    logging.info(f"Writing: {jitter_p}")
    jitter_p.write_text("\n".join(jitter_text))

    # I think we could re-do this for matplotlib.  Let's defer it for now.
    # Writing the data we will plot later in gnuplot
    data_p = image_location / (image_id + "_jitter_plot_py.txt")
    t1_shift = t1 - t0
    t2_shift = t2 - t0
    t3_shift = t3 - t0
    jittercheckx1_shift = min_jitter_check_x1 + np.real(x1[0])/2.0
    jitterchecky1_shift = min_jitter_check_y1 + np.real(y1[0])/2.0
    jittercheckx2_shift = min_jitter_check_x2 + np.real(x2[0])/2.0
    jitterchecky2_shift = min_jitter_check_y2 + np.real(y2[0])/2.0
    jittercheckx3_shift = min_jitter_check_x3 + np.real(x3[0])/2.0
    jitterchecky3_shift = min_jitter_check_y3 + np.real(y3[0])/2.0
    # ArrayXd* Data[] =
    #   {&ET_shift, &Sample, &Line,
    #    &t1_shift,
    #    &offx1, &xinterp1, &jittercheckx1_shift,
    #    &offy1, &yinterp1, &jitterchecky1_shift,
    #    &t2_shift,
    #    &offx2, &xinterp2, &jittercheckx2_shift,
    #    &offy2, &yinterp2, &jitterchecky2_shift,
    #    &t3_shift,
    #    &offx3, &xinterp3, &jittercheckx3_shift,
    #    &offy3, &yinterp3, &jitterchecky3_shift
    #   };

    # Note the comment before the first label
    string_labels = [
        "# ET_shift", "Sample", "Line",
        "t1_shift", "offx1", "xinterp1", "jittercheckx1_shift",
        "offy1", "yinterp1", "jitterchecky1_shift",
        "t2_shift", "offx2", "xinterp2", "jittercheckx2_shift",
        "offy2", "yinterp2", "jitterchecky2_shift",
        "t3_shift", "offx3", "xinterp3", "jittercheckx3_shift",
        "offy3", "yinterp3", "jitterchecky3_shift"
    ]

    # int numData = sizeof(Data)/sizeof(ArrayXd*);
    write_data_for_plotting(
        data_p, string_labels, et_shift, sample, line, t1_shift,
        offx1, xinterp1, jittercheckx1_shift, offy1, yinterp1,
        jitterchecky1_shift, t2_shift, offx2, xinterp2, jittercheckx2_shift,
        offy2, yinterp2, jitterchecky2_shift, t3_shift, offx3, xinterp3,
        jittercheckx3_shift, offy3, yinterp3, jitterchecky3_shift
    )

    write_csv(
        image_location / (image_id + "_jitter_plot_py.csv"),
        string_labels,
        et_shift, sample, line, t1_shift,
        offx1, xinterp1, jittercheckx1_shift, offy1, yinterp1,
        jitterchecky1_shift, t2_shift, offx2, xinterp2, jittercheckx2_shift,
        offy2, yinterp2, jitterchecky2_shift, t3_shift, offx3, xinterp3,
        jittercheckx3_shift, offy3, yinterp3, jitterchecky3_shift
    )

    gnuplot_p = image_location / (image_id + "_jitter_plot_py.plt")
    img_file_name = image_location / (image_id + "_jitter_plot_py.png")
    write_gnuplot_file(
        gnuplot_p, data_p, img_file_name, file_path1, file_path2, file_path3
    )

    return


def get_time_arrays(nfft, time_start, time_stop):
    # making time regularly spaced and interpolate offx and offy
    # into this time
    duration = time_stop - time_start
    nfft_arr = np.linspace(0, nfft - 1, nfft)
    et = nfft_arr * duration / (nfft - 1.0) + time_start
    et_shift = et - time_start
    return et, et_shift


def create_matrices(
    nfft: int, dt: float, time: np.array, tt: abc.Sequence,
    offx: np.array, offy: np.array, et_shift, t0, duration
):
    t_shift = time - t0

    fx = PchipInterpolator(t_shift, offx, extrapolate=False)
    fy = PchipInterpolator(t_shift, offy, extrapolate=False)
    xinterp = fx(et_shift)
    yinterp = fy(et_shift)

    np.nan_to_num(xinterp, copy=False, nan=np.mean(offx))
    np.nan_to_num(yinterp, copy=False, nan=np.mean(offy))

    # getting the frequencies of the Fourier transform
    freq = np.linspace(0, nfft / 2 - 1, int(nfft / 2))

    # taking the fourier transform of the offsets
    x = 2 * np.fft.fft(xinterp) / nfft
    y = 2 * np.fft.fft(yinterp) / nfft

    # separating sines and cosines
    xa = x[:int(nfft / 2)].real
    xb = -1 * x[:int(nfft / 2)].imag
    ya = y[:int(nfft / 2)].real
    yb = -1 * y[:int(nfft / 2)].imag

    # calculates the phase difference
    # ddt = mod(dt / duration * 2 * pi. * freq, 2 * pi);
    twopi = math.pi * 2
    ddt = (dt / duration * twopi) * freq
    for i, s in enumerate(ddt):
        ddt[i] -= twopi * math.floor(s / twopi)

    # the coeficients for the frequencies
    aaax = -0.5 * (-1 * xa * np.cos(ddt) + np.sin(ddt) * xb - xa) / np.sin(ddt)
    aaay = -0.5 * (-1 * ya * np.cos(ddt) + np.sin(ddt) * yb - ya) / np.sin(ddt)
    bbbx = -0.5 * (xb * np.cos(ddt) + np.sin(ddt) * xa + xb) / np.sin(ddt)
    bbby = -0.5 * (yb * np.cos(ddt) + np.sin(ddt) * ya + yb) / np.sin(ddt)

    # create series of sines and cosines
    overxx = np.empty([len(freq), len(tt)])
    overyy = np.empty([len(freq), len(tt)])
    for col, t in enumerate(tt):
        for row, f in enumerate(freq):
            ft = twopi * f * t
            sn = math.sin(ft)
            cn = math.cos(ft)
            overxx[row, col] = aaax[row] * sn + bbbx[row] * cn
            overyy[row, col] = aaay[row] * sn + bbby[row] * cn

    # Outputs
    # ArrayXd & tt, ArrayXd & ET, ArrayXd & ET_shift,
    # ArrayXd & ddt,
    # ArrayXd & xinterp, ArrayXd & yinterp,
    # ArrayXcd & X, ArrayXcd & Y, MatrixXd & overxx, MatrixXd & overyy
    #
    return ddt, overxx, overyy, xinterp, yinterp, x, y


def upper_power_of_two(value) -> int:
    """Returns the smallest power of 2 which is >= *value*."""
    result = 1
    while result < value:
        result <<= 1
    return result


def filter_data(nfft: int, c: float, data: np.array) -> np.array:
    # Apply a Gaussian filter to the data in the frequency domain.

    if len(data.shape) > 1:
        raise IndexError("The data array can only be 1D.")

    # Use padding so the data is not distorted.
    front_padding = math.floor(nfft - len(data)/ 2)
    back_padding = math.ceil(nfft - len(data)/ 2)

    # Apply the padding
    padded = np.concatenate((
        np.array([data[0]] * front_padding),
        data,
        np.array([data[-1]] * back_padding)
    ))

    freq = np.fft.fft(padded)

    # The exponential
    exp_vec = np.hstack(
        (np.linspace(0, nfft - 1, nfft), np.linspace(-1 * nfft, -1, nfft))
    )

    exponential = np.exp(-1 * c**2 * exp_vec**2)

    # The ifft of the product
    filtered = np.fft.ifft(freq * exponential)

    # Remove the padding and take the real part
    return filtered[front_padding:front_padding + len(data)].real


def make_null_problematic_frequencies(
    phasetol: float, ddt: np.array, x: np.array, y: np.array
):
    """Returns *x* and *y* as numpy masked arrays with 'problematic frequencies'
    masked.  It is assumed that *ddt* is a 1-D array which is the same
    size as axis 0 of *x* and *y*."""
    # Original function also returned a set of the indexes that were zeroed,
    # but I don't think we need that.  This could also use masked arrays
    # maybe?

    # TODO: x and y are 2d arrays.  Yes I want to use masked arrays, need to
    # figure out how to create the right kind of mask.

    if x.shape != y.shape:
        raise ValueError(
            f"The shape of x {x.shape} and y {y.shape} must be the same."
        )

    if x.shape[0] != ddt.size:
        raise ValueError(
            f"The size of ddt ({ddt.size}) must be the same as the first axis"
            f"of the x and y arrays {x.shape}"
        )

    # null the frequencies that cause a problem

    twopi = 2 * math.pi

    a = np.less(np.abs(ddt), phasetol)
    b = np.greater(np.abs(ddt), twopi - phasetol)
    null_positions = np.logical_or(a, b)

    # The assumption is that since ddt was a 1D array that
    # null_positions will be too, so we must reshape to a 2D
    # column, and then tile that across, so that each row
    # of the 2D matrix has the same value for all positions.
    null_2d = np.tile(null_positions.reshape(-1, 1), (1, x.shape[1]))

    x_masked = np.ma.array(x, mask=null_2d, fill_value=0)
    y_masked = np.ma.array(y, mask=null_2d, fill_value=0)

    return x_masked, y_masked


def overwrite_null_freq(
    array: np.ma.array, support1: np.ma.array, support2: np.ma.array
):
    """Returns a version of *array* with some masked values corrected
    by the average of *support1* and *support2*.

    For each masked value in *arr1*, an average of the non-null
    values in *arr2* and *arr3* are substituted.  If both *arr2*
    and *arr3* are also masked at that location, the mask will remain.
    """
    if not (array.shape == support1.shape == support2.shape):
        raise IndexError(
            f"The shape of the input array {array.shape}, and the shapes of"
            f"the support1 {support1.shape} and support2 {support2.shape}"
            f"must be the same."
        )

    stacked = np.ma.stack((support1, support2))
    smean = np.ma.mean(stacked, axis=0)
    # We put the original values into the smean array because if the input
    # array is int, then the averages do not carry through.  Instead we
    # put the "original" values into the smean array where they belong.
    smean[~array.mask] = array[~array.mask]
    return smean


def parse_file(
    file_path: os.PathLike, which: int, line_interval: float,
    window_size: int, window_width: int, et=None, et_shift=None,
    t0=None, duration=None
):
    """Parse the current text file, then extract and process the data into a
    matrix."""
    logging.info(f"Reading: {file_path}")

    flat = FlatFile(file_path)

    if flat["FROM"]["TdiMode"] == flat["MATCH"]["TdiMode"]:
        tdi = int(flat["FROM"]["TdiMode"])
    else:
        raise ValueError(
            f"The TdiMode is different for FROM ({flat['FROM']['TdiMode']}) "
            f"and MATCH ({flat['MATCH']['TdiMode']}) in {file_path}"
        )

    if flat["FROM"]["LineRate"] == flat["MATCH"]["LineRate"]:
        line_rate = float(flat["FROM"]["LineRate"].split()[0])
    else:
        raise ValueError(
            f"The LineRate is different for FROM ({flat['FROM']['LineRate']}) "
            f"and MATCH ({flat['MATCH']['LineRate']}) in {file_path}"
        )

    if which == 1:
        column = "MatchTime"
    else:
        column = "FromTime"

    # dt = which * (data[0][0] - data[0][3]);
    dt = which * (float(flat[0]["FromTime"]) - float(flat[0]["MatchTime"]))

    time = list()
    offset_x = list()
    offset_y = list()

    for row in flat:
        time.append(float(row[column]))
        offset_x.append(
            which * (float(row["RegSamp"]) - float(row["FromSamp"]))
        )
        offset_y.append(
            which * (float(row["RegLine"]) - float(row["FromLine"]))
        )

    offx_arr = np.array(offset_x)
    offy_arr = np.array(offset_y)

    magnitude = np.sqrt(offx_arr**2 + offy_arr**2)
    avemag = medfilt(magnitude, window_size)

    # Throw out the out-of-range values:
    good_data = list()
    for the_t, the_x, the_y, mag, ave in zip(
        time, offset_x, offset_y, magnitude, avemag
    ):
        high = ave + window_width
        low = ave - window_width
        if low < mag < high:
            good_data.append((the_t, the_x, the_y))

    # Some flat files have more than one measurement for the same timestamp.
    # For flat files that do not, this groupby loop is just a fancy copy.
    # However, for those that do, the multiple x and y offsets for the same
    # timestamp are averaged together.
    #
    # Also, the original MatLab code had no functions, but repeated code,
    # and there was a transcription error on the third pass of this averaging
    # such that the very first row was never included in the averaging.
    # The C++ code allowed for that broken behavior, but we won't here.
    t = list()
    offx = list()
    offy = list()
    for k, g in itertools.groupby(good_data, lambda val: val[0]):
        t.append(k)
        _, exes, whys = zip(*g)
        offx.append(statistics.mean(exes))
        offy.append(statistics.mean(whys))

    t_arr = np.array(t)

    nfft = upper_power_of_two(int(flat["FROM"]["Lines"]) / line_interval)
    filtered_x = filter_data(nfft, 2 / nfft, np.array(offx))
    filtered_y = filter_data(nfft, 2 / nfft, np.array(offy))

    tt = np.linspace(0, nfft - 1, nfft) / nfft

    if et is None:
        et, et_shift = get_time_arrays(nfft, t_arr[0], t_arr[-1])

    if t0 is None:
        t0 = t_arr[0]
        duration = t_arr[-1] - t0

    ddt, overxx, overyy, xinterp, yinterp, x, y, = create_matrices(
        nfft, dt, t_arr, tt, filtered_x, filtered_y, et_shift, t0, duration
    )

    return (tdi, nfft, line_rate,
            # t[0], t[-1] - t[0],
            tt, et, et_shift, t0, duration,
            dt, ddt,
            t_arr, filtered_x, filtered_y,
            xinterp, yinterp, x, y, overxx, overyy)
    # Outputs
    # TDI, nfft, linerate, t0, duration, tt, ET, ET_shift, dt, ddt,
    # t, offx, offy, xinterp, yinterp, X, Y, overxx, overyy
    # );


def pixel_smear(
    sample: np.array, line: np.array, t: np.array, linerate: float,
    tdi: int, image_location: os.PathLike, image_id: str
):
    # Inputs
    # ArrayXd const& Sample, ArrayXd const& Line, ArrayXd const& T,
    # double linerate, int TDI, string imageLocation, string imageId,
    # Outputs
    # double & maxSmearS, double & maxSmearL, double & maxSmearMag

    #  pixelSmear finds max smeared pixel amount in the image from derived
    #  jitter function.
    #
    #  Sample is the sample offsets from the jitter function
    #  Line is the line offsets from the jitter function
    #  T is the ephemeris time from the jitter function
    #  linerate is read from the flat file, equivalent to TDI
    #
    #  Pixel smear due to jitter is calculated by interpolating the jitter
    #  function at intervals equivalent to the linerate.  Then the
    #  difference is taken over that interval and multiplied by the TDI.
    #  This provides an estimated minimum for pixel smear. If the motion
    #  that caused the smear is not captured in the jitter derivation,
    #  then it cannot be plotted here.

    # The array T has large values. Use a shifted version of it when
    # interpolating to reduce the effect of those values on numerical
    # accuracy.
    shifted_t = t - t[0]

    # xi = T(1):linerate:T(end);
    n = math.floor((shifted_t[-1] - shifted_t[0])/linerate) + 1
    xi = np.linspace(0, n-1, n) * linerate + shifted_t[0]

    # Interpolate the jitter function at intervals equivalent to the linerate
    f_samp = PchipInterpolator(shifted_t, sample, extrapolate=False)
    f_line = PchipInterpolator(shifted_t, line, extrapolate=False)
    yis = f_samp(xi)
    yil = f_line(xi)

    np.nan_to_num(yis, copy=False, nan=0)
    np.nan_to_num(yil, copy=False, nan=0)

    # Undo the earlier shift
    xi += t[0]

    # Calculate the rate of change with respect to the linerate
    # in the sample direction
    dysdx = np.diff(yis) * tdi
    # in the line direction
    dyldx = np.diff(yil) * tdi

    # Calculate the magnitude of the smear
    mag_smear = np.sqrt(dysdx**2 + dyldx**2)

    # Find maxSmearS, the largest element by magnitude in dysdx
    max_smear_s = max(dysdx)

    # Find maxSmearL, the largest element by magnitude in dyldx
    max_smear_l = max(dyldx)

    # Find maxSmearMag, the largest element by magnitude in magSmear
    max_smear_mag = max(mag_smear)

    # Make a text file of the smear data
    smear_p = Path(image_location) / (image_id + "_smear_py.txt")
    smear_text = [
        f"""\
# Smear values are calculated from the derived jitter function for {image_id}.
# Maximum Cross-track pixel smear {max_smear_s}
# Maximum Down-track pixel smear {max_smear_l}
# Maximum Pixel Smear Magnitude {max_smear_mag}
# Sample                 Line                   EphemerisTime"""
    ]

    for ess, ell, exi in zip(dysdx, dyldx, xi):
        smear_text.append(f"{ess}     {ell}     {exi}")

    logging.info(f"Writing: {smear_p}")
    smear_p.write_text("\n".join(smear_text))

    # Outputs
    return max_smear_s, max_smear_l, max_smear_mag


def set_file_path(location: Path, file_path: Path):
    if file_path.parent == location:
        return file_path
    else:
        return location / file_path.name

# def uniform_interp(x, y, xi, extrapval):
#
#   I think this is just:
#       np.interp(xi, x, y, left=extrapval, right=extrapval)
#
#     # const ArrayXd & x, const ArrayXd & y, const ArrayXd & xi,
#     # double extrapval){
#
#     # Linear interpolation.
#     # We assume here that x is increasing and uniformly spaced.
#
#     # The value extrapval replaces the values outside of the interval
#     # spanned by X.
#
#     n  = len(x)
#     ni = len(xi)
#
#     if len(x) <= 1:
#         raise IndexError(
#             "Expecting at least two points in the input vector for "
#             "interpolation."
#         )
#
#     for x1, x2 in pairwise(x):
#         if x2 <= x1:
#             raise ValueError(
#             "Expecting an increasing vector in interpolation."
#             )
#
#     for (int s = 0; s < ni; s++){
#         if ( xi(s) < x(0) || xi(s) > x(n - 1) ){
#           yi(s) = extrapval;
#         }else{
#
#           int j  = (int)floor((xi(s) - x(0))/dx); if (j  <  0) j  = 0;
#           int jn = j + 1;                         if (jn >= n) jn = n - 1;
#
#           double slope = (y(jn) - y(j))/(x(jn) - x(j));
#           if ( x(jn) == x(j) ) slope = 0.0;
#           yi(s) = slope*(xi(s) - x(j)) + y(j);
#         }
#
#   return yi;


def write_data_for_plotting(path: os.PathLike, labels: list, *cols):
    # Given a collection of arrays of data, print those arrays as
    # columns in a text file. The columns need not have the same
    # number of elements.

    if len(labels) != len(cols):
        raise IndexError(
            "There is a different number of column labels than columns."
        )

    logging.info(f"Writing: {path}")
    with open(path, "w") as f:
        f.write("".join(map(lambda s: s.ljust(25), labels)) + "\n")

        # Print the data columns
        for zipped in itertools.zip_longest(*cols, fillvalue="nan"):
            f.write(
                "".join(map(lambda z: "{:>25.16}".format(z), zipped)) + "\n"
            )

        f.write("\n")
    return


def write_csv(path: os.PathLike, labels: list, *cols):
    logging.info(f"Writing: {path}")
    with open(path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(labels)
        for zipped in itertools.zip_longest(*cols, fillvalue="nan"):
            writer.writerow(zipped)


def write_gnuplot_file(
    gnuplot_path: Path, data_path: Path, img_path: Path,
    file_path1: Path, file_path2: Path, file_path3: Path
):
    logging.info(f"Writing: {gnuplot_path}")
    gnuplot_path.write_text(
        f"""\
dataFile  = '{data_path}'
imgFile   = '{img_path}'
filePath1 = '{file_path1}'
filePath2 = '{file_path2}'
filePath3 = '{file_path3}'

set terminal png size 1200, 900; set output imgFile
#set terminal pdfcairo;           set output 'fig.pdf'

set multiplot  # get into multiplot mode
set nokey      # no legend
set grid

set datafile missing 'nan'

w3 = 1.0/3.0; # will do 3 columns of plots
set size w3, w3

set title filePath1
set origin 0,    2*w3
plot dataFile using 4:5 with points pointtype 7 pointsize 0.6 lc rgb 'red', \
dataFile using 1:6 with lines lc rgb 'green', \
dataFile using 1:7 with lines lc rgb 'yellow'

set title filePath2
set origin w3,   2*w3
plot dataFile using 11:12 with points pointtype 7 pointsize 0.6 lc rgb 'red', \
dataFile using 1:13 with lines lc rgb 'green', \
dataFile using 1:14 with lines lc rgb 'yellow'

set title filePath3
set origin 2*w3, 2*w3
plot dataFile using 18:19 with points pointtype 7 pointsize 0.6 lc rgb 'red', \
dataFile using 1:20 with lines lc rgb 'green', \
dataFile using 1:21 with lines lc rgb 'yellow'

set title ''
set origin 0,     w3
plot dataFile using 4:8 with points pointtype 7 pointsize 0.6 lc rgb 'red', \
dataFile using 1:9 with lines lc rgb 'green', \
dataFile using 1:10 with lines lc rgb 'yellow'

set title ''
set origin w3,    w3
plot dataFile using 11:15 with points pointtype 7 pointsize 0.6 lc rgb 'red', \
dataFile using 1:16 with lines lc rgb 'green', \
dataFile using 1:17 with lines lc rgb 'yellow'

set title ''
set origin 2*w3,  w3
plot dataFile using 18:22 with points pointtype 7 pointsize 0.6 lc rgb 'red', \
dataFile using 1:23 with lines lc rgb 'green', \
dataFile using 1:24 with lines lc rgb 'yellow'

w2 = 0.5 # 1/2 of the plotting window
set size w2, w3

set title 'Cross-track Jitter'
set origin 0,     0
plot dataFile using 1:2 with lines lc rgb 'blue'

set title 'Down-track Jitter'
set origin w2,    0
plot dataFile using 1:3 with lines lc rgb 'blue'

unset multiplot # exit multiplot mode
"""
    )
    return
