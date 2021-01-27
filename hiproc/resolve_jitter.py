#!/usr/bin/env python
"""
Perform the jitter derivation for a HiRISE observation.

Outputs are the jitter in the x (sample) and y (line) directions in
pixels, and the time (ephemeris time at that translation), average error
between the derived jitter function and the original data, linerate
(line time read from flat files), and TDI. These outputs are
written to a text file. Another output is the pixel smear, also
written to a file.

The C++ version of this runs ~5x faster, FYI.
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
# Copyright 2020-2021, Ross A. Beyer (rbeyer@seti.org)
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
#
# The only write-up that I could find for the pre-cursor was:
#   S. Mattson, A. Boyd, R. L. Kirk, D. A. Cook, and E. Howington-Kraus,
#   HiJACK: Correcting spacecraft jitter in HiRISE images of Mars,
#   European Planetary Science Congress 2009, #604
#   https://ui.adsabs.harvard.edu/abs/2009epsc.conf..604M}
# However, that's really only a high-level description.

import argparse
import csv
import itertools
import logging
import math
import os
import sys
from collections import abc
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import medfilt
from scipy.interpolate import PchipInterpolator

import pvl

import hiproc.hirise as hirise
import hiproc.util as util
from hiproc.FlatFile import FlatFile

logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     parents=[util.parent_parser()])
    parser.add_argument(
        "-c",
        "--conf",
        required=False,
        default=Path(__file__).resolve().parent.parent /
                "data" / "ResolveJitter.conf",
        help="Path to a ResolveJitter.conf file.",
    )
    parser.add_argument(
        "-p",
        "--plot",
        required=False,
        action="store_true",
        help="Displays interactive plot.",
    )
    parser.add_argument(
        "--saveplot",
        required=False,
        nargs="?",
        default=False,
        const=True,
        help="Saves plot to a default filename in the output directory."
             "If a filename is provided it will be used to save the plot.",
    )
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

    util.set_logger(logger, args.log, args.logfile)

    if args.line_interval <= 0:
        raise ValueError("The parameter 'line_interval' must be positive.")

    start(
        args.file_path1, True if args.which1 != 1 else False,
        args.file_path2, True if args.which2 != 1 else False,
        args.file_path3, True if args.which3 != 1 else False,
        imgdir=args.image_location,
        obsid=args.image_id, lineint=args.line_interval, confpath=args.conf,
        plotshow=args.plot, plotsave=args.saveplot
    )
    return


def start(
    file_path1: Path, whichfrom1: bool,
    file_path2: Path, whichfrom2: bool,
    file_path3: Path, whichfrom3: bool,
    window_size=11, window_width=2,
    imgdir=None, obsid=None, lineint=None, confpath=None,
    plotshow=False, plotsave=False
):
    if imgdir is None:
        fp1 = file_path1
        fp2 = file_path2
        fp3 = file_path3
        image_location = fp1.parent
    else:
        fp1 = set_file_path(imgdir, file_path1)
        fp2 = set_file_path(imgdir, file_path2)
        fp3 = set_file_path(imgdir, file_path3)
        image_location = imgdir

    if lineint is None:
        if confpath is None:
            raise ValueError(
                f"lineint was None and so was confpath, so can't look up value."
            )
        else:
            line_interval = pvl.load(
                confpath
            )["AutoRegistration"]["ControlNet"]["Control_Lines"]
    else:
        line_interval = lineint

    oid1 = hirise.get_ObsID_fromfile(fp1)
    oid2 = hirise.get_ObsID_fromfile(fp2)
    oid3 = hirise.get_ObsID_fromfile(fp3)

    if oid1 == oid2 == oid3:
        oid = oid1
    else:
        raise ValueError(
            f"The observation IDs from the three file"
            f"paths ({file_path1}, {file_path2}, {file_path3}) do not match."
        )

    if obsid is None:
        image_id = str(oid)
    else:
        if obsid == str(oid):
            image_id = obsid
        else:
            raise ValueError(
                f"The provided obsid ({obsid}) does not match"
                f"the Observation IDs dervied from the images ({oid})."
            )

    if plotsave:
        try:
            plotsave = Path(plotsave)
        except TypeError:
            plotsave = image_location / (image_id + "_jitter_plot_py.pdf")

    t1, offx1, offy1, lines1, dt1, tdi1, linerate1 = parse_file(
        fp1, window_size, window_width, whichfrom1
    )
    t2, offx2, offy2, lines2, dt2, tdi2, linerate2 = parse_file(
        fp2, window_size, window_width, whichfrom2
    )
    t3, offx3, offy3, lines3, dt3, tdi3, linerate3 = parse_file(
        fp3, window_size, window_width, whichfrom3
    )

    # nfft is the number of "steps" that we will be using for our
    # fourier transforms.  These steps represent a certain number of
    # image lines, but also represent the number of "time steps" that
    # we will use.
    if lines1 == lines2 == lines3:
        nfft = upper_power_of_two(lines1 / line_interval)
    else:
        raise ValueError(
            "The number of lines in the three images is not identical."
        )

    if tdi1 == tdi2 == tdi3:
        tdi = tdi3
    else:
        raise ValueError("The values of tdi are not identical.")

    if linerate1 == linerate2 == linerate3:
        linerate = linerate3
    else:
        raise ValueError("The values of linerate are not identical.")

    offx1_filtered = filter_data(nfft, 2 / nfft, offx1)
    offy1_filtered = filter_data(nfft, 2 / nfft, offy1)
    offx2_filtered = filter_data(nfft, 2 / nfft, offx2)
    offy2_filtered = filter_data(nfft, 2 / nfft, offy2)
    offx3_filtered = filter_data(nfft, 2 / nfft, offx3)
    offy3_filtered = filter_data(nfft, 2 / nfft, offy3)

    # The values in tt are the fractional points in [0:1]
    # that correspond to the nfft number.
    tt = np.linspace(0, 1, nfft, endpoint=False)

    # The first file to be parsed sets nfftime, t0, and duration
    nfftime = np.linspace(t1[0], t1[-1], nfft)
    # et_shift = et - t1[0]
    t0 = t1[0]
    duration = t1[-1] - t0

    xinterp1, yinterp1, x1, y1, ddt1, overxx1, overyy1 = create_matrices(
        t1, offx1_filtered, offy1_filtered,
        dt1, duration, t0,
        nfft, nfftime, tt
    )
    xinterp2, yinterp2, x2, y2, ddt2, overxx2, overyy2 = create_matrices(
        t2, offx2_filtered, offy2_filtered,
        dt2, duration, t0,
        nfft, nfftime, tt
    )
    xinterp3, yinterp3, x3, y3, ddt3, overxx3, overyy3 = create_matrices(
        t3, offx3_filtered, offy3_filtered,
        dt3, duration, t0,
        nfft, nfftime, tt
    )

    logger.info("Searching for correct phasetol")
    # For the test data, the following while loop will *always* run the full
    # number of repetitions.  If it were guaranteed that there was only
    # one minima in the changing value of *error*, then we could exit
    # this loop early, once the error starts going back up, or adopt
    # strategies to divide and conquer the phasetol parameter space.
    # A better understanding of the error behavior could allow us to speed
    # these loops up.

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

    # No need to calculate these during every loop:
    rh0 = dict(
        x1=np.real(x1[0]) / 2.0,
        x2=np.real(x2[0]) / 2.0,
        x3=np.real(x3[0]) / 2.0,
        y1=np.real(y1[0]) / 2.0,
        y2=np.real(y2[0]) / 2.0,
        y3=np.real(y3[0]) / 2.0,
    )
    while error >= error_tol and k < repetitions:
        k += 1

        # setting the phase tolerance
        phasetol = k * tolcoef

        # null the frequencies that cause a problem (really zero them out)
        overxxx1, overyyy1 = mask_frequencies(
            phasetol, ddt1, overxx1, overyy1
        )
        overxxx2, overyyy2 = mask_frequencies(
            phasetol, ddt2, overxx2, overyy2
        )
        overxxx3, overyyy3 = mask_frequencies(
            phasetol, ddt3, overxx3, overyy3
        )

        # Adding all frequencies together
        stackedx = np.ma.stack((overxxx1, overxxx2, overxxx3))
        stackedy = np.ma.stack((overyyy1, overyyy2, overyyy3))

        overxxx = np.ma.mean(stackedx, axis=0)
        overyyy = np.ma.mean(stackedy, axis=0)

        # take the sum of each row
        overx = np.ma.sum(overxxx, axis=0)
        overy = np.ma.sum(overyyy, axis=0)

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

        error_vec = (
            np.abs(xinterp1 - (jittercheckx1 + rh0["x1"])) +
            np.abs(xinterp2 - (jittercheckx2 + rh0["x2"])) +
            np.abs(xinterp3 - (jittercheckx3 + rh0["x3"])) +
            np.abs(yinterp1 - (jitterchecky1 + rh0["y1"])) +
            np.abs(yinterp2 - (jitterchecky2 + rh0["y2"])) +
            np.abs(yinterp3 - (jitterchecky3 + rh0["y3"]))
        ) / 6.0

        error = error_vec.mean()
        logger.info(f"Error for phasetol {phasetol}: {error}")

        if error < min_avg_error:
            min_avg_error = error
            # min_k = k
            min_jitterx = jitterx
            min_jittery = jittery

    # end while
    logger.info(f"Minimum Error after phase filtering: {min_avg_error}")

    logger.info("Searching for correct filter size.")
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
            np.abs(xinterp1 - (jittercheckx1 + rh0["x1"])) +
            np.abs(xinterp2 - (jittercheckx2 + rh0["x2"])) +
            np.abs(xinterp3 - (jittercheckx3 + rh0["x3"])) +
            np.abs(yinterp1 - (jitterchecky1 + rh0["y1"])) +
            np.abs(yinterp2 - (jitterchecky2 + rh0["y2"])) +
            np.abs(yinterp3 - (jitterchecky3 + rh0["y3"]))
        )

        error = error_vec.mean()
        logger.info(f"Erorr for omega {omega}: {error}")

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
        nfftime, sample, line, linerate, tdi,
        path=(image_location / (image_id + "_smear_py.txt"))
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
    for s, l, e in zip(sample, line, nfftime):
        jitter_text.append(f"""{s}     {l}     {e}""")

    logger.info(f"Writing: {jitter_p}")
    jitter_p.write_text("\n".join(jitter_text))

    # I think we could re-do this for matplotlib.  Let's defer it for now.
    # Writing the data we will plot later in gnuplot
    data_p = image_location / (image_id + "_jitter_plot_py.txt")
    t1_shift = t1 - t0
    t2_shift = t2 - t0
    t3_shift = t3 - t0
    jittercheckx1_shift = min_jitter_check_x1 + rh0["x1"]
    jitterchecky1_shift = min_jitter_check_y1 + rh0["y1"]
    jittercheckx2_shift = min_jitter_check_x2 + rh0["x2"]
    jitterchecky2_shift = min_jitter_check_y2 + rh0["y2"]
    jittercheckx3_shift = min_jitter_check_x3 + rh0["x3"]
    jitterchecky3_shift = min_jitter_check_y3 + rh0["y3"]
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
        data_p, string_labels, nfftime - t0, sample, line, t1_shift,
        offx1_filtered, xinterp1, jittercheckx1_shift, offy1_filtered, yinterp1,
        jitterchecky1_shift, t2_shift, offx2_filtered, xinterp2,
        jittercheckx2_shift, offy2_filtered, yinterp2, jitterchecky2_shift,
        t3_shift, offx3_filtered, xinterp3,
        jittercheckx3_shift, offy3_filtered, yinterp3, jitterchecky3_shift
    )

    write_csv(
        image_location / (image_id + "_jitter_plot_py.csv"),
        string_labels,
        nfftime - t0, sample, line, t1_shift,
        offx1_filtered, xinterp1, jittercheckx1_shift, offy1_filtered, yinterp1,
        jitterchecky1_shift, t2_shift, offx2_filtered, xinterp2,
        jittercheckx2_shift, offy2_filtered, yinterp2, jitterchecky2_shift,
        t3_shift, offx3_filtered, xinterp3, jittercheckx3_shift,
        offy3_filtered, yinterp3, jitterchecky3_shift
    )

    gnuplot_p = image_location / (image_id + "_jitter_plot_py.plt")
    img_file_name = image_location / (image_id + "_jitter_plot_py.png")
    write_gnuplot_file(
        gnuplot_p, data_p, img_file_name, file_path1, file_path2, file_path3
    )

    if plotshow or plotsave:
        plot(
            t1_shift, offx1_filtered, offy1_filtered,
            xinterp1, yinterp1, jittercheckx1_shift, jitterchecky1_shift,
            file_path1.stem,
            t2_shift, offx2_filtered, offy2_filtered,
            xinterp2, yinterp2, jittercheckx2_shift, jitterchecky2_shift,
            file_path2.stem,
            t3_shift, offx3_filtered, offy3_filtered,
            xinterp3, yinterp3, jittercheckx3_shift, jitterchecky3_shift,
            file_path3.stem,
            nfftime - t0, sample, line, show=plotshow, save=plotsave
        )

    return


def create_matrices(
    time: np.array, offx: np.array, offy: np.array,
    dt: float, duration: float, t0: float,
    nfft: int, nfftime: np.array, tt: abc.Sequence
):
    """Returns a tuple of numpy arrays.

    :param time: Time values in seconds (numpy array).
    :param offx: Sample offsets (numpy array).
    :param offy: Line offsets (numpy array).
    :param dt: Time difference between the FROM and MATCH times (float).
    :param duration: The total time duration (float).
    :param t0: Zero time to start the Fourier Transform steps (float).
    :param nfft: Number of divisions to use for the Fourier Transform (int).
    :param nfftime: *nfft* time values starting at *t0* (numpy array).
    :param tt: Fractional values in [0:1] that correspond to the nfft number
        (numpy array).
    :return: There are seven elements in the tuple:
        0: Interpolation of *offx* at the values of *nfftime*
        1: Interpolation of *offy* at the values of *nfftime*
        2: Fourier transform of xinterp * (2 / *nfft*)
        3: Fourier transform of yinterp * (2 / *nfft*)
        4: phase difference
        5: overxx: ?
        6: overyy: ?
    """
    t_shift = time - t0

    fx = PchipInterpolator(t_shift, offx, extrapolate=False)
    fy = PchipInterpolator(t_shift, offy, extrapolate=False)
    xinterp = fx(nfftime - t0)
    yinterp = fy(nfftime - t0)

    np.nan_to_num(xinterp, copy=False, nan=np.mean(offx))
    np.nan_to_num(yinterp, copy=False, nan=np.mean(offy))

    # getting the frequencies of the Fourier transform
    freq = np.linspace(0, nfft / 2, int(nfft / 2), endpoint=False)

    # taking the fourier transform of the offsets
    x = 2 * np.fft.fft(xinterp) / nfft
    y = 2 * np.fft.fft(yinterp) / nfft

    # separating sines and cosines
    xa = x[:int(nfft / 2)].real
    xb = -1 * x[:int(nfft / 2)].imag
    ya = y[:int(nfft / 2)].real
    yb = -1 * y[:int(nfft / 2)].imag

    # calculates the phase difference
    twopi = math.pi * 2
    ddt_temp = (dt / duration * twopi) * freq
    ddt = ddt_temp - twopi * np.floor(ddt_temp / twopi)

    # the coeficients for the frequencies
    with np.errstate(divide="ignore"):
        aaax = -0.5 * (-1 * xa * np.cos(ddt) + np.sin(ddt) * xb - xa) / np.sin(ddt)
        aaay = -0.5 * (-1 * ya * np.cos(ddt) + np.sin(ddt) * yb - ya) / np.sin(ddt)

    with np.errstate(invalid="ignore"):
        bbbx = -0.5 * (xb * np.cos(ddt) + np.sin(ddt) * xa + xb) / np.sin(ddt)
        bbby = -0.5 * (yb * np.cos(ddt) + np.sin(ddt) * ya + yb) / np.sin(ddt)

    # create series of sines and cosines
    ft = freq.reshape(-1, 1) * tt.reshape(1, -1) * twopi
    sn = np.sin(ft)
    cn = np.cos(ft)
    aaax_rep = np.repeat(aaax.reshape(-1, 1), tt.size, axis=1)
    bbbx_rep = np.repeat(bbbx.reshape(-1, 1), tt.size, axis=1)
    aaay_rep = np.repeat(aaay.reshape(-1, 1), tt.size, axis=1)
    bbby_rep = np.repeat(bbby.reshape(-1, 1), tt.size, axis=1)
    with np.errstate(invalid="ignore"):
        overxx = aaax_rep * sn + bbbx_rep * cn
        overyy = aaay_rep * sn + bbby_rep * cn

    # Outputs
    # ArrayXd & tt, ArrayXd & ET, ArrayXd & ET_shift,
    # ArrayXd & ddt,
    # ArrayXd & xinterp, ArrayXd & yinterp,
    # ArrayXcd & X, ArrayXcd & Y, MatrixXd & overxx, MatrixXd & overyy
    #
    return xinterp, yinterp, x, y, ddt, overxx, overyy


def upper_power_of_two(value) -> int:
    """Returns the value of 2 raised to some power which is the smallest
     such value that is just >= *value*."""
    result = 1
    while result < value:
        result <<= 1
    return result


def filter_data(nfft: int, c: float, data: np.array) -> np.array:
    """Apply a Gaussian filter to the data in the frequency domain.

    :param nfft: The number of steps in the Fourier Transform (int).
    :param c: ?
    :param data: An array of data to be filtered (numpy array).
    """

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


def mask_frequencies(
    phasetol: float, ddt: np.array, x: np.array, y: np.array
):
    """Returns *x* and *y* as numpy masked arrays with 'problematic frequencies'
    masked.

    :param phasetol: Phase values from zero to *phasetol* and 2*pi - phasetol
        will be masked.
    :param ddt: phase difference
    :param x: overxx from create_matrices()
    :param y: overyy from create_matrices()

    It is assumed that *ddt* has the same size as axis 0 of *x* and *y*.
    """
    if x.shape != y.shape:
        raise ValueError(
            f"The shape of x {x.shape} and y {y.shape} must be the same."
        )

    if x.shape[0] != ddt.size:
        raise ValueError(
            f"The size of ddt ({ddt.size}) must be the same as the first axis"
            f"of the x and y arrays {x.shape}"
        )

    # mask the frequencies that cause a problem
    a = np.less(np.abs(ddt), phasetol)
    b = np.greater(np.abs(ddt), (2 * math.pi) - phasetol)
    null_positions = np.logical_or(a, b)

    # We must reshape to a 2D column, and then tile that across, so that each
    # row of the 2D matrix has the same value for all positions.
    null_2d = np.tile(null_positions.reshape(-1, 1), (1, x.shape[1]))

    x_masked = np.ma.array(x, mask=null_2d, fill_value=0)
    y_masked = np.ma.array(y, mask=null_2d, fill_value=0)

    return x_masked, y_masked


def parse_file(
    file_path: os.PathLike, window_size: int, window_width: int, whichfrom=True
):
    """Returns a tuple of information from the Flat file at *file_path*.

    There are seven elements in the tuple:
    0: unique time values (numpy array)
    1: sample offsets for each time (numpy array)
    2: line offsets for each time (numpy array)
    3: number of lines listed for the FROM file (int)
    4: seconds between the times in the FROM and MATCH files (float)
    5: The TDI for the FROM and MATCH files (int)
    6: The LineRate from the FROM and MATCH files (float)

    *window_size* is the kernel size for median filtering the offsets, should
    be an odd integer.

    *window_width* determines the boundaries above and below the filtered
    average beyond which to exclude outliers.

    The optional *whichfrom* parameter determines the sense of the offsets.
    The default value of True makes the offsets relative to the From cube,
    and False makes the offsets relative to the Match cube.
    """
    logger.info(f"Reading: {file_path}")

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

    if whichfrom == 1:
        column = "FromTime"
        which = -1
    else:
        column = "MatchTime"
        which = 1

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

    time_arr = np.array(time)
    offx_arr = np.array(offset_x)
    offy_arr = np.array(offset_y)

    magnitude = np.sqrt(offx_arr**2 + offy_arr**2)
    avemag = medfilt(magnitude, window_size)

    # Throw out the out-of-range values:
    high_window = avemag + window_width
    low_window = avemag - window_width
    good_idxs = np.nonzero(
        np.logical_and(low_window < magnitude, magnitude < high_window)
    )

    # Some flat files have more than one measurement for the same timestamp.
    # For those that do, the multiple x and y offsets for the same
    # timestamp are averaged together.
    #
    # Also, the original MatLab code had no functions, but repeated code,
    # and there was a transcription error on the third pass of this averaging
    # such that the very first row was never included in the averaging.
    # The C++ code allowed for that broken behavior, but we won't here.
    t_arr, unique_idxs = np.unique(time_arr[good_idxs], return_index=True)

    if unique_idxs.size == time_arr[good_idxs].size:
        offx = offx_arr[good_idxs]
        offy = offy_arr[good_idxs]
    else:
        x_means = list()
        for a in np.split(offx_arr[good_idxs], unique_idxs[1:]):
            x_means.append(np.mean(a))
        offx = np.array(x_means)

        y_means = list()
        for a in np.split(offy_arr[good_idxs], unique_idxs[1:]):
            y_means.append(np.mean(a))
        offy = np.array(y_means)

    return t_arr, offx, offy, int(flat["FROM"]["Lines"]), dt, tdi, line_rate


def pixel_smear(
    t: np.array,
    sample: np.array,
    line: np.array,
    linerate: float,
    tdi: int,
    path=None,
    image_id=None
):
    """Returns the smear values from the derived jitter function.

    Pixel smear due to jitter is calculated by interpolating the jitter
    function at intervals equivalent to the linerate.  Then the
    difference is taken over that interval and multiplied by the TDI.
    This provides an estimated minimum for pixel smear. If the motion
    that caused the smear is not captured in the jitter derivation,
    then it cannot be plotted here.

    :param t: Time values (numpy array).
    :param sample: Pixel offsets in the sample direction (numpy array).
    :param line: Pixel offsest in the line direction (numpy array).
    :param linerate: The image line rate (float).
    :param tdi: The image TDI value (int).
    :param path: Optional path to write out smear details to (Path).
    :param image_id: Optional Observation ID (string).
    :return:
    """
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
    if path is not None:
        id_str = f" for {image_id}"
        smear_text = [
            f"""\
# Smear values are calculated from the derived jitter function{id_str}.
# Maximum Cross-track pixel smear {max_smear_s}
# Maximum Down-track pixel smear {max_smear_l}
# Maximum Pixel Smear Magnitude {max_smear_mag}
# Sample                 Line                   EphemerisTime"""
        ]

        for ess, ell, exi in zip(dysdx, dyldx, xi):
            smear_text.append(f"{ess}     {ell}     {exi}")

        logger.info(f"Writing: {path}")
        path.write_text("\n".join(smear_text))

    # Outputs
    return max_smear_s, max_smear_l, max_smear_mag


def set_file_path(location: Path, file_path: Path):
    if file_path.parent == location:
        return file_path
    else:
        return location / file_path.name


def plot(
    t1, x1, y1, xinterp1, yinterp1, jittercheckx1, jitterchecky1, title1,
    t2, x2, y2, xinterp2, yinterp2, jittercheckx2, jitterchecky2, title2,
    t3, x3, y3, xinterp3, yinterp3, jittercheckx3, jitterchecky3, title3,
    et, sample, line, show=True, save=False
):
    plt.ioff()
    fig = plt.figure(constrained_layout=True)
    gs = fig.add_gridspec(3, 6)

    fig.suptitle("Resolve Jitter Results")

    ax00 = fig.add_subplot(gs[0, 0:2])
    ax00.set_title(title1)
    ax00.set_ylabel("Sample Offset")

    ax01 = fig.add_subplot(gs[0, 2:4])
    ax01.set_title(title2)

    ax02 = fig.add_subplot(gs[0, 4:])
    ax02.set_title(title3)

    ax10 = fig.add_subplot(gs[1, 0:2])
    ax10.set_ylabel("Line Offset")
    ax10.set_xlabel("Seconds")
    ax11 = fig.add_subplot(gs[1, 2:4])
    ax11.set_xlabel("Seconds")
    ax12 = fig.add_subplot(gs[1, 4:])
    ax12.set_xlabel("Seconds")

    ax20 = fig.add_subplot(gs[2, 0:3])
    ax20.set_title("Cross-Track Jitter")
    ax20.set_ylabel("Sample Offset")
    ax20.set_xlabel("Seconds")

    ax21 = fig.add_subplot(gs[2, 3:])
    ax21.set_title("Down-Track Jitter")
    ax21.set_ylabel("Line Offset")
    ax21.set_xlabel("Seconds")

    ax00.plot(t1, x1, "o", c="red")
    ax00.plot(et, xinterp1, c="green")
    ax00.plot(et, jittercheckx1, c="yellow")

    ax01.plot(t2, x2, "o", c="red")
    ax01.plot(et, xinterp2, c="green")
    ax01.plot(et, jittercheckx2, c="yellow")

    ax02.plot(t3, x3, "o", c="red")
    ax02.plot(et, xinterp3, c="green")
    ax02.plot(et, jittercheckx3, c="yellow")

    ax10.plot(t1, y1, "o", c="red")
    ax10.plot(et, yinterp1, c="green")
    ax10.plot(et, jitterchecky1, c="yellow")

    ax11.plot(t2, y2, "o", c="red")
    ax11.plot(et, yinterp2, c="green")
    ax11.plot(et, jitterchecky2, c="yellow")

    ax12.plot(t3, y3, "o", c="red")
    ax12.plot(et, yinterp3, c="green")
    ax12.plot(et, jitterchecky3, c="yellow")

    ax20.plot(et, sample, c="blue")
    ax21.plot(et, line, c="blue")

    if save:
        logger.info(f"Writing: {save}")
        plt.savefig(save)

    if show:
        plt.show()


def write_csv(path: os.PathLike, labels: list, *cols, fillvalue="nan"):
    """Identical to write_data_for_plotting(), but writes a CSV file
    instead of a fixed-width text file.
    """
    logger.info(f"Writing: {path}")
    with open(path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(labels)
        for zipped in itertools.zip_longest(*cols, fillvalue=fillvalue):
            writer.writerow(zipped)


def write_data_for_plotting(
    path: os.PathLike, labels: list, *cols, fillvalue="nan"
):
    """Given a collection of arrays of data (*cols*), write those arrays as
    fixed-width columns (25 characters wide) in a text file at *path*, prefixed
    by the *labels* in the first row. The columns need not have the same number
    of elements as each other, and any short columns will be filled with
    *fillvalue*.

    This is mostly historic to provide a text file that gnuplot can use.
    """

    if len(labels) != len(cols):
        raise IndexError(
            "There is a different number of column labels than columns."
        )

    logger.info(f"Writing: {path}")
    with open(path, "w") as f:
        f.write("".join(map(lambda s: s.ljust(25), labels)) + "\n")

        # Print the data columns
        for zipped in itertools.zip_longest(*cols, fillvalue=fillvalue):
            f.write(
                "".join(map(lambda z: "{:>25.16}".format(z), zipped)) + "\n"
            )

        f.write("\n")


def write_gnuplot_file(
    gnuplot_path: Path, data_path: Path, img_path: Path,
    file_path1: Path, file_path2: Path, file_path3: Path
):
    """Writes a gnuplot file that will plot the contents of a file
    written by write_data_for_plotting().

    The file will be written to *gnuplot_path*.  *data_path* should be the
    file written by write_data_for_plotting().  The *img_path* is the png
    file that will be created when the file at *gnuplot_path* is run by
    gnuplot.  *file_path1*, *file_path2*, and *file_path3* are just used
    to provide titles to the plots.
    """
    logger.info(f"Writing: {gnuplot_path}")
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
