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

WARNING: This program and its modules may not be fully functional and
may return poor results.
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
import pkg_resources
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import medfilt
from scipy.interpolate import PchipInterpolator
from scipy.optimize import minimize_scalar

import pvl

import hiproc.hirise as hirise
import hiproc.util as util
from hiproc.FlatFile import FlatFile

logger = logging.getLogger(__name__)


def arg_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, parents=[util.parent_parser()]
    )
    parser.add_argument(
        "-c",
        "--conf",
        type=argparse.FileType('r'),
        default=pkg_resources.resource_stream(
            __name__,
            'data/ResolveJitter.conf'
        ),
        help="Path to a ResolveJitter.conf file, only needed if "
        "--lineinterval isn't given.",
    )
    parser.add_argument(
        "--csv",
        action="store_true",
        help="This program writes out a fixed-width data file with extra "
        "information for plotting.  This also writes out a "
        "comma-separated version.",
    )
    parser.add_argument(
        "--lineinterval",
        type=float,
        help="The number of lines to use to set the number of Fourier "
        "transform intervals, defaults to Control_Lines in the "
        "ResolveJitter.conf file.",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        help="Output directory.  Defaults to the directory of the first "
        "input file.",
    )
    parser.add_argument(
        "--outprefix",
        help="Prefix string for output files.  If not given, will default "
        "to the Observation ID of the images.",
    )
    parser.add_argument(
        "-p", "--plot", action="store_true", help="Displays interactive plot.",
    )
    parser.add_argument(
        "--saveplot",
        nargs="?",
        default=False,
        const=True,
        help="Saves plot to a default filename in the output directory. "
        "If a filename is provided it will be used to save the plot.",
    )
    parser.add_argument(
        "--whichmatch1",
        action="store_false",
        dest="which1",
        help="If specified, the sense of the offsets for the first "
        "file will be relative to the MATCH cube, rather than the FROM "
        "cube.",
    )
    parser.add_argument(
        "--whichmatch2",
        action="store_false",
        dest="which2",
        help="If specified, the sense of the offsets for the second "
        "file will be relative to the MATCH cube, rather than the FROM "
        "cube.",
    )
    parser.add_argument(
        "--whichmatch3",
        action="store_false",
        dest="which3",
        help="If specified, the sense of the offsets for the third "
        "file will be relative to the MATCH cube, rather than the FROM "
        "cube.",
    )
    # parser.add_argument(
    #     "--optimize",
    #     action="store_true",
    #     help="Will run experimental optimizer."
    # )
    parser.add_argument(
        "files",
        nargs="*",
        help="Three flat.txt files that are the output of ISIS hijitreg.",
    )
    return parser


def main():
    parser = arg_parser()
    args = parser.parse_args()

    util.set_logger(args.verbose, args.logfile, args.log)

    if len(args.files) == 3:
        # With just three arguments, these are the expected flat files.
        fp1, fp2, fp3 = map(Path, args.files)
        which1 = args.which1
        which2 = args.which2
        which3 = args.which3
    elif len(args.files) == 9:
        # This is the old-style positional calling
        oldparser = argparse.ArgumentParser()
        oldparser.add_argument("outdir", type=Path)
        oldparser.add_argument("outprefix", type=str)
        oldparser.add_argument("lineinterval", type=float)
        oldparser.add_argument("file_path1", type=Path)
        oldparser.add_argument("which1", type=int, choices=[-1, 1])
        oldparser.add_argument("file_path2", type=Path)
        oldparser.add_argument("which2", type=int, choices=[-1, 1])
        oldparser.add_argument("file_path3", type=Path)
        oldparser.add_argument("which3", type=int, choices=[-1, 1])
        args = oldparser.parse_args(args.files, namespace=args)
        fp1 = args.file_path1
        fp2 = args.file_path2
        fp3 = args.file_path3
        which1 = True if args.which1 != 1 else False
        which2 = True if args.which2 != 1 else False
        which3 = True if args.which2 != 1 else False
    else:
        parser.error("Only takes 3 or 9 positional arguments.")

    if args.lineinterval is None:
        args.lineinterval = pvl.load(args.conf)["AutoRegistration"][
            "ControlNet"
        ]["Control_Lines"]
    elif args.lineinterval <= 0:
        raise ValueError("--lineinterval must be positive.")

    if args.outdir is None:
        outdir = fp1.parent
    else:
        outdir = args.outdir
        fp1 = set_file_path(outdir, fp1)
        fp2 = set_file_path(outdir, fp2)
        fp3 = set_file_path(outdir, fp3)

    with util.main_exceptions(args.verbose):
        start(
            fp1,
            which1,
            fp2,
            which2,
            fp3,
            which3,
            line_interval=args.lineinterval,
            outdir=outdir,
            outprefix=args.outprefix,
            plotshow=args.plot,
            plotsave=args.saveplot,
            writecsv=args.csv,
            opt=args.optimize
        )
    return


def start(
    file_path1: Path,
    whichfrom1: bool,
    file_path2: Path,
    whichfrom2: bool,
    file_path3: Path,
    whichfrom3: bool,
    line_interval: float,
    outdir: Path,
    outprefix=None,
    plotshow=False,
    plotsave=False,
    writecsv=False,
    opt=False
):
    oid1 = hirise.get_ObsID_fromfile(file_path1)
    oid2 = hirise.get_ObsID_fromfile(file_path2)
    oid3 = hirise.get_ObsID_fromfile(file_path3)

    if oid1 == oid2 == oid3:
        oid = oid1
    else:
        raise ValueError(
            f"The observation IDs from the three file"
            f"paths ({file_path1}, {file_path2}, {file_path3}) do not match."
        )

    (
        nfftime,
        sample,
        line,
        linerate,
        tdi,
        t0,
        t1,
        t2,
        t3,
        offx_filtered,
        offy_filtered,
        xinterp,
        yinterp,
        min_avg_error,
        min_k,
        jittercheckx,
        jitterchecky,
        rh0,
    ) = resolve_jitter(
        file_path1,
        whichfrom1,
        file_path2,
        whichfrom2,
        file_path3,
        whichfrom3,
        line_interval,
        optimize=opt
    )

    # The outputs
    logger.info(f"Average error is {min_avg_error} at min index {min_k}")
    logger.info(f"linerate is {linerate}")
    logger.info(f"TDI = {tdi}")

    # To do: Remove the py suffix from output filenames.

    if outprefix is None:
        outprefix = str(oid)
    else:
        outprefix = str(outprefix)

    # Characterize the smear:
    (
        max_smear_sample,
        max_smear_line,
        max_smear_mag,
        dysdx,
        dyldx,
        xi,
    ) = pixel_smear(nfftime, sample, line, linerate, tdi)

    write_smear_data(
        (outdir / (outprefix + "_smear_py.txt")),
        max_smear_sample,
        max_smear_line,
        max_smear_mag,
        dysdx,
        dyldx,
        xi,
        oid,
    )

    # Make a text file of the jitter data
    jitter_p = outdir / (outprefix + "_jitter_py.txt")
    jitter_text = [
        f"""# Using image {oid} the jitter was found with an
# Average Error of {min_avg_error}
# Maximum Cross-track pixel smear {max_smear_sample}
# Maximum Down-track pixel smear {max_smear_line}
# Maximum Pixel Smear Magnitude {max_smear_mag}
#
# Sample                 Line                   ET"""
    ]
    for s, l, e in zip(sample, line, nfftime):
        jitter_text.append(f"""{s}     {l}     {e}""")

    logger.info(f"Writing: {jitter_p}")
    jitter_p.write_text("\n".join(jitter_text))

    # Create data for plotting
    data_p = outdir / (outprefix + "_jitter_plot_py.txt")
    t_shift = [t1 - t0, t2 - t0, t3 - t0]
    jittercheckx_shift = [
        jittercheckx[0] + rh0["x1"],
        jittercheckx[1] + rh0["x2"],
        jittercheckx[2] + rh0["x3"],
    ]
    jitterchecky_shift = [
        jitterchecky[0] + rh0["y1"],
        jitterchecky[1] + rh0["y2"],
        jitterchecky[2] + rh0["y3"],
    ]

    # Note the comment before the first label
    # This ordering and naming is historic to the original file output.
    string_labels = [
        "# ET_shift",
        "Sample",
        "Line",
        "t1_shift",
        "offx1",
        "xinterp1",
        "jittercheckx1_shift",
        "offy1",
        "yinterp1",
        "jitterchecky1_shift",
        "t2_shift",
        "offx2",
        "xinterp2",
        "jittercheckx2_shift",
        "offy2",
        "yinterp2",
        "jitterchecky2_shift",
        "t3_shift",
        "offx3",
        "xinterp3",
        "jittercheckx3_shift",
        "offy3",
        "yinterp3",
        "jitterchecky3_shift",
    ]
    data_to_plot = (
        string_labels,
        nfftime - t0,
        sample,
        line,
        t_shift[0],
        offx_filtered[0],
        xinterp[0],
        jittercheckx_shift[0],
        offy_filtered[0],
        yinterp[0],
        jitterchecky_shift[0],
        t_shift[1],
        offx_filtered[1],
        xinterp[1],
        jittercheckx_shift[1],
        offy_filtered[1],
        yinterp[1],
        jitterchecky_shift[1],
        t_shift[2],
        offx_filtered[2],
        xinterp[2],
        jittercheckx_shift[2],
        offy_filtered[2],
        yinterp[2],
        jitterchecky_shift[2],
    )

    write_data_for_plotting(data_p, *data_to_plot)

    if writecsv:
        write_csv(outdir / (outprefix + "_jitter_plot_py.csv"), *data_to_plot)

    gnuplot_p = outdir / (outprefix + "_jitter_plot_py.plt")
    img_file_name = outdir / (outprefix + "_jitter_plot_py.png")
    write_gnuplot_file(
        gnuplot_p, data_p, img_file_name, file_path1, file_path2, file_path3
    )

    if plotsave:
        try:
            plotsave = Path(plotsave)
        except TypeError:
            plotsave = outdir / (outprefix + "_jitter_plot_py.pdf")

    if plotshow or plotsave:
        plot(
            t_shift,
            offx_filtered,
            offy_filtered,
            xinterp,
            yinterp,
            jittercheckx_shift,
            jitterchecky_shift,
            [file_path1.stem, file_path2.stem, file_path3.stem],
            nfftime - t0,
            sample,
            line,
            show=plotshow,
            save=plotsave,
        )

    return


def resolve_jitter(
    file_path1: Path,
    whichfrom1: bool,
    file_path2: Path,
    whichfrom2: bool,
    file_path3: Path,
    whichfrom3: bool,
    line_interval: float,
    window_size=11,
    window_width=2,
    optimize=False
):
    """
    Returns a large tuple of information that is the result of solving for
    the jitter based on the three input files.

    The first file path sets some variables for all of the runs.  The
    whichfrom booleans determine determines the sense of the offsets.
    A value of True makes the offsets relative to the From cube specified
    in the flat.tab, and False makes the offsets relative to the Match cube.

    :param file_path1: File Path for first flat.tab file.
    :param whichfrom1: Offsets relative to FROM for first flat.tab file.
    :param file_path2: File Path for second flat.tab file.
    :param whichfrom2: Offsets relative to FROM for second flat.tab file.
    :param file_path3: File Path for third flat.tab file.
    :param whichfrom3: Offsets relative to FROM for third flat.tab file.
    :param line_interval: The number of lines to use to set the number of
        Fourier transform intervals.
    :param window_size: The kernel size for median filtering the offsets,
        should be an odd integer.
    :param window_width: Sets the boundaries above and below the filtered
        average beyond which to exclude outliers.

    :return: A gnarly 18-tuple

    The values of the 18-tuple are:

    0
        Time values starting at *t0* (numpy array).
    1
        Sample direction jitter for each time (numpy array).
    2
        Line direction jitter for each time (numpy array).
    3
        The LineRate from the FROM and MATCH files (float).
    4
        The TDI for the FROM and MATCH files (int).
    5
        Zero time at which the Fourier Transform steps start (float).
    6
        Unique time values from file one (numpy array).
    7
        Unique time values from file two (numpy array).
    8
        Unique time values from file three (numpy array).
    9
        List of three numpy arrays representing the Gaussian filtered
        offsets in the sample direction from the three input files.
    10
        List of three numpy arrays representing the Gaussian filtered
        offsets in the line direction from the three input files.
    11
        List of three numpy arrays representing the interpolation of
        sample offsets for each input file at each time.
    12
        List of three numpy arrays representing the interpolation of
        line offsets for each input file at each time.
    13
        The average error of the jitter solution with respect to the
        measured offsets.
    14
        min_k
    15
        [min_jitter_check_x1, min_jitter_check_x2, min_jitter_check_x3],
    16
        [min_jitter_check_y1, min_jitter_check_y2, min_jitter_check_y3],
    17
        The first value of the sample and line from each file (dict).

    """
    t1, offx1, offy1, lines1, dt1, tdi1, linerate1 = parse_file(
        file_path1, window_size, window_width, whichfrom1
    )
    t2, offx2, offy2, lines2, dt2, tdi2, linerate2 = parse_file(
        file_path2, window_size, window_width, whichfrom2
    )
    t3, offx3, offy3, lines3, dt3, tdi3, linerate3 = parse_file(
        file_path3, window_size, window_width, whichfrom3
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

    offx_filtered = list(
        map(
            filter_data,
            itertools.repeat(nfft),
            itertools.repeat(2 / nfft),
            (offx1, offx2, offx3),
        )
    )
    offy_filtered = list(
        map(
            filter_data,
            itertools.repeat(nfft),
            itertools.repeat(2 / nfft),
            (offy1, offy2, offy3),
        )
    )

    # The values in tt are the fractional points in [0:1]
    # that correspond to the nfft number.
    tt = np.linspace(0, 1, nfft, endpoint=False)

    # The first file to be parsed sets nfftime, t0, and duration
    nfftime = np.linspace(t1[0], t1[-1], nfft)
    # et_shift = et - t1[0]
    t0 = t1[0]
    duration = t1[-1] - t0

    xinterp = [None] * 3
    yinterp = [None] * 3
    xinterp[0], yinterp[0], x1, y1, ddt1, overxx1, overyy1 = create_matrices(
        t1,
        offx_filtered[0],
        offy_filtered[0],
        dt1,
        duration,
        t0,
        nfft,
        nfftime,
        tt,
    )
    xinterp[1], yinterp[1], x2, y2, ddt2, overxx2, overyy2 = create_matrices(
        t2,
        offx_filtered[1],
        offy_filtered[1],
        dt2,
        duration,
        t0,
        nfft,
        nfftime,
        tt,
    )
    xinterp[2], yinterp[2], x3, y3, ddt3, overxx3, overyy3 = create_matrices(
        t3,
        offx_filtered[2],
        offy_filtered[2],
        dt3,
        duration,
        t0,
        nfft,
        nfftime,
        tt,
    )

    logger.info("Searching for correct phasetol")
    # For the test data, the following while loop will *always* run the full
    # number of repetitions.  If it were guaranteed that there was only
    # one minima in the changing value of *error*, then we could exit
    # this loop early, once the error starts going back up, or adopt
    # strategies to divide and conquer the phasetol parameter space.
    # A better understanding of the error behavior could allow us to speed
    # these loops up.
    #
    # Tests with scipy.optimize.minimize_scalar() indicate that this
    # data is not convex, and local minima confuse the minimizers.

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

    if optimize:
        # the scipy.optimize.minimize_scaler() methods got distracted by
        # a local minima
        opt_res = minimize_scalar(
            jitter_error,
            args=(
                tolcoef,
                tt,
                duration,
                rh0,
                (dt1, dt2, dt3),
                xinterp,
                yinterp,
                (ddt1, ddt2, ddt3),
                (overxx1, overxx2, overxx3),
                (overyy1, overyy2, overyy3)
            ),
            method="Brent",
            bracket=(0, 50),
            # method="bounded",
            # bounds=(0, 50),
            tol=error_tol,
            options=dict(maxiter=repetitions, disp=True)
        )
        logger.info(opt_res)
        min_k = opt_res.x
    else:
        while error >= error_tol and k < repetitions:
            k += 1

            error = jitter_error(
                k,
                tolcoef,
                tt,
                duration,
                rh0,
                (dt1, dt2, dt3),
                xinterp,
                yinterp,
                (ddt1, ddt2, ddt3),
                (overxx1, overxx2, overxx3),
                (overyy1, overyy2, overyy3)
            )

            if error < min_avg_error:
                min_avg_error = error
                min_k = k
        logger.info(f"Minimum Error after phase filtering: {min_avg_error}")

        # end while

    min_jitterx, min_jittery = jitterxy(
        min_k * tolcoef,
        (ddt1, ddt2, ddt3),
        (overxx1, overxx2, overxx3),
        (overyy1, overyy2, overyy3),
    )

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
        c = omega / (2.0 * nfft)

        jitterxx = filter_data(nfft, c, min_jitterx)
        jitteryy = filter_data(nfft, c, min_jittery)

        jitterxx = jitterxx - jitterxx[0]
        jitteryy = jitteryy - jitteryy[0]

        # checking
        jittercheckx1 = (
            np.interp(tt + dt1 / duration, tt, jitterxx, left=0, right=0)
            - jitterxx
        )
        jitterchecky1 = (
            np.interp(tt + dt1 / duration, tt, jitteryy, left=0, right=0)
            - jitteryy
        )

        jittercheckx2 = (
            np.interp(tt + dt2 / duration, tt, jitterxx, left=0, right=0)
            - jitterxx
        )
        jitterchecky2 = (
            np.interp(tt + dt2 / duration, tt, jitteryy, left=0, right=0)
            - jitteryy
        )

        jittercheckx3 = (
            np.interp(tt + dt3 / duration, tt, jitterxx, left=0, right=0)
            - jitterxx
        )
        jitterchecky3 = (
            np.interp(tt + dt3 / duration, tt, jitteryy, left=0, right=0)
            - jitteryy
        )

        error_vec = (
            1.0
            / 6.0
            * (
                np.abs(xinterp[0] - (jittercheckx1 + rh0["x1"]))
                + np.abs(xinterp[1] - (jittercheckx2 + rh0["x2"]))
                + np.abs(xinterp[2] - (jittercheckx3 + rh0["x3"]))
                + np.abs(yinterp[0] - (jitterchecky1 + rh0["y1"]))
                + np.abs(yinterp[1] - (jitterchecky2 + rh0["y2"]))
                + np.abs(yinterp[2] - (jitterchecky3 + rh0["y3"]))
            )
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

    return (
        nfftime,
        sample,
        line,
        linerate,
        tdi,
        t0,
        t1,
        t2,
        t3,
        offx_filtered,
        offy_filtered,
        xinterp,
        yinterp,
        min_avg_error,
        min_k,
        [min_jitter_check_x1, min_jitter_check_x2, min_jitter_check_x3],
        [min_jitter_check_y1, min_jitter_check_y2, min_jitter_check_y3],
        rh0,
    )


def jitterxy(phasetol, ddt, overxx, overyy):
    # it = iter((ddt, overxx, overyy))
    # the_len = len(next(it))
    # if not all(len(l) == the_len for l in it):
    #     raise ValueError('Not all lists have same length!')

    # null the frequencies that cause a problem (really zero them out)
    masked_overxx = list()
    masked_overyy = list()
    for (d, x, y) in zip(ddt, overxx, overyy):
        xxx, yyy = mask_frequencies(phasetol, d, x, y)
        masked_overxx.append(xxx)
        masked_overyy.append(yyy)

    # Adding all frequencies together
    stackedx = np.ma.stack(masked_overxx)
    stackedy = np.ma.stack(masked_overyy)

    overxxx = np.ma.mean(stackedx, axis=0)
    overyyy = np.ma.mean(stackedy, axis=0)

    # take the sum of each row
    overx = np.ma.sum(overxxx, axis=0)
    overy = np.ma.sum(overyyy, axis=0)

    jitterx = overx - overx[0]
    jittery = overy - overy[0]

    return jitterx, jittery


def jitter_error(
    k, tolcoef, tt, duration, rh0, dt, xinterp, yinterp, ddt, overxx, overyy
):
    # setting the phase tolerance
    phasetol = k * tolcoef

    jitterx, jittery = jitterxy(
        phasetol,
        ddt,
        overxx,
        overyy,
    )

    # checking
    jittercheckx1 = (
        np.interp(tt + dt[0] / duration, tt, jitterx, left=0, right=0)
        - jitterx
    )
    jitterchecky1 = (
        np.interp(tt + dt[0] / duration, tt, jittery, left=0, right=0)
        - jittery
    )

    jittercheckx2 = (
        np.interp(tt + dt[1] / duration, tt, jitterx, left=0, right=0)
        - jitterx
    )
    jitterchecky2 = (
        np.interp(tt + dt[1] / duration, tt, jittery, left=0, right=0)
        - jittery
    )

    jittercheckx3 = (
        np.interp(tt + dt[2] / duration, tt, jitterx, left=0, right=0)
        - jitterx
    )
    jitterchecky3 = (
        np.interp(tt + dt[2] / duration, tt, jittery, left=0, right=0)
        - jittery
    )

    error_vec = (
                    np.abs(xinterp[0] - (jittercheckx1 + rh0["x1"]))
                    + np.abs(xinterp[1] - (jittercheckx2 + rh0["x2"]))
                    + np.abs(xinterp[2] - (jittercheckx3 + rh0["x3"]))
                    + np.abs(yinterp[0] - (jitterchecky1 + rh0["y1"]))
                    + np.abs(yinterp[1] - (jitterchecky2 + rh0["y2"]))
                    + np.abs(yinterp[2] - (jitterchecky3 + rh0["y3"]))
                ) / 6.0

    error = error_vec.mean()
    logger.info(f"Error for phasetol {phasetol}: {error}")
    return error


def create_matrices(
    time: np.array,
    offx: np.array,
    offy: np.array,
    dt: float,
    duration: float,
    t0: float,
    nfft: int,
    nfftime: np.array,
    tt: np.array,
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
    xa = x[: int(nfft / 2)].real
    xb = -1 * x[: int(nfft / 2)].imag
    ya = y[: int(nfft / 2)].real
    yb = -1 * y[: int(nfft / 2)].imag

    # calculates the phase difference
    twopi = math.pi * 2
    ddt_temp = (dt / duration * twopi) * freq
    ddt = ddt_temp - twopi * np.floor(ddt_temp / twopi)

    # the coeficients for the frequencies
    with np.errstate(divide="ignore"):
        aaax = (
            -0.5
            * (-1 * xa * np.cos(ddt) + np.sin(ddt) * xb - xa)
            / np.sin(ddt)
        )
        aaay = (
            -0.5
            * (-1 * ya * np.cos(ddt) + np.sin(ddt) * yb - ya)
            / np.sin(ddt)
        )

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
    front_padding = math.floor(nfft - len(data) / 2)
    back_padding = math.ceil(nfft - len(data) / 2)

    # Apply the padding
    padded = np.concatenate(
        (
            np.array([data[0]] * front_padding),
            data,
            np.array([data[-1]] * back_padding),
        )
    )

    freq = np.fft.fft(padded)

    # The exponential
    exp_vec = np.hstack(
        (np.linspace(0, nfft - 1, nfft), np.linspace(-1 * nfft, -1, nfft))
    )

    exponential = np.exp(-1 * c ** 2 * exp_vec ** 2)

    # The ifft of the product
    filtered = np.fft.ifft(freq * exponential)

    # Remove the padding and take the real part
    return filtered[front_padding : front_padding + len(data)].real


def mask_frequencies(phasetol: float, ddt: np.array, x: np.array, y: np.array):
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

    magnitude = np.sqrt(offx_arr ** 2 + offy_arr ** 2)
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
    t: np.array, sample: np.array, line: np.array, linerate: float, tdi: int,
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
    :return: six-tuple which contains:
        0. sample maximum smear value (float)
        1. line maximum smear value (float)
        2. Maximum smear magnitude (float)
        3. Sample smear values (numpy array)
        4. Line smear values (numpy array)
        5. Ephemeris times at linerate intervals (numpy array)
    """
    # The output here differs from the C++ output in a very, very minor
    # way.  The absolute difference between the line and sample smear
    # is 0.002 pixels.  I could not track this down.  It may be due to
    # precision differences in the double values that C++ uses, but I'm
    # just not sure.

    # The array T has large values. Use a shifted version of it when
    # interpolating to reduce the effect of those values on numerical
    # accuracy.
    shifted_t = t - t[0]

    # xi = T(1):linerate:T(end);
    n = math.floor((shifted_t[-1] - shifted_t[0]) / linerate) + 1
    # This is from the original code, but by definition shifted_t[0] is zero.
    # xi = np.linspace(0, n - 1, n) * linerate + shifted_t[0]
    xi = np.linspace(0, n - 1, n) * linerate

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
    mag_smear = np.sqrt(dysdx ** 2 + dyldx ** 2)

    # Find maxSmearS, the largest element by magnitude in dysdx
    msi = np.argmax(np.abs(dysdx))
    max_smear_s = dysdx[msi]

    # Find maxSmearL, the largest element by magnitude in dyldx
    msi = np.argmax(np.abs(dyldx))
    max_smear_l = dyldx[msi]

    # Find maxSmearMag, the largest element by magnitude in magSmear
    max_smear_mag = np.max(mag_smear)

    # Outputs
    return max_smear_s, max_smear_l, max_smear_mag, dysdx, dyldx, xi


def set_file_path(location: Path, file_path: Path):
    if file_path.parent == location:
        return file_path
    else:
        return location / file_path.name


def plot(
    t,
    x,
    y,
    xinterp,
    yinterp,
    jittercheckx,
    jitterchecky,
    title,
    et,
    sample,
    line,
    show=True,
    save=False,
):
    plt.ioff()
    fig = plt.figure(constrained_layout=True)
    gs = fig.add_gridspec(3, 6)

    fig.suptitle("Resolve Jitter Results")

    ax00 = fig.add_subplot(gs[0, 0:2])
    ax00.set_title(title[0])
    ax00.set_ylabel("Sample Offset")

    ax01 = fig.add_subplot(gs[0, 2:4])
    ax01.set_title(title[1])

    ax02 = fig.add_subplot(gs[0, 4:])
    ax02.set_title(title[2])

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

    ax00.plot(t[0], x[0], "o", c="red")
    ax00.plot(et, xinterp[0], c="green")
    ax00.plot(et, jittercheckx[0], c="yellow")

    ax01.plot(t[1], x[1], "o", c="red")
    ax01.plot(et, xinterp[1], c="green")
    ax01.plot(et, jittercheckx[1], c="yellow")

    ax02.plot(t[2], x[2], "o", c="red")
    ax02.plot(et, xinterp[2], c="green")
    ax02.plot(et, jittercheckx[2], c="yellow")

    ax10.plot(t[0], y[0], "o", c="red")
    ax10.plot(et, yinterp[0], c="green")
    ax10.plot(et, jitterchecky[0], c="yellow")

    ax11.plot(t[1], y[1], "o", c="red")
    ax11.plot(et, yinterp[1], c="green")
    ax11.plot(et, jitterchecky[1], c="yellow")

    ax12.plot(t[2], y[2], "o", c="red")
    ax12.plot(et, yinterp[2], c="green")
    ax12.plot(et, jitterchecky[2], c="yellow")

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
    with open(path, "w", newline="") as csvfile:
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


def write_smear_data(
    path: Path,
    max_smear_s,
    max_smear_l,
    max_smear_mag,
    dysdx,
    dyldx,
    et,
    image_id=None,
):
    # Make a text file of the smear data
    if image_id is None:
        id_str = ""
    else:
        id_str = f" for {image_id}"

    smear_text = [
        f"""\
# Smear values are calculated from the derived jitter function{id_str}.
# Maximum Cross-track pixel smear {max_smear_s}
# Maximum Down-track pixel smear {max_smear_l}
# Maximum Pixel Smear Magnitude {max_smear_mag}
# Sample                 Line                   EphemerisTime"""
    ]

    for ess, ell, exi in zip(dysdx, dyldx, et):
        smear_text.append(f"{ess}     {ell}     {exi}")

    logger.info(f"Writing: {path}")
    path.write_text("\n".join(smear_text))

    return


def write_gnuplot_file(
    gnuplot_path: Path,
    data_path: Path,
    img_path: Path,
    file_path1: Path,
    file_path2: Path,
    file_path3: Path,
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
