#!/usr/bin/env python
"""Perform fast-fourier-transform cleaning on HiRISE Channels."""

# Copyright 2021, Ross A. Beyer (rbeyer@seti.org)
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
# This program is based on procedures and code written by
# Audrie Fennema and Sarah Mattson circa April 2013.

import argparse
import logging
import os
import shutil
import sys
from datetime import datetime
from pathlib import Path

import pvl
import kalasiris as isis

import hiproc.util as util

logger = logging.getLogger(__name__)


def arg_parser():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        parents=[util.parent_parser()],
    )
    parser.add_argument(
        "-o", "--output",
        required=False,
        default=".fftclean.cub",
        help="The name of the output .cub file to write.  Optionally, if "
             "it starts with a '.' it is considered a suffix "
             "and the input file will have its suffix removed, and this will "
             "be added to input file name. Default: %(default)s",
    )
    parser.add_argument("file", help="ISIS Cube file to clean.")
    return parser


def main():
    args = arg_parser().parse_args()

    util.set_logger(args.verbose, args.logfile, args.log)

    out_p = util.path_w_suffix(args.output, args.file)

    with util.main_exceptions(args.verbose):
        clean(
            args.file,
            out_p,
            keep=args.keep,
        )
    sys.exit(0)


def clean(
    in_path: os.PathLike,
    out_path: os.PathLike,
    keep=False,
):
    """The file at *out_path* will be the result of running fft
    cleaning of the file at *in-path*.

    If *keep* is True, then all intermediate files will be preserved,
    otherwise, this function will clean up any intermediary files
    it creates.
    """
    in_p = Path(in_path)
    out_p = Path(out_path)

    temp_token = datetime.now().strftime("fftclean-%y%m%d%H%M%S")
    to_delete = isis.PathSet()

    # Check input channel cube for gaps. If they exist, interpolate and run
    # fft_clean on the interpolated gaps cube. Replace gaps in final cleaned
    # image.
    gaplog = to_delete.add(in_p.with_suffix(f".{temp_token}.gaplog.pvl"))
    isis.findgaps(in_p, log=gaplog)
    gappvl = pvl.load(gaplog)
    next_cube = in_p
    if "Gap" in gappvl:
        logger.info(f"{in_p.name} has gaps.")
        next_cube = to_delete.add(
            in_p.with_suffix(f".{temp_token}.fillgap.cub")
        )
        isis.fillgap(in_p, to=next_cube, interp="linear", direction="sample")
    else:
        logger.info(f"{in_p.name} has no gaps.")

    # Take the FFT of the input channel cube
    magcube = to_delete.add(in_p.with_suffix(f".{temp_token}.mag.cub"))
    phasecube = to_delete.add(in_p.with_suffix(f".{temp_token}.pha.cub"))
    isis.fft(next_cube, magnitude=magcube, phase=phasecube)

    samps = int(isis.getkey_k(magcube, "Dimensions", "Samples"))
    lines = int(isis.getkey_k(magcube, "Dimensions", "Lines"))

    # Make the filter for the magnitude image
    # Make circle mask based on butterworth filter idea (changed /45 to /25)
    circ = to_delete.add(in_p.with_suffix(f".{temp_token}.circ.cub"))
    isis.fx(
        to=circ,
        equation="1 - (1/{--pi*[1+((line - 51)/45)^2]} + "
                 "1/{--pi*[1 + ((sample - 51)/45)^2]})",
        mode="outputonly",
        lines=101,
        samples=101
    )
    circ_str = to_delete.add(
        in_p.with_suffix(f".{temp_token}.circ.stretch.cub")
    )
    isis.stretch(
        circ,
        to=circ_str,
        usepercentages=True,
        pairs="0.0:1 0.01:1 75:.99 80:0.9 100:0"
    )
    circ_enl = to_delete.add(
        in_p.with_suffix(f".{temp_token}.circ.enl.cub")
    )
    isis.enlarge(
        circ_str,
        to=circ_enl,
        interp="cubic",
        mode="total",
        ons=samps,
        onl=lines
    )

    # Make full size cross mask
    crossfull = to_delete.add(
        in_p.with_suffix(f".{temp_token}.cross_full.cub")
    )
    isis.fx(
        to=crossfull,
        equation=f"[1 - e^(--pi*[(line - (( {lines} / 2 ) +1))^6])] * "
                 f"[1 - e^(--pi*[(sample - (( {samps} / 2) +1))^6])]",
        mode="outputonly",
        lines=lines,
        samples=samps
    )

    maskcube = to_delete.add(
        in_p.with_suffix(f".{temp_token}.mask.cub")
    )
    isis.fx(
        to=maskcube,
        f1=crossfull,
        f2=circ_enl,
        equation="f1 * f2"
    )

    # Filter the mag image before detecting spikes
    # (filter out low frequency information)
    magfilt = in_p.with_suffix(f".{temp_token}.magfilt.cub")
    isis.fx(f1=magcube, f2=maskcube, to=magfilt, equation="f2 * abs(f1)")

    # Normalize the filtered magnitude image
    magnorm = float(pvl.loads(isis.stats(magfilt).stdout)["Results"]["Median"])

    magfilt_cn = in_p.with_suffix(f".{temp_token}.magfilt.cubenorm.cub")
    isis.fx(f1=magfilt, to=magfilt_cn, equation=f"f1/{magnorm}")

    # Run a Sobel (gradient) filter on the filtered magnitude image to
    # emphasize spikes.
    magfilt_sobel = in_p.with_suffix(f".{temp_token}.magfilt.sobel.cub")
    isis.gradient(magfilt_cn, to=magfilt_sobel, gradtype="sobel")

    # Threshold the normalized, filtered mag image to find spikes.
    # N times average empirically derived value. Typically n=5
    threshold = 5
    threshcube = in_p.with_suffix(f".{temp_token}.threshold.cub")
    isis.fx(
        f1=magfilt_sobel,
        to=threshcube,
        equation=f"f1 > ( {threshold} * cubeavg(f1) )"
    )

    # Use a Gaussian filter to broaden the spikes.
    gaus_cube = in_p.with_suffix(f".{temp_token}.threshold.gaus.cub")
    isis.gauss(threshcube, to=gaus_cube, size=5, stddev=6.0)

    # Threshold the Gaussian filtered cube and re-apply the cross mask to it.
    # This is to re-exclude the cross which the Gaussian filter might have
    # spread into.
    thresh2 = in_p.with_suffix(f".{temp_token}.threshold2.cub")
    isis.fx(f1=gaus_cube, to=thresh2, equation=" f1 > 0 ")

    wcross = in_p.with_suffix(f".{temp_token}.threshold2.wcross.cub")
    isis.fx(f1=thresh2, f2=crossfull, to=wcross, equation=" f1 * f2 ")

    # Multiply the binary threshold image by the magnitude maximum, plus 1.
    magmax = float(pvl.loads(isis.stats(magcube).stdout)["Results"]["Maximum"])
    magfilt_divby = in_p.with_suffix(f".{temp_token}.magfilt.divideby.cub")
    isis.fx(f1=wcross, to=magfilt_divby, equation=f"1 + (f1 * {magmax})")

    # Divide the unfiltered magnitude image by the thresholded image to
    # minimize the noise spikes.
    magclean = in_p.with_suffix(f".{temp_token}.magclean.cub")
    isis.fx(f1=magcube, f2=magfilt_divby, equation="f1 / f2", to=magclean)

    # Perform the inverse transform with the cleaned magnitude image and the
    # original phase image.
    cleaned = in_p.with_suffix(f".{temp_token}.clean.cub")
    isis.ifft(magnitude=magclean, phase=phasecube, to=cleaned)

    if "Gap" in gappvl:
        gapclean = in_p.with_suffix(f".{temp_token}.cleangapped.cub")
        isis.fx(f1=cleaned, f2=in_p, to=gapclean, equation=" f1 + (0 * f2 )")
        cleaned = gapclean

    # Subtract the cleaned image from the cubenormed image to make a
    # difference cube. Gather stats.
    diff = in_p.with_suffix(f".{temp_token}.diff.cub")
    isis.algebra(in_p, from2=cleaned, to=diff, operator="subtract")
    # Command("stats from=$base.cubenorm.cub");
    # Command("stats from=$base.clean.cub");
    # Command("stats from=$base.diff.cub > $base.diff.stats.txt");

    # Command("/bin/rm -f $base.magfilt.divideby.cub $base.magfilt.cubenorm.cub");

    # Make a crop of the difference cubes for validation.
    crop_line_start = int((lines / 2) - (samps / 2))
    diff_crop = in_p.with_suffix(f".{temp_token}.diff.crop.cub")
    isis.crop(
        diff,
        to=diff_crop,
        sample=1,
        nsamples=samps,
        line=crop_line_start,
        nline=samps
    )

    shutil.copy(cleaned, out_p)

    if not keep:
        to_delete.unlink()

    return


if __name__ == "__main__":
    main()
