#!/usr/bin/env python
"""Deal with stuck bits in HiRISE pixels via unflipping.

This program will attempt to unflip the bit-flipped pixels.
Only the image area is affected. This is also an older approach
and is less supported than bitflips.py.

Please see bitflips.py for a more detailed explanation of
the HiRISE bit-flip problem.

This file contains functions to deal with 'unflipping' the bits.
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
# Parts of this were inspired by clean_bit_flips.pro by Alan Delamere,
# Oct 2019, but the implementation here was written from scratch.

import argparse
import itertools
import logging
import math
import shutil
import subprocess
import sys
from pathlib import Path

import kalasiris as isis

import hiproc.util as util

logger = logging.getLogger(__name__)


def main():
    try:
        parser = argparse.ArgumentParser(
            description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter,
            parents=[util.parent_parser()],
        )
        parser.add_argument(
            "-o", "--output", required=False, default=".bitunflip.cub"
        )
        parser.add_argument("cube", help="ISIS Cube file.")

        args = parser.parse_args()

        util.set_logger(args.verbose, args.logfile, args.log)

        out_p = util.path_w_suffix(args.output, args.cube)

        unflip(Path(args.cube), out_p, keep=args.keep)

    except subprocess.CalledProcessError as err:
        print("Had an ISIS error:", file=sys.stderr)
        print(" ".join(err.cmd), file=sys.stderr)
        print(err.stdout, file=sys.stderr)
        print(err.stderr, file=sys.stderr)


def get_range(start, stop, step=None) -> range:
    """Returns a range object given floating point start and
    stop values (the optional step must be an int).  If
    step is none, it will be set to 1 or -1 depending on the
    relative values of start and stop.
    """
    if start < stop:
        start = math.floor(start)
        stop = math.ceil(stop)
        if step is None:
            step = 1
    else:
        start = math.ceil(start)
        stop = math.floor(stop)
        if step is None:
            step = -1

    return range(start, stop, step)


def find_gap(hist: list, start, stop, step=None, findfirst=False) -> int:
    """Returns a DN between *start* and *stop* which is determined to be a gap.

    Given the range specified by *start*, *stop* and *step*, the
    DNs are extracted from the list of namedtuples (which must
    conform to the HistRow namedtuple from kalasiris.Histogram).
    The 'gaps' between continuous ranges of DN are found (DNs where
    the histogram has no pixels).  The largest gap is found, and
    then the DN value closest to start from that largest gap is
    returned.

    If step is not specified or None, it will be set to 1 or -1
    depending on the relative values of start and stop.

    If *findfirst* is `True`, then rather than finding the 'biggest'
    gap, it will return the DN from the gap that is closest to
    *start*.
    """
    dn_window = get_range(start, stop, step)

    hist_set = set(filter(lambda d: d in dn_window, (int(x.DN) for x in hist)))

    missing = sorted(
        set(dn_window).difference(hist_set), reverse=(dn_window.step < 0)
    )
    # logging.info(f'{bs} missing: ' + str(missing))

    if len(missing) > 0:
        sequences = list()
        for k, g in itertools.groupby(
            missing,
            (lambda x, c=itertools.count(): next(c) - (dn_window.step * x)),
        ):
            sequences.append(list(g))

        if findfirst:
            return sequences[0][0]
        else:
            # find biggest gap
            return max(sequences, key=len)[0]
    else:
        # There was no gap in DN.
        raise ValueError(
            "There was no gap in the DN window from " f"{start} to {stop}."
        )


def find_min_dn(hist: list, start, stop, step=None) -> int:
    """Returns the DN from *hist* with the lowest pixel count.

    Where *hist* is a list of namedtuples (which must conform to
    the HistRow namedtuple from kalasiris.Histogram).  The *start*,
    *stop*, and *step* parameters are DN values in that *hist* list.

    Given how the min() mechanics work, this should return the first
    DN value of the entry with the lowest pixel count, in range order
    (if there is more than one).
    """
    r = get_range(start, stop, step)
    hist_list = sorted(
        filter(lambda x: int(x.DN) in r, hist),
        key=lambda x: int(x.DN),
        reverse=(r.step < 0),
    )

    h = min(hist_list, key=lambda x: int(x.Pixels))
    return int(h.DN)


def subtract_over_thresh(
    in_path: Path, out_path: Path, thresh: int, delta: int, keep=False
):
    """This is a convenience function that runs ISIS programs to add or
    subtract a value to DN values for pixels that are above or below
    a threshold.

    For all pixels in the *in_path* ISIS cube, if *delta* is positive,
    then pixels with a value greater than *thresh* will have *delta*
    subtracted from them.  If *delta* is negative, then all pixels
    less than *thresh* will have *delta* added to them.
    """

    # Originally, I wanted to just do this simply with fx:
    # eqn = "\(f1 + ((f1{glt}={thresh}) * {(-1 * delta)}))"
    # However, fx writes out floating point pixel values, and we really
    # need to keep DNs as ints as long as possible.  Sigh.

    shutil.copyfile(in_path, out_path)

    mask_p = in_path.with_suffix(".threshmask.cub")
    mask_args = {"from": in_path, "to": mask_p}
    if delta > 0:
        mask_args["min"] = thresh
    else:
        mask_args["max"] = thresh
    isis.mask(**mask_args)

    delta_p = in_path.with_suffix(".delta.cub")
    isis.algebra(
        mask_p, from2=in_path, to=delta_p, op="add", a=0, c=(-1 * delta)
    )

    isis.handmos(delta_p, mosaic=out_path)

    if not keep:
        mask_p.unlink()
        delta_p.unlink()

    return


def mask_gap(in_path: Path, out_path: Path):
    """Attempt to mask out pixels beyond the central DNs of the median
    based on gaps in the histogram.

    This approach worked well based on 'ideal' reverse-clocked data
    or 'dark' images, but in 'real' HiRISE images of Mars, the reality
    is that there are 'gaps' everywhere along the DN range, and this
    approach ends up being too 'dumb'.
    """

    hist = isis.Histogram(in_path)

    median = math.trunc(float(hist["Median"]))
    # std = math.trunc(float(hist['Std Deviation']))

    high = find_gap(
        hist, median, math.trunc(float(hist["Maximum"])), findfirst=True
    )

    low = find_gap(
        hist, median, math.trunc(float(hist["Minimum"])), findfirst=True
    )

    highdist = high - median
    lowdist = median - low

    maskmax = find_gap(hist, median, (median + (2 * highdist)))
    maskmin = find_gap(hist, median, (median - (2 * lowdist)))

    isis.mask(in_path, to=out_path, minimum=maskmin, maximum=maskmax)

    return


def get_unflip_thresh(hist: list, far, near, lowlimit) -> int:
    """Provides a robust threshhold for the unflipping process."""
    # try:
    thresh = find_gap(hist, far, near)
    # except ValueError:
    #     logging.info('No zero threshold found.')
    #     # d is not defined here, and I think this was the victim of
    #     # some cut'n'paste work, but would need to go analyze the
    #     # repo history to track it down.  For the moment, this will
    #     # just error!
    #     if d >= lowlimit:
    #         thresh = find_min_dn(hist, far, near)
    #         logging.info(f'Found a minimum threshold: {thresh}')
    #         if thresh == far or thresh == near:
    #             raise ValueError(f'Minimum, {thresh}, at edge of range.')
    #     else:
    #         raise ValueError(f'Below {lowlimit}, skipping.')
    return thresh


def unflip(in_p: Path, out_p: Path, keep=False):
    """Attempt to indentify DNs whose bits have been flipped in the
    ISIS cube indicated by *in_p*, and unflip them.

    Only operates on image-area.
    """
    to_del = isis.PathSet()

    # This full suite of deltas works well for the reverse-clock area
    # and even 'dark' images, but not 'real' images.
    # deltas = (8192, 4096, 2048, 1024, 512, 256, 128, 64)

    # Restricting the number of deltas might help, but this seems
    # arbitrary.
    deltas = (8192, 4096, 2048, 1024)

    count = 0
    suffix = ".bf{}-{}{}.cub"
    this_p = to_del.add(in_p.with_suffix(suffix.format(count, 0, 0)))
    this_p.symlink_to(in_p)

    median = math.trunc(float(isis.stats_k(in_p)["Median"]))

    for (sign, pm, extrema) in ((+1, "m", "Maximum"), (-1, "p", "Minimum")):
        # logging.info(pm)
        for delt in deltas:
            d = sign * delt
            far = median + d
            near = median + (d / 2)

            hist = isis.Histogram(this_p)

            logger.info(
                f"bitflip position {pm}{delt}, near: {near} "
                f"far: {far}, extrema: {hist[extrema]}"
            )
            if (sign > 0 and far < float(hist[extrema])) or (
                sign < 0 and far > float(hist[extrema])
            ):
                count += 1
                s = suffix.format(count, pm, delt)
                next_p = to_del.add(this_p.with_suffix("").with_suffix(s))
                try:
                    thresh = get_unflip_thresh(hist, far, near, d)
                except ValueError as err:
                    logger.info(err)
                    count -= 1
                    break
                subtract_over_thresh(this_p, next_p, thresh, d, keep=keep)
                this_p = next_p
            else:
                logger.info(
                    "The far value was beyond the extrema. " "Didn't bother."
                )

    shutil.move(this_p, out_p)

    if not keep:
        to_del.remove(this_p)
        to_del.unlink()

    return


if __name__ == "__main__":
    main()
