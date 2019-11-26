#!/usr/bin/env python
"""Deal with stuck bits in HiRISE pixels.

If there is some process which is flipping one of the HiRISE
bits in a 14-bit pixel the wrong way, it will inadvertently add
or subtract a certain DN value to a pixel.  Naturally, it goes
from 2^13 (8192) all the way down to 2.

In practice, these 'bitflipped' pixels show up as spikes in the
image histogram at DN values approximately centered on the image
Median plus or minus the bitflip value (8192, 4096, etc.).

These are not narrow spikes and have variable width, so care
must be taken when dealing with them.  This program will do one
of two things, it will either attempt to un-flip the bits on
the bit-flipped pixels or it will attempt to mask them out.
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
# This program was originally inspired by clean_bit_flips.pro by Alan Delamere,
# Oct 2019, but was written from scratch.

import argparse
import itertools
import logging
import math
import os
import shutil
import subprocess
from pathlib import Path

import numpy as np
from scipy.signal import find_peaks

import pvl

import pyrise.hirise as hirise
import pyrise.util as util
import kalasiris as isis


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     parents=[util.parent_parser()])
    parser.add_argument('-o', '--output',
                        required=False, default='.bitflip.cub')
    parser.add_argument('-m', '--mask', required=False, action='store_true',
                        help=('If set, the program will mask the '
                              'bit-flipped pixels, otherwise will '
                              'try to collapse them.'))
    parser.add_argument('cube', metavar="some.cub-file", nargs='+',
                        help='More than one can be listed here.')

    args = parser.parse_args()

    util.set_logging(args.log, args.logfile)

    if(len(args.cube) > 1 and
       not args.output.startswith('.')):
        logging.critical('With more than one input cube file, the --output '
                         'must start with a period, and it '
                         f'does not: {args.output}')
        sys.exit()

    for i in args.cube:
        out_p = util.path_w_suffix(args.output, i)

        try:
            if args.mask:
                mask(Path(i), out_p, keep=args.keep)
            else:
                unflip(Path(i), out_p, keep=args.keep)
        except subprocess.CalledProcessError as err:
            print('Had an ISIS error:')
            print(' '.join(err.cmd))
            print(err.stdout)
            print(err.stderr)
            raise err
    return


def get_range(start, stop, step=None):
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
    """Given the range specified by start, stop and step, the DNs are
    extracted from the list of namedtuples (from kalasiris.Histogram)
    the 'gaps' between continuous ranges of DN are found (DNs where
    the histogram has no pixels).  The largest gap is found, and
    then the DN value closest to start from that largest gap is
    returned.

    If step is not specified or None, it will be set to 1 or -1
    depending on the relative values of start and stop.

    If findfirst is True, then rather than finding the 'biggest'
    gap, it will return the DN from the gap that is closest to
    start.
    """
    dn_window = get_range(start, stop, step)

    hist_set = set(filter(lambda d: d in dn_window,
                          (int(x.DN) for x in hist)))

    missing = sorted(set(dn_window).difference(hist_set),
                     reverse=(dn_window.step < 0))
    # logging.info(f'{bs} missing: ' + str(missing))

    if len(missing) > 0:
        sequences = list()
        for k, g in itertools.groupby(missing,
                                      (lambda x,
                                       c=itertools.count(): next(c) -
                                       (dn_window.step * x))):
            sequences.append(list(g))

        if findfirst:
            return sequences[0][0]
        else:
            # find biggest gap
            return max(sequences, key=len)[0]
    else:
        # There was no gap in DN.
        raise ValueError('There was no gap in the DN window from '
                         f'{start} to {stop}.')


def find_min_dn(hist: list, start, stop, step=None):
    """Given how the min() mechanics work, this should return the first
    DN value of the entry with the lowest pixel count, in range order
    (if there is more than one).
    """
    r = get_range(start, stop, step)
    hist_list = sorted(filter(lambda x: int(x.DN) in r, hist),
                       key=lambda x: int(x.DN),
                       reverse=(r.step < 0))

    h = min(hist_list, key=lambda x: int(x.Pixels))
    return int(h.DN)


def subtract_over_thresh(in_path: os.PathLike, out_path: os.PathLike,
                         thresh: int, delta: int, keep=False):
    """For all pixels in the in_path ISIS cube, if delta is positive, then
    pixels with a value greater than thresh will have delta subtracted from
    them.  If delta is negative, then all pixels less than thresh
    will have delta added to them.
    """

    # Originally, I wanted to just do this simply with fx:
    # eqn = "\(f1 + ((f1{glt}={thresh}) * {(-1 * delta)}))"
    # However, fx writes out floating point pixel values, and we really
    # need to keep DNs as ints as long as possible.  Sigh.

    shutil.copyfile(in_path, out_path)

    # stats = isis.stats_k(in_path)
    # dnmin = math.trunc(float(stats['Minimum']))
    # dnmax = math.trunc(float(stats['Maximum']))

    # suffix = f'+SignedWord+{dnmin}:{dnmax}'
    mask_p = in_path.with_suffix('.threshmask.cub')
    # mask_args = {'from': in_path, 'to': str(mask_p) + suffix}
    mask_args = {'from': in_path, 'to': mask_p}
    if delta > 0:
        mask_args['min'] = thresh
    else:
        mask_args['max'] = thresh
    util.log(isis.mask(**mask_args).args)

    delta_p = in_path.with_suffix('.delta.cub')
    # util.log(isis.algebra(mask_p, to=str(delta_p) + suffix, op='unary',
    #                      a=0, c=delta).args)
    util.log(isis.algebra(mask_p, from2=in_path, to=delta_p,
                          op='add', a=0, c=(-1 * delta)).args)

    util.log(isis.handmos(delta_p, mosaic=out_path).args)

    if not keep:
        mask_p.unlink()
        delta_p.unlink()

    return


def histogram(in_path: Path, hist_path: Path):
    """This is just a convenience function to facilitate logging."""
    util.log(isis.hist(in_path, to=hist_path).args)
    return isis.Histogram(hist_path)


def mask(in_path: Path, out_path: Path, keep=False):
    """Attempt to mask out pixels beyond the central DNs of the median
    based on minima in the histogram."""
    from PyRISE.HiCal import analyze_cubenorm_stats2

    to_del = isis.PathSet()

    hist_p = in_path.with_suffix('.hist')
    hist = histogram(in_path, hist_p)

    median = math.trunc(float(hist['Median']))

    cubenorm_stats_file = to_del.add(in_path.with_suffix('.cn.stats'))
    util.log(isis.cubenorm(in_path, stats=cubenorm_stats_file).args)
    (mindn, maxdn) = analyze_cubenorm_stats2(cubenorm_stats_file, median,
                                             5)
    (maskmin, maskmax) = find_smart_window(hist, mindn, maxdn, median,
                                           plot=True)

    # util.log(isis.mask(in_path, to=out_path, minimum=maskmin,
    #                    maximum=maskmax).args)

    if not keep:
        hist_p.unlink()
    return


def find_minima_index(central_idx: int, limit_idx: int,
                      minima_idxs: np.ndarray, pixel_counts: np.ndarray) -> int:
    '''Searches the pixel_counts at the minima_idxs positions between
       central_idx and limit_idx to find the one with the lowest pixel_count
       value (if multiple, choose the closest to the limit_idx position).

       If no minima_idx values are found in the range, the next value in
       minima_idx past the limit_idx will be returned.
    '''
    info_str = 'Looking for {} {} range ...'
    try:
        if limit_idx < central_idx:
            logging.info(info_str.format('min', 'inside'))
            inrange_iter = filter(lambda i: i < central_idx and i >= limit_idx,
                                  minima_idxs)
            select_idx = 0
        else:
            logging.info(info_str.format('max', 'inside'))
            inrange_iter = filter(lambda i: i > central_idx and i <= limit_idx,
                                  minima_idxs)
            select_idx = -1

        inrange_i = np.fromiter(inrange_iter, int)
        logging.info(str(inrange_i) + ' are the indexes inside the range.')

        value = min(np.take(pixel_counts, inrange_i))
        logging.info(f'{value} is the minimum Pixel count amongst those '
                     'indexes.')

        idx = inrange_i[np.asarray(pixel_counts[inrange_i] ==
                                   value).nonzero()][select_idx]
    except ValueError:
        if limit_idx < central_idx:
            logging.info(info_str.format('min', 'outside'))
            idx = max(filter(lambda i: i < limit_idx, minima_idxs))
        else:
            logging.info(info_str.format('max', 'outside'))
            idx = min(filter(lambda i: i > limit_idx, minima_idxs))

    logging.info(f'{idx} is the minimum index.')
    return idx


def find_smart_window(hist: list, mindn: int, maxdn: int,
                      centraldn: int, plot=False) -> tuple:
    '''Returns a minimum and maximum DN value from hist which are
       based on using the find_minima_index() function with the
       given mindn, maxdn, and centraldn values.

       If plot is True, this function will display a plot describing
       its work.  The curve represents the hist values, the shaded
       area marks the window between the given mindn and maxdn.  The
       'x'es mark all the detected minima, and the red dots indicate
       the minimum and maximum DN values that this function will
       return.
    '''
    hist_list = sorted(hist, key=lambda x: int(x.DN))
    pixel_counts = np.fromiter((int(x.Pixels) for x in hist_list), int)
    dn = np.fromiter((int(x.DN) for x in hist_list), int)

    mindn_i = (np.abs(dn - mindn)).argmin()
    maxdn_i = (np.abs(dn - maxdn)).argmin()
    central_i = np.where(dn == centraldn)

    minima_i, _ = find_peaks(np.negative(pixel_counts))

    min_i = find_minima_index(central_i, mindn_i, minima_i, pixel_counts)
    max_i = find_minima_index(central_i, maxdn_i, minima_i, pixel_counts)

    if plot:
        import matplotlib.pyplot as plt

        indices = np.arange(0, len(pixel_counts))
        dn_window = np.fromiter(map((lambda i: i >= mindn_i and i <= maxdn_i),
                                    (x for x in range(len(pixel_counts)))),
                                dtype=bool)
        plt.yscale('log')
        plt.fill_between(indices, pixel_counts, where=dn_window,
                         color='lightgray')
        plt.axvline(x=central_i, c='gray')
        plt.axvline(x=np.argmax(pixel_counts), c='lime', ls='--')
        plt.plot(pixel_counts)
        plt.plot(minima_i, pixel_counts[minima_i], "x")
        plt.plot(min_i, pixel_counts[min_i], "o", c='red')
        plt.plot(max_i, pixel_counts[max_i], "o", c='red')
        plt.show()

    return (dn[min_i], dn[max_i])


def mask_gap(in_path: Path, out_path: Path, keep=False):
    """Attempt to mask out pixels beyond the central DNs of the median
    based on gaps in the histogram."""

    # This approach worked well based on 'ideal' reverse-clocked data
    # or 'dark' images, but in 'real' HiRISE images of Mars, the reality
    # is that there are 'gaps' everywhere along the DN range, and this
    # approach ends up being too 'dumb'.

    to_del = isis.PathSet()

    hist_p = in_path.with_suffix('.hist')
    hist = histogram(in_path, hist_p)

    median = math.trunc(float(hist['Median']))
    # std = math.trunc(float(hist['Std Deviation']))

    high = find_gap(hist, median,
                    math.trunc(float(hist['Maximum'])),
                    findfirst=True)

    low = find_gap(hist, median,
                   math.trunc(float(hist['Minimum'])),
                   findfirst=True)

    highdist = high - median
    lowdist = median - low

    maskmax = find_gap(hist, median, (median + (2 * highdist)))
    maskmin = find_gap(hist, median, (median - (2 * lowdist)))

    util.log(isis.mask(in_path, to=out_path, minimum=maskmin,
                       maximum=maskmax).args)

    if not keep:
        hist_p.unlink()

    return


def get_unflip_thresh(hist, far, near, lowlimit):
    try:
        thresh = find_gap(hist, far, near)
    except ValueError:
        logging.info('No zero threshold found.')
        if d >= lowlimit:
            thresh = find_min_dn(hist, far, near)
            logging.info(f'Found a minimum threshold: {thresh}')
            if(thresh == far or thresh == near):
                raise ValueError(f'Minimum, {thresh}, at edge of range.')
        else:
            raise ValueError(f'Below {lowlimit}, skipping.')
    return thresh


def unflip(in_p: Path, out_p: Path, keep=False):
    """Attempt to indentify DNs whose bits have been flipped, and
    unflip them.
    """
    to_del = isis.PathSet()

    # This full suite of deltas works well for the reverse-clock area
    # and even 'dark' images, but not 'real' images.
    # deltas = (8192, 4096, 2048, 1024, 512, 256, 128, 64)

    # Restricting the number of deltas might help, but this seems
    # arbitrary.
    deltas = (8192, 4096, 2048, 1024)

    count = 0
    suffix = '.bf{}-{}{}.cub'
    this_p = to_del.add(in_p.with_suffix(suffix.format(count, 0, 0)))
    this_p.symlink_to(in_p)

    median = math.trunc(float(isis.stats_k(in_p)['Median']))

    for (sign, pm, extrema) in ((+1, 'm', 'Maximum'),
                                (-1, 'p', 'Minimum')):
        # logging.info(pm)
        for delt in deltas:
            d = sign * delt
            far = median + d
            near = median + (d / 2)

            try:
                hist_p = to_del.add(this_p.with_suffix('.hist'))
                hist = histogram(this_p, hist_p)
            except ValueError:
                # Already have this .hist, don't need to remake.
                pass

            logging.info(f'bitflip position {pm}{delt}, near: {near} '
                         f'far: {far}, extrema: {hist[extrema]}')
            if((sign > 0 and far < float(hist[extrema])) or
               (sign < 0 and far > float(hist[extrema]))):
                count += 1
                s = suffix.format(count, pm, delt)
                next_p = to_del.add(this_p.with_suffix('').with_suffix(s))
                try:
                    thresh = get_unflip_thresh(hist, far, near, d)
                except ValueError as err:
                    logging.info(err)
                    count -= 1
                    break
                subtract_over_thresh(this_p, next_p, thresh, d, keep=keep)
                this_p = next_p
            else:
                logging.info("The far value was beyond the extrema. "
                             "Didn't bother.")

    shutil.move(this_p, out_p)

    if not keep:
        to_del.remove(this_p)
        to_del.unlink()

    return
