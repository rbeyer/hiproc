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

# Copyright 2019, Ross A. Beyer (rbeyer@seti.org)
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
# This program is based on clean_bit_flips.pro by Alan Delamere, Oct 2019.


import argparse
import itertools
import logging
import math
import os
import shutil
import subprocess
from pathlib import Path

import pvl

import PyRISE.hirise as hirise
import PyRISE.util as util
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


def find_thresh(hist: list, start, stop, step=None,
                findfirst=False) -> int:
    '''Given the range specified by start, stop and step, the DNs are
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
    '''
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

    dn_window = range(start, stop, step)

    hist_set = set(filter(lambda d: d in dn_window,
                          (int(x.DN) for x in hist)))

    missing = sorted(set(dn_window).difference(hist_set), reverse=(step < 0))
    # logging.info(f'{bs} missing: ' + str(missing))

    if len(missing) > 0:
        sequences = list()
        for k, g in itertools.groupby(missing,
                                      (lambda x,
                                       c=itertools.count(): next(c) - (step * x))):
            sequences.append(list(g))

        if findfirst:
            return sequences[0][0]
        else:
            # find biggest gap
            return max(sequences, key=len)[0]
    else:
        # There was no gap in DN.
        # Maybe look into the 'Pixel' values of each DN, and find a minimum?
        raise ValueError('There was no gap in the DN window from '
                         f'{start} to {stop}.')


def subtract_over_thresh(in_path: os.PathLike, out_path: os.PathLike,
                         thresh: int, delta: int, keep=False):
    '''For all pixels in the in_path ISIS cube, if delta is positive, then
       pixels with a value greater than thresh will have delta subtracted from
       them.  If delta is negative, then all pixels less than thresh
       will have delta added to them.
    '''

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
    '''This is just a convenience function to facilitate logging.'''
    util.log(isis.hist(in_path, to=hist_path).args)
    return isis.Histogram(hist_path)


def mask(in_path: Path, out_path: Path, keep=False):
    '''Attempt to mask out pixels beyond the central DNs of the median.'''

    to_del = isis.PathSet()

    hist_p = in_path.with_suffix('.hist')
    hist = histogram(in_path, hist_p)

    median = math.trunc(float(hist['Median']))
    # std = math.trunc(float(hist['Std Deviation']))

    high = find_thresh(hist, median,
                       math.trunc(float(hist['Maximum'])),
                       findfirst=True)

    low = find_thresh(hist, median,
                      math.trunc(float(hist['Minimum'])),
                      findfirst=True)

    highdist = high - median
    lowdist = median - low

    maskmax = find_thresh(hist, median, (median + (2 * highdist)))
    maskmin = find_thresh(hist, median, (median - (2 * lowdist)))

    util.log(isis.mask(in_path, to=out_path, minimum=maskmin,
                       maximum=maskmax).args)

    if not keep:
        hist_p.unlink()

    return


def unflip(in_p: Path, out_p: Path, keep=False):
    '''Attempt to indentify DNs whose bits have been flipped, and
       unflip them.
    '''
    to_del = isis.PathSet()

    deltas = (8192, 4096, 2048, 1024, 512, 256, 128, 64)

    count = 0
    suffix = '.bf{}-{}{}.cub'
    this_p = to_del.add(in_p.with_suffix(suffix.format(count, 0, 0)))
    this_p.symlink_to(in_p)

    median = math.trunc(float(isis.stats_k(in_p)['Median']))

    for (sign, pm, extrema) in ((+1, 'm', 'Maximum'),
                                (-1, 'p', 'Minimum')):
        logging.info(pm)
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
                    thresh = find_thresh(hist, far, near)
                except ValueError:
                    logging.info('No threshold found.')
                    count -= 1
                    break
                subtract_over_thresh(this_p, next_p, thresh, d, keep=keep)
                this_p = next_p
            else:
                logging.info("The far value was beyond the extrema. "
                             "Didn't bother.")

    if not keep:
        to_del.unlink()

    return
