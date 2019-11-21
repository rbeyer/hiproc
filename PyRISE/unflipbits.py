#!/usr/bin/env python
"""Unflip stuck bits in HiRISE pixels."""

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
                        required=False, default='.unflip.cub')
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
        print(out_p)

        try:
            unflip(i, out_p, keep=args.keep)
        except subprocess.CalledProcessError as err:
            print('Had an ISIS error:')
            print(' '.join(err.cmd))
            print(err.stdout)
            print(err.stderr)
            raise err
    return


def find_thresh(hist: isis.Histogram, far: float, near: float) -> int:
    if far > near:
        start = math.ceil(far)
        stop = math.floor(near)
        step = -1
        rev = True
    else:
        start = math.floor(far)
        stop = math.ceil(near)
        step = 1
        rev = False

    dn_window = range(start, stop, step)

    hist_set = set(filter(lambda d: d in dn_window,
                          (int(x.DN) for x in hist)))

    missing = sorted(set(dn_window).difference(hist_set), reverse=rev)
    # logging.info(f'{bs} missing: ' + str(missing))

    if len(missing) > 0:
        sequences = list()
        for k, g in itertools.groupby(missing,
                                      (lambda x,
                                       c=itertools.count(): next(c) - (step * x))):
            sequences.append(list(g))

        return max(sequences, key=len)[0]
    else:
        # There was no gap in DN.
        # Maybe look into the 'Pixel' values of each DN, and find a minimum?
        raise ValueError('There was no gap in the DN window from '
                         f'{start} to {stop}.')


def subtract_over_thresh(in_path: os.PathLike, out_path: os.PathLike,
                         thresh: int, delta: int, keep=False):
    '''For all pixels in in_path, if delta is positive, then pixels
       with a value greater than thresh will have delta subtracted from
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


def unflip(cube: os.PathLike, out_path: os.PathLike, keep=False):
    in_p = Path(cube)
    out_p = Path(out_path)

    to_del = isis.PathSet()

    # If there is some process which is flipping one of the HiRISE bits
    # in a 14-bit pixel the wrong way, it will inadvertently add or subtract
    # these numbers of DN to a pixel.  Naturally, it goes all the way down
    # the powers of two to 2 itself, but it is less certain that such small
    # values are erroneous.  These are good guesses, but maybe some exclusion
    # should happen based on standard deviations or something, so we don't dip
    # too far down?
    deltas = (8192, 4096, 2048, 1024, 512, 256, 128, 64)

    count = 0
    suffix = '.uf{}-{}{}.cub'
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
                util.log(isis.hist(this_p, to=hist_p).args)
                hist = isis.Histogram(hist_p)
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
