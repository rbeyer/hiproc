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
import logging
import os
import subprocess
from pathlib import Path

import pvl

import PyRISE.hirise as hirise
import PyRISE.util as util
import kalasiris as isis


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     parents=[util.parent_parser()])
    parser.add_argument('-o', '--output', required=False, default='.unflip.cub')
    parser.add_argument('-n', '--noise', required=False, default=20)
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
            unflip(i, out_p, args.noise, keep=args.keep)
        except subprocess.CalledProcessError as err:
            print('Had an ISIS error:')
            print(' '.join(err.cmd))
            print(err.stdout)
            print(err.stderr)
            raise err
    return


def unflip(cube: os.PathLike, out_path: os.PathLike, noise_margin: int,
           keep=False):
    in_p = Path(cube)
    out_p = Path(out_path)

    to_del = isis.PathSet()

    stats = isis.stats(in_p)
    util.log(stats.args)
    median = int(pvl.loads(stats.stdout)['Results']['Median'])

    # If there is some process which is flipping one of the HiRISE bits
    # in a 14-bit pixel the wrong way, it will inadvertently add or subtract
    # these numbers of DN to a pixel.  Naturally, it goes all the way down
    # the powers of two to 2 itself, but it is less certain that such small
    # values are erroneous.  These are good guesses, but maybe some exclusion
    # should happen based on standard deviations or something, so we don't dip
    # too far down?
    steps = (-4096, -2048, -1024, -512, -256, -128, -64,
             1024, 512, 256, 128, 64)

    paths = list()
    for s in steps:
        paths.append(to_del.add(in_p.with_suffix('.fx{}.cub'.format(s))))

    in_paths = [in_p] + paths[:-1]
    out_paths = paths[:-1] + [out_p]

    eqn = "\(f1 + (f1{}={})*{})"

    for (i, o, s) in zip(in_paths, out_paths, steps):
        if s < 0:
            glt = '>'
            thresh = (median - s - noise_margin)
        else:
            glt = '<'
            thresh = (median + s + noise_margin)

        util.log(isis.fx(f1=i,
                         to=o,
                         equation=eqn.format(glt, thresh, s)).args)

    if not keep:
        to_del.unlink()

    return
