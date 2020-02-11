#!/usr/bin/env python
"""Cleans anomalous pixels in the reverse-clock area of HiRISE EDRs.

These anomalous pixels result from bit-flip noise and are detailed
in bitflips.py.

This program is intended to operate on the reverse-clock area of a
HiRISE EDR, to determine the valid DN window, and correct it, by
replacing any DN in the reverse-clock area that are not within the
window with ISIS NULL values.
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
# This program is based on clean_bit_flips.pro by Alan Delamere, Oct 2019.


import argparse
import csv
import logging
import math
import subprocess
from datetime import datetime
from pathlib import Path
from shutil import copyfile

import kalasiris as isis

import PyRISE.util as util
from PyRISE.bitflips import mask


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     parents=[util.parent_parser()])
    parser.add_argument('-o', '--output',
                        required=False, default='.revclean.cub')
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
            start(Path(i), out_p, keep=args.keep)
        except subprocess.CalledProcessError as err:
            print('Had an ISIS error:')
            print(' '.join(err.cmd))
            print(err.stdout)
            print(err.stderr)
            # raise err
    return


def start(cube_path: Path, out_path: Path, keep=False):
    # This string will get placed in the filename for all of our
    # temporary files. It will (hopefully) prevent collisions with
    # existing files and also allow for easy clean-up if keep=True
    temp_token = datetime.now().strftime('RevClean-%y%m%d%H%M%S')

    to_del = isis.PathSet()

    # Make the higlob
    higlob = to_del.add(cube_path.with_suffix(f'.{temp_token}.higlob.cub'))
    util.log(isis.higlob(cube_path, to=higlob).args)

    # Get just the reverse-clock lines
    reverse = to_del.add(cube_path.with_suffix(f'.{temp_token}.reverse.cub'))
    util.log(isis.crop(higlob, to=reverse, nline=20).args)

    # run mask() to get DN range
    masked = to_del.add(cube_path.with_suffix(f'.{temp_token}.mask.cub'))
    (mindn, maxdn) = mask(reverse, masked, line=True, plot=False, keep=keep)

    # run tabledump to get CSV file
    rev_csv = to_del.add(cube_path.with_suffix(f'.{temp_token}.reverse.csv'))
    util.log(isis.tabledump(cube_path, to=rev_csv,
                            name='HiRISE Calibration Image').args)
    rev_rows = list()
    with open(rev_csv) as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            out_row = list()
            for elem in row:
                try:
                    if(int(elem) < 0 or
                       (int(elem) >= mindn and int(elem) <= maxdn)):
                        out_row.append(elem)
                    else:
                        # Replace with the integer value of the ISIS NULL pixel
                        # out_row.append(-32768')
                        out_row.append('0')
                except ValueError:
                    out_row.append(elem)
            rev_rows.append(out_row)

    # open CSV file, read, filter, write out new CSV file
    rev_clean_csv = to_del.add(
        cube_path.with_suffix(f'.{temp_token}.rev-clean.csv'))

    with open(rev_clean_csv, 'w') as cleaned:
        writer = csv.writer(cleaned)
        writer.writerows(rev_rows)

    # copy input to output
    copyfile(cube_path, out_path)

    # run csv2table
    util.log(isis.csv2table(csv=rev_clean_csv, to=out_path,
                            tablename='"HiRISE Calibration Image"').args)

    if not keep:
        to_del.unlink()

    return
