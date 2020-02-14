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
import os
import struct
import subprocess
from datetime import datetime
from pathlib import Path
from shutil import copyfile


import kalasiris as isis
import pvl

import pyrise.util as util
from pyrise.bitflips import mask


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
    # masked = to_del.add(cube_path.with_suffix(f'.{temp_token}.mask.cub'))
    (mindn, maxdn) = mask(reverse, 'dummy', line=True, plot=False, keep=keep)

    # # run tabledump to get CSV file
    # rev_csv = to_del.add(cube_path.with_suffix(f'.{temp_token}.reverse.csv'))
    # util.log(isis.tabledump(cube_path, to=rev_csv,
    #                         name='HiRISE Calibration Image').args)
    # rev_rows = list()
    # with open(rev_csv) as csvfile:
    #     reader = csv.reader(csvfile)
    #     for row in reader:
    #         out_row = list()
    #         for elem in row:
    #             try:
    #                 if(int(elem) < 0 or
    #                    (int(elem) >= mindn and int(elem) <= maxdn)):
    #                     out_row.append(elem)
    #                 else:
    #                     # Replace with the integer value of the ISIS NULL pixel
    #                     # out_row.append(-32768')
    #                     out_row.append('0')
    #             except ValueError:
    #                 out_row.append(elem)
    #         rev_rows.append(out_row)

    # # open CSV file, read, filter, write out new CSV file
    # rev_clean_csv = to_del.add(
    #     cube_path.with_suffix(f'.{temp_token}.rev-clean.csv'))

    # with open(rev_clean_csv, 'w') as cleaned:
    #     writer = csv.writer(cleaned)
    #     writer.writerows(rev_rows)

    # # copy input to output
    # copyfile(cube_path, out_path)

    # # run csv2table
    # util.log(isis.csv2table(csv=rev_clean_csv, to=out_path,
    #                         tablename='"HiRISE Calibration Image"').args)

    filter_reverse(cube_path, out_path, 1450, 1500)

    if not keep:
        to_del.unlink()

    return


def filter_reverse(cube, out_cube, mindn, maxdn, replacement=None):
    """Creates a new file at *out_cube* (or overwrites it) based
    on filtering the contents of the reverse-clocked area.

    The contents of the reverse-clocked area (not its buffer or
    dark pixel areas) are examined.  If the DN value is in the
    allowable range, it is checked against *mindn* and *maxdn*, if
    it is between them, it is left unchanged.  If it is outside of
    the range, the value of *replacement* is written into that
    location in *out_path*.  If *replacement* is ``None`` then the
    appropriate value for ISIS NULL will be used.

    No other areas of the image are affected, and the new file at
    *out_path* should otherwise be identical to the file given by
    *img_path*.
    """
    label = pvl.load(str(cube))
    cal_start, samples, data_fmt, null, width = get_cal_table_info(label)
    lines = 20
    # width = None
    # b_order = None
    # b_signed = None
    # data_fmt = None
    # pixels_d = label['IsisCube']['Core']['Pixels']
    # if pixels_d['Type'].casefold() == 'SignedWord'.casefold():
    #     width = 4
    #     b_signed = True
    #     if replacement is None:
    #         replacement = -32768
    #     else:
    #         raise ValueError(
    #             "Don't know how to deal with Pixel Type" + pixels_d['Type'])
    if replacement is None:
        replacement = null

    # if pixels_d['ByteOrder'].casefold() == 'Lsb'.casefold():
    #     b_order = 'little'
    # else:
    #     raise ValueError("Don't know how to deal with "
    #                      "ByteOrder '{}'".format(pixels_d['ByteOrder']))

    # print(f"min: {mindn}, max: {maxdn}")
    with open(cube, mode='rb') as (f):
        with open(out_cube, mode='wb') as (of):
            of.write(f.read(cal_start))  # Copy first part of file
            # Filter the Calibration table
            for lineno in range(lines):
                for sampno in range(samples):
                    b_pixel = f.read(width)
                    dn = struct.unpack(data_fmt, b_pixel)[0]
                    # dn = int.from_bytes(b_pixel, byteorder=b_order,
                    #                     signed=b_signed)
                    print(f"dn: {dn}")
                    if 0 < dn < mindn or maxdn < dn < 16383:
                        # logging.info(f"foudn dn to replace: {dn}")
                        # print(f'repl: {replacement}
                        # of type {type(replacement)}')
                        b_pixel = struct.pack(data_fmt, replacement)
                        # of.write(replacement.to_bytes(width,
                        #                               byteorder=b_order,
                        #                               signed=b_signed))
                    # else:
                    #    of.write(b_pixel)
                    of.write(b_pixel)

            of.write(f.read())  # Copy last part of file


def get_cal_table_info(label: dict):
    start = None
    samples = None
    data_fmt = None
    size = None

    data_sizes = {'Integer': 4,
                  'Double': 8,
                  'Real': 4,
                  'Text': 1}
    data_formats = {'Integer': 'i',
                    'Double': 'd',
                    'Real': 'f'}

    null_val = {'Integer': -32768}

    for t in label.getlist('Table'):
        if t['Name'] == 'HiRISE Calibration Image':
            start = t['StartByte'] - 1
            samples = t['Field']['Size']
            data_fmt = data_formats[t['Field']['Type']]
            null = null_val[t['Field']['Type']]
            size = data_sizes[t['Field']['Type']]
            break

    return (start, samples, data_fmt, null, size)


# This function comes from
# https://github.com/USGS-Astrogeology/ale/blob/master/ale/base/data_isis.py
# commit 5454875e75e900444abf77a9b01ae2e3fd600beb (tag: 0.7.3)
# Make no changes that are not pushed upstream to this file.
# Perhaps one day we will just be able to simply import.
def read_table_data(table_label, cube):
    """
    Helper function to read all of the binary table data

    Parameters
    ----------
    table_label : PVLModule
                  The ISIS table label
    cube : file
           The ISIS cube file

    Returns
    -------
    bytes :
        The binary portion of the table data
    """
    cubehandle = open(cube, "rb")
    cubehandle.seek(table_label['StartByte'] - 1)
    return cubehandle.read(table_label['Bytes'])


# This function comes from
# https://github.com/USGS-Astrogeology/ale/blob/master/ale/base/data_isis.py
# commit 5454875e75e900444abf77a9b01ae2e3fd600beb (tag: 0.7.3)
# Make no changes that are not pushed upstream to this file.
# Perhaps one day we will just be able to simply import.
def parse_table(table_label, data):
    """
    Parse an ISIS table into a dictionary.

    Parameters
    ----------
    table_label : PVLModule
                  The ISIS table label
    data : bytes
           The binary component of the ISIS table

    Returns
    -------
    dict :
           The table as a dictionary with the keywords from the label and the
           binary data
    """
    data_sizes = {'Integer': 4,
                  'Double': 8,
                  'Real': 4,
                  'Text': 1}
    data_formats = {'Integer': 'i',
                    'Double': 'd',
                    'Real': 'f'}

    # Parse the binary data
    fields = table_label.getlist('Field')
    results = {field['Name']: [] for field in fields}
    offset = 0
    for record in range(table_label['Records']):
        for field in fields:
            if field['Type'] == 'Text':
                field_data = data[
                    offset:offset + field['Size']].decode(encoding='latin_1')
            else:
                data_format = data_formats[field['Type']] * field['Size']
                field_data = struct.unpack_from(data_format, data[offset:])
                if len(field_data) == 1:
                    field_data = field_data[0]

            results[field['Name']].append(field_data)
            offset += data_sizes[field['Type']] * field['Size']

    # Parse the keywords from the label
    results.update({key: value for key,
                    value in table_label.items() if not isinstance(
                        value, pvl._collections.PVLGroup)})

    return results


if __name__ == '__main__':
    main()
