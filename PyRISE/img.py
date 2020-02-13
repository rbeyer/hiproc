#!/usr/bin/env python
"""Functions to work with HiRISE .img files.
"""

# Copyright 2020, Ross A. Beyer (rbeyer@seti.org)
#
# Reuse is permitted under the terms of the license.
# The LICENSE file is at the top level of this library.


import argparse
import csv
# import logging
from datetime import datetime
import os
import statistics
import subprocess
from pathlib import Path

import pvl
import kalasiris as isis

import pyrise.util as util
from pyrise.bitflips import mask


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     parents=[util.parent_parser()])
    parser.add_argument('-o', '--output',
                        required=False, default='.clean.IMG')
    parser.add_argument('-c', '--cube', required=False)
    parser.add_argument('img', metavar="some.img-file",
                        help='A HiRISE EDR.img file.')

    args = parser.parse_args()

    util.set_logging(args.log, args.logfile)

    out_p = util.path_w_suffix(args.output, args.img)

    if args.cube is not None:
        filter_reverse_cube(args.img, args.cube, out_p, keep=args.keep)
    else:
        filter_reverse(args.img, out_p, 500, 1500)

    # rev = get_reverse(args.img)

    # try:
    #     rev_cube = Path(args.img).with_suffix('.test-reverse.cub')
    #     to_isis(rev, rev_cube)
    # except subprocess.CalledProcessError as err:
    #         print('Had an ISIS error:')
    #         print(' '.join(err.cmd))
    #         print(err.stdout)
    #         print(err.stderr)

    return


def filter_reverse_cube(img_path: os.PathLike, cube_path: os.PathLike,
                        out_path: Path, keep=False):
    temp_token = datetime.now().strftime('ImgRevClean-%y%m%d%H%M%S')
    to_del = isis.PathSet()

    cube_p = Path(cube_path)

    higlob = to_del.add(cube_p.with_suffix(f".{temp_token}.higlob.cub"))
    util.log(isis.higlob(cube_p, to=higlob).args)

    reverse = to_del.add(cube_p.with_suffix(f".{temp_token}.reverse.cub"))
    util.log(isis.crop(higlob, to=reverse, nline=20).args)

    mindn, maxdn = mask(reverse, 'dummy', line=True, plot=False, keep=keep)
    # filter_reverse(cube_path, out_path, mindn, maxdn)
    filter_reverse(img_path, out_path, mindn, maxdn)

    if not keep:
        to_del.unlink()

    return


def get_params(label: dict) -> tuple:
    """Returns a tuple of specific parameters from the *label*
    dict-like that are used in reading and writing the binary
    contents of a .IMG file.
    """

    width = None
    if label['CALIBRATION_IMAGE']['SAMPLE_BITS'] == 16:
        width = 2
    elif label['CALIBRATION_IMAGE']['SAMPLE_BITS'] == 8:
        width = 1
    else:
        raise ValueError(
            "Don't know how to handle '{}'".format(
                label['CALIBRATION_IMAGE']['SAMPLE_BITS']) +
            "as a value of SAMPLE_BITS, should be 8 or 16.")

    b_order = 'little'
    if label['CALIBRATION_IMAGE']['SAMPLE_TYPE'].startswith('MSB'):
        b_order = 'big'

    b_signed = True
    if 'UNSIGNED' in label['CALIBRATION_IMAGE']['SAMPLE_TYPE']:
        b_signed = False

    return (
        label['INSTRUMENT_SETTING_PARAMETERS']['MRO:LOOKUP_CONVERSION_TABLE'],
        int(label['^CALIBRATION_IMAGE'].value),
        width, b_order, b_signed,
        label['CALIBRATION_IMAGE']['LINE_PREFIX_BYTES'],
        label['CALIBRATION_IMAGE']['LINE_SAMPLES'],
        label['CALIBRATION_IMAGE']['LINE_SUFFIX_BYTES'])


def get_reverse(img_path: os.PathLike) -> list:
    # FYI: comparing the list of lists with what is seen from ISIS
    # higlob with FLIP set is a little misleading.  higlob only flips
    # the 'image area' pixels, but the buffer and dark stay right where
    # they are.
    rev_list = list()

    (lut, cal_start, width, b_order, b_signed,
     prefix, samp, suffix) = get_params(pvl.load(img_path))
    samples = prefix + samp + suffix
    lines = 20

    print(f'{b_order} {b_signed}')

    with open(img_path, mode='rb') as f:
        f.seek(cal_start)
        for lineno in range(lines):
            line_list = list()
            for sampno in range(samples):
                line_list.append(unlut(lut, int.from_bytes(f.read(width),
                                                           byteorder=b_order,
                                                           signed=b_signed)))
            # The first six bytes are prefix
            # 12 buffer
            # 16 dark
            # rev_list.append(line_list[(6 + 12):-16])
            rev_list.append(line_list[6:])

    return rev_list


def filter_reverse(img_path: os.PathLike, out_path: os.PathLike,
                   mindn: int, maxdn: int, replacement=None):
    """Creates a new file at *out_path* (or overwrites it) based
    on filtering the contents of the reverse-clocked area.

    The contents of the reverse-clocked area (not its buffer or dark
    pixel areas) are examined and unlutted to 14-bit HiRISE DN values.
    If the DN value isn't a special value (that would end up being NULL,
    HIS, or LIS in the ISIS image) it is checked against *mindn* and *maxdn*,
    if it is between them, it is left unchanged.  If it is outside of the
    range, the value of *replacement* is converted to bytes and written
    into that spot of *out_path*.  If *replacement* is ``None`` then
    the appropriate value for ISIS NULL will be used.

    No other areas of the image are affected, and the new file at
    *out_path* should otherwise be identical to the file given by
    *img_path*.
    """

    (lut, cal_start, width, b_order, b_signed,
     prefix, samples, suffix) = get_params(pvl.load(img_path))
    lines = 20

    if width == 1:  # 8 bit
        if replacement is None:
            replacement = 255
    elif width == 2:  # 16 bit
        if replacement is None:
            replacement = -1

    # print(f'mindn: {mindn}, maxdn: {maxdn}')
    with open(img_path, mode='rb') as f:
        with open(out_path, mode='wb') as of:
            # Copy the data from beginning to the beginning of the
            # Calibration Image
            of.write(f.read(cal_start))
            for lineno in range(lines):
                # Copy the prefix from each line
                of.write(f.read(prefix))
                # This is the meat:
                for sampno in range(samples):
                    byte = f.read(width)
                    dn = unlut(lut, int.from_bytes(byte,
                                                   byteorder=b_order,
                                                   signed=b_signed))
                    # print(f'dn: {dn}')
                    if((0 < dn < mindn) or (maxdn < dn < 16383)):
                        of.write((replacement).to_bytes(width,
                                                        byteorder=b_order,
                                                        signed=b_signed))
                    else:
                        of.write(byte)
                # Copy the suffix from each line
                of.write(f.read(suffix))
            # Copy the rest.
            of.write(f.read())


def unlut(lut_table: list, pixel: int) -> int:
    if len(lut_table) == 1:
        # That means there's a single pair "(0, 0)" which indicates
        # that no lut has been applied, so this function should no-op.
        return pixel
    if pixel == 0:
        return 0
    return int(statistics.mean(lut_table[pixel]))


# def relut(lut_table: list, real: float) -> int:
#     for i, pair in enumerate(lut_table):
#         if real >= pair[0] and real <= pair[1]:
#             return i
#
#     return 255


def to_isis(values, out_pathlike: os.PathLike):
    out_path = Path(out_pathlike)
    csv_path = out_path.with_suffix('.csv')

    with open(csv_path, 'w') as csvfile:
        writer = csv.writer(csvfile)
        for line in values:
            writer.writerow(reversed(line))

    util.log(isis.ascii2isis(csv_path, to=out_path, lines=len(values),
                             samples=len(values[0])).args)


if __name__ == "__main__":
    main()
