#!/usr/bin/env python
"""Functions to work with HiRISE .img files.

   At the moment, if just given a .img file, it will
   perform some basic, canned, filtering of the reverse-clock
   area.

   If also given an ISIS cube, then it will assume that the
   ISIS cube is a result of running hi2isis on the given .img
   file.  It will analyze that ISIS cube's reverse-clock area
   and come up with a masking window and the resulting .img file
   will have its reverse-clock area filtered with that window.
"""

# Copyright 2020, Ross A. Beyer (rbeyer@seti.org)
#
# Reuse is permitted under the terms of the license.
# The LICENSE file is at the top level of this library.


import argparse
import csv
import logging
import os
import statistics
import subprocess

from datetime import datetime
from itertools import repeat
from pathlib import Path

import pvl
import kalasiris as isis

import pyrise.util as util
from pyrise.bitflips import mask


class decoder:

    def __init__(self, lut, b_order, b_signed, min_dn, max_dn,
                 width, replacement=None):
        self.b_order = b_order
        self.b_signed = b_signed
        self.min_dn = min_dn
        self.max_dn = max_dn
        self.width = width

        if replacement is None:
            if self.width == 1:  # 8 bit
                replacement = 255
            elif self.width == 2:  # 16 bit
                replacement = -1
            else:
                raise ValueError("Don't know how to set NULL value for "
                                 f"bit-width {width}")
        self.replacement = replacement.to_bytes(self.width,
                                                byteorder=self.b_order,
                                                signed=self.b_signed)

        # Pre-compute the lookups:
        self.lut_table = list()
        for pair in lut:
            self.lut_table.append(int(statistics.mean(pair)))

        # Ensure the first element is zero.
        if lut[0][0] == 0:
            self.lut_table[0] = 0
        else:
            ValueError(f"The first LUT pair doesn't have a zero: {lut[0]}")

    def unlut(self, pixel: int) -> int:
        if len(self.lut_table) == 1:
            # That means there's a single pair "(0, 0)" which indicates
            # that no lut has been applied, so this function should no-op.
            return pixel

        return self.lut_table[pixel]

    def process(self, byte, area):
        dn = self.unlut(int.from_bytes(byte,
                                       byteorder=self.b_order,
                                       signed=self.b_signed))

        if(area is not None and ((self.min_dn < dn < area[0])
                                 or (area[1] < dn < self.max_dn))):
            return self.replacement
        else:
            return byte


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     parents=[util.parent_parser()])
    parser.add_argument('-o', '--output',
                        required=False, default='.bitclean.IMG')
    parser.add_argument('-c', '--cube', required=False)
    parser.add_argument('img', metavar="some.img-file",
                        help='A HiRISE EDR.img file.')

    args = parser.parse_args()

    util.set_logging(args.log, args.logfile)

    out_p = util.path_w_suffix(args.output, args.img)

    if args.cube is not None:
        clean_from_cube(args.img, args.cube, out_p,
                        rev=True, img=True,
                        keep=args.keep)
    else:
        clean(args.img, out_p, (500, 1500))

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


def clean_from_cube(img_path: os.PathLike, cube_path: os.PathLike,
                    out_path: Path, rev=False, buf=False, img=False,
                    dark=False, keep=False):
    temp_token = datetime.now().strftime('ImgRevClean-%y%m%d%H%M%S')
    to_del = isis.PathSet()

    cube_p = Path(cube_path)

    higlob = to_del.add(cube_p.with_suffix(f".{temp_token}.higlob.cub"))
    util.log(isis.higlob(cube_p, to=higlob).args)

    if rev:
        logging.info('clean_from_cube Reverse-Clock DN Window')
        revcub = to_del.add(cube_p.with_suffix(f".{temp_token}.reverse.cub"))
        util.log(isis.crop(higlob, to=revcub, nline=20).args)

        rev_tup = mask(revcub, 'dummy', line=True, plot=False, keep=keep)
    else:
        rev_tup = None

    if buf:
        logging.info('clean_from_cube Buffer DN Window')
        bufcub = to_del.add(cube_p.with_suffix(f".{temp_token}.buffer.cub"))
        util.log(isis.crop(higlob, to=bufcub, samp=1, nsamp=12).args)

        buf_tup = mask(bufcub, 'dummy', line=False, plot=False, keep=keep)
    else:
        buf_tup = None

    if dark:
        logging.info('clean_from_cube Dark DN Window')
        darkcub = to_del.add(cube_p.with_suffix(f".{temp_token}.dark.cub"))
        samples = isis.getkey_k(higlob, 'Dimensions', 'Samples')
        startsamp = samples - 15
        util.log(isis.crop(higlob, to=darkcub, samp=startsamp).args)

        dark_tup = mask(darkcub, 'dummy', line=False, plot=False, keep=keep)
    else:
        dark_tup = None

    if img:
        logging.info('clean_from_cube Image Area DN Window')
        img_tup = mask(cube_p, 'dummy', line=False, plot=False, keep=keep)
    else:
        img_tup = None

    clean(img_path, out_path,
          reverse=rev_tup, ramp=None, buf=buf_tup, image=img_tup,
          dark=dark_tup, replacement=None)

    if not keep:
        to_del.unlink()

    return


def get_info(name: str, label: dict) -> dict:
    """Returns a dict of parameters from the *label*
    dict-like that are used in reading and writing the binary
    contents of a .IMG file based on *name*.
    """
    info = dict()

    info['offset'] = int(label[f'^{name}'].value)

    if 'SAMPLE_BITS' in label[name]:
        if label[name]['SAMPLE_BITS'] == 16:
            info['width'] = 2
        elif label[name]['SAMPLE_BITS'] == 8:
            info['width'] = 1
        else:
            raise ValueError(
                "Don't know how to handle '{}'".format(
                    label[name]['SAMPLE_BITS']) +
                f"as a value of SAMPLE_BITS in {name}, should be 8 or 16.")
    # elif 'COLUMN' in label[name]:
    #     for t in label[name].getlist('COLUMN'):
    #         if 'Pixels' in t['NAME']:
    #             info['width'] = t['ITEM_BYTES']
    else:
        raise ValueError("Can't determine how to set the width.")

    if label['CALIBRATION_IMAGE']['SAMPLE_TYPE'].startswith('MSB'):
        info['b_order'] = 'big'
    else:
        info['b_order'] = 'little'

    if 'UNSIGNED' in label['CALIBRATION_IMAGE']['SAMPLE_TYPE']:
        info['b_signed'] = False
    else:
        info['b_signed'] = True

    info['lines'] = label[name]['LINES']
    info['samples'] = label[name]['LINE_SAMPLES']

    if label[name]['LINE_PREFIX_BYTES'] == 18:
        info['prefix'] = 12
    else:
        raise Exception()

    if label[name]['LINE_SUFFIX_BYTES'] == 16:
        info['suffix'] = 16
    else:
        raise Exception()

    return info


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


def clean_check(reverse, ramp, buf, image, dark,
                rr_info, img_info) -> tuple:
    """Checks various values for consistency."""
    width = None

    if not any(map(isinstance, (reverse, ramp, buf, image, dark),
                   repeat(tuple))):
        raise ValueError("No two-tuples were given for any area.")

    w_set = set(i['width'] for i in (rr_info, img_info))
    if len(w_set) == 1:
        width = w_set.pop()
    else:
        raise ValueError("The different areas have different bit-widths: "
                         f"{w_set}")

    image_offset = img_info['offset']

    o_set = set(i['b_order'] for i in (rr_info, img_info))
    if len(o_set) == 1:
        b_order = o_set.pop()
    else:
        raise ValueError("The different areas have different byte orders: "
                         f"{o_set}")

    s_set = set(i['b_signed'] for i in (rr_info, img_info))
    if len(s_set) == 1:
        b_signed = s_set.pop()
    else:
        raise ValueError("The different areas have different signed bytes: "
                         f"{s_set}")

    return width, b_order, b_signed


def clean(img_path: os.PathLike, out_path: os.PathLike,
          reverse=None, ramp=None, buf=None, image=None, dark=None,
          replacement=None):
    """Creates a new file at *out_path* (or overwrites it) based
    on filtering the selected contents of the image.

    The five 'areas' are indicated by the *reverse*, *ramp*, *buf*,
    *image*, and *dark* arguments.  If they are ``None`` then that
    area will not be altered.  Otherwise, they should be a two-tuple
    of ``int`` values that define a low and a high DN value.

    For each area that is not ``None``, the pixels are examined and
    unlutted to 14-bit HiRISE DN values.  If the DN value isn't a
    special value (that would end up being NULL, HIS, or LIS in the
    ISIS image) it is checked against the first and second items
    of the two-tuple. If the value is between them, the pixel
    is left unchanged.  If it is outside of the range, the value of
    *replacement* is converted to bytes and written into that spot
    of *out_path*.  If *replacement* is ``None`` then the appropriate
    value for ISIS NULL will be used.

    No other areas of the image are affected, and the new file at
    *out_path* should otherwise be identical to the file given by
    *img_path*.
    """
    label = pvl.load(img_path)
    rr_info = get_info('CALIBRATION_IMAGE', label)
    img_info = get_info('IMAGE', label)

    width, b_order, b_signed = clean_check(reverse, ramp, buf, image, dark,
                                           rr_info, img_info)
    d = decoder(
        label['INSTRUMENT_SETTING_PARAMETERS']['MRO:LOOKUP_CONVERSION_TABLE'],
        b_order, b_signed,
        0, 16383,  # min and max of the allowable 14-bit range
        width, replacement)

    # print(f'mindn: {mindn}, maxdn: {maxdn}')
    with open(img_path, mode='rb') as f:
        with open(out_path, mode='wb') as of:
            # Copy the data from beginning to the beginning of the
            # Calibration Image
            of.write(f.read(rr_info['offset']))
            for lineno in range(20):
                of.write(f.read(6))
                for sampno in range(rr_info['prefix']
                                    + rr_info['samples']
                                    + rr_info['suffix']):
                    of.write(d.process(f.read(width), reverse))

            for lineno in range(rr_info['lines'] - 20):
                of.write(f.read(6))
                for sampno in range(rr_info['prefix']
                                    + rr_info['samples']
                                    + rr_info['suffix']):
                    of.write(d.process(f.read(width), ramp))

            for lineno in range(img_info['lines']):
                of.write(f.read(6))
                for sampno in range(img_info['prefix']):
                    of.write(d.process(f.read(width), buf))

                for sampno in range(img_info['samples']):
                    of.write(d.process(f.read(width), image))

                for sampno in range(img_info['suffix']):
                    of.write(d.process(f.read(width), dark))

            # Write the rest
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
