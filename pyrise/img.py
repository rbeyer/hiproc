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
    elif 'COLUMN' in label[name]:
        for t in label[name].getlist('COLUMN'):
            if 'Pixels' in t['NAME']:
                info['width'] = t['ITEM_BYTES']
    else:
        raise ValueError("Can't determine how to set the width.")

    if label['CALIBRATION_IMAGE']['SAMPLE_TYPE'].startswith('MSB'):
        info['b_order'] = 'big'

    if 'UNSIGNED' in label['CALIBRATION_IMAGE']['SAMPLE_TYPE']:
        info['b_signed'] = False

    info['lines'] = label[name]['LINES']
    info['samples'] = label[name]['LINE_SAMPLES']
    info['prefix'] = label[name]['LINE_PREFIX_BYTES']
    info['suffix'] = label[name]['LINE_SUFFIX_BYTES']

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


def filter(img_path: os.PathLike, out_path: os.PathLike,
           reverse=None, buf=None, image=None, dark=None,
           replacement=None):
    """Creates a new file at *out_path* (or overwrites it) based
    on filtering the selected contents of the image.

    The four 'areas' are indicated by the *reverse*, *buf*, *image*,
    and *dark* arguments.  If they are ``None`` then that area will
    not be altered.  Otherwise, they should be a two-tuple of ``int``
    values that define a low and a high DN value.

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
    if not all(map(isinstance, (reverse, buf, image, dark), repeat(tuple))):
        raise ValueError("No two-tuples were given for any area.")

    label = pvl.load(img_path)
    lut = label['INSTRUMENT_SETTING_PARAMETERS']['MRO:LOOKUP_CONVERSION_TABLE']

    rev_info = get_info('CALIBRATION_IMAGE', label)
    rev_info['lines'] = 20

    buf_info = get_info('LINE_PREFIX_TABLE', label)
    img_info = get_info('IMAGE', label)
    dark_info = get_info('LINE_SUFFIX_TABLE', label)

    w_set = set(i['width'] for i in (rev_info, buf_info, img_info, dark_info))
    if len(w_set) == 1:
        width = w_set[0]
    else:
        raise ValueError("The different areas have different bit-widths: "
                         f"{w_set}")

    if replacement is None:
        if width == 1:  # 8 bit
            replacement = 255
        elif width == 2:  # 16 bit
            replacement = -1
        else:
            raise ValueError("Don't know how to set NULL value for bit-width "
                             f"{width}")

    if(buf_info['offset'] == img_info['offset']
       and img_info['offset'] == dark_info['offset']):
        image_offset = img_info['offset']
    else:
        raise ValueError("The offsets for the buffer, image, and dark areas "
                         "should all be the same, but they're not.")

    # print(f'mindn: {mindn}, maxdn: {maxdn}')
    with open(img_path, mode='rb') as f:
        with open(out_path, mode='wb') as of:
            # Copy the data from beginning to the beginning of the
            # Calibration Image
            of.write(f.read(rev_info['offset']))
            if reverse is None:
                of.write(f.read(image_offset - rev_info['offset']))
            else:
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
