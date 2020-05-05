#!/usr/bin/env python
"""Process .dat files from Bob King's processing workflow.

This program is just an adaptor to get Bob's data run through
the algorithms in bitflips.py.
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

import argparse
import collections
import logging
import os
import sys
from pathlib import Path

import numpy as np

import pyrise.bitflips as bf
import pyrise.util as util

SpecialPixels = collections.namedtuple('SpecialPixels',
                                       ['Min', 'Null',
                                        'Lrs', 'Lis',
                                        'His', 'Hrs',
                                        'Max'])
# Could also have done:
# import kalasiris.specialpixels.SpecialPixels as SpecialPixels
# but it seemed a little silly to do an import for one namedtuple definition


def main():
    # try:
        parser = argparse.ArgumentParser(
            description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter,
            parents=[util.parent_parser()])
        parser.add_argument('-o', '--output',
                            required=False, default='.bitflip.dat')
        parser.add_argument('-w', '--width', required=False, default=5,
                            help="The number of medstd widths for bit-flip "
                            "cleaning.")
        parser.add_argument('-r', '--replacement', required=False, type=int,
                            help="By default, the program will replace "
                            "identified pixels with an appropriate NULL data "
                            "value, but if provided this value will be used "
                            "instead.")
        parser.add_argument('-p', '--plot', required=False,
                            action='store_true',
                            help="Displays plot for each area.")
        parser.add_argument('-n', '--dryrun', required=False,
                            action='store_true',
                            help="Does not produce a cleaned output file.")
        parser.add_argument('file', help='A .dat file to clean.')

        args = parser.parse_args()

        util.set_logging(args.log, args.logfile)

        out_p = util.path_w_suffix(args.output, args.file)

        clean(args.file, out_p, width=args.width,
              replacement=args.replacement,
              plot=args.plot, dryrun=args.dryrun, keep=args.keep)

        sys.exit(0)
    # except Exception as err:
    #     print(err, file=sys.stderr)
    #     sys.exit(1)


def clean(in_path: os.PathLike, out_path: os.PathLike, width=5,
          replacement=None, plot=False, dryrun=False, keep=False):
    """The file at *out_path* will be the result of running bit-flip
    cleaning of the file at *in-path*.

    **WARNING**: at this time, in_path is checked to conform to being
    a series of 20-line observations, and if so, then reverse-clock
    cleaning will be applied.

    If *replacement* is not specified, pixels that are identified
    as being beyond the allowable DN window (which may be defined
    differently for each of the image areas), will be set to the
    appropriate value from the LUT table, otherwise they will be
    set to this value.

    If *keep* is True, then all intermediate files will be preserved,
    otherwise, this function will clean up any intermediary files
    it creates.
    """

    in_p = Path(in_path)
    out_p = Path(out_path)

    header_count = 5
    header = np.fromfile(in_p, dtype=np.uint32, count=header_count, sep='')
    # Looks like: [   2  256   20    4 5120]
    # cols = header[1]
    # rows = header[2]

    if header[2] != 20:
        raise ValueError("The file did not contain 20-line observations, "
                         "and is not a set of reverse-clock data.")

    a = np.fromfile(in_p, dtype=np.float32, count=-1, sep='', offset=20)

    if a.shape[0] % (header[1] * header[2]) > 0:
        raise ValueError(f"The file had {a.shape} float values that were "
                         f"not evenly divisible by {header[1]} x {header[2]}.")

    obs_num = int(a.shape[0] / header[2] / header[1])
    rev_clocks = a.reshape((obs_num, header[2], header[1]))
    # print(rev_clocks)
    # print(rev_clocks.shape)

    # Special Pixels
    # ISIS has 5 special pixel 'values', not all of which are needed for these
    # calculations, but the called functions expect the structure, so we must
    # create one to use.
    #
    # For reasons that are unknown to me, there is one Special Pixel situation
    # for images that have no LUT applied, and a different one for those that
    # have.  To understand the differences, please see the specialpixels()
    # function of the LUT_Table class in the img.py file of this distro.
    #
    # This data does not (to my knowledge) have a way to figure out which
    # way to go, so we will use the unlutted arrangement:

    specialpix = SpecialPixels(Min=1, Null=65535, Lrs=0, Lis=0,
                               His=16383, Hrs=65535, Max=16382)
    # I haven't examined *all* LUT tables, but I think the only difference
    # between the unlutted and LUT table 'specialpixels' above, would be
    # to set Null=-9998.

    if replacement is not None:
        specialpix.Null = replacement

    for obs in range(rev_clocks.shape[0]):
        this_slice = np.s_[obs, :, :]
        r = rev_clocks[this_slice]

        if r[r == 0].size:
            percent = (r[r == 0].size / r.size) * 100
            desc = f"Observation {obs} has {r[r == 0].size} zeros, {percent}%."

            if percent > 40:
                print(f"{desc} Too many zeros, will not process.")
                continue

            print(desc)
            r_masked = np.ma.masked_outside(r, specialpix.Min, specialpix.Max)
            r_clean = bf.clean_array(
                r_masked, width=width, axis=1,
                plot=(f"Observation {obs} Reverse-Clock" if plot else False))
            rev_clocks[this_slice] = bf.apply_special_pixels(r_clean,
                                                             specialpix)
    if not dryrun:
        output = np.append(header, rev_clocks.flat)
        with open(out_p, mode='wb') as f:
            f.write(header.tobytes())
            f.write(rev_clocks.tobytes())

    return


if __name__ == "__main__":
    main()
