#!/usr/bin/env python
"""Examines slither files to determine whether HiJACK processing
is needed, or if HiNoProj is sufficient.

If the average of the absolute difference between the slither fit
and the measurements in the slither.txt file is greater than the
threshhold in the configuration file, this program will tell you
that.
"""

# Copyright 2008-2020, Arizona Board of Regents on behalf of the Lunar and
# Planetary Laboratory at the University of Arizona.
#   - Orignal Perl program.
#
# Copyright 2020, Ross A. Beyer (rbeyer@seti.org)
#   - Elements of this Python program are are based on the original Perl
#     but the logic here is rewritten from scratch to emulate functionality.
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
# This program is based on HiPrecision version 3.2.3 (2020/06/23),
# and on the Perl HiPrecisionInit program ($Revision: 1.28 $
#                                          $Date: 2020/02/15 00:36:20 $)
# by Audrie Fennema and Sarah Mattson as employees of the University of
# Arizona.

# This program is *substantially* less capable than its
# Perl predecessor.  This is because most of the processing
# tasks were moved into the Python ResolveJitter / HiJACK
# programs.

import argparse
import logging
import os
import pkg_resources

import pvl

import hiproc.SlitherStats as sstats
import hiproc.util as util

logger = logging.getLogger(__name__)


def arg_parser():
    parser = argparse.ArgumentParser(
        description=__doc__,
        parents=[util.parent_parser()],
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-c",
        "--conf",
        required=False,
        type=argparse.FileType('r'),
        default=pkg_resources.resource_stream(
            __name__,
            'data/HiPrecisionInit.conf'
        ),
        help="Path to the HiPrecisionInit config file.  Defaults to "
             "HiPrecisionInit.conf distributed with the library."
    )
    parser.add_argument(
        "slither_text",
        metavar="slither.txt-files",
        nargs="+",
        help="The text files created by HiSlither."
    )
    return parser


def main():
    args = arg_parser().parse_args()

    # Ignore args.log to always print info when run from the command line.
    util.set_logger("info", args.logfile, args.log)

    with util.main_exceptions(args.verbose):
        hijack, averages, thresh = check(
            args.slither_text, pvl.load(args.conf)
        )

    print("\n".join(message(hijack, averages, thresh, args.slither_text)))

    return


def check(slither_paths: list, conf: dict):

    yes_HiJACK = list()
    avediffs = list()
    thresh = float(conf["HiPrecisionInit"]["Mean_Jitter_Magnitude_Threshold"])
    for s in slither_paths:
        (need, avediff) = needs_HiJACK(s, thresh)
        yes_HiJACK.append(need)
        avediffs.append(avediff)

    return yes_HiJACK, avediffs, thresh


def needs_HiJACK(slither_path: os.PathLike, threshold: float):
    (_, avediff, _) = sstats.Polyfit(slither_path)
    need = None
    if avediff > threshold:
        need = True
    else:
        need = False

    return need, avediff


def message(yes_hijack: list, averages: list, thresh: float, paths: list):
    lines = list()
    lines.append(f"Mean_Jitter_Magnitude_Threshold: {thresh}")
    lines.append("Average\tProcess \tFile Name")
    for a, j, s in zip(averages, yes_hijack, paths):
        lines.append("{:.2}\t{}\t{}".format(
            a, "HiJACK  " if j else "HiNoProj", s
        ))

    if any(yes_hijack):
        lines.append("Probably should run HiJACK.")
    else:
        lines.append("Can safely run NoProj.")

    return lines
