#!/usr/bin/env python
"""Prepares an observation for HiPrecision processing.

   This program is *substantially* less capable than its
   Perl predecessor.  This is because most of the processing
   tasks were moved into ResolveJitter / HiJACK.
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

import argparse
import logging
import os
from pathlib import Path

import pvl

import hiproc.SlitherStats as sstats
import hiproc.util as util

logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, parents=[util.parent_parser()]
    )
    parser.add_argument(
        "-c",
        "--conf",
        required=False,
        default=Path(__file__).resolve().parent.parent
        / "data"
        / "HiPrecisionInit.conf",
    )
    parser.add_argument("slither_text", metavar="slither.txt-files", nargs="+")

    args = parser.parse_args()

    # Ignore args.log to always print info when run from the command line.
    util.set_logger("info", args.logfile, args.log)

    start(args.slither_text, args.conf)

    return


def start(slither_paths: list, conf_path: os.PathLike):
    conf = pvl.load(str(conf_path))

    yes_HiJACK = list()
    thresh = float(conf["HiPrecisionInit"]["Mean_Jitter_Magnitude_Threshold"])
    logger.info(f"Mean_Jitter_Magnitude_Threshold: {thresh}")
    logger.info(f"Average\tProcess \tFile Name")
    for s in slither_paths:
        (need, avediff) = needs_HiJACK(s, thresh)
        if need:
            terminus = "HiJACK  "
        else:
            terminus = "HiNoProj"

        logger.info("{:.2}\t{}\t{}".format(avediff, terminus, s))
        yes_HiJACK.append(need)

    return yes_HiJACK


def needs_HiJACK(slither_path: os.PathLike, threshold: float):
    (_, avediff, _) = sstats.Polyfit(slither_path)
    need = None
    if avediff > threshold:
        need = True
    else:
        need = False

    return need, avediff
