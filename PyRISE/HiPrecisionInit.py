#!/usr/bin/env python
"""Prepares an observation for HiPrecision processing.

   This program is *substantially* less capable than its
   Perl predecessor.  This is because most of the processing
   tasks were moved into ResolveJitter / HiJACK.
"""

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


# This program is based on HiPrecision version 1.20 2017/07/19
# and on the Perl HiPrecisionInit program: ($Revision: 1.27 $
#                                           $Date: 2017/08/30 # 21:29:23 $)
# by Audrie Fennema and Sarah Mattson
# which is Copyright(C) 2008 Arizona Board of Regents, under the GNU GPL.
#
# Since that suite of software is under the GPL, none of it can be directly
# incorporated in this program, since I wish to distribute this software
# under the Apache 2 license.  Elements of this software (written in an entirely
# different language) are based on that software but rewritten from scratch to
# emulate functionality.

import argparse
import logging
import os
from pathlib import Path

import pvl

import PyRISE.SlitherStats as sstats
import PyRISE.util as util


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     parents=[util.parent_parser()])
    parser.add_argument('-c', '--conf',    required=False,
                        default=Path(__file__).resolve().parent.parent /
                        'data' / 'HiPrecisionInit.conf')
    parser.add_argument('slither_text', metavar="slither.txt-files", nargs='+')

    args = parser.parse_args()

    # Ignore args.log to always print info when run from the command line.
    util.set_logging('info')

    start(args.slither_text, args.conf)

    return


def start(slither_paths: list, conf_path: os.PathLike):
    conf = pvl.load(str(conf_path))

    yes_HiJACK = list()
    thresh = float(conf['HiPrecisionInit']['Mean_Jitter_Magnitude_Threshold'])
    logging.info(f'Mean_Jitter_Magnitude_Threshold: {thresh}')
    logging.info(f'Average\tProcess \tFile Name')
    for s in slither_paths:
        (need, avediff) = needs_HiJACK(s, thresh)
        if need:
            terminus = 'HiJACK  '
        else:
            terminus = 'HiNoProj'

        logging.info('{:.2}\t{}\t{}'.format(avediff, terminus, s))
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
