#!/usr/bin/env python
"""Convert one of Alan's MDR files (in DN) to a cube file (in I/F)."""

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


# This program is based on equations in make_IOF.pro, by Alan Delamere,
# 2019-05-22.

import argparse
import collections
import logging
import statistics
import subprocess
from pathlib import Path

import PyRISE.hirise as hirise
import PyRISE.util as util
import kalasiris as isis


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-o', '--output',  required=False, default='.iof.cub')
    parser.add_argument('mdr', metavar="MDR_file")
    parser.add_argument('-l', '--log',  required=False, default='WARNING',
                        help="The log level to show for this program, can be a "
                        "named log level or a numerical level.")
    parser.add_argument('-k', '--keep', required=False, default=False,
                        action='store_true',
                        help="Normally, the program will clean up any "
                        "intermediary files, but if this option is given, it "
                        "won't.")

    args = parser.parse_args()

    util.set_logging(args.log)

    mdr_path = Path(args.mdr)

    h2i_path = mdr_path.with_suffix('.MDR.cub')
    out_path = util.path_w_suffix(args.output, h2i_path)

    # The first thing Alan's program did was to crop the image down to only the
    # 'imaging' parts.  We're not doing that so the resultant file has a
    # geometry similar # to what comes out of ISIS hical.

    logging.info(isis.hi2isis(mdr_path, to=h2i_path).args)

    fpa_t = statistics.mean([float(isis.getkey_k(h2i_path, 'Instrument',
                                                 'FpaPositiveYTemperature')),
                             float(isis.getkey_k(h2i_path, 'Instrument',
                                                 'FpaNegativeYTemperature'))])

    cid = hirise.get_CCDID_fromfile(h2i_path)
    tdg = t_dep_gain(cid.ccdname, fpa_t)
    suncorr = solar_correction()
    sed = float(isis.getkey_k(h2i_path, 'Instrument', 'LineExposureDuration'))
    zbin = 1  # Ideally, this comes from the hical conf file GainUnitConversionBinFactor

    # The 'ziof' name is from the ISIS HiCal/GainUnitConversion.h, it is a
    # divisor in the calibration equation.
    ziof = (zbin * tdg * sed * 1e-6 * suncorr)
    eqn = f"\(F1 / {ziof})"
    try:
        logging.info(isis.fx(f1=h2i_path, to=out_path, equ=eqn).args)
    except subprocess.CalledProcessError as err:
        print(err.stderr)

    if not args.keep:
        h2i_path.unlink()


def solar_correction() -> float:
    # Ideally, get the value of the distance to the Sun for this observation.
    au = 1.498  # Placeholder from Alan
    s = 1.5 / au  # Not sure about this, comes from ISIS
    return s * s


def t_dep_gain(ccd: str, t: float) -> float:
    '''Given the type of CCD, and the FPA temperature in C,
       calculate the temperature dependent gain.'''
    # Equivalent to getTempDepGain() in ISIS HiCal/GainUnitConversion.h
    # But numbers all come from Alan.
    t_factor = 1 + (t - 18.9)
    if ccd == 'RED':
        return 157709797 * t_factor * 0.0005704 * 6.37688
    elif ccd == 'IR':
        return 56467454 * t_factor * 0.002696 * 6.99017
    elif ccd == 'BG':
        return 115121269 * t_factor * 0.00002295 * 7.00042
    else:
        raise Exception(f'{ccd} is not a valid CCD type.')
