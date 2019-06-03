#!/usr/bin/env python
"""Convert one of Alan's MDR files (in DN or I/F) to a cube file (in I/F)."""

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
# and on equations in the ISIS program hical.
# 2019-05-22.

import argparse
import collections
import logging
import shutil
import statistics
import subprocess
import sys
from pathlib import Path
from osgeo import gdal

import pvl

import PyRISE.hirise as hirise
import PyRISE.util as util
import kalasiris as isis


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-o', '--output',  required=False, default='.mdr.iof.cub')
    parser.add_argument('-e', '--edr',  required=True)
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

    edr_path = Path(args.edr)
    mdr_path = Path(args.mdr)

    to_del = isis.PathSet()

    h2i_path = to_del.add(edr_path.with_suffix('.hi2isis.cub'))
    out_path = util.path_w_suffix(args.output, edr_path)

    # The first thing Alan's program did was to crop the image down to only the
    # 'imaging' parts.  We're not doing that so the resultant file has a
    # geometry similar to what comes out of ISIS hical.

    # Convert the EDR to a cube file
    logging.info(isis.hi2isis(edr_path, to=h2i_path).args)

    # Convert Alan's MDR to a cube file
    mdr_cub_path = to_del.add(mdr_path.with_suffix('.alan.cub'))
    logging.info(f'Running gdal_translate {mdr_path} -of ISIS3 {mdr_cub_path}')
    gdal.Translate(str(mdr_cub_path), str(mdr_path), format='ISIS3')

    h2i_s = isis.getkey_k(h2i_path, 'Dimensions', 'Samples')
    h2i_l = isis.getkey_k(h2i_path, 'Dimensions', 'Lines')
    mdr_s = isis.getkey_k(mdr_cub_path, 'Dimensions', 'Samples')
    mdr_l = isis.getkey_k(mdr_cub_path, 'Dimensions', 'Lines')

    if h2i_s != mdr_s:
        logging.critical(f'The number of samples in {h2i_path} ({h2i_s}) '
                         f'and {mdr_cub_path} ({mdr_s}) are different. Exiting.')
        sys.exit()

    if h2i_l != mdr_l:
        logging.critical(f'The number of lines in {h2i_path} ({h2i_l}) '
                         f'and {mdr_cub_path} ({mdr_l}) are different. Exiting.')
        sys.exit()

    # Convert the EDR to the right bit type for post-HiCal Pipeline:
    h2i_16b_p = to_del.add(h2i_path.with_suffix('.16bit.cub'))
    logging.info(isis.bit2bit(h2i_path, to=h2i_16b_p, bit='16bit',
                              clip='minmax', minval=0, maxval=1.5).args)

    shutil.copyfile(h2i_16b_p, out_path)

    # Is the MDR in DN or I/F?
    maximum_pxl = float(pvl.loads(isis.stats(mdr_cub_path).stdout)['Results']['Maximum'])
    if maximum_pxl < 1.5:
        logging.info('MDR is already in I/F units.')
        mdr_16b_p = to_del.add(mdr_cub_path.with_suffix('.16bit.cub'))
        logging.info(isis.bit2bit(mdr_cub_path, to=mdr_16b_p, bit='16bit',
                                  clip='minmax', minval=0, maxval=1.5).args)
        logging.info(isis.handmos(mdr_16b_p, mosaic=out_path).args)
    else:
        logging.info('MDR is in DN units and will be converted to I/F.')

        fpa_t = statistics.mean([float(isis.getkey_k(h2i_16b_p, 'Instrument',
                                                     'FpaPositiveYTemperature')),
                                 float(isis.getkey_k(h2i_16b_p, 'Instrument',
                                                     'FpaNegativeYTemperature'))])
        print(f'fpa_t {fpa_t}')

        cid = hirise.get_CCDID_fromfile(h2i_16b_p)
        tdg = t_dep_gain(cid.ccdname, fpa_t)
        suncorr = solar_correction()
        sed = float(isis.getkey_k(h2i_16b_p, 'Instrument', 'LineExposureDuration'))
        zbin = 1  # Ideally, this comes from the hical conf file GainUnitConversionBinFactor

        # The 'ziof' name is from the ISIS HiCal/GainUnitConversion.h, it is a
        # divisor in the calibration equation.
        print(f'zbin {zbin}')
        print(f'tdg {tdg}')
        print(f'sed {sed}')
        print(f'suncorr {suncorr}')
        ziof = (zbin * tdg * sed * 1e-6 * suncorr)
        eqn = f"\(F1 / {ziof})"

        mdriof_p = to_del.add(mdr_cub_path.with_suffix('.iof.cub'))
        to_s = '{}+SignedWord+{}:{}'.format(mdriof_p, 0, 1.5)
        logging.info(isis.fx(f1=mdr_cub_path, to=to_s, equ=eqn).args)

        logging.info(isis.handmos(mdriof_p, mosaic=out_path).args)

    if not args.keep:
        to_del.unlink()


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
    # These equations are really g * (1 + t - baseT) * Q * absgainTDI
    # where these variables are from the various hical .conf files.
    #   g = FilterGainCorrection (DN/s) - Alan's numbers don't match
    #   baseT = IoverFbasetemperature (deg C)  - Alan's numbers are the same
    #   Q = QEpercentincreaseperC (1/deg C) - Alan's numbers are the same
    #   absgainTDI = AbsGain_TDI128 (? units) - Alan's numbers don't match
    zgain = dict(RED=157709797, IR=56467454, BG=115121269)
    baseT = dict(RED=18.9, IR=18.9, BG=18.9)
    QEpcntC = dict(RED=0.0005704, IR=0.002696, BG=0.00002295)
    absgainTDI = dict(RED=6.37688, IR=6.99017, BG=7.00042)

    return zgain[ccd] * (1 + (t - baseT[ccd]) * QEpcntC[ccd] * absgainTDI[ccd])
