#!/usr/bin/env python
"""Stitch together two HiRISE channel files to create a single CCD image file."""

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


# This program is based on HiStitch version 1.32 2016/08/04
# and on the Perl HiStitch program: ($Revision: 1.24 $ $Date: 2016/08/05 18:05:28 $)
# by Eric Eliason
# which is Copyright(C) 2016 Arizona Board of Regents, under the GNU GPL.
#
# Since that suite of software is under the GPL, none of it can be directly
# incorporated in this program, since I wish to distribute this software
# under the Apache 2 license.  Elements of this software (written in an entirely
# different language) are based on that software but rewritten from scratch to
# emulate functionality.

import argparse
import os
import re
import shutil
import sys
from pathlib import Path

import kalasiris as isis
import hirise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--db',          required=False, default='HiCat.db')
    parser.add_argument('-o', '--output', required=False, default='.HiStitch.cub')
    parser.add_argument('-c', '--conf',    required=False, default='HiStitch.conf')
    parser.add_argument('-k', '--keep',   required=False, default=False)
    parser.add_argument('cube0', metavar=".cub-file")
    parser.add_argument('cube1', metavar=".cub-file", required=False)

    # The original Perl needed the specific output of cubenorm from a
    # particular step in the HiCal pipeline before here.  However, rather
    # than keep track of those files, open them, and read them again, we
    # can (and did) perform the relevant check in HiCal.py, and then save
    # that in the db that we can check now.

    args = parser.parse_args()

    (pid0, pid1) = get_pids(args.cube0, args.cube1)

    outcub_path = set_outcube(args.output, pid0)

    # GetConfigurationParameters()
    conf = pvl.load(args.conf)

    # GetProductFiles()
    channel_cub = {'0': None, '1': None}
    if not args.cube1:
        channel_cub['0'] = Path(args.cube0)
    else:
        channel_cub = sort_input_cubes(args.cube0, args.cube1)

    # PrepareDBStatements()
    #   select from HiCat.EDR_Products
    db0
    db1

    (truthchannel, balanceratio) = HiStitch(channel_cub['0'], channel_cub['1'],
                                            output_path,
                                            conf['HiStitch'], db0, db1,
                                            int(pid0.ccdnumber),
                                            keep=args.keep)

    # insert ObsID, pid0.ccdname+pid0.ccdnumber, truthchannel, balanceratio
    # into HiCat.CCD_Processing_Statistics


def get_pids(cub_a: Path, cub_b=None) -> tuple:
    pid_a = hirise.get_ProdID_fromfile(cub_a)
    pid_b = None
    if cub_b:
        pid_b = hirise.get_ProdID_fromfile(cub_b)
        if pid_a != pid_b:
            raise ValueError(f'These do not appear to be channels '
                             'from the same Observation: {args.cube0}: '
                             '{pid_a} and {args.cube1}: {pid_b}')
    return (pid_a, pid_b)


def set_outcube(output: os.PathLike, pid: hirise.ProductID) -> Path:
    if args.output.startswith('.'):
        stitch_prod = '{}_{}_{}_{}{}'.format(pid.phase,
                                             pid.orbit_number,
                                             pid.latesque,
                                             pid.ccdname,
                                             pid.ccdnumber)
        return Path(output).with_name(stitch_prod + output)
    else:
        return Path(output)


def sort_input_cubes(a_cub: os.PathLike, b_cub: os.PathLike) -> dict:
    '''Figures out which one is Channel 0 and which is Channel 1.'''

    d = {'0': None, '1': None}

    a_chan = isis.getkey_k(a_cub, 'Instrument', 'ChannelNumber')
    b_chan = isis.getkey_k(b_cub, 'Instrument', 'ChannelNumber')

    try:
        d[a_chan] = Path(a_cub)
        d[b_chan] = Path(b_cub)
        for k, v in d.items():
            if v is None:
                raise RuntimeError(f'These have the same channel: {k}.')
    except KeyError as err:
        raise RuntimeError(f'This file had a channel other than 0 or 1: {err}')

    return d


def HiStitch(cub0, cub1, out_cub, conf, db0, db1, ccd_number: int, keep=False) -> tuple:
    # Allows for indexing in lists ordered by bin value.
    b = 1, 2, 4, 8, 16

    # This string will get placed in the filename for all of our
    # temporary files. It will (hopefully) prevent collisions with
    # existing files and also allow for easy clean-up if keep=True
    temp_token = datetime.now().strftime('HiStitch-%y%m%d%H%M%S')

    # ProcessingStep() - mostly sets up stuff
    max_mean = max(db0.IMAGE_MEAN, db1.IMAGE_MEAN)
    flags = set_flags(conf, db0, db1, ccd_number, b.index(db0.bin), max_mean)

    # HiStitchStep() - runs HiStitch, and inserts to db
    HiStitchStep(cub0, cub1, out_cub, db0.bin, ccd_number,
                 flags.balance, flags.equalize)

    if cub1 and balance:
        # run getkey from ?HiStitch_output?
        truthchannel = isis.getkey_k(stitch_cub, 'HiStitch', 'TruthChannel')
        balanceratio = isis.getkey_k(stitch_cub, 'HiStitch', 'BalanceRatio')

    if flags.furrow:
        furrow_file = out_cube.with_suffix(f'.{temp_token}.temp.cub')
        HiFurrow_Fix(stitch_cub, furrow_file, max_mean, keep=keep)
        furrow_file.rename(out_cube)

    if cub1 and balance:
        return(truthchannel, balanceratio)
    else:
        return(None, None)


def set_flags(conf, db0, db1, ccdnum: int, bindex: int, max_mean) -> collections.namedtuple:
    '''Set various processing flags based on various configuration
    parameters.'''
    HiStitchFlags = collections.namedtuple('HiStitchFlags',
                                           ['furrow', 'balance', 'equalize'])

    max_lis = max(db0.LOW_SATURATED_PIXELS, db1.LOW_SATURATED_PIXELS)
    max_std = max(db0.CAL_MASK_STANDARD_DEVIATION,
                  db1.CAL_MASK_STANDARD_DEVIATION)
    max_dstd = max(db0.IMAGE_DARK_STANDARD_DEVIATION,
                   db1.IMAGE_DARK_STANDARD_DEVIATION)
    max_gper = max(db0.GAP_PIXELS_PERCENT, db1.GAP_PIXELS_PERCENT)

    # We don't run FurrowCheck() here, but do it earlier in HiCal.py
    # and store the result to get now.

    # Determine if the Furrow Normalization is to occur
    furrow = False
    if(db0.zap or db1.zap
       and max_mean >= conf['HiStitch_Furrow_Normalization']
       and (db0.bin == 2 or db0.bin == 4)):
        furrow = True

    # Determine if the balance or equalize options are to occur
    balance = False
    equalize = False

    if (conf['HiStitch_Balance'] or conf['HiStitch_Equalize'] or
        (conf['HiStitch_Balance_Processing'][ccdnum] == 0 or
         (max_lis < conf['HiStitch_LIS_Pixels'] and
          # How its written below is how the original had it, but I think
          # it needs to be indexed by bindex:
          # max_std < conf['HiStitch_Balance_Bin_Mask_STD'][bindex] and
          # max_dstd < conf['HiStitch_Balance_Bin_DarkPixel_STD'][bindex]))):
          max_std < conf['HiStitch_Balance_Bin_Mask_STD'][db0.bin] and
          max_dstd < conf['HiStitch_Balance_Bin_DarkPixel_STD'][db0.bin]))):
        equalize = conf['HiStitch_Equalize'] and max_gper < conf['HiStitch_Gap_Percent']
        balance = conf['HiStitch_Balance']

    return HiStitchFlags(furrow, balance, equalize)


def HiStitchStep(in_cub0, in_cub1, out_cub, binning, ccdnum, balance, equalize):
    zpad_bin = '{:0>2}'.format(binning)
    histitch_args = {'from1': in_cub0, 'to': out_cub}
    if in_cub1:
        histitch_args['from2'] = in_cub1
    if not balance and not equalize:
        histitch_args['balance'] = False
    if equalize or balance:
        histitch_args['skip'] = conf[f'HiStitch_Bin{zpad_bin}_Area'][0]
        histitch_args['seamsize'] = conf[f'HiStitch_Bin{zpad_bin}_Area'][1]
        histitch_args['channel'] = conf['HiStitch_Control_Channel'][ccdnum]
        if balance:
            histitch_args['balance'] = 'TRUE'
        if equalize:
            histitch_args['balance'] = 'EQUALIZE'
            histitch_args['width'] = ''
            histitch_args['operator'] = ''

    return isis.histitch(**histitch_args)


def HiFurrow_Fix(in_cub, out_cub, max_mean, keep=False):
    '''Perform a normalization of the furrow region of bin 2 or 4
       HiRISE images. The input to this script is a HiRISE stitch
       product containing both channels of a CCD.'''

    binning = isis.getkey_k(in_cub, 'Instrument', 'Summing')
    lines = isis.getkey_k(in_cub, 'Dimensions', 'Lines')
    samps = isis.getkey_k(in_cub, 'Dimensions', 'Samples')

    if binning != '2' and binning != '4':
        raise Exception('HiFurrow_Fix only supports correction for '
                        'bin 2 or 4 data.')
    if binning == '2' and samps != '1024':
        raise Exception('Improper number of samples, {} for a stitch '
                        'product with bin 2.'.format(samps))

    # This string will get placed in the filename for all of our
    # temporary files. It will (hopefully) prevent collisions with
    # existing files and also allow for easy clean-up if keep=True
    temp_token = datetime.now().strftime('HFF-%y%m%d%H%M%S')
    to_del = isis.PathSet()

    # For bin2 and bin4 imaging, specify width of furrow based on
    # image average DN range
    range_low['2'] = (512, 513)  # 2 pixel furrow width
    range_mid['2'] = (511, 514)  # 4 pixel furrow width
    range_hgh['2'] = (511, 514)  # 4 pixel furrow width
    range_max['2'] = (510, 515)  # 6 pixel furrow width

    # Original code had low/mid/hgh for bin2 and bin4, but they
    # were hard-coded to be identical.
    dn_range_low = 9000
    dn_range_mid = 10000
    dn_range_hgh = 12000

    range_low['4'] = (256, 257)  # 2 pixel furrow width
    range_mid['4'] = (255, 258)  # 4 pixel furrow width
    range_hgh['4'] = (255, 258)  # 4 pixel furrow width
    range_max['4'] = (254, 259)  # 6 pixel furrow width

    if max_avg > dn_range_hgh:
        dn_range = range_max[binning]
    elif max_avg > dn_range_mid:
        dn_range = range_hgh[binning]
    elif max_avg > dn_range_low:
        dn_range = range_mid[binning]
    else:
        dn_range = range_low[binning]

    lpf_samp = int((dn_range[1] - dn_range[0] + 1) / 2) * 4 + 1
    lpf_line = int(lpf_samp / 2) * 20 + 1

    # Create a mask file
    # DN=1 for non-furrow area
    # DN=0 for furrow area
    fx_cub = to_del.add(in_cub.with_suffix(f'.{temp_token}.fx.cub'))
    isis.fx(to=fx_cub, mode='OUTPUTONLY', lines=lines, samples=samps,
            equation=f'(1*(sample<{dnrange[0]})+ 1*(sample>{dn_range[1]}) + 0)')

    # Create a file where the furrow area is set to null
    mask1_cub = to_del.add(in_cub.with_suffix(f'.{temp_token}.mask1.cub'))
    isis.mask(in_cub, mask=fx_cub, to=mask1_cub, min_=1, max_=1,
              preserve='INSIDE', spixels='NULL')

    # Lowpass filter to fill in the null pixel area
    lpf_cub = to_del.add(in_cub.with_suffix(f'.{temp_token}.lpf.cub'))
    isis.lowpass(mask1_cub, to=lpf_cub, sample=lpf_samp, line=lpf_line,
                 null=True, hrs=False, his=False, lis=False)

    # Create a file where non-furrow columns are set to null
    mask2_cub = to_del.add(in_cub.with_suffix(f'.{temp_token}.mask2.cub'))
    isis.mask(in_cub, mask=fx_cub, to=mask2_cub, min_=0, max_=0,
              preserve='INSIDE', spixels='NULL')

    # Highpass filter the furrow region
    hpf_cub = to_del.add(in_cub.with_suffix(f'.{temp_token}.hpf.cub'))
    isis.highpass(mask2_cub, to=hpf_cub, sample=1, line=lpf_line)

    # Add lowpass and highpass together to achieve desired result
    alg_cub = to_del.add(in_cub.with_suffix(f'.{temp_token}.alg.cub'))
    isis.algebra(from_=lpf_cub, from2=hpf_cub, to=alg_cub,
                 operator='ADD', A=1.0, B=1.0)

    # copy the input file to the output file then mosaic the
    # furrow area as needed.
    shutil.copyfile(in_cub, out_cube)
    isis.handmos(alg_cub, mosaic=out_cube, outsample=1, outline=1, outband=1,
                 insample=1, inline=1, inband=1, create='NO')

    if not keep:
        to_del.unlink()

    return


if __name__ == "__main__":
    main()
