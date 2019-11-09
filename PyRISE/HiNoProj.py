#!/usr/bin/env python
"""Resamples images from individual CCDs to ideal camera geometry and
mosaicks them into a single cube."""

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
# and on the Perl HiNoProj program: ($Revision: 1.9 $ $Date: 2017/08/30 21:26:39 $)
# by Audrie Fennema
# which is Copyright(C) 2008 Arizona Board of Regents, under the GNU GPL.
#
# Since that suite of software is under the GPL, none of it can be directly
# incorporated in this program, since I wish to distribute this software
# under the Apache 2 license.  Elements of this software (written in an entirely
# different language) are based on that software but rewritten from scratch to
# emulate functionality.

import argparse
import itertools
import logging
import os
import re
import shutil
import sys
from datetime import datetime
from pathlib import Path

import pvl

import kalasiris as isis
import PyRISE.hirise as hirise
import PyRISE.util as util
import PyRISE.HiColorInit as hci
import PyRISE.HiColorNorm as hcn


class Cube(hci.HiColorCube):

    def __init__(self, pathlike):
        super().__init__(pathlike)
        self.next_path = None
        self.line_offset = 0
        self.samp_offset = 0


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     parents=[util.parent_parser()])
    parser.add_argument('-o', '--output', required=False, default='_RED.NOPROJ.cub')
    parser.add_argument('-c', '--conf',   required=False,
                        default=Path(__file__).resolve().parent.parent /
                        'resources' / 'HiNoProj.conf')
    parser.add_argument('-b', '--base_ccd_number', required=False, default=5)
    parser.add_argument('cubes', metavar="balance.cub-files", nargs='+')

    args = parser.parse_args()

    util.set_logging(args.log)

    cubes = list(map(Cube, args.cubes))
    cubes.sort()

    if not all(c.ccdname == 'RED' for c in cubes):
        logging.critical('Not all of the input files are RED CCD files.')
        sys.exit()

    sequences = list()
    for k, g in itertools.groupby((int(c.ccdnumber) for c in cubes),
                                  lambda x, c=itertools.count(): next(c) - x):
        sequences.append(list(g))

    if len(sequences) != 1:
        logging.critical('The given cubes are not a single run of sequential '
                         'HiRISE CCDs, instead there are '
                         '{} groups with these CCD numbers: {}.'.format(len(sequences),
                                                                        sequences))
        sys.exit()

    base_ccd = list(filter(lambda x: x.ccdnumber == str(args.base_ccd_number),
                           cubes))
    if len(base_ccd) != 1:
        logging.critical(f'The base ccd, number {args.base_ccd_number}, is not '
                         'one of the given cubes.')
        sys.exit()

    conf = pvl.load(str(args.conf))
    conf_check(conf)

    outcub_path = hcn.set_outpath(args.output, cubes)

    import subprocess
    try:
        HiNoProj(cubes, base_ccd[0], outcub_path, conf['HiNoProj'], keep=args.keep)
    except subprocess.CalledProcessError as err:
        print('Had an ISIS error:')
        print(err.cmd)
        print(err.stdout)
        print(err.stderr)
        raise err
    return


def conf_check(conf: dict) -> None:
    '''Various checks on parameters in the configuration.'''

    t = 'Shape'
    util.conf_check_strings(t, ('ELLIPSOID', 'SYSTEM', 'USER'), conf['HiNoProj'][t])

    return


def is_polar(cubes, pole_tolerance: float, temp_token: str) -> bool:

    abs_lats = list()
    for c in cubes:
        temp_p = c.path.with_suffix(f'.{temp_token}.spice.cub')
        shutil.copyfile(c.path, temp_p)

        util.log(isis.spiceinit(temp_p).args)

        cam_cp = isis.camrange(temp_p)
        util.log(cam_cp.args)

        cpvl = pvl.loads(cam_cp.stdout)

        abs_lats.append(abs(float(cpvl['UniversalGroundRange']['MinimumLatitude'])))
        abs_lats.append(abs(float(cpvl['UniversalGroundRange']['MaximumLatitude'])))

        temp_p.unlink()

    if max(abs_lats) > float(pole_tolerance):
        return True
    else:
        return False


def handmos_side(cubes, base_cube, out_p: os.PathLike, left=True):
    '''Runs handmos to add cubes which are to one side or the other
       of the *base_cube*.'''
    # handmos left side
    side = 1
    priority = 'beneath'
    if left:
        side = -1
        priority = 'ontop'
    ssm = 1
    slm = 1
    for c in cubes[(cubes.index(base_cube) + side)::side]:
        slm -= side * round(c.line_offset)
        ssm -= side * round(c.samp_offset)
        util.log(isis.handmos(c.next_path, mosaic=out_p, priority=priority,
                              outline=slm, outsample=ssm, outband=1).args)
    return


def copy_and_spice(inpath: os.PathLike, outpath: os.PathLike,
                   conf: dict, polar=False):
    shutil.copyfile(inpath, outpath)

    spiceinit_args = {'from': outpath, 'shape': conf['Shape'],
                      'cksmithed': True, 'spksmithed': True,
                      'spkrecon': True, 'spkpredicted': False}

    if conf['Shape'] == 'USER':
        if polar:
            spiceinit_args['model'] = conf['Polar_Shape_Model_Path']
        else:
            spiceinit_args['model'] = conf['Shape_Model_Path']

    util.log(isis.spiceinit(**spiceinit_args).args)

    return


def get_offsets(cube: os.PathLike, match: os.PathLike, flat: os.PathLike) -> tuple:
    util.log(isis.hijitreg(cube, match=match, flat=flat).args)

    with open(flat, 'r') as f:
        flat_text = f.read()

        match = re.search(r'#\s+Average Line Offset:\s+(\S+)', flat_text)
        avg_line_offset = float(match.group(1))

        match = re.search(r'#\s+Average Sample Offset:\s+(\S+)', flat_text)
        avg_samp_offset = float(match.group(1))

    return (avg_line_offset, avg_samp_offset)


def add_offsets(cubes: list, base_ccdnumber: int, temp_token: str, keep=False) -> tuple:
    flats = isis.PathSet()
    for i, c in enumerate(cubes[:-1]):
        pair = '{}-{}'.format(cubes[i].get_ccd(), cubes[i + 1].get_ccd())
        flat_p = flats.add(c.path.with_suffix(f'.{temp_token}.{pair}.flat.tab'))

        (avg_line_offset, avg_samp_offset) = get_offsets(cubes[i].next_path,
                                                         cubes[i + 1].next_path,
                                                         flat_p)

        j = i
        if int(c.ccdnumber) >= base_ccdnumber:
            j = i + 1

        cubes[j].line_offset = avg_line_offset
        cubes[j].samp_offset = avg_samp_offset

    if not keep:
        flats.unlink()
        return (cubes, None)
    else:
        return (cubes, flats)


def fix_labels(cubes: list, path: os.PathLike, matched_cube: str, prodid: str) -> None:
    util.log(isis.editlab(path, option='modkey', grpname='Archive',
                          keyword='ProductId',
                          value=prodid).args)
    util.log(isis.editlab(path, option='modkey', grpname='Instrument',
                          keyword='MatchedCube',
                          value=str(matched_cube)).args)

    # Fix ck kernel in InstrumentPointing in RED label
    # This doesn't seem to be needed, maybe it was HiROC-specific.

    #  Add SourceProductIds to Archive group in label
    logging.info('Original Perl just assumes that both channels are included '
                 'in the balance cube.')
    source_ids = list()
    for c in cubes:
        source_ids.append(isis.getkey_k(c.path, 'Instrument', 'StitchedProductIds'))

    util.log(isis.editlab(path, option='ADDKEY', grpname='Archive',
                          keyword='SourceProductId',
                          value='({})'.format(', '.join(source_ids))).args)
    return


def HiNoProj(cubes: list, base_cube, outcub_path: os.PathLike, conf: dict, keep=False):

    temp_token = datetime.now().strftime('HiNoProj-%y%m%d%H%M%S')
    to_del = isis.PathSet()
    out_p = Path(outcub_path)

    polar = False
    if conf['Shape'] == 'USER':
        polar = is_polar(cubes, conf['Pole_Tolerance'], temp_token)

    for c in cubes:

        temp_p = to_del.add(c.path.with_suffix(f'.{temp_token}.spiced.cub'))
        copy_and_spice(c.path, temp_p, conf, polar)

        util.log(isis.spicefit(temp_p).args)

        c.next_path = to_del.add(c.path.with_suffix(f'.{temp_token}.noproj.cub'))
        c.path = temp_p

    for c in cubes:
        util.log(isis.noproj(c.path, match=base_cube.path, to=c.next_path,
                             source='frommatch').args)

    # Run hijitreg on adjacent noproj'ed ccds to get average line/sample offset
    (cubes, _) = add_offsets(cubes, int(base_cube.ccdnumber), temp_token, keep=keep)

    # Mosaic noproj'ed ccds using average line/sample offset
    shutil.copyfile(base_cube.next_path, out_p)
    logging.info('Original Perl hard codes this file copy from RED5, even if '
                 'another cube is selected as the base_ccd.')

    handmos_side(cubes, base_cube, out_p, left=True)
    handmos_side(cubes, base_cube, out_p, left=False)

    util.log(isis.editlab(out_p, option='addkey', grpname='Instrument',
                          keyword='ImageJitterCorrected', value=0).args)
    fix_labels(cubes, out_p, base_cube,
               '{}_{}'.format(str(cubes[0].get_obsid()),
                              cubes[0].ccdname))

    if not keep:
        to_del.unlink()

    return
