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
# on the Perl HiPrecisionInit program: (Revision: 1.27 $ $Date: 2017/08/30)
# by Audrie Fennema and Sarah Mattson
# which is Copyright(C) 2008 Arizona Board of Regents, under the GNU GPL.
# And on the Perl ResolveJitter program: (Revision: 1.24 $ $Date: 2017/11/08)
# by Sarah Mattson and Audrie Fennema
# And on the Perl HiJACK program: (Revision: 1.22 $ $Date: 2018/09/27)
# by Audrie Fennema and Sarah Mattson
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
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import pvl

import kalasiris as isis
import PyRISE.hirise as hirise
import PyRISE.util as util
import PyRISE.HiColorInit as hci
import PyRISE.HiColorNorm as hcn
import PyRISE.HiJitReg as hjr
import PyRISE.HiNoProj as hnp


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     parents=[util.parent_parser()])
    parser.add_argument('-o', '--outdir', required=False, default='./HiJACK')
    parser.add_argument('-c', '--conf_dir',   required=False,
                        default=Path(__file__).resolve().parent.parent /
                        'resources', help='Directory where ResolveJitter.conf '
                        'and HiJACK.conf can be found.')
    parser.add_argument('-b', '--base_ccd_number', required=False, default=5)
    parser.add_argument('cubes', metavar="balance.cub-files", nargs='+')

    args = parser.parse_args()

    util.set_logging(args.log)

    cubes = list(map(hnp.Cube, args.cubes))
    cubes.sort()

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

    # Need both red4 and red5
    red45 = list(filter(lambda x: x == 4 or x == 5,
                        (int(c.ccdnumber) for c in cubes)))
    if len(red45) != 2:
        logging.critical('The RED3 and RED5 cubes are not in the given'
                         'list of files: {}'.format(*cubes))
        sys.exit()

    import subprocess
    try:
        HiJACK(cubes, base_ccd[0], args.outdir, args.conf_dir, keep=args.keep)

    except subprocess.CalledProcessError as err:
        print('Had an ISIS error:')
        print(err.cmd)
        print(err.stdout)
        print(err.stderr)
        raise err
    return


def match_red(cubes: list, base_cube, flat_path, elargement_ratio=1.0006):
    '''Prepare cubs for hijitreg by matching size and binning to center red
       ccds.
    '''
    # Where does the 1.0006 value for OPTICAL_ENLARGEMENT_RATIO come from?
    # Not sure.

    red4 = list(filter(lambda x: int(x.ccdnumber) == 4, cubes))[0]
    red5 = list(filter(lambda x: int(x.ccdnumber) == 5, cubes))[0]

    if red4.bin != red5.bin:
        raise ValueError('RED4 and RED5 binning not equal.')

    for c in cubes:
        bin_ratio = c.bin / base_cube.bin

        # Differences in field of view between color and red are corrected for by
        # multiplying bin ratio by pre-determined constant for BG and
        # dividing bin ratio by that constant for IR

        mag = bin_ratio
        if c.ccdname == 'BG':
            mag = bin_ratio * elargement_ratio
        elif c.ccdname == 'IR':
            mag = bin_ratio / elargement_ratio

        # scale cubes as needed
        if mag > 1:
            util.log(isis.enlarge(c.path, to=c.next_path, sscale=mag,
                                  lscale=bin_ratio, interp='BILINEAR').args)
        elif mag < 1:
            util.log(isis.reduce(c.path, to=c.next_path, sscale=mag,
                                 lscale=bin_ratio, mode='SCALE').args)
        else:
            shutil.copy(c.path, c.next_path)

        if c.ccdname != 'RED':
            offset = (200 * (c.bin - base_cube.bin) + c.tdi - base_cube.tdi) / base_cube.bin
            mos_path = c.next_path.with_suffix('.mosaiced.cub')
            util.log(isis.handmos(c.next_path, mosaic=mos_path, create='Y',
                                  nlines=base_cube.lines,
                                  nsamples=base_cube.samps,
                                  nbands=1, outline=offset,
                                  outsamp=1, outband=1).args)
            shutil.move(mos_path, c.next_path)

        if bin_ratio != 1:
            util.log(isis.editlab(c.next_path, options='MODKEY',
                                  grpname='INSTRUMENT', KEYWORD='SUMMING',
                                  value=base_cube.bin).args)

    # This section creates flat.tabs for RED4-RED5 pair only
    rows = int(red5.lines / 50)
    util.log(isis.hijitreg(red4.next_path, match=red5.next_path, row=rows,
                           flat=flat_path).args)

    return


def find_common(ccds: set) -> str:
    overlap_red4 = frozenset(('RED3', 'RED4', 'RED5', 'IR10', 'BG12'))
    overlap_red5 = frozenset(('RED4', 'RED5', 'RED6', 'IR11', 'BG13'))

    if ccds.issubset(overlap_red5):
        return 'RED5'
    elif ccds.issubset(overlap_red4):
        return 'RED4'
    else:
        raise ValueError('Cannot determine common ccd.')


def ccd_set_generator(conf: dict):
    set_names = ('CCDs', 'Alt_1_CCDs', 'Alt_2_CCDs', 'Alt_3_CCDs')
    for n in set_names:
        yield frozenset(conf[n].split(','))


def determine_resjit_set(cubes: list, set_gen) -> tuple:
    cube_d = dict()
    for c in cubes:
        cube_d[c.get_ccd()] = c
    cube_s = frozenset(cube_d.keys())

    for s in set_gen:
        if s.issubset(cube_s):
            resjit_cubes = list()
            for x in s:
                resjit_cubes.append(cube_d[x])
            return (resjit_cubes, cube_d[find_common(s)])

    raise ValueError('No suitable ccd sets found.')


def analyze_flat(cube: hjr.JitterCube, fraction: float) -> int:
    # Need to be concerned about the derived CanSlither value?
    ret = hjr.Analyze_Flat(cube, 0, (fraction * 2))

    # Need to override some of the checking in hjr.Analyze_Flat()
    if ret == -1 and cube['SuspectCount'] <= 3:
        ret = 1

    return ret


def make_flats(cubes, common_cube, conf, temp_token, keep=False):
    # If the flat files already exist, don't remake them.
    # "$OBS_ID"."_".$ccd."-".$common.".flat.tab"
    jitter_cubes = list()
    n_row = common_cube.lines / conf['AutoRegistration']['ControlNet']['Control_Lines']

    for c in cubes:
        if c == common_cube:
            continue
        jitter_cubes.append(hjr.JitterCube(c.next_path,
                                           matchccd=common_cube.get_ccd(),
                                           config=conf))

    successful_flats = list()

    if not all(x.flattab_path.exists() for x in jitter_cubes):
        confauto = conf['AutoRegistration']
        for c in jitter_cubes:
            params = {'ROWS': n_row}
            params['TOLERANCE'] = confauto['Algorithm']['Tolerance']
            min_fraction_good = confauto['AnaylyzeFlat']['Minimum_Good']
            if c.ccdname == 'RED':
                redcolor = 'Red'
            else:
                redcolor = 'Color'

            params['COLS'] = confauto['ControlNet']['Control_Cols_' + redcolor]
            params['PATTERN_SAMPLES'] = confauto['PatternChip' + redcolor]['Samples']
            params['PATTERN_LINES'] = confauto['PatternChip' + redcolor]['Lines']
            params['SEARCH_SAMPLES'] = confauto['SearchChip' + redcolor]['Samples']
            params['SEARCH_LINES'] = confauto['SearchChip' + redcolor]['Lines']

            if c.ccdname != 'RED':
                channels = isis.getkey_k(c.path, 'Instrument', 'StitchedProductIds')
                if len(channels) < 2:
                    logging.info(f'Increasing columns because {c.path} is '
                                 'missing a channel.')
                    params['COLS'] += 1
                    min_fraction_good *= 0.5
            logging.info('The minimum allowable Fraction Good '
                         f'Matches = {min_fraction_good}')
            step = 0
            while step <= confauto['Algorithm']['Steps']:
                logging.info(f'Step {step} begin')

                hjr.run_HiJitReg(common_cube.next_path, c, params, 'ResolveJitter',
                                 'no cnetbin2pvl?')

                ret = analyze_flat(c, min_fraction_good)

                if ret == 1:
                    successful_flats.append(c.flattab_path)
                    break
                else:
                    step += 1
                    c.regdef_path.unlink()
                    c.flattab_path.unlink()
                    c.cnet_path.unlink()
                    params['TOLERANCE'] -= (confauto['Algorithm']['Increment'] * step)
            else:
                raise RuntimeError('Flat file for {c} is not within tolerances.')

    else:
        successful_flats = list(x.flattab_path for x in jitter_cubes)

    return successful_flats


def ResolveJitter(cubes: list, conf: dict, jitter_path: os.PathLike,
                  temp_token: str, keep=False):

    set_gen = ccd_set_generator(conf['ResolveJitter'])

    while True:
        (resjit_cubes, common_cube) = determine_resjit_set(cubes,
                                                           set_gen)
        flats = make_flats(resjit_cubes, common_cube, conf,
                           temp_token, keep=keep)
        if len(flats) + 1 == len(resjit_cubes):
            break

    # keep_regdefs() writes 'KEEP' into the files status of the regdef files.
    # Not sure that this is required, so skipping.

    # run resolveJitter3HiJACK.m or resolveJitter4HiJACK.cc
    rj_args = ['resolveJitter',
               str(jitter_path.parent),
               str(common_cube.get_obsid()),
               conf['AutoRegistration']['ControlNet']['Control_Lines']]

    for f in flats:
        rj_args.append(f)
        rj_args.append('-1')

    subprocess.run(rj_args, check=True)

    # Just double-check that the file we expect was created:
    if not jitter_path.exists():
        raise FileNotFoundError(jitter_path)

    return


def mosaic_dejittered(cubes: list, out_p: os.PathLike, prodid: str):
    if len(cubes) != 2:
        raise ValueError('Input cube list did not have 2 cubes.')

    hnp.handmos_side(cubes, cubes[1], out_p, left=True)
    hnp.fix_labels(cubes, out_p, str(cubes[1]), prodid)


def flat_reader(flat_path: os.PathLike):
    dialect = csv.Dialect
    dialect.delimiter = ' '
    dialect.skipinitialspace = True
    dialect.quoting = csv.QUOTE_NONE
    dialect.lineterminator = '\n'

    with open(_flat, 'r') as f:
        flat_data = f.read()

    fromtime = list()
    fromsamp = list()
    fromline = list()
    regsamp = list()
    regline = list()
    reader = csv.DictReader(itertools.filterfalse(lambda x:
                                                  x.startswith('#') or
                                                  x.isspace() or
                                                  len(x) == 0,
                                                  flat_data.splitlines()),
                            dialect=dialect)
    for row in reader:
        fromtime.append(float(row['FromTime']))
        fromsamp.append(int(row['FromSamp']))
        fromline.append(int(row['FromLine']))
        regsamp.append(float(row['RegSamp']))
        regline.append(float(row['RegLine']))

    return (np.array(fromtime),
            np.array(fromsamp),
            np.array(fromline),
            np.array(regsamp),
            np.array(regline))


def plot_flats(pre_flat, dejit_flat):
    (pre_time, pre_samp, pre_line, pre_regs, pre_regl) = flat_reader(pre_flat)
    (dej_time, dej_samp, deg_line, deg_regs, deg_regl) = flat_reader(dejit_flat)

    plt.ioff()
    fig, (ax0, ax1) = plt.subplots(2, 1, sharex=True)
    fig.suptitle(f'Jitter Corrected Offsets')

    ax0.set_ylabel('Sample Offsets (pixels)')
    ax0.scatter(pre_time, (pre_samp - pre_regs), label=pre_flat)
    ax0.scatter(dej_time, (dej_samp - dej_regs), label=dejit_flat)
    ax0.legend()

    ax1.set_ylabel('Line Offsets (pixels)')
    ax1.set_xlabel('Image Time (s)')
    ax1.scatter(pre_time, (pre_line - pre_regl), label=pre_flat)
    ax1.scatter(dej_time, (dej_line - dej_regl), label=dejit_flat)

    plt.show()

    return


def HiJACK(cubes: list, base_cube, outdir: os.PathLike,
           conf_dir: os.PathLike, keep=False):

    temp_token = datetime.now().strftime('HiJACK-%y%m%d%H%M%S')
    outdir_p = Path(outdir)
    outdir_p.mkdir(exist_ok=True)

    to_del = isis.PathSet()

    # Pre-HiJACK
    for c in cubes:
        c.next_path = to_del.add((outdir_p / c.path.name).with_suffix(f'.{temp_token}.prehijack.cub'))

    flat_path = (outdir_p /
                 '{}.{}.{}'.format(cubes[0].get_obsid(),
                                   temp_token,
                                   'RED5-RED4_prehijack.flat.tab'))
    match_red(cubes, base_cube, flat_path, keep=keep)

    jitter_path = outdir_p / (str(cubes[0].get_obsid()) + '_jitter.txt')
    if not jitter_path.exists():
        resolvejitter_conf = Path(conf_dir) / 'ResolveJitter.conf'
        ResolveJitter(cubes, pvl.load(str(resolvejitter_conf)),
                      jitter_path, temp_token, keep=keep)

    # SmearStats just updates the db, not needed.

    # HiJACK proper:
    hijack_conf = pvl.load(str(Path(conf_dir) / 'HiJACK.conf'))
    hnp.conf_check(hijack_conf)

    polar = False
    if hijack_conf['Shape'] == 'USER':
        polar = hnp.is_polar(cubes, hijack_conf['Pole_Tolerance'], temp_token)

    for c in cubes:
        hnp.copy_and_spice(c.next_path,
                           to_del.add(c.path.with_suffix(f'.{temp_token}.spiced.cub')),
                           hijack_conf, polar)
    # Run hijitter
    inlist = list()
    outlist = list()
    for c in cubes:
        c.next_path = (outdir_p / c.path.name).with_suffix('').with_suffix('.dejittered.cub')
        inlist.append(c.path)
        outlist.append(c.next_path)

    jitterck_p = outdir_p / (str(hirise.ObservationID(cubes[0])) + '.jittery.bc')

    hijitregdef_p = Path(os.environ['ISIS3DATA']) / 'mro/calibration/hijitreg.p1745.s3070.def'

    with isis.fromlist.temp(inlist) as inlist_f:
        with isis.fromlist.temp(outlist) as outlist_f:
            util.log(isis.hijitter(fromlist=inlist_f,
                                   jitter=jitter_path,
                                   regdef=hijitregdef_p,
                                   tolist=outlist_f,
                                   jitterck=jitterck_p).args)
    # Mosaic dejittered cubs
    # Remember hijitter makes all the individual cubes the size of the entire image.
    #  with the image data in the appropriate space for the ccd.

    out_root = outdir_p / str(hirise.ObservationID(cubes[0]))

    # All REDs
    red_p = out_root.with_name(out_root.name + '_RED.NOPROJ.cub')
    shutil.copyfile(base_cube.next_path, red_p)
    logging.info('Original Perl hard codes this file copy from RED5, even if '
                 'another cube is selected as the base_ccd.')
    (red_cubes,
     red_flat_files) = hnp.add_offsets(filter(lambda x: x.ccdname == 'RED', cubes),
                                       5, temp_token, keep=True)
    hnp.handmos_side(red_cubes, base_cube, red_p, left=True)
    hnp.handmos_side(red_cubes, base_cube, red_p, left=False)
    hnp.fix_labels(red_cubes, red_p, base_cube,
                   '{}_RED'.format(str(red_cubes[0].get_obsid())))

    # Center RED for color
    center_red_p = out_root.with_name(out_root.name + '_RED4-5.NOPROJ.cub')
    center_red_cubes = list(filter(lambda x: x.ccdnumber == '4' or x.ccdnumber == '5',
                                   red_cubes))
    shutil.copyfile(center_red_cubes[1].next_path, center_red_p)
    mosaic_dejittered(center_red_cubes, center_red_p,
                      '{}_RED4-5'.format(str(center_red_cubes[0].get_obsid())))

    # for IR and BG we handmos, because I'm not sure if the dejittered cubes
    # have the same line and sample numbers as the base_cube.  Need to
    # experiment.  If not, may also need to do additional label fixing?
    # IR next
    ir_p = out_root.with_name(out_root.name + '_IR.NOPROJ.cub')
    (ir_cubes, _) = hnp.add_offsets(filter(lambda x: x.ccdname == 'IR', cubes),
                                    5, temp_token, keep=keep)
    util.log(isis.handmos(ir_cubes[1].next_path, mosaic=ir_p,
                          create='yes', nlines=base_cube.lines,
                          nsamples=base_cube.samps, nbands=1).args)
    mosaic_dejittered(ir_cubes, ir_p,
                      '{}_IR'.format(str(ir_cubes[0].get_obsid())))

    # BG next
    bg_p = out_root.with_name(out_root.name + '_BG.NOPROJ.cub')
    (bg_cubes, _) = hnp.add_offsets(filter(lambda x: x.ccdname == 'BG', cubes),
                                    5, temp_token, keep=keep)
    util.log(isis.handmos(bg_cubes[1].next_path, mosaic=bg_p,
                          create='yes', nlines=base_cube.lines,
                          nsamples=base_cube.samps, nbands=1).args)
    mosaic_dejittered(bg_cubes, bg_p,
                      '{}_BG'.format(str(bg_cubes[0].get_obsid())))

    #  Create color product
    irb_p = out_root.with_name(out_root.name + '_IRB.NOPROJ.cub')
    with isis.fromlist.temp(ir_p, center_red_p, bg_p) as f:
        util.log(isis.cubeit(fromlist=f, to=irb_p, proplab=center_red_p).args)

    # Make plot of before and after flat.tab results
    dejit_flat = list(filter(lambda x: 'RED5-RED4' in x, red_flat_files))[0]
    plot_flats(flat_path, dejit_flat)

    if not keep:
        flat_path.unlink()
        red_flat_files.unlink()
        to_del.unlink()

    return