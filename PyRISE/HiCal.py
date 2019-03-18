#!/usr/bin/env python
"""Generate radiometrically corrected HiRISE Channel products."""

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


# This program is based on HiCal version 1.61 (2016/12/05),
# and on the Perl HiCal program: ($Revision: 1.52 $ $Date: 2016/12/05 19:14:24 $)
# by Eric Eliason and Richard Leis
# which is Copyright(C) 2004 Arizona Board of Regents, under the GNU GPL.
#
# Since that suite of software is under the GPL, none of it can be directly
# incorporated in this program, since I wish to distribute this software
# under the Apache 2 license.  Elements of this software (written in an entirely
# different language) are based on that software but rewritten from scratch to
# emulate functionality.

import argparse
import collections
import csv
import re
import statistics
import sys
from datetime import datetime
from pathlib import Path

import pvl

import hirise
import kalasiris as isis

from .util import PathSet


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--db',            required=False, default='HiCat.db')
    parser.add_argument('-o', '--output',  required=False, default='.HiCal.cub')
    parser.add_argument('-c', '--conf',    required=False, default='HiCal.conf')
    parser.add_argument('-k', '--keep',    required=False, default=False)
    parser.add_argument('cube', metavar=".cub-file")

    args = parser.parse_args()

    in_cube = Path(args.cube)

    if args.output.startswith('.'):
        out_cube = in_cube.with_suffix(args.output)
    else:
        out_cube = Path(args.output)

    # Get Configuration Parameters
    conf = pvl.load(args.conf)

    # Setup00 builds data structures, seems long and painful, need?
    # Setup01 - don't need - about output data routeing
    # Setup02 - db
    #   select from HiCat.EDR_Products, written by EDR_Stats
    # Setup03 - gets binning information to find binning for ObsID ... hmmm, alternate?
    #   select from HiCat.Planned_Observations
    bin2  # boolean
    bin4  # boolean

    # These come from selections on HiCat.EDR_Products
    binning  # BINNING
    mask_std  # CAL_MASK_STANDARD_DEVIATION
    dark_std  # IMAGE_DARK_STANDARD_DEVIATION
    lis_pixels  # LOW_SATURATED_PIXELS
    image_lines  # IMAGE_LINES
    line_samples  # LINE_SAMPLES
    image_mean  # select IMAGE_MEAN from HiCat.EDR_Products
    image_buffer_mean  # IMAGE_BUFFER_MEAN

    # setup done

    ccdchan = hirise.getccdchannel(in_cube)

    lis_per = db.lis_pixels / (db.image_lines * db.line_samples) * 100.0

    if (ccdchan == ('IR10', '1') and
       lis_per > conf['HiCal']['HiCal_bypass']):
        print('Bypassing IR10_1')
        sys.exit()

    (std, diff_std, zapped) = HiCal(in_cube, out_cube, ccdchan,
                                    conf, db_stuff, keep=args.keep)

    # insert these into HiCat.EDR_Products:
    # std as HIGH_PASS_FILTER_CORRECTION_STANDARD_DEVIATION
    # diff_std as DESTRIPED_DIFFERENCE_STANDARD_DEVIATION
    # The zapped flag is new, I want to put it in the DB to use later.


def HiCal(in_cube, out_cube, ccdchan, conf, db, keep=False):
    # Allows for indexing in lists ordered by bin value.
    b = 1, 2, 4, 8, 16

    # Keep from having to write out ['HiCal'] all the time.
    hconf = conf['HiCal']

    # This string will get placed in the filename for all of our
    # temporary files. It will (hopefully) prevent collisions with
    # existing files and also allow for easy clean-up if keep=True
    temp_token = datetime.now().strftime('HiCal-%y%m%d%H%M%S')

    flags = set_flags(hconf, db, ccdchan, b.index(db.bin))

    # Start processing cube files
    to_delete = PathSet()
    next_cube = in_cube.with_suffix(f'.{temp_token}.cub')
    if(db.bin > 1 and db.image_mean > 7000.0):
        furrow_cube = to_delete.add(in_cube.with_suffix(f'.{temp_token}.ffix.cub'))
        furrows_found = furrow_nulling(in_cube, furrow_cube, db.bin, ccdchan)
        next_cube = furrow_cube

    # Run hical
    lis_per = db.lis_pixels / (db.image_lines * db.line_samples) * 100.0
    hical_file = to_delete.add(next_cube.with_suffix('.hical.cub'))
    run_hical(next_cube, hical_file, hconf,
              lis_per, db.image_buffer_mean, db.bin,
              flags.noise_filter)
    next_cube = hical_file

    if furrows_found:
        lpfz_file = to_delete.add(next_cube.with_suffix('.lpfz.cub'))
        isis.lowpass(next_cube, to=lpfz_file, lines=3, samples=3,
                     minopt='COUNT', minimum=5, filter_='OUTSIDE')
        next_cube = lpfz_file

    # Perform gain-drift correction
    if(db.bin != 8):  # There is no gain fix for bin8 imaging
        higain_file = to_delete.add(next_cube.with_suffix('.fx.cub'))
        HiGainFx(next_cube, higain_file,
                 hconf['HiGainFx_Coefficient_Path'],
                 hconf['HiGainFx_Version'])
        next_cube = higain_file

    # Perform the high-pass filter cubenorm steps
    (sl, nl) = set_lines(hconf, db)
    crop_file = to_delete.add(next_cube.with_suffix('.crop.cub'))
    isis.crop(next_cube, to=crop_file, line=sl, nlines=nl)

    # This file needs to be kept for the next step, HiStitch.
    stats_file = next_cube.with_suffix('.cubenorm.tab')
    stats_fix_file = to_delete.add(next_cube.with_suffix('.cubenorm_fix.tab'))
    isis.cubenorm(crop_file, stats=stats_file, format_='TABLE')

    (std_final, zapped) = Cubenorm_Filter(stats_file, stats_fix_file,
                                          boxfilter=5, pause=True,
                                          divide=flags.divide,
                                          chan=ccdchan[1])

    # Now perform the cubnorm_plus correction
    div_or_sub = {True: 'DIVIDE', False: 'SUBTRACT'}
    cubenorm_args = {'direction': 'COLUMN',
                     'statsource': 'TABLE',
                     'normalize': 'AVERAGE',
                     'preserve': 'FALSE',
                     'mode': div_or_sub[flags.divide]}
    cubenorm_file = to_delete.add(next_cube.with_suffix('.cn.cub'))
    to_s = '{}+SignedWord+{}:{}'.format(cubenorm_file,
                                        hconf['HiCal_Normalization_Minimum'],
                                        hconf['HiCal_Normalization_Maximum'])
    isis.cubenorm(next_cube, fromstats=stats_fix_file, to=to_s, **cubenorm_args)
    next_cube = cubenorm_file

    # NoiseFilter() - external Noise_Filter
    if flags.noise_filter:
        noisefilter_file = to_delete.add(next_cube.with_suffix('.nf.cub'))
        NoiseFilter(next_cube, output=noisefilter_file,
                    conf=conf['NoiseFilter'],
                    minimum=hconf['HiCal_Normalization_Minimum'],
                    maximum=hconf['HiCal_Normalization_Maximum'], zapc=flags.zapcols)
        next_cube = noisefilter_file

    # Hidestripe() - isis.[hidestripe,hipass,lowpass,algebra]
    diff_std_dev = None
    if flags.destripe_filter:
        hidestripe_file = to_delete.add(next_cube.with_suffix('.hd.cub'))
        diff_std_dev = Hidestripe(next_cube, hidestripe_file, db.bin,
                                  hconf['HiCal_Normalization_Minimum'],
                                  hconf['HiCal_Normalization_Maximum'],
                                  hconf['HiCal_Hidestripe_Correction'],
                                  db.line_samples, keep=keep)
        next_cube = hidestripe_file

    # Create final output file
    to_delete.remove(next_cube)
    next_cube.rename(out_cub)

    if not keep:
        to_delete.unlink()

    return(std_final, diff_std_dev, zapped)


def FurrowCheck(vpnts: list, channel: int) -> int:
    # This function was brought forward from HiStitch, because
    # it is more appropriate to perform the check here, and then
    # store the result for later.
    zap = False
    i = -1 * channel

    if vpnts[i] != max(vpnts):
        zap = True

    return zap


def set_flags(conf, db, ccdchan: tuple, bindex: int) -> collections.namedtuple:
    '''Set various processing flags based on various configuration
    parameters.'''
    HiCalFlags = collections.namedtuple('HiCalFlags', ['destripe_filter',
                                                       'noise_filter',
                                                       'zapcols',
                                                       'divide'])
    destripe_filter = False
    if(1 == db.bin and (db.bin2 or db.bin4)):
        destripe_filter = True
    elif(2 == db.bin and db.bin4):
        destripe_filter = True

    noise_filter = False
    if ((process_this(ccdchan, conf['HiCal_Noise_Processing'])) and
        ((db.dark_std >=
          conf['HiCal_Noise_Bin_DarkPixel_STD'][bindex]) or
         (db.mask_std >=
          conf['HiCal_Noise_Bin_Mask_STD'][bindex]) or
         (db.lis_pixels >= conf['HiCal_Noise_LIS_Count']))):
        noise_filter = True

    zapcols = False
    if(lis_pixels >= conf['HiCal_Noise_LIS_Count']):
        zapcols = True

    divide = False
    if 'DIVIDE' == conf['HiCal_HPF_Cubenorm']:
        divide = True

    return HiCalFlags(destripe_filter, noise_filter, zapcols, divide)


def set_lines(conf, db) -> collections.namedtuple:
    '''Determine the right values for the start lines and number of lines.'''
    Lines = collections.namedtuple('Lines', ['start', 'number'])

    sl = conf['HiCal_Bin1_Skip_Top_Lines'] / db.bin + 1
    nl = db.image_lines - (hconf['HiCal_Bin1_Skip_Top_Lines'] +
                           hconf['HiCal_Bin1_Skip_Bot_Lines']) / db.bin
    if nl < 2000 / db.bin:
        sl = 1
        nl = db.image_lines

    return Lines(sl, nl)


def run_hical(in_cube, hical_cub, hconf,
              lis_per, image_buffer_mean, binning, noise_filter):

    to_s = '{}+SignedWord+{}:{}'.format(hical_cub,
                                        hconf['HiCal_Normalization_Minimum'],
                                        hconf['HiCal_Normalization_Maximum'])
    hical_args = {'to': to_s, 'units': 'IOF'}
    if(hconf['HiCal_ISIS_Conf'] != 'DEFAULT'):
        if(lis_per < 5 and image_buffer_mean > 0):
            hical_args['conf'] = hconf['HiCal_ISIS_Conf']
        else:
            hical_args['conf'] = hconf['HiCal_ISIS_Conf_Noise']

    if noise_filter:
        mask_cube = in_cube.with_sufffix('.mask.cub')
        mask(mask_cube, out_cube,
             hconf['NoiseFilter_Raw_Min'],
             hconf['NoiseFilter_Raw_Max'], binning)
        isis.hical(mask_cube, **hical_args)
        mask_cube.unlink()
    else:
        isis.hical(in_cube, **hical_args)

    return


def process_this(ccdchan: tuple, flag_list: list) -> int:
    if len(flag_list) != 28:
        raise IndexError('the list must have 28 elements')
    ccd_number = int(re.match(r"RED|IR|BG(\d+)", ccdchan[0]).group())
    i = (2 * ccd_number) + int(ccdchan[1])
    return flag_list[i]


def furrow_nulling(cube, out_cube, binning, ccdchan, keep=False):
    furrows_found = False

    to_del = PathSet()
    fcrop_file = to_del.add(out_cube.with_sufffix('.fcrop.cub'))

    isis.mask(cube, mask=cube, to=out_cube, minimum=1000000, maximum=1000000)

    furrow_values = furrow_setup(ccdchan[0], binning)
    chan_samp = chan_samp_setup(ccdchan[0], binning)

    # Crop out the portion of the image that will not be furrow checked.
    # ^- that's what the original said, but this is really cropping out
    #    the part that will be furrow corrected.
    isis.crop(cube, to=fcrop_file,
              samp=chan_samp.ssamp, nsamp=chan_samp.nsamp)

    isis.handmos(fcrop_file, mosaic=out_cube,
                 insamp=1, outsamp=chan_samp.ssamp, create='no')

    # for each column subject to furrowing, crop out each column,
    # run mask to null pixels above the furrow threshold, and
    # mosaic into the output image.
    for (s, furrow_v) in zip(chan_samp.samp, furrow_values):
        fscrop_file = to_del.add(out_cube.with_sufffix(f'.crop{s}.cub'))
        isis.crop(cube, to=fscrop_file, samp=s, nsamp=1)

        fsmask_file = to_del.add(out_cube.with_sufffix(f'.mask{s}.cub'))
        isis.mask(fscrop_file, mask=fscrop_file, to=fsmask_file,
                  minimum=0, maximum=furrow_v)

        isis.handmos(fsmask_file, mosaic=out_cube, insamp=1, outsamp=s)

    else:
        # finally, check to see if any furrow pixels were zapped.
        # The original code checked on the fist iteration, this
        # mechanism checks on the last, not sure if there was anything
        # magical about the first, or maybe we should do it each time?
        fscrop_stats = pvl.loads(isis.stats(fscrop_file).stdout)
        fsmask_stats = pvl.loads(isis.stats(fsmask_file).stdout)

        if fsmask_stats['Results']['NullPixels'] > fscrop_stats['Results']['NullPixels']:
            furrows_found = True
            # Perform simple 3x3 LPFZ filter to clean up the edges along the furrow
            trim_file = out_cube.with_sufffix('.trim.cub')
            isis.trimfilter(out_cube, to=trim_file,
                            lines=3, samples=3, minopt='COUNT')
            trim_file.rename(out_cube)

    if not keep:
        to_del.unlink()

    return furrows_found


def mask(in_cube, out_cube, noisefilter_min, noisefilter_max, binning):
    '''mask out unwanted pixels'''
    to_del = PathSet()
    temp_cube = to_del.add(out_cube.with_suffix('.mask_temp.cub'))

    isis.mask(in_cube, mask=in_cube, to=temp_cube,
              minimum=noisefilter_min, maximum=noisefilter_max,
              preserve='INSIDE', spixels='NONE')

    cubenorm_stats_file = to_del.add(temp_cube.with_suffix('.cn.stats'))
    isis.cubenorm(temp_cube, stats=cubenorm_stats_file)
    (mindn, maxdn) = analyze_cubenorm_stats(cubenorm_stats_file, binning)

    isis.mask(temp_cube, mask=temp_cube, to=out_cube,
              minimum=mindn, maximum=maxdn,
              preserve='INSIDE', spixels='NONE')

    if not keep:
        to_del.unlink()

    return


def analyze_cubenorm_stats(statsfile, binning):
    with open(statsfile) as csvfile:
        valid_points = list()
        std_devs = list()
        mins = list()
        maxs = list()
        reader = csv.DictReader(csvfile, dialect=isis.cubenormDialect)
        for row in reader:
            valid_points.append(row['ValidPoints'])
            std_devs.append(row['StdDev'])
            mins.append(row['Minimum'])
            maxs.append(row['Maximum'])

    maxvp = max(valid_points)

    # Original note:
    # # Get the median standard deviation value for all columns that have
    # # the maximum valid pixel count
    # #
    # # 2016-12-02 Note: this may not pick the median standard
    # # deviation value but the value from an index 0.95 times the number
    # # of entries in the sorted cubenorm statistics file.
    #
    # That seems to be exactly what it does.  I think the term 'median'
    # in the original (apparently pre-2016) comment is wrong.
    # The variable is called 'facstd' and when it is used below, it refers
    # to this being 'medstd + tol' which indicates that it really is
    # meant to be a factor above the median, which it is.
    std_w_maxvp = list()
    for (vp, std) in zip(valid_points, std_devs):
        if vp == maxvp:
            std_w_maxvp.append(std)

    std_w_maxvp.sort()
    facstd = std_w_maxvp[int((len(std_w_maxvp) - 1) * 0.95)]

    # Original note:
    # # find the minimum of minimums and the maximum of maximums for any
    # # column whose std is less than or equal to $medstd + $tol;
    min_w_maxvp = list()
    max_w_maxvp = list()
    for (vp, std, mi, ma) in zip(valid_points, std_devs, mins, maxs):
        if(vp >= (maxvp * 0.9) and std < facstd):
            min_w_maxvp.append(mi)
            max_w_maxvp.append(ma)

    # The original code sorts these min and max_w_maxvp arrays,
    # but then ignores this sorting.  To get that original
    # behavior, comment out these two lines:
    min_w_maxvp.sort()
    max_w_maxvp.sort()

    mindn = min_w_maxvp[int((len(min_w_maxvp) - 1) * 0.05)]
    maxdn = max_w_maxvp[int((len(max_w_maxvp) - 1) * 0.95)]

    if(1 == binning):
        mindn *= 0.7
        maxdn *= 1.3
    elif(2 == binning):
        mindn *= 0.6
        maxdn *= 1.4
    else:
        mindn *= 0.5
        maxdn *= 1.5

    return(mindn, maxdn)


def HiGainFx(cube, outcube, coef_path, version, keep=False):
    '''Perform a Gain-Drift correction on an HiIRSE Channel Image.'''
    binning = isis.getkey(cube, 'Instrument', 'Summing')
    ccd = isis.getkey(cube, 'Instrument', 'CcdId')
    chan = isis.getkey(cube, 'Instrument', 'ChannelNumber')

    coef_f_path = coef_path / f'HiRISE_Gain_Drift_Correction_Bin{binning}.{version}.csv'

    with open(coef_path) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            if hirise.getccdchannel(row['CCD CH']) == (ccd, chan):
                max_line = row['Max line']
                a_coef = (row['R(0)'], row['R(1)'], row['R(2)'])

    eqn = "'((F1/({0}+({1}*line)+({2}*line*line)))*(line<{3}) + (F1*(line>={3})))'".format(a_coef, max_line)

    tfile = outcube.with_suffix('.tempfx.cub')
    isis.fx(f1=cube, to=tfile, mode='CUBES', equation=eqn)
    isis.specadd(tfile, to=outcube, match=cube)
    if not keep:
        tfile.unlink()
    return


def Cubenorm_Filter(cubenorm_tab, outfile, pause=False,
                    boxfilter=5, divide=False, chan=None):
    '''Perform a highpass filter on the cubenorm table output of the columnar average and median values.'''
    if boxfilter < 3:
        raise ValueError(f'boxfilter={boxfilter} is less than 3')
    if not chan:
        chan = hirise.getccdchannel(cubenorm_tab)[1]

    # Make a list to receive each column
    valid_points = list()
    averages = list()
    medians = list()
    other_cols = list()
    header = list()
    with open(cubenorm_tab) as csvfile:
        reader = csv.DictReader(csvfile, dialect=isis.cubenormDialect)
        header = reader.fieldnames
        for row in reader:
            valid_points.append(row.pop('ValidPoints'))
            averages.append(row.pop('Average'))
            medians.append(row.pop('Median'))
            other_cols.append(row)

    zapped = FurrowCheck(valid_points, chan)

    avgflt = Cubenorm_Filter_filter(averages,
                                    boxfilter=boxfilter, iterations=50,
                                    chan=chan, pause=pause, vpoints=valid_points,
                                    divide=divide)
    medflt = Cubenorm_Filter_filter(medians,
                                    boxfilter=boxfilter, iterations=50,
                                    chan=chan, pause=pause, vpoints=valid_points,
                                    divide=divide)

    with open(outfile) as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=header, dialect=isis.cubenormDialect)
        writer.writeheader()
        for (d, vp, av, md) in zip(other_cols, valid_points, avgflt, medflt):
            d['ValidPoints'] = vp
            d['Average'] = av
            d['Median'] = md
            writer.writerow(d)

    # Calculate the standard deviation value for the filtered average:
    return (statistics.stdev(avgflt), zapped)


def cut_size(chan: int, length: int) -> collections.namedtuple:
    '''Determine the appropriate values for the cut sizes.'''
    Cut = collections.namedtuple('Cut', ['left', 'right'])
    c = Cut(6, 6)

    if 0 == chan:
        if 511 == length:
            c.left = 40
        elif 255 == length:
            c.left = 50
    elif 1 == chan:
        if 511 == length:
            c.right = 40
        elif 255 == length:
            c.right = 50
    else:
        raise ValueError(f'chan={chan} is not 0 or 1')

    return c


def Cubenorm_Filter_filter_boxfilter(inlist: list, origlist: list, boxfilter: int) -> list:
    x = inlist.copy()
    hwidth = int(boxfilter / 2)
    frac = 0.25
    for step in range(3):
        for it in range(1, iterations):
            xflt = x.copy()
            for i, _ in enumerate(x):
                if x[i] != 0.0:
                    xflt[i] = statistics.mean(x[i - hwidth:i + hwidth])
            x = xflt.copy()
        # Zap any columns that are different from the average by more then 25%
        if step == 2:
            frac = 0.125
        for i, (orig, new) in enumerate(zip(origlist, x)):
            if orig != 0 and new != 0 and abs(orig - new) / new > frac:
                x[i] = 0
    return x


def Cubenorm_Filter_filter(inlist, boxfilter, iterations, chan, pause, vpoints, divide):
    '''This performs highpass filtering on the passed list.'''

    x = inlist.copy()
    cut = cut_size(chan, len(x))

    # zap the left edge
    x[:cut.left] = [0] * cut.left

    # zap the right edge
    x[(len(x) - cut.right):] = [0] * cut.right

    # zap the pause point pixels
    if pause and 1023 == len(x):
        # 1st pixel = index 1
        ch_pause[0] = 252, 515, 778  # Channel 0 pause point sample locations
        ch_pause[1] = 247, 510, 773  # Channel 1 pause point sample locations
        ch_width[0] = 17, 17, 17  # Number of pixels to cut from pause point
        ch_width[1] = -17, -17, -17

        for samp, width in zip(ch_pause[chan], ch_width[chan]):
            zap_slice = pause_slicer(samp, width)
            x[zap_slice] = [0] * abs(width)

    # boxfilter
    x = Cubenorm_Filter_filter_boxfilter(x, inlist, boxfilter)

    # Perform the highpass difference of divide the original from lowpass
    maxvp = max(vpoints)
    for i, (orig, vp) in enumerate(zip(inlist, vpoints)):
        if orig != 0 and x[i] != 0 and vp == maxvp:
            if divide:
                x[i] = orig / x[i]
            else:
                x[i] = orig - x[i]
        else:
            x[i] = None

    # Need to patch up any of those None values with neighboring values:
    if 0 == chan:
        first_i = -1
        e = reversed(list(enumerate(x)))
        patch = 1
    else:
        first_i = 0
        e = list(enumerate(x))
        patch = -1

    if x[first_i] is None:
        x[first_i] = int(divide)
    for (i, x_val) in e:
        if x_val is None:
            x[i] = x[i + patch]

    return(x)


def pause_slicer(samp: int, width: int) -> slice:
    '''Returns a slice object which satisfies the range of indexes for a pause
       point.

       The incoming numbers for samp are 1-based pixel numbers, so must
       subtract 1 to get a list index.
       The width values are the number of pixels to affect, including the
       pause point pixel.  If positive they start with the pause point pixel
       and count 'up.'  If negative, they start with the pause point pixel
       and count 'down.'
    '''
    # We don't need to protect for indices less than zero or greater than the
    # length of the list, because slice objects can take values that would not
    # be valid for item access.
    s_start = None
    s_stop = None
    if width > 0:
        s_start = samp - 1
        s_stop = s_start + width
    else:
        s_start = samp + width
        s_stop = samp
    return slice(s_start, s_stop)


def highlow_destripe(in_cube, out_cube, conf, isisnorm='',
                     lnull=True, lhrs=True, lhis=True, llrs=True, llis=True,
                     keep=False):
    # Perform highpass/lowpass filter vertical destripping
    to_delete = PathSet()
    lpf_cub = to_delete.add(out_cube.with_suffix('.lpf.cub'))
    isis.lowpass(in_cub, to=lpf_cub,
                 line=conf['NoiseFilter_LPF_Line'],
                 samp=conf['NoiseFilter_LPF_Samp'],
                 minopt='PERCENT', replace='NULL',
                 minimum=['NoiseFilter_LPF_Minper'],
                 null=lnull, hrs=lhrs, his=lhis, lrs=llrs, lis=llis)

    hpf_cub = to_delete.add(out_cube.with_suffix('.hpf.cub'))
    isis.highpass(in_cub, to=hpg_cub, minopt='PERCENT',
                  line=conf['NoiseFilter_HPF_Line'],
                  samp=conf['NoiseFilter_HPF_Samp'],
                  minimum=conf['NoiseFilter_HPF_Minper'])

    isis.algebra(from_=lpf_cub, from2=hpf_cub,
                 to=out_cube.with_suffix('.cub' + isisnorm),
                 operator='ADD')
    if not keep:
        to_delete.unlink()

    return


def getHistVal(histogram, conf):
    '''Return information about the histogram'''
    lisper = 0
    if int(histogram['Total Pixels']) - int(histogram['Null Pixels']) > 0:
        lisper = int(histogram['Lis Pixels']) / (int(histogram['Total Pixels'])
                                                 - int(histogram['Null Pixels'])) * 100
    cumper = conf['NoiseFilter_HighEnd_Percent']
    if cumper < 99.0:
        cumper = 99.0

    if lisper > conf['NoiseFilter_Hard_Tolmax']:
        hard_high_end = conf['NoiseFilter_Hard_HighEnd_Percent']
        if hard_high_end < 99.9:
            hard_high_end = 99.9
        cumper = hard_high_end

    for row in histogram:
        if row.CumulativePercent > cumper:
            maxval = row.DN
            break

    return(lisper, maxval)


def NoiseFilter_noisefilter(from_cube, to_cube, flattol, conf,
                            maxval, tolmin, tolmax):
    '''Convenience function for the repeated noisefiltering.'''
    isis.noisefilter(from_cub, to=to_cube, flattol=flattol,
                     low=conf['NoiseFilter_Minimum_Value'],
                     high=maxval,
                     tolmin=tolmin, tolmax=tolmax,
                     sample=conf['NoiseFilter_Noise_Samp'],
                     line=conf['NoiseFilter_Noise_Line'],
                     toldef='STDDEV', replace='NULL',
                     lisisnoise=True, lrsisnoise=True)


def NoiseFilter_cubenorm_edit(in_tab: Path, out_tab: Path):
    '''This function zaps the relevent pixels in the cubenorm output and
       creates an edited cubenorm file.'''
    # Slightly different values from other function, not entirely sure why.
    # Pause point locations are 1-based pixel numbers, so -1 to get list index.
    # With values are the number of pixels to affect, including the pause point pixel
    ch_pause[0] = 1, 252, 515, 778  # Channel 0 pause point sample locations
    ch_pause[1] = 247, 510, 773, 1024  # Channel 1 pause point sample locations
    ch_width[0] = 3, 6, 6, 6  # Number of pixels to cut from pause point
    ch_width[1] = -8, -7, -6, -3  # sign indicates direction of cut from pause point

    vpnts = list()
    other_cols = list()
    header = list()
    with open(in_tab) as csvfile:
        reader = csv.DictReader(csvfile, dialect=isis.cubenormDialect)
        header = reader.fieldnames
        for row in reader:
            vpnts.append(row.pop['ValidPoints'])
            other_cols.append(row)

    max_vpnts = max(vpnts)
    if max_vpnts <= 0:
        max_vpnts = 1

    # Create a 'unity array' (original code has two, but they're identical).
    # Zap any columns with less then NoiseFilter_Zap_Fraction
    norm = [1] * len(vpnts)
    for i, v in enumerate(vpnts):
        if(zapc and v / max_vpnts < conf['NoiseFilter_Zap_Fraction']):
            norm[i] = 0

    # Determine if the pause point pixels need to be zapped
    for samp, width in zip(ch_pause[chan], ch_width[chan]):
        zap_slice = pause_slicer(samp, width)
        for i in range(zap_slice):
            if vpnts[i] / max_vpnts < conf['NoiseFilter_Nonvalid_Fraction']:
                norm[i] = 0

    with open(out_tab) as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=header, dialect=isis.cubenormDialect)
        writer.writeheader()
        for (d, vp, n) in zip(other_cols, vpnts, norm):
            d['ValidPoints'] = vp
            d['Average'] = n
            d['Median'] = n
            d['StdDev'] = n
            d['Minimum'] = n
            d['Maximum'] = n
            writer.writerow(d)
    return


def NoiseFilter(in_cube, output, conf, minimum=None, maximum=None, zapc=False, keep=False):
    '''Perform salt/pepper noise removal.'''
    binning = isis.getkey_k(in_cube, 'Instrument', 'Summing')
    (ccd, chan) = hirise.getccdchannel(isis.getkey_k(in_cube, 'Archive', 'ProductId'))
    isisnorm = ''
    if minimum is not None and maximum is not None:
        isisnorm = '+SignedWord+{}:{}'.format(minimum, maximum)

    to_delete = PathSet()

    h = isis.Histogram(in_cube)
    (LisP, MaxVal) = getHistVal(h, conf)

    cn_tab = to_delete.add(output.with_suffix('.cn.tab'))
    isis.cubenorm(in_cube, stats=cn_tab, format_='TABLE', direction='COLUMN')

    cn2_tab = to_delete.add(output.with_suffix('.cn2.tab'))
    NoiseFilter_cubenorm_edit(cn_tab, cn2_tab)

    # Zap the bad colmns for the highpass and lowpass filter
    zap_cub = to_delete.add(output.with_suffix('.zap.cub'))
    isis.cubenorm(in_cube, to=zap_cub, fromstats=cn2_tab, statsource='TABLE',
                  mode='DIVIDE', norm='AVE', preserve='FALSE')

    # Perform highpass/lowpass filter vertical destripping
    add_cub = to_delete.add(output.with_suffix('.add.cub'))
    highlow_destripe(zap_cub, add_cub, conf, isisnorm, llis=False, keep=keep)

    # Perform the 1st noise filter
    tolmin = conf['NoiseFilter_Tolmin']
    tolmax = conf['NoiseFilter_Tolmax']
    if LisP >= conf['NoiseFilter_Hard_Filtering']:
        tolmin = conf['NoiseFilter_Hard_Tolmin']
        tolmax = conf['NoiseFilter_Hard_Tolmax']
    flattol = h['Std Deviation'] * conf['NoiseFilter_Flattol']
    if flattol < 0.00001:
        flattol = 0.00001

    nf1_cub = to_delete.add(output.with_suffix('.nf1.cub'))
    NoiseFilter_noisefilter(add_cub, nf1_cub.with_suffix('.cub' + isisnorm),
                            flattol, conf, MaxVal, tolmin, tolmax)

    # Perform the 2nd noise filter
    nf2_cub = to_delete.add(output.with_suffix('.nf2.cub'))
    NoiseFilter_noisefilter(nf1_cub, nf2_cub.with_suffix('.cub' + isisnorm),
                            flattol, conf, MaxVal, tolmin, tolmax)

    # Perform the 3rd noise filter
    nf3_cub = to_delete.add(output.with_suffix('.nf3.cub'))
    NoiseFilter_noisefilter(nf2_cub, nf3_cub.with_suffix('.cub' + isisnorm),
                            flattol, conf, MaxVal, tolmin, tolmax)

    # Perform another highpass/lowpass filter now that the
    # data are much cleaner
    add2_cub = to_delete.add(output.with_suffix('.add2.cub'))
    highlow_destripe(nf3_cub, add2_cub, conf, isisnorm,
                     lnull=False, lhrs=False, lhis=False, llrs=False,
                     llis=False, keep=keep)

    if ccd.startswith('RED'):
        # Perform LPFZ  filters if we have a RED filter image.
        # For IR and BG filter data, assume that the HiColorNorm pipeline
        # step will interpolate using the BG/RED and IR/RED ratio data.
        lowmin = int(conf['NoiseFilter_LPFZ_Line'] *
                     conf['NoiseFilter_LPFZ_Samp'] / 3)
        lpfz_cub = to_delete.add(output.with_suffix('.lpfz.cub'))
        isis.lowpass(add2_cub, to=lpfz_cub.with_suffix('.cub' + isisnorm),
                     sample=3, line=3, minopt='COUNT', minimum=1,
                     filter_='OUTSIDE',
                     null=True, hrs=False, his=True, lrs=True, lis=True)
        isis.lowpass(lpfz_cub, to=output.with_suffix('.cub' + isisnorm),
                     sample=conf['NoiseFilter_LPFZ_Samp'],
                     line=conf['NoiseFilter_LPFZ_Line'],
                     minopt='COUNT', minimum=lowmin, filter_='OUTSIDE',
                     null=True, hrs=False, his=True, lrs=True, lis=True)
    else:
        to_delete.remove(add2_cub)
        add2_cub.rename(output)

    if not keep:
        to_delete.unlink()
    return


def Hidestripe(in_cube, out_cube, binning, minimum, maximum, hidcorr,
               line_samples, keep=False) -> str:
    # SignedWord+$HiCal_Normalization_Minimum:$HiCal_Normalization_Maximum
    to_s = '+SignedWord+{}:{}'.format(minimum, maximum)
    to_del = PathSet()
    if 1 == binning:
        temp_cub = to_del.add(out_cube.with_suffix('.hd.cub'))
        isis.hidestripe(in_cube, to=temp_cub + to_s, parity='EVEN', correction=hidcorr)
        isis.hidestrip(temp_cub, to=out_cube + to_s, parity='ODD', correction=hidcorr)
    else:
        lpf_cube = to_del.add(out_cube.with_suffix('.lpf.cub'))
        hpf_cube = to_del.add(out_cube.with_suffix('.hpf.cub'))
        boxsamp = (2 * line_samples) - 1

        isis.lowpass(in_cube, to=lpf_cube, samples=boxsamp, lines=3,
                     null=False, hrs=False, his=False, lrs=False, lis=False)
        isis.highpass(in_cube, to=hpf_cube, samples=boxsamp, lines=1,
                      propagate=True)
        isis.algebra(from_=lpf_cube, from2=hpf_cube, to=out_cube + to_s,
                     operator='ADD', a='1.0', b='1.0')

    # Standard deviation of difference between cube after noise
    # filter and cube after hidestripe
    diff_cube = to_del.add(out_cube.with_suffix('.diff.cub'))
    isis.algebra(from_=out_cube, from2=in_cube, to=diff_cube, operator='subtract')
    stddev = pvl.loads(isis.stats(diff_cube).stdout)['Results']['StandardDeviation']

    if not keep:
        to_delete.unlink()

    return stddev


def chan_samp_setup(channel, binning):
    '''Returns a named tuple which contains a list and two numbers.'''
    samp = dict()
    # samp[chan][binning]
    samp['0']['2'] = 1,  2,  3,  4,  5,  6,  7,  8,  9, 10
    samp['1']['2'] = 512, 511, 510, 509, 508, 507, 506, 505, 504, 503
    samp['0']['4'] = 1,  2,  3,  4,  5,  6
    samp['1']['4'] = 256, 255, 254, 253, 252, 251

    ssamp = dict()
    # ssamp[chan][binning]
    ssamp['0']['2'] = 11
    ssamp['1']['2'] = 1
    ssamp['0']['4'] = 7
    ssamp['1']['4'] = 1

    nsamp = dict()
    # nsamp[chan][binning]
    nsamp['0']['2'] = 502
    nsamp['1']['2'] = 502
    nsamp['0']['4'] = 250
    nsamp['1']['4'] = 250

    ChanSamp = collections.namedtuple('ChanSamp', ['samp', 'ssamp', 'nsamp'])

    c = ChanSamp(samp[channel][binning], ssamp[channel][binning], nsamp[channel][binning])
    return c


def furrow_setup(ccd, binning):
    '''Returns the right tuple of furrow.'''
    # Expect to call the returned dict like this: d[ccd][binning]
    # Assume each channel for each CCD will have the same threshold
    d = dict()
    d['RED0']['2'] = 8000, 8100,  8700,  9200,  9600, 10000, 12000, 12000, 12000, 12000
    d['RED0']['4'] = 8000, 9000,  9500,  9900,  9900, 10000
    d['RED1']['2'] = 7200, 7200,  7800,  8400,  9000,  9500, 12000, 12000, 12000, 12000
    d['RED1']['4'] = 8000, 8100,  9200,  9600,  9800, 10000
    d['RED2']['2'] = 7800, 7800,  8400,  9000,  9600, 10000, 12000, 12000, 12000, 12000
    d['RED2']['4'] = 8000, 8700,  9500,  9800,  9900, 10000
    d['RED3']['2'] = 7800, 8100,  8300,  9200,  9600, 10000, 12000, 12000, 12000, 12000
    d['RED3']['4'] = 7900, 9200,  9700,  9900, 10000, 10500
    d['RED4']['2'] = 7800, 7800,  8300,  9000,  9500,  9900, 12000, 12000, 12000, 12000
    d['RED4']['4'] = 8000, 8700,  9700, 10000, 10300, 10600
    d['RED5']['2'] = 7900, 8200,  8600,  9200,  9600, 10000, 12000, 12000, 12000, 12000
    d['RED5']['4'] = 8000, 9300,  9700,  9900, 10200, 10700
    d['RED6']['2'] = 7500, 7500,  8100,  8500,  9200, 10000, 12000, 12000, 12000, 12000
    d['RED6']['4'] = 8000, 8400,  9700, 10000, 10500, 10700
    d['RED7']['2'] = 7600, 8300,  8900,  9400,  9900, 11000, 12000, 12000, 12000, 12000
    d['RED7']['4'] = 7700, 9600, 10000, 10200, 11000, 12000
    d['RED8']['2'] = 7200, 7200,  7900,  8500,  9000,  9400, 12000, 12000, 12000, 12000
    d['RED8']['4'] = d['RED7']['4']
    d['RED9']['2'] = 7600, 8300,  8600,  9200,  9600, 10000, 12000, 12000, 12000, 12000
    d['RED9']['4'] = 8000, 8800,  9200,  9400,  9800, 10500
    d['IR10']['2'] = d['RED0']['2']
    d['IR10']['4'] = 7600, 8300,  9000, 10000, 10500, 12000
    d['IR11']['2'] = d['RED0']['2']
    d['IR11']['4'] = d['IR10']['4']
    d['BG12']['2'] = d['RED0']['2']
    d['BG12']['4'] = d['IR10']['4']
    d['BG13']['2'] = d['RED0']['2']
    d['BG13']['4'] = d['IR10']['4']

    return d[ccd][binning]


if __name__ == "__main__":
    main()
