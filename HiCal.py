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

import argparse, collections, csv
import isis, pvl

def main():
    parser = argparse.ArgumentParser( description=__doc__ )
    parser.add_argument('--db',          required=False, default='HiCat.db' )
    parser.add_argument('-o','--output', required=False, default='.HiCal.cub')
    parser.add_argument('-k','--keep',   required=False, default=False )
    parser.add_argument('cube', metavar=".cub-file" )

    args = parser.parse_args()

    ofile_cub = ''
    if args.output.startswith('.'): 
        ofile_cub = os.path.splitext( args.cube )[0] + args.output
    else: ofile_cub = args.output

    # GetConfigurationParameters() - Read some stuff from HiCal.conf, could be done with pvl
    HiCal_Noise_Processing # a set of bools to indicate which channel gets noise processed
    HiCal_Noise_Bin_DarkPixel_STD # a list of numbers, indexable by binning
    HiCal_Noise_Bin_Mask_STD # a list of numbers, indexable by binning
    HiCal_Noise_LIS_Count # single value
    NoiseFilter_Raw_Min # single value
    NoiseFilter_Raw_Max # single value
    HiCal_Normalization_Minimum # single value
    HiCal_Normalization_Maximum # single value
    HiCal_ISIS_Conf # single value
    HiCal_ISIS_Conf_Noise # single value
    HiGainFx_Coefficient_Path # single value
    HiGainFx_Version # single value
    HiCal_Bin1_Skip_Top_Lines # single value
    HiCal_Bin1_Skip_Bot_Lines # single value
    HiCal_HPF_Cubenorm # single value

    # Setup00 builds data structures, seems long and painful, need?
    # Setup01 - don't need - about output data routeing
    # Setup02 - db
    #   select from HiCat.EDR_Products, written by EDR_Stats
    # Setup03 - gets binning information to find binning for ObsID ... hmmm, alternate?
    #   select from HiCat.Planned_Observations
    # ProcessingSwitches() sets a variety of booleans
        
    binning # select BINNING from HiCat.EDR_Products
    mask_std # CAL_MASK_STANDARD_DEVIATION
    dark_std # IMAGE_DARK_STANDARD_DEVIATION
    lis_pixels # LOW_SATURATED_PIXELS
    image_lines # IMAGE_LINES
    line_samples # LINE_SAMPLES
    image_mean # select IMAGE_MEAN from HiCat.EDR_Products 
    image_buffer_mean # IMAGE_BUFFER_MEAN

    lis_per = lis_pixels/( image_lines * line_samples)*100.0
    
    ### setup done ###

    if( binning > 1 and image_mean > 7000.0):
        (furrows_found, furrow_cube) = furrow_nulling( args.cube, binning )

    noise_filter = False
    if( HiCal_Noise_Processing[ getccdchannel(args.cube) ] ):
        if dark_std   >= HiCal_Noise_Bin_DarkPixel_STD[binning] : noise_filter = True
        if mask_std   >= HiCal_Noise_Bin_Mask_STD[binning]      : noise_filter = True
        if lis_pixels >= HiCal_Noise_LIS_Count)                 : noise_filter = True

    # Run hical
    hical_file = os.path.splitext( args.cube )[0] + '.hical.cub'
    hical_args = {'to': f'{hical_file}+SignedWord+{HiCal_Normalization_Minimum}:{HiCal_Normalization_Maximum}',
                  'units':'IOF'}
    if( HiCal_ISIS_Conf != 'DEFAULT' ):
        if( lis_per < 5 and image_buffer_mean > 0): 
              hical_args['conf']=HiCal_ISIS_Conf
        else: hical_args['conf']=HiCal_ISIS_Conf_Noise
            
    if noise_filter: 
        isis.hical( mask(furrow_cube, NoiseFilter_Raw_Min, NoiseFilter_Raw_Max, binning), **hical_args )
    else:
        isis.hical( furrow_cube, **hical_args ) 

    next_file = hical_file

    if furrows_found:
        lpfz_file = os.path.splitext( args.cube )[0] + '.lpfz.cub'
        isis.lowpass( next_file, to=lpfz_file, lines=3, samples=3, minopt='COUNT', minimum=5, filter='OUTSIDE'
        next_file = lpfz_file

    # Perform gain-drift correction
    if( 8 != binning ): # There is no gain fix for bin8 imaging
        higain_file = os.path.splitext( args.cube )[0] + '.fx.cub'
        HiGainFx( next_file, higain_file, HiGainFx_Coefficient_Path, HiGainFx_Version )
        next_file = higain_file
    
    # Perform the high-pass filter cubenorm step
    sl = HiCal_Bin1_Skip_Top_Lines/binning + 1
    nl = image_lines - (HiCal_Bin1_Skip_Top_Lines + HiCal_Bin1_Skip_Bot_Lines)/binning
    if nl < 2000/binning: 
        sl = 1
        nl = image_lines
    crop_file = os.path.splitext( args.cube )[0] + '.crop.cub'
    isis.crop( next_file, to=crop_file, line=sl, nlines=nl )
    stats_file = os.path.splitext( args.cube )[0] + '.cubenorm.tab'
    isis.cubenorm( crop_file, stats=stats_file, format='TABLE' )

    Cubenorm_Filter( stats_file, outfile, boxfilter=5, pause=True, divide=HiCal_HPF_Cubenorm )

    #   insert into HiCat.EDR_Products
    #   isis.cubenorm again
    # NoiseFilter() - external Noise_Filter
    # Hidestripe() - isis.[hidestripe,hipass,lowpass,algebra] 
    #   insert into HiCat.EDR_Products

def furrow_nulling( cube, binning ):
    furrows_found = False

    f = os.path.splitext( cube )[0]
    furrow_fix_file = f + '.furrowfix.cub'
    fcrop_file      = f + '.fcrop.cub'
    trim_file       = f + '.trim.cub'

    isis.mask( cube, mask=cube, to=furrow_fix_file, minimum=1000000, maximum=1000000 )

    ccdchan = getccdchannel(cube)

    furrow_values = furrow_setup( ccdchan[0], binning )
    chan_samp = chan_samp_setup( ccdchan[0], binning )

    # Crop out the portion of the image that will not be furrow checked.
    # ^- that's what the original said, but this is really cropping out
    #    the part that will be furrow corrected.
    isis.crop( cube, to=fcrop_file, 
               samp=chan_samp.ssamp, nsamp=chan_samp.nsamp )

    isis.handmos( fcrop_file, mosaic=furrow_fix_file, 
                  insamp=1, outsamp=chan_samp.ssamp, create='no' )

    # for each column subject to furrowing, crop out each column,
    # run mask to null pixels above the furrow threshold, and
    # mosaic into the output image. 
    for (s, furrow_v) in zip(chan_samp.samp, furrow_values):
        fscrop_file = f + f'.crop{s}.cub'
        isis.crop( cube, to=fscrop_file, samp=s, nsamp=1 )

        fsmask_file = f + f'.mask{s}.cub'
        isis.mask( fscrop_file, mask= fscrop_file, to= fsmask_file, 
                   minimum=0, maximum=furrow_v )

        isis.handmos( fsmask_file, mosaic=furrow_fix_file, insamp=1, outsamp=s )

    else: # finally, check to see if any furrow pixels were zapped.
          # The original code checked on the fist iteration, this 
          # mechanism checkson the last, not sure if there was anything
          # magical about the first, or maybe we should do it each time?
        fscrop_stats = pvl.loads( isis.stats( fscrop_file ).stdout )
        fsmask_stats = pvl.loads( isis.stats( fsmask_file ).stdout )

        if fsmask_stats['Results']['NullPixels'] > fscrop_stats['Results']['NullPixels']:
            furrows_found = True
    # clean up these files we made in this previous loop?

    if furrows_found:
        # Perform simple 3x3 LPFZ filter to clean up the edges along the furrow
        isis.trimfilter( furrow_fix_file, to=trim_file, 
                         lines=3, samples=3, minopt='COUNT' )
        return (True, trim_file)
    else: 
        return (False, furrow_fix_file)


def mask( cube, noisefilter_min, noisefilter_max, binning ):
    '''mask out unwanted pixels'''
    f = os.path.splitext( cube )[0]
    mask_file = f+'.mask.cub'
    tmask_file = f+'.tmask.cub'

    isis.mask( cube, mask=cube, to=tmask_file, 
               minimum=noisefilter_min, maximum=noisefilter_max, 
               preserve='INSIDE', spixels='NONE' )

    cubenorm_stats_file = os.path.splitext( tmask_file )[0]+'.cubenorm.stats'
    isis.cubenorm( tmask_file, stats= cubenorm_stats_file )
    (mindn, maxdn) = analyze_cubenorm_stats( cubenorm_stats_file, binning )

    isis.mask( tmask_file, mask=tmask_file, to=mask_file, 
               minimum=mindn, maximum=maxdn,
               preserve='INSIDE', spixels='NONE' )

    return( mask_file )
    

def analyze_cubenorm_stats( statsfile, binning ):
    with open( statsfile ) as csvfile:
        valid_points = list()
        std_devs     = list()
        mins         = list()
        maxs         = list()
        reader = csv.DictReader( csvfile, dialect=isis.cubenormDialect )
        for row in reader:
            valid_points.append( row['ValidPoints'] )
            std_devs.append(     row['StdDev']      )
            mins.append(         row['Minimum']     )
            maxs.append(         row['Maximum']     )

    maxvp = max( valid_points )

    # Original note:
    ## Get the median standard deviation value for all columns that have
    ## the maximum valid pixel count
    ## 
    ## 2016-12-02 Note: this may not pick the median standard
    ## deviation value but the value from an index 0.95 times the number
    ## of entries in the sorted cubenorm statistics file.
    #
    # That seems to be exactly what it does.  I think the term 'median'
    # in the original (apparently pre-2016) comment is wrong.
    # The variable is called 'facstd' and when it is used below, it refers
    # to this being 'medstd + tol' which indicates that it really is 
    # meant to be a factor above the median, which it is.
    std_w_maxvp = list()
    for (vp, std) in zip(valid_points, std_devs):
        if vp == maxvp: std_w_maxvp.append( std )

    std_w_maxvp.sort()
    facstd = std_w_maxvp[ int( (len(std_w_maxvp)-1)*0.95 ) ]
    
    # Original note:
    ## find the minimum of minimums and the maximum of maximums for any
    ## column whose std is less than or equal to $medstd + $tol;
    min_w_maxvp = list()
    max_w_maxvp = list()
    for (vp,std,mi,ma) = zip( valid_points,std_devs,mins,maxs ):
        if( vp >= (maxvp * 0.9) and std < facstd ):
            min_w_maxvp.append( mi )
            max_w_maxvp.append( ma )

    # The original code sorts these min and max_w_maxvp arrays,
    # but then ignores this sorting.  To get that original 
    # behavior, comment out these two lines:
    min_w_maxvp.sort()
    max_w_maxvp.sort()

    mindn = min_w_maxvp[ int( (len(min_w_maxvp)-1)*0.05 ) ]
    maxdn = max_w_maxvp[ int( (len(max_w_maxvp)-1)*0.95 ) ]

    if( 1 == binning ):
        mindn *= 0.7
        maxdn *= 1.3
    elif( 2 == binning ):
        mindn *= 0.6
        maxdn *= 1.4
    else:
        mindn *= 0.5
        maxdn *= 1.5

    return( mindn, maxdn )


def HiGainFx( cube, outcube, coef_path, version ):
    '''Perform a Gain-Drift correction on an HiIRSE Channel Image.'''
    binning = isis.getkey( cube, 'Instrument', 'Summing' )
    ccd =     isis.getkey( cube, 'Instrument', 'CcdId' )
    chan =    isis.getkey( cube, 'Instrument', 'ChannelNumber' )

    coef_f_path = os.path.join( coef_path, 
                                f'HiRISE_Gain_Drift_Correction_Bin{binning}.{version}.ccv' )

    with open( coef_path ) as csvfile:
        reader = csv.DictReader(csvfile):
            for row in reader:
                if hirise.getccdchannel( row['CCD CH'] ) == (ccd, chan):
                    max_line = row['Max line']
                    a_coef = ( row['R(0)'], row['R(1)'], row['R(2)'] )

    eqn = "'((F1/({0}+({1}*line)+({2}*line*line)))*(line<{3}) + (F1*(line>={3})))'".format( a_coef, max_line )

    tfile = os.path.splitext( cube )[0]+'.tempfx.cub'
    isis.fx( f1=cube, to=tfile, mode='CUBES', equation=eqn )
    isis.specadd( tfile, to=outcube, match=cube )
    return

def Cubenorm_Filter( cubenorm_tab, outfile, pause=False, boxfilter=5, divide=False, chan=None)
    '''Perform a highpass filter on the cubenorm table output of the columnar average and median values.'''
    if boxfilter < 3: raise ValueError(f'boxfilter={boxfilter} is less than 3')
    if not chan: chan = hirise.getccdchannel( cubenorm_tab )[1]

    # Make a list to receive each column
    valid_points = list()
    averages     = list()
    medians      = list()
    other_cols   = list()
    header       = list()
    with open( cubenorm_tab ) as csvfile:
        reader = csv.DictReader( csvfile, dialect=isis.cubenormDialect )
        header = reader.fieldnames
        for row in reader:
            valid_points.append( row.pop('ValidPoints') )
            averages.append(     row.pop('Average')      )
            medians.append(      row.pop('Median')     )
            other_cols.append( row )
    maxvp = max( valid_points )

    avgflt = Cubenorm_Filter_filter( averages, boxfilter=boxfilter, iterations=50, chan=chan, pause=pause )
    # do stuff

    with open( outfile ) as csvfile:
        writer = csv.DictWriter( csvfile, fieldnames=header, dialect=isis.cubenormDialect )
        writer.writeheader()
        for (d, vp, av, md) in zip(other_cols, valid_points, averages, medians):
            d['ValidPoints'] = vp
            d['Average'] = av
            d['Median'] = md
            writer.writerow( d )

def Cubenorm_Filter_filter( inlist, boxfilter, iterations, chan, pause )
    '''This performs highpass filtering on the passed list.'''

    x = inlist.copy()

    left_cut = 6
    right_cut = 6
    ch_pause[0] = 252,515,778 # Channel 0 pause point sample locations (1st pixel = index 1)
    ch_pause[1] = 247,510,773 # Channel 1 pause point sample locations
    ch_width[0] = 17,17,17 # Number of pixels to cut from pause point
    ch_width[1] = 17,17,17
    ch_direc[0] = 'right' # Direction of cut
    ch_direc[1] = 'left'

    if 0 == chan:
        if 511 == len(x): left_cut = 40
        elif 255 == len(x): left_cut = 50
    elif 1 == chan: 
        if 511 == len(x): right_cut = 40
        elif 255 == len(x): right_cut = 50
    else:
        raise ValueError(f'chan={chan} is not 0 or 1')

    # zap the left edge
    x[:(left_cut-1)] = 0

    # zap the right edge
    x[(len(x)-right_cut):] = 0

    # zap the pause point pixels
    if pause and 1023 == len(x):
        for samp, width in zip(ch_pause[chan], ch_width[chan]):
            if 'left' == ch_direc[chan]:
                i1 = samp - width
                i2 = samp - 1
            else:
                i1 = samp - 1
                i2 = samp + width - 2

            if i1 < 0        : i1 = 0
            if i2 > len(x)-1 : i2 = len(x) - 1
            x[i1:i2] = 0

    # boxfilter
    hwidth = int( boxfilter/2 )
    for step in range(1,3):
        for it in range(1,iterations):
            xflt = list()
            for i in x:
                i1 = i - hwidth
                if i1 < 0: i1 = 0
                i2 = i + hwidth
                if i2 > len(x)-1: i2 = len(x) - 1
                xflt.append( statistics.mean( x[i1:i2] ) )
            x = xflt.copy()
        if step < 3:
            frac = 0.25
            if step >= 2: frac = 0.125

    

def chan_samp_setup( channel, binning )
    '''Returns a named tuple which contains a list and two numbers.'''
    samp = dict()
    #samp[chan][binning]
    samp['0']['2'] =   1,  2,  3,  4,  5,  6,  7,  8,  9, 10
    samp['1']['2'] = 512,511,510,509,508,507,506,505,504,503
    samp['0']['4'] =   1,  2,  3,  4,  5,  6
    samp['1']['4'] = 256,255,254,253,252,251

    ssamp = dict()
    #ssamp[chan][binning]
    ssamp['0']['2'] = 11
    ssamp['1']['2'] = 1
    ssamp['0']['4'] = 7
    ssamp['1']['4'] = 1

    nsamp = dict()
    #nsamp[chan][binning]
    nsamp['0']['2'] = 502
    nsamp['1']['2'] = 502
    nsamp['0']['4'] = 250
    nsamp['1']['4'] = 250

    ChanSamp = collections.namedtuple( 'ChanSamp', ['samp', 'ssamp', 'nsamp'] )

    c = ChanSamp( samp[channel][binning], ssamp[channel][binning], nsamp[channel][binning] )
    return c



def furrow_setup( ccd, binning )
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
