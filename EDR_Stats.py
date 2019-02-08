#!/usr/bin/env python
"""Create an ISIS cube from a HiRISE EDR .img file and place HiRISE EDR image statistics into HiCat's EDR_Products table."""

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


# This program is based on EDR_Stats version 2.16.1 (2016/06/16)
# by Eric Eliason and Audrie Fennema 
# which is Copyright(C) 2004 Arizona Board of Regents, under the GNU GPL.
#
# Since that suite of software is under the GPL, none of it can be directly
# incorporated in this program, since I wish to distribute this software 
# under the Apache 2 license.  Elements of this software (written in an entirely
# different language) are based on that software but rewritten from scratch to
# emulate functionality.

import argparse, csv, math, os
import hirise, pvl

def main():
    parser = argparse.ArgumentParser( description=__doc__ )
    parser.add_argument('--db',          required=False, default='HiCat.db' )
    parser.add_argument('-t','--table',  required=False, default='EDR_Products')
    parser.add_argument('-o','--output', required=False, default='.EDR_Stats.cub')
    parser.add_argument('--histmax',     required=False, default=0.01 )
    parser.add_argument('--histmin',     required=False, default=99.99 )
    parser.add_argument('-k','--keep',   required=False, default=False )
    parser.add_argument('img', metavar=".img-file" )

    args = parser.parse_args()

    # Need to parse a .conf file?  Not sure.

    ofile_cub = ''
    if args.output.startswith('.'): 
        ofile_cub = os.path.splitext( args.img )[0] + args.output
    else: ofile_cub = args.output

    # Convert to .cub
    isis.hi2isis( args.img, to=ofile_cub )

    # Get some info from the new cube:
    product_id    = isis.getkey( ofile_cub, 'Archive','ProductId')
    image_lines   = isis.getkey( ofile_cub, 'Dimensions','Lines')
    image_samples = isis.getkey( ofile_cub, 'Dimensions','Samples')

    histats = parse_histat( isis.histat( ofile_cub, useoffsets=True,
                                            leftimage=0,     rightimage=1, 
                                            leftcalbuffer=3, rightcalbuffer=1,
                                            leftcaldark=3,   rightcaldark=1,
                                            leftbuffer=3,    rightbuffer=1,
                                            leftdark=3,      rightdark=1 ) )


    dncnt = get_dncnt( ofile_cub, args.histmin, args.histmax, keep=args.keep )
    snr = snr( ofile_cub, gainspvl, histats ) 
    gapp = ( histats['IMAGE_MAXIMUM'] / (image_lines * image_samples) )*100.0

    # DB stuff
    # add a bunch of stuff from the histats call, gapp, snr, dncnt


def parse_histat( pvltext ):
    '''Parse the output of histat into a dictionary'''

    p = pvl.loads( pvltext )
    d = {}

    # Image Area Statistics  
    d['IMAGE_MEAN']               = p['IMAGE']['Average']
    d['IMAGE_STANDARD_DEVIATION'] = p['IMAGE']['StandardDeviation']
    d['IMAGE_MINIMUM']            = p['IMAGE']['Minimum']
    d['IMAGE_MAXIMUM']            = p['IMAGE']['Maximum']
    d['GAP_PIXELS']               = p['IMAGE']['NullPixels']
    d['LOW_SATURATED_PIXELS']     = p['IMAGE']['LisPixels']
    d['HIGH_SATURATED_PIXELS']    = p['IMAGE']['HisPixels']

    # Calibration Reverse-readout Statistics   
    d['CAL_REVERSE_MEAN']               = p['CAL_REVERSE']['Average']
    d['CAL_REVERSE_STANDARD_DEVIATION'] = p['CAL_REVERSE']['StandardDeviation']
    d['CAL_REVERSE_MINIMUM']            = p['CAL_REVERSE']['Minimum']
    d['CAL_REVERSE_MAXIMUM']            = p['CAL_REVERSE']['Maximum']
   
    # Calibration Mask Statistics   
    d['CAL_MASK_MEAN']               = p['CAL_MASK']['Average']
    d['CAL_MASK_STANDARD_DEVIATION'] = p['CAL_MASK']['StandardDeviation']
    d['CAL_MASK_MINIMUM']            = p['CAL_MASK']['Minimum']
    d['CAL_MASK_MAXIMUM']            = p['CAL_MASK']['Maximum']
      
    # Calibration Ramp Statistics   
    d['CAL_RAMP_MEAN']               = getkey( pvl, 'CAL_RAMP','Average' )
    d['CAL_RAMP_STANDARD_DEVIATION'] = getkey( pvl, 'CAL_RAMP','StandardDeviation' )
    d['CAL_RAMP_MINIMUM']            = getkey( pvl, 'CAL_RAMP','Minimum' )
    d['CAL_RAMP_MAXIMUM']            = getkey( pvl, 'CAL_RAMP','Maximum' )
   
    # Image Dark Reference Statistics
    d['IMAGE_DARK_MEAN']               = getkey( pvl, 'IMAGE_DARK','Average' )
    d['IMAGE_DARK_STANDARD_DEVIATION'] = getkey( pvl, 'IMAGE_DARK','StandardDeviation' )
    d['IMAGE_DARK_MINIMUM']            = getkey( pvl, 'IMAGE_DARK','Minimum' )
    d['IMAGE_DARK_MAXIMUM']            = getkey( pvl, 'IMAGE_DARK','Maximum' )

    # Image Buffer Area
    d['IMAGE_BUFFER_MEAN']               = getkey( pvl, 'IMAGE_BUFFER','Average' )
    d['IMAGE_BUFFER_STANDARD_DEVIATION'] = getkey( pvl, 'IMAGE_BUFFER','StandardDeviation' )
    d['IMAGE_BUFFER_MINIMUM']            = getkey( pvl, 'IMAGE_BUFFER','Minimum' )
    d['IMAGE_BUFFER_MAXIMUM']            = getkey( pvl, 'IMAGE_BUFFER','Maximum' )

    # Calibration Image Dark Reference
    d['CAL_DARK_MEAN']               = getkey( pvl, 'CAL_DARK','Average' )
    d['CAL_DARK_STANDARD_DEVIATION'] = getkey( pvl, 'CAL_DARK','StandardDeviation' )
    d['CAL_DARK_MINIMUM']            = getkey( pvl, 'CAL_DARK','Minimum' )
    d['CAL_DARK_MAXIMUM']            = getkey( pvl, 'CAL_DARK','Maximum' )
   
    # Calibration Image Buffer Area  
    d['CAL_BUFFER_MEAN']               = getkey( pvl, 'CAL_BUFFER','Average' )
    d['CAL_BUFFER_STANDARD_DEVIATION'] = getkey( pvl, 'CAL_BUFFER','StandardDeviation' )
    d['CAL_BUFFER_MINIMUM']            = getkey( pvl, 'CAL_BUFFER','Minimum' )
    d['CAL_BUFFER_MAXIMUM']            = getkey( pvl, 'CAL_BUFFER','Maximum' )
   
    # Calibration Dark Ramp Area
    d['CAL_DARK_RAMP_MEAN']               = getkey( pvl,'CAL_DARK_RAMP','Average' )
    d['CAL_DARK_RAMP_STANDARD_DEVIATION'] = getkey( pvl,'CAL_DARK_RAMP','StandardDeviation' )
    d['CAL_DARK_RAMP_MINIMUM']            = getkey( pvl,'CAL_DARK_RAMP','Minimum' )
    d['CAL_DARK_RAMP_MAXIMUM']            = getkey( pvl,'CAL_DARK_RAMP','Maximum' )
   
    # Image Post Ramp Area
    d['IMAGE_POST_RAMP_MEAN']               = getkey( pvl,'IMAGE_POSTRAMP','Average' )
    d['IMAGE_POST_RAMP_STANDARD_DEVIATION'] = getkey( pvl,'IMAGE_POSTRAMP','StandardDeviation' )
    d['IMAGE_POST_RAMP_MINIMUM']            = getkey( pvl,'IMAGE_POSTRAMP','Minimum' )
    d['IMAGE_POST_RAMP_MAXIMUM']            = getkey( pvl,'IMAGE_POSTRAMP','Maximum' )

    return d


def get_dncnt( cub, hmin, hmax, keep=False ):
    '''Extract DN count from the histogram of a cub file'''

    histfile = os.path.splitext( ofile_cub )[0] + '.hist'
    if not os.path.isfile( histfile ): isis.hist( cub, histfile )

    c = 0
    with open( histfile ) as csvfile:
        headers[]
        for line in csvfile:
            if 'DN,' in line: headers = line.split(',')
        reader = csv.DictReader( csvfile, fieldnames= headers )
        for row in reader:
            if( row['Percent'] >= hmin and row['Percent'] <= hmax ): c += 1

    if not keep: os.remove( histfile )
    return c


def snr( cub, gainsfile, histats ):
    '''Calculate the signal to noise ratio.'''

    summing = isis.getkey( cub, 'Instrument','Summing')
    ccdchan = '{}_{}'.format( hirise.getccdchannel(cub) )

    gainspvl = pvl.load( gainsfile )
    gain = float( gainspvl['Gains'][ccdchan]['Bin'+summing] )

    img_mean   = float( histats['IMAGE_MEAN'] )
    lis_pixels = float( histats['LOW_SATURATED_PIXELS'] )
    buf_mean   = float( histats['IMAGE_BUFFER_MEAN'] )

    snr = -9999
    r = 90 # Note from original file: 
           # 150 e *Changed value to 90 e- 1/31/2012 to bring closer 
           # to HIPHOP value for read noise. SM

    if( 0 == lis_pixels and img_mean > 0.0 and buf_mean > 0 ):
        s = (img_mean - buf_mean) * gain
        snr = s/math.sqrt(s + r*r)
        print( '\nCalculation of Signal/Noise Ratio:' )
        print( '\tIMAGE_MEAN: {}'.format(img_mean) )
        print( '\tIMAGE_BUFFER_MEAN: {}'.format(buf_mean) )
        print( '\tR (electrons/DN): {}'.format(r) )
        print( '\tGain: {}'.format(gain) )
        print( 'Signal/Noise ratio: {}'.format(snr) )

    return snr


if __name__ == "__main__":
    main()
