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

def main():
    parser = argparse.ArgumentParser( description=__doc__ )
    parser.add_argument('--db',          required=False, default='HiCat.db' )
    parser.add_argument('-o','--output', required=False, default='.HiCal.cub')
    parser.add_argument('-k','--keep',   required=False, default=False )
    parser.add_argument('cube', metavar=".cub-file" )

    args = parser.parse_args()

    ofile_cub = ''
    if args.output.startswith('.'): 
        ofile_cub = os.path.splitext( args.img )[0] + args.output
    else: ofile_cub = args.output

    # GetConfigurationParameters() - Read some stuff from HiCal.conf, could be done with pvl
    # Setup00 builds data structures, seems long and painful, need?
    # Setup01 - don't need - about output data routeing
    # Setup02 - db
    #   select from HiCat.EDR_Products, written by EDR_Stats
    # Setup03 - gets binning information to find binning for ObsID ... hmmm, alternate?
    #   select from HiCat.Planned_Observations
    # ProcessingSwitches() sets a variety of booleans

    image_mean # select IMAGE_MEAN from HiCat.EDR_Products 
    binning # select BINNING from HiCat.EDR_Products
    
    ### setup done ###

    if( binning > 1 and image_mean > 7000.0):
        furrow_nulling( args.cube )

    # Mask()
    #   isis.mask
    #   isis.cubenorm
    #   isis.mask
    #   AnalyzeCubenormStats()
    # HiCal() runs isis.hical
    # Lpfz() runs isis.lowpass
    # GainDrift() - runs external HiGainFx program
    # CubeNorm() - isis.crop, isis.cubenorm, external Cubenorm_Filter
    #   insert into HiCat.EDR_Products
    #   isis.cubenorm again
    # NoiseFilter() - external Noise_Filter
    # Hidestripe() - isis.[hidestripe,hipass,lowpass,algebra] 
    #   insert into HiCat.EDR_Products

def furrow_nulling( cube ):
    # Furrow() - sets filenames, run isis.mask

    furrow_fix_file = os.path.splitext( cube )[0] + '.furrowfix.cub'
    isis.mask( cube, mask=cube, to=furrow_fix_file, min=1000000, max=1000000 )

    furrow_values = furrow_setup()

    #   FurrowSetup() - contains and then sets furrow_values
    #   isis.crop
    #   isis.handmos
    #   isis.crop and isis.mask
    #   isis.stats, isis.handmos, isis.trimfilter


def furrow_setup()
    '''Just builds this dict of dicts which holds tuples.'''
    # Expect to call it like this: d[ccd][binning]
    # Assume each channel for each CCD will have the same threshold
    d['RED0']['2'] = 8000, 8100, 8700,  9200,  9600, 10000, 12000, 12000, 12000, 12000
    d['RED0']['4'] = 8000, 9000, 9500,  9900,  9900, 10000
    d['RED1']['2'] = 7200, 7200, 7800,  8400,  9000,  9500, 12000, 12000, 12000, 12000
    d['RED1']['4'] = 8000, 8100, 9200,  9600,  9800, 10000
    d['RED2']['2'] = 7800, 7800, 8400,  9000,  9600, 10000, 12000, 12000, 12000, 12000
    d['RED2']['4'] = 8000, 8700, 9500,  9800,  9900, 10000
    d['RED3']['2'] = 7800, 8100, 8300,  9200,  9600, 10000, 12000, 12000, 12000, 12000
    d['RED3']['4'] = 7900, 9200, 9700,  9900, 10000, 10500
    d['RED4']['2'] = 7800, 7800, 8300,  9000,  9500,  9900, 12000, 12000, 12000, 12000
    d['RED4']['4'] = 8000, 8700, 9700, 10000, 10300, 10600
    d['RED5']['2'] = 7900, 8200, 8600,  9200,  9600, 10000, 12000, 12000, 12000, 12000
    d['RED5']['4'] = 8000, 9300, 9700,  9900, 10200, 10700
    d['RED6']['2'] = 7500, 7500, 8100,  8500,  9200, 10000, 12000, 12000, 12000, 12000
    d['RED6']['4'] = 8000, 8400, 9700, 10000, 10500, 10700
    d['RED7']['2'] = 7600, 8300,  8900,  9400, 9900, 11000, 12000, 12000, 12000, 12000
    d['RED7']['4'] = 7700, 9600, 10000, 10200, 11000, 12000
    d['RED8']['2'] = 7200, 7200, 7900,  8500,  9000,  9400, 12000, 12000, 12000, 12000
    d['RED8']['4'] = 7700, 9600, 10000, 10200, 11000, 12000
    d['RED9']['2'] = 7600, 8300, 8600,  9200,  9600, 10000, 12000, 12000, 12000, 12000
    d['RED9']['4'] = 8000, 8800, 9200,  9400,  9800, 10500
    d['IR10']['2'] = 8000, 8100, 8700,  9200,  9600, 10000, 12000, 12000, 12000, 12000
    d['IR10']['4'] = 7600, 8300, 9000, 10000, 10500, 12000
    d['IR11']['2'] = 8000, 8100, 8700,  9200,  9600, 10000, 12000, 12000, 12000, 12000
    d['IR11']['4'] = 7600, 8300, 9000, 10000, 10500, 12000
    d['BG12']['2'] = 8000, 8100, 8700,  9200,  9600, 10000, 12000, 12000, 12000, 12000
    d['BG12']['4'] = 7600, 8300, 9000, 10000, 10500, 12000
    d['BG13']['2'] = 8000, 8100, 8700,  9200,  9600, 10000, 12000, 12000, 12000, 12000
    d['BG13']['4'] = 7600, 8300, 9000, 10000, 10500, 12000
    return d

if __name__ == "__main__":
    main()
