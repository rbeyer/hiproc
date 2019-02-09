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

    # Read some stuff from HiCal.conf, could be done with pvl
    # Setup00 builds data structures, seems long and painful, need?
    # Setup01 - don't need - about output data routeing
    # Setup02 - pulls from db EDR_Products, written by EDR_Stats
    # Setup03 - gets binning information from Planned_Observations to find binning for ObsID
    # ProcessingSwitches() sets a variety of booleans
    #
    # above is setup
    #
    # SetHiCalVersion() we don't need to record the version info of isis.hical
    # Furrow() - yeah, need to do this
    # Mask() - yeah
    # HiCal()
    # Lpfz()
    # GainDrift()
    # CubeNorm() - writes stuff back to EDR_Products
    # NoiseFilter()
    # Hidestripe() - writes to EDR_Products



if __name__ == "__main__":
    main()
