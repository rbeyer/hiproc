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


# This program is based on HiCal version 1.32 2016/08/04,
# and on the Perl HiCal program: ($Revision: 1.24 $ $Date: 2016/08/05 18:05:28 $)
# by Eric Eliason 
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

    # GetConfigurationParameters()
    # GetProductFiles() - just gets the names of two products, could easily be explicit inputs
    # PrepareDBStatements()
    #   select from HiCat.EDR_Products
    #   insert HiCat.CCD_Processing_Statistics
    # ProcessingStep() - mostly sets up stuff
    #   FurrowCheck() - reads info from .cubenorm.tab
    # HiStitchStep() - runs HiStitch, and inserts to db
    # HiFurrowStep() - runs external HiFurrow_Fix
    # skip ThumbBrowseStep()



if __name__ == "__main__":
    main()
