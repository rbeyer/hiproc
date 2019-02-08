#!/usr/bin/env python
"""This module contains HiRISE utility functions."""

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

import re

def getobsid( s ):
    '''Extracts a HiRISE Observation ID from the string'''
    match = re.search( r"\w{3}_\d{5,6}_\d{4}", s )
    if( match ):
        return match.group()
    else: raise Exception( "Couldn't find an Obs ID in "+s )

def getccdchannel( s ):
    '''Extracts a HiRISE CCD name and channel number from the string'''
    ccd_pattern = r"(RED\d|IR1[0-1]|BG1[2-3])"
    match = re.search( ccd_pattern+r"_([0-1])", s )
    if( match ):
        return match.group(1,2)
    else:
        match = re.search( ccd_pattern, s )
        if( match ):
            return ( match.group(), None )

# def orbit_parser(text):
#     pattern = re.compile(r"([a-zA-Z]{3})_(\d{6})_(\d{4})")
#     match = pattern.search( text )
#     if match:
#         (phase,orbit_number,latesque) = match.groups()
#         ObsID = '_'.join([phase,orbit_number,latesque])
#         return (ObsID, orbit_number, phase)
#     else:
#         raise Exception( "No HiRISE Observation ID identified in "+text )
