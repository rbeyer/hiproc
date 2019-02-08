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

# Create some compiled regex Patterns to use in this module.
phase_re = re.compile(r"\w{3}")
orbit_re = re.compile(r"\d{5,6}")
lat_re   = re.compile(r"\d{4}")

#obsid_re = re.compile(r"("+phase_re.pattern+r")_("+orbit_re.pattern+r")_("+lat_re.pattern+r")")
#obsid_re = re.compile( r"({})_({})_({})".format(phase_re.pattern, orbit_re.pattern, lat_re.pattern) )
obsid_re = re.compile( fr"(?<!\w)({phase_re.pattern})_({orbit_re.pattern})_({lat_re.pattern})(?!\d)" )

ccd_re = re.compile( r"(RED\d|IR1[0-1]|BG1[2-3])" )
ccdchan_re = re.compile( fr"{ccd_re.pattern}_([0-1])" )


class ObservationID:
    """A class for HiRISE Observation IDs."""

    def __init__(self, *args):
        if len(args) < 1 or len(args) > 3: raise IndexError('accepts 1 to 3 arguments')
        if len(args) == 1:
            (self.phase, self.orbit_number, self.latesque) = parseObsID(args[0])
        else: 
          if phase_re.fullmatch(args[0]): self.phase = args[0]
          else: raise ValueError(f'phase argument, {args[0]}, did not match regex: {phase_re.pattern}')

          if orbit_re.fullmatch(args[1]): self.orbit_number = args[1]
          else: raise ValueError(f'orbit argument, {args[0]}, did not match regex: {orbit_re.pattern}')

          if lat_re.fullmatch(args[2]): self.latesque = args[2]
          else: raise ValueError(f'latitude argument, {args[2]}, did not match regex: {lat_re.pattern}')

    def __str__(self):
        return f'{self.phase}_{self.orbit_number}_{self.latesque}'

    def __repr__(self):
        return (f'{self.__class__.__name__}('
                f'\'{self.phase}_{self.orbit_number}_{self.latesque}\')')

def parseObsID(s):
    '''Parses a string to find the first occurence of a valid ObsID'''
    match = obsid_re.search( str(s) )
    if(match):
        return match.groups()
    else: raise ValueError('{} did not match regex: {}'.format(s,obsid_re.pattern))

def getObsID( s ):
    '''Extracts a HiRISE Observation ID from the string'''
    return '_'.join( parseObsID(s) )

def getccd( s ):
    '''Extracts a HiRISE CCD name from the string'''
    match = ccd_re.search( s )
    if match: return match.group()
    else:     raise ValueError(f'{s} did not match regex: {ccd_re.pattern}')

def getccdchannel( s ):
    '''Extracts a HiRISE CCD name and channel number from the string'''
    match = ccdchan_re.search( s )
    if match: return match.groups()
    else:     raise ValueError(f'{s} did not match regex: {ccdchan_re.pattern}')
