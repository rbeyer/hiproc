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
from itertools import repeat

phase_names = ('INT', 'CAL', 'ATL', 'KSC', 'LAU', 'CRU', 'APR', 'AEB',
               'TRA', 'PSP',  # 'REL'
               'ESP')


phase_max_orbit = dict(zip(phase_names, repeat(0)),
                       TRA=1000, PSP=11248, ESP=None)

# Create some compiled regex Patterns to use in this module.
phase_re = re.compile("|".join(phase_names))
orbit_re = re.compile(r"\d{1,6}")
lat_re = re.compile(r"\d{1,4}")

# obsid_re = re.compile(r"("+phase_re.pattern+r")_("+orbit_re.pattern+r")_("+lat_re.pattern+r")")
# obsid_re = re.compile( r"({})_({})_({})".format(phase_re.pattern, orbit_re.pattern, lat_re.pattern) )
obsid_re = re.compile(fr"(?<!\w)(?P<phase>{phase_re.pattern})?_?(?P<orbit>{orbit_re.pattern})_(?P<latesque>{lat_re.pattern})(?!\d)")

ccd_name_re = re.compile(r"RED|IR|BG")
ccd_re = re.compile(r"RED\d|IR1[0-1]|BG1[2-3]")
ccdchan_re = re.compile(fr"({ccd_re.pattern})_([0-1])")
prodid_re = re.compile(fr"{obsid_re.pattern}_(?P<ccd>{ccd_re.pattern})?_(?P<channel>[0-1])?(?!\d)")


class ObservationID:
    """A class for HiRISE Observation IDs."""

    def __init__(self, *args):
        if len(args) == 1:
            parsed = obsid_re.search(str(s)).groupdict()
            phase = parsed['phase']
            orbit = parsed['orbit']
            lat = parsed['latesque']
        elif len(args) == 2:
            (phase, orbit, lat) = None, args
        elif len(args) == 3:
            (phase, orbit, lat) = args
        else:
            raise IndexError('accepts 1 to 3 arguments')

        if phase:
            if phase in phase_names:
                self.phase = phase
            else:
                raise ValueError(f'phase argument, {phase}, did not match '
                                 'any known HiRISE phase:' + str(phase_names))
        else:
            self.phase = getphase(orbit)

        if orbit_re.fullmatch(orbit):
            self.orbit_number = '{:0>6}'.format(orbit)
        else:
            raise ValueError(f'orbit argument, {orbit}, did not match '
                             'regex: {orbit_re.pattern}')

        if lat_re.fullmatch(lat):
            self.latesque = '{:0>4}'.format(lat)
        else:
            raise ValueError(f'latitude argument, {lat}, did not match '
                             'regex: {lat_re.pattern}')

    def __str__(self):
        return f'{self.phase}_{self.orbit_number}_{self.latesque}'

    def __repr__(self):
        return (f'{self.__class__.__name__}('
                f'\'{self.phase}_{self.orbit_number}_{self.latesque}\')')


class ProductID(ObservationID):
    """A class for HiRISE Product IDs."""

    def __init__(self, *args):
        self.ccdname = None
        self.ccdnumber = None
        self.channel = None
        if len(args) == 1:
            parsed = prodid_re.search(str(s)).groupdict()
            phase = parsed['phase']
            orbit = parsed['orbit']
            lat = parsed['latesque']
            ccd = parsed['ccd']
            chan = parsed['channel']
        elif len(args) == 2:
            super().__init__(args)
        elif len(args) == 3 or len(args) == 4:
            if args[0] in phase_names:
                try:
                    ccd = getccd(args[-1])
                    args.pop()
                except ValueError:
                    pass
            last = args.pop()

            try:
                ccd = getccd(last)
            except ValueError:
                ccd = getccd(args.pop())

            super().__init__(args)

            ccd = getccd(args[2])
            self.ccdname = getccdname(ccd)
            self.ccdnumber = getccdnumber(ccd)

            if len(args) == 4:
                pass
        elif len(args) == 5:
            pass
        elif len(args) == 6:
            pass
        else:
            raise IndexError('accepts 1 to 6 arguments')

    def __str__(self):
        if self.ccdname and self.ccdnumber:
            ccd = self.ccdname + self.ccdnumber
            return '_'.join(self.phase, self.orbit_number, self.latesque,
                            ccd, self.channel)
        else:
            return super().__str__()

    def __repr__(self):
        return (f'{self.__class__.__name__}(\'{self.__str__()}\')')


def parseObsID(s: str) -> tuple:
    '''Parses a string to find the first occurence of a valid ObsID'''
    match = obsid_re.search(str(s))
    if(match):
        return match.groups()
    else:
        raise ValueError('{} did not match regex: {}'.format(s, obsid_re.pattern))


def parseProdID(s: str) -> tuple:
    '''Parses a string to find the first occurence of a valid ProductID'''
    match = prodid_re.search(str(s))
    if(match):
        if len(match) == 3:
            return match.groups()
        elif len(match) == 4:
            if match.group(1) in phase_names:
                raise ValueError('{} did not match regex: {}'.format(s, prodid_re.pattern))
            else:
                (name, number) = re.match(fr"({ccd_name_re.pattern})(\d+)", match.group(3))
                return (match.group[1:2], name, number, group(4))
        elif len(match) == 5:
            (name, number) = re.match(fr"({ccd_name_re.pattern})(\d+)", match.group(4))
            return (match.group[1:3], name, number, group(5))
        else:
            raise RunTimeError('Did not expect this many match elements.')
    else:
        raise ValueError('{} did not match regex: {}'.format(s, prodid_re.pattern))


def getObsID(s: str) -> str:
    '''Extracts a HiRISE Observation ID from the string'''
    return str(ObservationID(s))


def getProdID(s: str) -> str:
    '''Extracts a HiRISE Product ID from the string.'''
    return str(ProductID(s))


def getphase(orbit: int) -> str:
    '''Returns the name of the HiRISE mission phase based on the MRO orbit
    number.'''
    last_phase = None
    for (k, v) in reversed(list(phase_max_orbit.items())):
        if v is None:
            last_phase = k
            continue
        if orbit > v:
            if last_phase is not None:
                return last_phase
            else:
                raise IndexError(f'The orbit number {orbit} is greater than '
                                 'the last HiRISE orbit.')
    else:
        raise IndexError(f'The orbit number {orbit} is not a HiRISE orbit '
                         'number.')


def getccd(s: str) -> str:
    '''Extracts a HiRISE CCD from the string (would get
    'RED5' from 'PSP_010502_2090_RED5_0.EDR_Stats.cub').'''
    match = ccd_re.search(s)
    if match:
        return match.group()
    else:
        raise ValueError(f'{s} did not match regex: {ccd_re.pattern}')


def getccdname(s: str) -> str:
    '''Extracts the name of a HiRISE CCD from the string (would get
    'RED' from 'PSP_010502_2090_RED5_0.EDR_Stats.cub').'''
    match = ccd_name_re.search(s)
    if match:
        return match.group()
    else:
        raise ValueError(f'{s} did not match regex: {ccd_name_re.pattern}')


def getccdnumber(s: str) -> str:
    '''Extracts the number of a HiRISE CCD from the string (would get
    '5' from 'PSP_010502_2090_RED5_0.EDR_Stats.cub').'''
    match = ccd_re.search(s)
    if match:
        return re.match(fr"{ccd_name_re.pattern}(\d+)", match.group())
    else:
        raise ValueError(f'{s} did not match regex: {ccd_re.pattern}')


def getccdnamenumber(s: str) -> tuple:
    '''Extracts the name and number of a HiRISE CCD from the string (would get
    '('RED', '5')' from 'PSP_010502_2090_RED5_0.EDR_Stats.cub').'''
    nn = getccd(s)
    return(getccdname(nn), getccdnumber(nn))


def getccdchannel(s: str) -> tuple:
    '''Extracts a HiRISE CCD and channel number from the string (would get
    ('RED5', '1') from 'PSP_010502_2090_RED5_0.EDR_Stats.cub').'''
    match = ccdchan_re.search(s)
    if match:
        return match.groups()
    else:
        raise ValueError(f'{s} did not match regex: {ccdchan_re.pattern}')


def get_ObsID_fromfile(path: os.PathLike) -> ObservationID:
    '''Reads the file to get the ObservationID, if an ISIS cube,
       otherwise parses the filepath.'''
    p = Path(path)

    try:
        return getObsID(isis.getkey_k(p, 'Archive', 'ObservationId'))
    except CalledProcessError:
        for part in reversed(p.parts):
            try:
                return getObsID(part)
            except ValueError:
                continue
        raise ValueError('Could not extract a HiRISE Observation ID from ' + path)
