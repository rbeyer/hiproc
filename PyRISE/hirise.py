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

import os
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
lat_re = re.compile(r"[0-3]?\d?\d?[05]")

# obsid_re = re.compile(r"("+phase_re.pattern+r")_("+orbit_re.pattern+r")_("+lat_re.pattern+r")")
# obsid_re = re.compile( r"({})_({})_({})".format(phase_re.pattern, orbit_re.pattern, lat_re.pattern) )
obsid_core_re = re.compile(fr"(?P<phase>{phase_re.pattern})?_?(?P<orbit>{orbit_re.pattern})_(?P<latesque>{lat_re.pattern})")
obsid_re = re.compile(fr"(?<!\w){obsid_core_re.pattern}(?!\d)")

ccd_name_re = re.compile(r"RED|IR|BG")
ccd_re = re.compile(r"RED\d|IR1[0-1]|BG1[2-3]")
chan_re = re.compile(r"(?P<channel>[01])")
ccdchan_re = re.compile(fr"(?P<ccd>{ccd_re.pattern})_(?P<channel>[0-1])")
prodid_re = re.compile(fr"(?<!\w){obsid_core_re.pattern}_?(?P<ccd>{ccd_re.pattern})?_?(?P<channel>[0-1])?(?!\d)")


class ObservationID:
    """A class for HiRISE Observation IDs."""

    def __init__(self, *args):
        if len(args) == 1:
            match = obsid_re.search(str(args[0]))
            if match:
                parsed = match.groupdict()
                phase = parsed['phase']
                orbit = parsed['orbit']
                lat = parsed['latesque']
            else:
                raise ValueError('{} did not match regex: {}'.format(args[0], obsid_re.pattern))
        elif len(args) == 2:
            (phase, orbit, lat) = None, *args
        elif len(args) == 3:
            (phase, orbit, lat) = args
        else:
            raise IndexError('accepts 1 to 3 arguments')

        if orbit_re.fullmatch(orbit):
            self.orbit_number = '{:0>6}'.format(orbit)
        else:
            raise ValueError(f'orbit argument, {orbit}, did not match '
                             'regex: {orbit_re.pattern}')

        if phase:
            if phase in phase_names:
                if orbit_in_phase(self.orbit_number, phase):
                    self.phase = phase
                else:
                    raise ValueError('The orbit, {}, is outside the allowed '
                                     'range ({}, {}) in HiRISE mission phase '
                                     '{}'.format(self.orbit_number,
                                                 phase_max_orbit[prev_phase(phase)],
                                                 phase_max_orbit[phase],
                                                 phase))
            else:
                raise ValueError(f'phase argument, {phase}, did not match '
                                 'any known HiRISE phase:' + str(phase_names))
        else:
            self.phase = getphase(orbit)

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

        ccdname = None
        ccdnumber = None
        chan = None

        items = list(args)

        if len(items) == 1:
            match = prodid_re.search(str(items[0]))
            if match:
                parsed = match.groupdict()
                items = (parsed['phase'], parsed['orbit'], parsed['latesque'])
                if parsed['ccd']:
                    (ccdname, ccdnumber) = getccdnamenumber(parsed['ccd'])
                chan = parsed['channel']
            else:
                raise ValueError('{} did not match regex: '
                                 '{}'.format(args[0], prodid_re.pattern))
        elif len(items) == 2:
            pass
        elif len(items) == 3:
            if items[0] not in phase_names:
                (ccdname, ccdnumber) = getccdnamenumber(items.pop())
        elif len(items) == 4:
            if items[0] not in phase_names:
                (ccdname, ccdnumber, chan) = _namenumchan(items[-2], items[-1])
                items.pop()
                items.pop()
            else:
                (ccdname, ccdnumber) = getccdnamenumber(items.pop())
        elif len(items) == 5:
            if items[0] not in phase_names:
                chan = _match_chan(items.pop())
                ccdnumber = _match_num(items.pop())
                ccdname = getccdname(items.pop())
            else:
                (ccdname, ccdnumber, chan) = _namenumchan(items[-2], items[-1])
                items.pop()
                items.pop()
        elif len(items) == 6:
            chan = _match_chan(items.pop())
            ccdnumber = _match_num(items.pop())
            ccdname = getccdname(items.pop())
        else:
            raise IndexError('accepts 1 to 6 arguments')
        super().__init__(*items)

        self.ccdname = ccdname
        self.ccdnumber = ccdnumber
        self.channel = chan

    def __str__(self):
        if self.ccdname and self.ccdnumber:
            ccd = self.ccdname + self.ccdnumber
            if self.channel:
                ccd += '_' + self.channel
            return '_'.join([self.phase, self.orbit_number, self.latesque,
                            ccd])
        else:
            return super().__str__()

    def __repr__(self):
        return (f'{self.__class__.__name__}(\'{self.__str__()}\')')


def _match_chan(s: str):
    matched = chan_re.fullmatch(s)
    if matched:
        return matched.group()
    else:
        return None


def _match_num(s: str):
    matched = re.fullmatch(r"1?\d", s)
    if matched:
        return matched.group()
    else:
        return None


def _namenumchan(one: str, two: str) -> tuple:
    ccdname = None
    ccdnumber = None
    chan = None
    try:
        (ccdname, ccdnumber) = getccdnamenumber(one)
        chan = _match_chan(two)
    except ValueError:
        ccdname = getccdname(one)
        ccdnumber = _match_num(two)
    return (ccdname, ccdnumber, chan)


def parseObsID(s: str) -> tuple:
    '''Parses a string to find the first occurence of a valid ObsID'''
    match = obsid_re.search(str(s))
    if(match):
        return match.groups()
    else:
        raise ValueError('{} did not match regex: {}'.format(s, obsid_re.pattern))


def parseProdID(s: str) -> tuple:
    '''Parses a string to find the first occurence of a valid ProductID'''
    pid = ProductID(s)
    items = [pid.phase, pid.orbit_number, pid.latesque,
             pid.ccdname, pid.ccdnumber, pid.channel]
    return tuple(filter(lambda x: x is not None, items))


def prev_phase(s: str) -> str:
    '''Returns the name of the previous HiRISE mission phase.'''
    prev_index = phase_names.index(s) - 1
    if prev_index < 0:
        raise IndexError(f'HiRISE mission phase, {s}, is the first one, '
                         'it has no previous phase.')
    else:
        return phase_names[prev_index]


def orbit_in_phase(orbit, phase: str):
    '''Determines whether the given orbit is within the given HiRISE
    mission phase.'''
    previous = prev_phase(phase)
    if int(orbit) <= phase_max_orbit[previous]:
        if(phase_max_orbit[phase] is not None and
           phase_max_orbit[phase] == phase_max_orbit[previous] and
           phase_max_orbit[phase] == int(orbit)):
            return True
        return False
    if phase_max_orbit[phase] is not None and int(orbit) > phase_max_orbit[phase]:
        return False
    return True


def getObsID(s: str) -> str:
    '''Extracts a HiRISE Observation ID from the string'''
    return str(ObservationID(s))


def getProdID(s: str) -> str:
    '''Extracts a HiRISE Product ID from the string.'''
    return str(ProductID(s))


def getphase(orbit) -> str:
    '''Returns the name of the HiRISE mission phase based on the MRO orbit
    number.'''
    last_phase = None
    for (k, v) in reversed(list(phase_max_orbit.items())):
        if v is not None and int(orbit) > v:
            if last_phase is not None:
                return last_phase
            else:
                raise IndexError(f'The orbit number {orbit} is greater than '
                                 'the last HiRISE orbit.')
        last_phase = k
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
        return re.search(r"\d{1,2}", match.group()).group()
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
