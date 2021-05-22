#!/usr/bin/env python
"""This module contains HiRISE utility functions."""

# Copyright 2020, Ross A. Beyer (rbeyer@seti.org)
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
import subprocess
from itertools import repeat
from pathlib import Path

import kalasiris as isis

phase_names = (
    "INT",
    "CAL",
    "ATL",
    "KSC",
    "LAU",
    "CRU",
    "APR",
    "AEB",
    "TRA",
    "PSP",  # 'REL'
    "ESP",
)

phase_max_orbit = dict(
    zip(phase_names, repeat(0)), TRA=1000, PSP=11248, ESP=None
)

ccd_numbers = dict(
    RED=(0, 1, 2, 3, 4, 5, 6, 7, 8, 9), IR=(10, 11), BG=(12, 13)
)

# Create some compiled regex Patterns to use in this module.
phase_re = re.compile("|".join(phase_names))
orbit_re = re.compile(r"\d{1,6}")
lat_re = re.compile(r"[0-3]?\d?\d?[05]")
target_re = re.compile(fr"{lat_re.pattern}|9(?:[0-2]\d\d|30[0-3])")

obsid_core_re = re.compile(
    fr"(?P<phase>{phase_re.pattern})?_?"
    fr"(?P<orbit>{orbit_re.pattern})_"
    fr"(?P<target>{target_re.pattern})"
)
obsid_re = re.compile(fr"(?<!\w){obsid_core_re.pattern}(?!\d)")

ccd_name_re = re.compile(r"RED|IR|BG")
ccd_re = re.compile(r"RED\d|IR1[0-1]|BG1[2-3]")
ccdid_core_re = re.compile(
    fr"{obsid_core_re.pattern}_(?P<ccd>{ccd_re.pattern})"
)
ccdid_re = re.compile(fr"(?<!\w){ccdid_core_re.pattern}(?!\d)")

chan_re = re.compile(r"(?P<channel>[01])")
ccdchan_re = re.compile(fr"(?P<ccd>{ccd_re.pattern})_(?P<channel>[0-1])")
chanid_core_re = re.compile(fr"{ccdid_core_re.pattern}_(?P<channel>[0-1])")
chanid_re = re.compile(fr"(?<!\w){chanid_core_re.pattern}(?!\d)")

prodid_re = re.compile(
    fr"(?<!\w){obsid_core_re.pattern}_?"
    fr"(?P<ccd>{ccd_re.pattern})?_?"
    fr"(?P<channel>[0-1])?(?!\d)"
)


class ObservationID:
    """A class for HiRISE Observation IDs.

    :ivar phase: The three-letter string indicating the mission phase
        of this Observation ID.
    :ivar orbit_number: The six-digit orbit number (with leading zeros)
        as a string.
    :ivar target: The four digit code (with leading zeroes) is defined as
        follows: Values between 0000 and 3595, inclusive, reflect
        the central latitude of the observation to within a
        half-degree, multiplied by 10.  (0905, for example, is thus
        a latitude of 90.5 degrees.) 0.0 degrees is the night-side
        equator, which is where the orbit number increments. 90.0
        south pole on the ascending pass. 180.0 is the day-side
        equator on the ascending pass. 270.0 is the north pole on
        the descending pass. Text values between 9000 and 9303,
        inclusive, represent off-planet targets such as Phobos,
        Deimos, or stars.
    """

    def __init__(self, *args):
        if len(args) == 1:
            if isinstance(args[0], ObservationID):
                phase = args[0].phase
                orbit = args[0].orbit_number
                target = args[0].target
            else:
                match = obsid_re.search(str(args[0]))
                if match:
                    parsed = match.groupdict()
                    phase = parsed["phase"]
                    orbit = parsed["orbit"]
                    target = parsed["target"]
                else:
                    raise ValueError(
                        f"{args[0]} did not match regex: {obsid_re.pattern}"
                    )
        elif len(args) == 2:
            (phase, orbit, target) = None, *args
        elif len(args) == 3:
            (phase, orbit, target) = args
        else:
            raise IndexError("accepts 1 to 3 arguments")

        self.orbit_number = self.format_orbit(orbit)
        self.target = self.format_target(target)

        if phase is not None:
            if phase in phase_names:
                if is_orbit_in_phase(self.orbit_number, phase):
                    self.phase = phase
                else:
                    raise ValueError(
                        "The orbit, {}, is outside the allowed "
                        "range ({}, {}) in HiRISE mission phase "
                        "{}".format(
                            self.orbit_number,
                            phase_max_orbit[prev_phase(phase)],
                            phase_max_orbit[phase],
                            phase,
                        )
                    )
            else:
                raise ValueError(
                    f"phase argument, {phase}, did not match "
                    f"any known HiRISE phase: {phase_names}"
                )
        else:
            self.phase = get_phase(orbit)

    def __str__(self):
        return f"{self.phase}_{self.orbit_number}_{self.target}"

    def __repr__(self):
        return (
            f"{self.__class__.__name__}("
            f"'{self.phase}_{self.orbit_number}_{self.target}')"
        )

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return (
                self.phase == other.phase
                and self.orbit_number == other.orbit_number
                and self.target == other.target
            )
        return False

    def __lt__(self, other):
        if isinstance(other, self.__class__):
            return (
                phase_names.index(self.phase),
                int(self.orbit_number),
                int(self.target),
            ) < (
                phase_names.index(other.phase),
                int(other.orbit_number),
                int(other.target),
            )
        else:
            return NotImplemented

    def __hash__(self):
        return hash((self.phase, self.orbit_number, self.target))

    @staticmethod
    def format_orbit(orbit) -> str:
        if orbit_re.fullmatch(str(orbit)):
            return "{:0>6}".format(orbit)
        else:
            raise ValueError(
                f"orbit argument, {orbit}, did not match "
                f"regex: {orbit_re.pattern}"
            )

    @staticmethod
    def format_target(val) -> str:
        if target_re.fullmatch(str(val)):
            return "{:0>4}".format(val)
        else:
            raise ValueError(
                f"target argument, {val}, did not match "
                f"regex: {target_re.pattern}"
            )


class CCDID(ObservationID):
    """A class for HiRISE CCD IDs.

    :ivar ccdname: The CCD name: 'RED', 'IR', or 'BG'.
    :ivar ccdnumber: The CCD number as a string: 0 through 13.
    """

    def __init__(self, *args):
        items = list(args)

        if len(items) == 1:
            if isinstance(args[0], CCDID):
                items = (args[0].phase, args[0].orbit_number, args[0].target)
                ccdname = args[0].ccdname
                ccdnumber = args[0].ccdnumber
            else:
                match = ccdid_re.search(str(items[0]))
                if match:
                    parsed = match.groupdict()
                    items = (
                        parsed["phase"],
                        parsed["orbit"],
                        parsed["target"],
                    )
                    (ccdname, ccdnumber) = get_ccdnamenumber(parsed["ccd"])
                else:
                    raise ValueError(
                        f"Could not construct a CCDID. {args[0]} did not match "
                        f"regex: {ccdid_re.pattern}"
                    )
        elif len(items) == 3:
            try:
                (ccdname, ccdnumber) = get_ccdnamenumber(items.pop())
            except ValueError as err:
                raise ValueError(
                    "Could not construct a CCDID. " + str(err)
                ) from err
        elif len(items) == 4:
            if items[0] in phase_names:
                (ccdname, ccdnumber) = get_ccdnamenumber(items.pop())
            else:
                ccdnumber = _match_num(items.pop())
                ccdname = get_ccdname(items.pop())
        elif len(items) == 5:
            ccdnumber = _match_num(items.pop())
            ccdname = get_ccdname(items.pop())
        else:
            raise IndexError(
                "CCDID accepts 1, 3, 4, or 5 arguments, "
                f"{len(items)} were provided."
            )
        super().__init__(*items)

        self.ccdname = ccdname
        self.ccdnumber = ccdnumber

    def __str__(self):
        return "_".join([super().__str__(), self.get_ccd()])

    def __repr__(self):
        return f"{self.__class__.__name__}('{self.__str__()}')"

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return (
                super().__eq__(other)
                and self.ccdname == other.ccdname
                and self.ccdnumber == other.ccdnumber
            )
        return False

    def __lt__(self, other):
        if isinstance(other, self.__class__):
            if super().__eq__(other):
                return int(self.ccdnumber) < int(other.ccdnumber)
            else:
                return super().__lt__(other)
        else:
            return NotImplemented

    def __hash__(self):
        return hash((super().__hash__(), self.ccdname, self.ccdnumber))

    def get_ccd(self) -> str:
        return self.ccdname + self.ccdnumber

    def get_obsid(self) -> ObservationID:
        return ObservationID(self.phase, self.orbit_number, self.target)


class ChannelID(CCDID):
    """A class for HiRISE Channel IDs.

    :ivar channel: The Channel number (0 or 1) as a string.
    """

    def __init__(self, *args):
        items = list(args)

        if len(items) == 1:
            match = chanid_re.search(str(items[0]))
            if match:
                parsed = match.groupdict()
                items = (parsed["phase"], parsed["orbit"], parsed["target"])
                items += get_ccdnamenumber(parsed["ccd"])
                chan = parsed["channel"]
            else:
                raise ValueError(
                    "{} did not match regex: "
                    "{}".format(args[0], chanid_re.pattern)
                )
        elif 2 <= len(items) <= 6:
            maybechan = items.pop()
            try:
                chan = self.format_chan(maybechan)
            except ValueError as err:
                raise ValueError("The last item, " + str(err))

            if len(items) == 1:
                if isinstance(items[0], CCDID):
                    items = [
                        items[0].phase,
                        items[0].orbit_number,
                        items[0].target,
                        items[0].ccdname,
                        items[0].ccdnumber,
                    ]
                else:
                    raise ValueError(
                        "{} was not a CCDID object".format(items[0])
                    )

            elif len(items) == 2:
                (ccdname, ccdnumber) = get_ccdnamenumber(items.pop())
                if isinstance(items[0], ObservationID):
                    items = [
                        items[0].phase,
                        items[0].orbit_number,
                        items[0].target,
                        ccdname,
                        ccdnumber,
                    ]
                else:
                    raise ValueError(
                        "{} was not a ObservationID object".format(items[0])
                    )

        else:
            raise IndexError("accepts 1 to 6 arguments")
        super().__init__(*items)

        self.channel = str(chan)

    def __str__(self):
        return "_".join([super().__str__(), self.channel])

    def __repr__(self):
        return f"{self.__class__.__name__}('{self.__str__()}')"

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return super().__eq__(other) and self.channel == other.channel
        return False

    def __lt__(self, other):
        if isinstance(other, self.__class__):
            if super().__eq__(other):
                return int(self.channel) < int(other.channel)
            else:
                return super().__lt__(other)
        else:
            return NotImplemented

    def __hash__(self):
        return hash((super().__hash__(), self.channel))

    @staticmethod
    def format_chan(chan) -> str:
        if isinstance(chan, int):
            if chan == 0 or chan == 1:
                return chan
            else:
                raise ValueError(f"{chan}, is not zero or one.")
        else:
            matched = chan_re.fullmatch(chan)
            if matched:
                return matched.group()
            else:
                raise ValueError(
                    f"{chan}, did not match a " "channel {chan_re.pattern}"
                )


# def _match_chan(s: str):
#     matched = chan_re.fullmatch(s)
#     if matched:
#         return matched.group()
#     else:
#         return None


def _match_num(s: str):
    matched = re.fullmatch(r"1?\d", s)
    if matched:
        return matched.group()
    else:
        return None


# def _namenumchan(one: str, two: str) -> tuple:
#     ccdname = None
#     ccdnumber = None
#     chan = None
#     try:
#         (ccdname, ccdnumber) = get_ccdnamenumber(one)
#         chan = _match_chan(two)
#     except ValueError:
#         ccdname = get_ccdname(one)
#         ccdnumber = _match_num(two)
#     return (ccdname, ccdnumber, chan)


def parseObsID(s: str) -> tuple:
    """Parses a string to find the first occurence of a valid ObsID.

    Returns a tuple of strings that were matched.
    """
    match = obsid_re.search(str(s))
    if match:
        return match.groups()
    else:
        raise ValueError(
            "{} did not match regex: {}".format(s, obsid_re.pattern)
        )


def prev_phase(s: str) -> str:
    """Returns the name of the previous HiRISE mission phase."""
    prev_index = phase_names.index(s) - 1
    if prev_index < 0:
        raise IndexError(
            f"HiRISE mission phase, {s}, is the first one, "
            "it has no previous phase."
        )
    else:
        return phase_names[prev_index]


def is_orbit_in_phase(orbit, phase: str):
    """Determines whether the given orbit is within the given HiRISE
    mission phase."""
    previous = prev_phase(phase)
    if int(orbit) <= phase_max_orbit[previous]:
        if (
            phase_max_orbit[phase] is not None
            and phase_max_orbit[phase] == phase_max_orbit[previous]
            and phase_max_orbit[phase] == int(orbit)
        ):
            return True
        return False
    if (
        phase_max_orbit[phase] is not None
        and int(orbit) > phase_max_orbit[phase]
    ):
        return False
    return True


def ObsIDstr(s: str) -> str:
    """Extracts a HiRISE Observation ID string from the given string."""
    return str(ObservationID(s))


def get_phase(orbit) -> str:
    """Returns the name of the HiRISE mission phase based on the MRO orbit
    number.

    The given orbit must be an int or something that can be converted
    to an int via int().
    """
    last_phase = None
    for (k, v) in reversed(list(phase_max_orbit.items())):
        if v is not None and int(orbit) > v:
            if last_phase is not None:
                return last_phase
            else:
                raise IndexError(
                    f"The orbit number {orbit} is greater than "
                    "the last HiRISE orbit."
                )
        last_phase = k
    else:
        raise IndexError(
            f"The orbit number {orbit} is not a HiRISE orbit " "number."
        )


def get_ccd(s: str) -> str:
    """Extracts a HiRISE CCD from the string (would get
    'RED5' from 'PSP_010502_2090_RED5_0.EDR_Stats.cub')."""
    match = ccd_re.search(s)
    if match:
        return match.group()
    else:
        raise ValueError(f"{s} did not match regex: {ccd_re.pattern}")


def _getccdname_fromint(ccdnum: int) -> str:
    for k, v in ccd_numbers.items():
        if ccdnum in v:
            return k
    else:
        raise ValueError(f"The value {ccdnum} must be a value from 0 to 13.")


def get_ccdname(item) -> str:
    """Extracts the name of a HiRISE CCD from the string (would get
    'RED' from 'PSP_010502_2090_RED5_0.EDR_Stats.cub').  If an integer
    is given instead of a string, it will provide the corresponding name.
    """
    try:
        match = ccd_name_re.search(item)
        if match:
            return match.group()
        else:
            try:
                return _getccdname_fromint(int(item))
            except ValueError:
                raise ValueError(
                    f"{item} did not match regex: {ccd_name_re.pattern}"
                )
    except TypeError:
        try:
            if isinstance(item, int):
                return _getccdname_fromint(int(item))
            else:
                raise TypeError("Expected string or integer object.")
        except ValueError:
            raise


def get_ccdnumber(s: str) -> str:
    """Extracts the number of a HiRISE CCD from the string (would get
    '5' from 'PSP_010502_2090_RED5_0.EDR_Stats.cub')."""
    match = ccd_re.search(s)
    if match:
        return re.search(r"\d{1,2}", match.group()).group()
    else:
        raise ValueError(f"{s} did not match regex: {ccd_re.pattern}")


def get_ccdnamenumber(item) -> tuple:
    """Extracts the name and number of a HiRISE CCD from the input (would get
    '('RED', '5')' from 'PSP_010502_2090_RED5_0.EDR_Stats.cub')."""
    try:
        nn = get_ccd(item)
        return get_ccdname(nn), get_ccdnumber(nn)
    except TypeError:
        return get_ccdname(item), str(item)
    except ValueError:
        return get_ccdname(int(item)), item


def get_ccdchannel(s: str) -> tuple:
    """Extracts a HiRISE CCD and channel number from the string (would get
    ('RED5', '1') from 'PSP_010502_2090_RED5_0.EDR_Stats.cub')."""
    match = ccdchan_re.search(s)
    if match:
        return match.groups()
    else:
        raise ValueError(f"{s} did not match regex: {ccdchan_re.pattern}")


def _reverse_recurse(path: os.PathLike, IDclass, name):
    p = Path(path)
    for part in reversed(p.parts):
        try:
            return IDclass(part)
        except ValueError:
            continue
    raise ValueError(f"Could not extract a {name} from {path}")


def _get_fromfile(path: os.PathLike, IDclass, name, archivekey):
    p = Path(path)

    try:
        return IDclass(isis.getkey_k(p, "Archive", archivekey))
    except (subprocess.CalledProcessError, ValueError):
        # The CalledProcessError is if there is some problem with running
        # getkey, the ValueError is if a CCDID can't be extracted from
        # the labels.  This allows this function to be called on any kind of
        # file, as it will just try and read the CCDID from the filename
        # or could also reverse recurse up to directory names.
        return _reverse_recurse(p, IDclass, name)


def get_ObsID_fromfile(path: os.PathLike) -> ObservationID:
    """Reads the file to get the ObservationID, if an ISIS cube,
       otherwise parses the filepath."""

    return _get_fromfile(
        path, ObservationID, "HiRISE Observation ID", "ObservationId"
    )


def get_CCDID_fromfile(path: os.PathLike) -> CCDID:
    """Reads the file to get the CCDID, if an ISIS cube,
       otherwise parses the filepath."""

    return _get_fromfile(path, CCDID, "HiRISE CCD ID", "ProductId")


def get_ChannelID_fromfile(path: os.PathLike) -> ChannelID:
    """Reads the file to get the ChannelID, if an ISIS cube,
       otherwise parses the filepath."""

    return _get_fromfile(path, ChannelID, "HiRISE Channel ID", "ProductId")
