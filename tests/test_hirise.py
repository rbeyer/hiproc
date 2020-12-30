#!/usr/bin/env python
"""This module has tests for the HiRISE utility functions."""

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

import contextlib
import itertools
import shutil
import unittest
from pathlib import Path

import kalasiris as isis

from hiproc import hirise
from .utils import resource_check as rc

# Hardcoding this, but I sure would like a better solution.
HiRISE_img = Path("test-resources") / "PSP_010502_2090_RED5_0.img"
img = HiRISE_img


class TestResources(unittest.TestCase):
    """Establishes that the test image exists."""

    def test_resources(self):
        (truth, test) = rc(img)
        self.assertEqual(truth, test)


class TestInternal(unittest.TestCase):

    # def test_match_chan_yes(self):
    #     self.assertEqual('1', hirise._match_chan('1'))

    # def test_match_chan_no(self):
    #     self.assertIsNone(hirise._match_chan('5'))

    def test_match_num_yes(self):
        strings = ("1", "9", "12")
        for s in strings:
            with self.subTest(s=s):
                self.assertEqual(s, hirise._match_num(s))

    def test_match_num_no(self):
        self.assertIsNone(hirise._match_num("22"))


class TestParse(unittest.TestCase):
    def test_parseObsID_good_embed(self):
        truth = "ESP", "123456", "1235"
        s = "This is an ObservationID: " + "_".join(truth)
        tokens = hirise.parseObsID(s)
        self.assertTupleEqual(truth, tokens)

    def test_parseObsID_no_ObsID(self):
        s = "There is no Observation ID here."
        self.assertRaises(ValueError, hirise.parseObsID, s)

    def test_ObsIDstr(self):
        truth = "ESP_123456_1235"
        s = "This is an ObservationID: " + truth
        self.assertEqual(truth, hirise.ObsIDstr(s))


class TestPhase(unittest.TestCase):
    def setUp(self):
        self.tuples = (
            ("0", "AEB", True),
            ("0", "TRA", False),
            ("500", "TRA", True),
            ("1100", "TRA", False),
            ("1100", "PSP", True),
            ("12000", "PSP", False),
            ("12000", "ESP", True),
            ("999999999999", "ESP", True),
        )

    def test_prev_phase_first(self):
        self.assertRaises(IndexError, hirise.prev_phase, hirise.phase_names[0])

    def test_prev_phase(self):
        these_phases = hirise.phase_names[1:]
        for prev, this in zip(hirise.phase_names, these_phases):
            with self.subTest():
                self.assertEqual(prev, hirise.prev_phase(this))

    def test_orbit_in_phase(self):
        for t in self.tuples:
            with self.subTest():
                self.assertEqual(hirise.is_orbit_in_phase(t[0], t[1]), t[2])

    def test_get_phase(self):
        for t in filter(lambda x: int(x[0]) and x[-1], self.tuples):
            with self.subTest():
                self.assertEqual(t[1], hirise.get_phase(t[0]))

    def test_get_phase_bad(self):
        for t in filter(
            lambda y: int(y[0]),
            itertools.filterfalse(lambda x: x[-1], self.tuples),
        ):
            with self.subTest():
                self.assertNotEqual(t[1], hirise.get_phase(t[0]))

    def test_get_phase_Error(self):
        for t in filter(lambda x: int(x[0]) == 0, self.tuples):
            with self.subTest():
                self.assertRaises(IndexError, hirise.get_phase, t[0])


class TestObsID(unittest.TestCase):
    def test_init_tuple(self):
        tuples = (
            (("PSP", "005632", "1225"), ("PSP", "005632", "1225")),
            (("PSP", "005632", "1225"), ("5632", "1225")),
            (("TRA", "000001", "0005"), ("1", "5")),
            (("ESP", "057866", "1670"), ("ESP", "57866", "1670")),
            (("ESP", "034783", "1850"), ("ESP", "034783", "1850")),
            (("AEB", "000000", "0000"), ("AEB", "000000", "0000")),
        )
        for truth, t in tuples:
            with self.subTest():
                oid = hirise.ObservationID(*t)
                self.assertTupleEqual(
                    truth, (oid.phase, oid.orbit_number, oid.target)
                )

    def test_init_string(self):
        string_tuples = (
            ("PSP_005632_1225", "PSP_005632_1225"),
            ("PSP_005632_1225", "5632_1225"),
            ("PSP_005632_0005", "5632_5"),
            ("TRA_000001_0005", "1_5"),
            ("ESP_057866_1670", "ESP_057866_1670"),
            ("ESP_057866_1670", "ESP_57866_1670"),
            ("ESP_034783_1855", "ESP_034783_1855"),
            ("ESP_034783_1850", "ESP_034783_1850"),
        )
        for s in string_tuples:
            with self.subTest():
                oid = hirise.ObservationID(s[1])
                self.assertEqual(s[0], str(oid))

    def test_init_bad_tuples(self):
        tuples = (
            ("ABCD", "123456", "1235"),
            ("AB", "123456", "1230"),
            ("ESP", "1234", "1235"),
            ("ESP", "1234567", "1230"),
            ("ESP", "12345", "12345"),
            ("ESP", "12345", "123"),
        )
        for t in tuples:
            with self.subTest(t):
                self.assertRaises(ValueError, hirise.ObservationID, *t)

    def test_init_bad_strings(self):
        strings = (
            "ABC-123456-1235",
            "AB_123456_1235",
            "ESPA_123456_1235",
            "PSP_1234_1234",
            "ESP_1234567_1235",
            "ESP_123456_12345",
            "ESP_123456",
            "foobar",
        )
        for s in strings:
            with self.subTest(s):
                self.assertRaises(ValueError, hirise.ObservationID, s)

    def test_init_too_many_args(self):
        self.assertRaises(
            IndexError, hirise.ObservationID, "ABC", "123456", "1234", "extra"
        )

    def test_init_zero_args(self):
        self.assertRaises(IndexError, hirise.ObservationID)

    def test_repr(self):
        s = "This is an Observation ID: ESP_057866_1670"
        oid = hirise.ObservationID(s)
        self.assertEqual("ObservationID('ESP_057866_1670')", repr(oid))

    def test_lt(self):
        oid1 = hirise.ObservationID("PSP_005632_1225")
        oid2 = hirise.ObservationID("ESP_057866_1670")
        oid3 = hirise.ObservationID("ESP_057867_1670")
        oid4 = hirise.ObservationID("ESP_057867_1675")
        self.assertTrue(oid1 < oid2)
        self.assertTrue(oid2 < oid3)
        self.assertTrue(oid3 < oid4)


class TestCCDID(unittest.TestCase):
    def test_init_tuple(self):
        tuples = (
            (
                ("ESP", "034783", "1850", "RED", "5"),
                ("ESP", "034783", "1850", "RED", "5"),
            ),
            (
                ("ESP", "034783", "1850", "RED", "5"),
                ("ESP", "034783", "1850", "RED5"),
            ),
            (
                ("ESP", "034783", "1850", "RED", "5"),
                ("034783", "1850", "RED", "5"),
            ),
            (
                ("ESP", "034783", "1850", "RED", "5"),
                ("034783", "1850", "RED5"),
            ),
            (("ESP", "034783", "1850", "RED", "5"), ("34783", "1850", "RED5")),
            (
                ("ESP", "034783", "1850", "RED", "5"),
                ("ESP", "034783", "1850", "RED", "5"),
            ),
            (
                ("ESP", "034783", "1850", "RED", "5"),
                ("034783", "1850", "RED", "5"),
            ),
            (
                ("ESP", "034783", "1850", "RED", "5"),
                ("ESP", "034783", "1850", "RED5"),
            ),
            (
                ("ESP", "034783", "1850", "RED", "5"),
                ("ESP", "034783", "1850", "RED", "5"),
            ),
        )
        for truth, t in tuples:
            with self.subTest(truth=truth, test=t):
                cid = hirise.CCDID(*t)
                self.assertTupleEqual(
                    truth,
                    (
                        cid.phase,
                        cid.orbit_number,
                        cid.target,
                        cid.ccdname,
                        cid.ccdnumber,
                    ),
                )

    def test_init_string(self):
        string_tuples = (
            ("ESP_034783_1850_RED5", "ESP_034783_1850_RED5"),
            ("ESP_034783_1850_RED5", "34783_1850_RED5"),
            ("ESP_034783_1850_RED5", "ESP_034783_1850_RED5_0"),
            ("ESP_034783_1850_RED5", "034783_1850_RED5_0"),
        )
        for s in string_tuples:
            cid = hirise.CCDID(s[1])
            with self.subTest(truth=s[0], ccid=cid):
                self.assertEqual(s[0], str(cid))

    def test_init_bad_tuples(self):
        tuples = (
            ("ESP", "034783", "1850"),
            ("34783", "1850", "BLU1"),
            ("34783", "1850", "RED"),
        )
        for t in tuples:
            with self.subTest(test=t):
                self.assertRaises(ValueError, hirise.CCDID, *t)

    def test_init_bad_strings(self):
        strings = (
            "ABC-123456-1235",
            "AB_123456_1235",
            "ESPA_123456_1235",
            "PSP_1234_1234",
            "ESP_1234567_1235",
            "ESP_123456_12345",
            "ESP_123456",
            "ESP_034783_1850",
            "ESP_034783_1850_IR13",
            "foobar",
        )
        for s in strings:
            with self.subTest(test=s):
                self.assertRaises(ValueError, hirise.CCDID, s)

    def test_init_wrong_arg_count(self):
        self.assertRaises(
            IndexError, hirise.CCDID, "ESP", "123456", "1235", "RED", "5", "0"
        )
        self.assertRaises(IndexError, hirise.CCDID, "123456", "1235")
        self.assertRaises(IndexError, hirise.CCDID)

    def test_repr(self):
        s = "This is a CCDID: ESP_034783_1850_RED5"
        cid = hirise.CCDID(s)
        self.assertEqual("CCDID('ESP_034783_1850_RED5')", repr(cid))

    def test_lt(self):
        cid1 = hirise.CCDID("ESP_034783_1850_RED4")
        cid2 = hirise.CCDID("ESP_034783_1850_RED5")
        cid3 = hirise.CCDID("ESP_034783_1850_IR10")
        cid4 = hirise.CCDID("ESP_034783_1850_BG13")
        self.assertTrue(cid1 < cid2)
        self.assertTrue(cid2 < cid3)
        self.assertTrue(cid3 < cid4)

        ccds = [cid3, cid4, cid1, cid2]
        self.assertEqual(sorted(ccds), [cid1, cid2, cid3, cid4])

    def test_get_ccd(self):
        s = "This is a CCDID: ESP_034783_1850_RED5"
        cid = hirise.CCDID(s)
        self.assertEqual(cid.get_ccd(), "RED5")

    def test_get_obsid(self):
        s = "This is a CCDID: ESP_034783_1850_RED5"
        oid = hirise.ObservationID(s)
        cid = hirise.CCDID(s)
        self.assertEqual(oid, cid.get_obsid())


class TestChannelID(unittest.TestCase):
    def test_init_tuple(self):
        tuples = (
            (
                ("ESP", "034783", "1850", "RED", "5", "0"),
                ("034783", "1850", "RED5", "0"),
            ),
            (
                ("ESP", "034783", "1850", "RED", "5", "0"),
                ("034783", "1850", "RED", "5", "0"),
            ),
            (
                ("ESP", "034783", "1850", "RED", "5", "0"),
                ("ESP", "034783", "1850", "RED5", "0"),
            ),
            (
                ("ESP", "034783", "1850", "RED", "5", "0"),
                ("ESP", "034783", "1850", "RED", "5", "0"),
            ),
        )

        for truth, t in tuples:
            with self.subTest(truth=truth, test=t):
                cid = hirise.ChannelID(*t)
                self.assertTupleEqual(
                    truth,
                    (
                        cid.phase,
                        cid.orbit_number,
                        cid.target,
                        cid.ccdname,
                        cid.ccdnumber,
                        cid.channel,
                    ),
                )

    def test_init_string(self):
        string_tuples = (
            ("ESP_034783_1850_RED5_0", "ESP_034783_1850_RED5_0"),
            ("ESP_034783_1850_RED5_1", "34783_1850_RED5_1"),
        )
        for s in string_tuples:
            cid = hirise.ChannelID(s[1])
            with self.subTest(truth=s[0], channelid=cid):
                self.assertEqual(s[0], str(cid))

    def test_init_CCDID(self):
        obsid = hirise.ObservationID("ESP_034783_1850")
        ccdid = hirise.CCDID("ESP_034783_1850_RED5")
        chanid = hirise.ChannelID("ESP_034783_1850_RED5_0")
        cid1 = hirise.ChannelID(ccdid, "0")
        cid2 = hirise.ChannelID(obsid, 5, 0)
        self.assertEqual(cid1, chanid)
        self.assertEqual(cid2, chanid)

    def test_init_bad_tuples(self):
        tuples = (
            ("ESP", "034783", "1850", "RED5"),
            ("ESP", "034783", "1850", "RED5", "A"),
            ("ESP", "034783", "1850", "RED5", "2"),
            ("123456", "1235", "123"),
            ("123456", "1235"),
            ("ESP", "034783", "1850", "RED5", "not a channel"),
        )
        for t in tuples:
            with self.subTest(test=t):
                self.assertRaises(ValueError, hirise.ChannelID, *t)

    def test_init_bad_strings(self):
        strings = (
            "ABC-123456-1235",
            "AB_123456_1235",
            "ESPA_123456_1235",
            "PSP_1234_1234",
            "ESP_1234567_1235",
            "ESP_123456_12345",
            "ESP_123456",
            "ESP_034783_1850",
            "ESP_034783_1850_IR13",
            "ESP_034783_1850_RED5",
            "ESP_034783_1850_RED5_A",
            "ESP_034783_1850_RED5_2",
            "ESP_034783_1850_RED5_10" "foobar",
        )
        for s in strings:
            with self.subTest(test=s):
                self.assertRaises(ValueError, hirise.ChannelID, s)

    def test_init_wrong_arg_count(self):
        self.assertRaises(
            IndexError,
            hirise.ChannelID,
            "ESP",
            "123456",
            "1235",
            "RED",
            "5",
            "0",
            "extra",
        )
        self.assertRaises(IndexError, hirise.ChannelID)

    def test_repr(self):
        s = "This is a CCDID: ESP_034783_1850_RED5_0"
        cid = hirise.ChannelID(s)
        self.assertEqual("ChannelID('ESP_034783_1850_RED5_0')", repr(cid))

    def test_lt(self):
        cid1 = hirise.ChannelID("ESP_034783_1850_RED5_0")
        cid2 = hirise.ChannelID("ESP_034783_1850_RED5_1")
        cid3 = hirise.ChannelID("ESP_034783_1850_IR10_1")
        self.assertTrue(cid1 < cid2)
        self.assertTrue(cid2 < cid3)

    def test_get_ccd(self):
        s = "This is a ChannelID: ESP_034783_1850_RED5_0"
        cid = hirise.ChannelID(s)
        self.assertEqual(cid.get_ccd(), "RED5")

    def test_get_obsid(self):
        s = "This is a ChannelID: ESP_034783_1850_RED5_0"
        oid = hirise.ObservationID(s)
        cid = hirise.ChannelID(s)
        self.assertEqual(oid, cid.get_obsid())


class TestStrings(unittest.TestCase):
    def test_str_(self):
        test = "This is an Obs Id: ESP_057866_1670"
        truth = "ESP_057866_1670"
        self.assertEqual(truth, hirise.ObservationID(test).__str__())


class TestGetters(unittest.TestCase):
    def setUp(self):
        self.ccds = (
            ("RED", "0"),
            ("RED", "5"),
            ("IR", "10"),
            ("IR", "11"),
            ("BG", "12"),
            ("BG", "13"),
        )
        self.channels = ("0", "1")

    def test_get_ccd_good(self):
        s = "This is a good CCD name: RED4"
        truth = "RED4"
        self.assertIn(hirise.get_ccd(s), truth)

    def test_get_ccd_bulk(self):
        for c in self.ccds:
            with self.subTest(c):
                s = "".join(c)
                self.assertEqual(hirise.get_ccd(s), s)

    def test_get_ccd_bad(self):
        s = "There is no CCD name in here."
        self.assertRaises(ValueError, hirise.get_ccd, s)

    def test_get_ccdname(self):
        for c in self.ccds:
            with self.subTest(c):
                s = "".join(c)
                self.assertEqual(hirise.get_ccdname(s), c[0])

    def test_get_ccdname_int(self):
        self.assertEqual(hirise.get_ccdname(5), "RED")
        self.assertEqual(hirise.get_ccdname(str(5)), "RED")

    def test_get_ccdname_bad(self):
        s = "There is no CCD name in here: ESP_057866_1670_YEL5_0"
        self.assertRaises(ValueError, hirise.get_ccd, s)
        self.assertRaises(TypeError, hirise.get_ccdname, 5.44)

    def test_get_ccdnumber(self):
        for c in self.ccds:
            with self.subTest(c=c):
                s = "".join(c)
                self.assertEqual(hirise.get_ccdnumber(s), c[1])

    def test_get_ccdnumber_bad(self):
        s = "There is no CCD number in here: ESP_057866_1670_RED_0"
        self.assertRaises(ValueError, hirise.get_ccdnumber, s)

    def test_get_ccdnamenumber(self):
        for c in self.ccds:
            with self.subTest(c):
                s = "".join(c)
                self.assertEqual(hirise.get_ccdnamenumber(s), c)

    def test_get_ccdnamenumber_just_number(self):
        for c in self.ccds:
            with self.subTest(c):
                self.assertEqual(hirise.get_ccdnamenumber(c[1]), c)

    def test_get_ccdnamenumber_int(self):
        for c in self.ccds:
            with self.subTest(c):
                self.assertEqual(hirise.get_ccdnamenumber(int(c[1])), c)

    def test_get_ccdnamenumber_bad(self):
        self.assertRaises(ValueError, hirise.get_ccdnamenumber, "RED")

    def test_get_ccdchannel_good(self):
        s = "IR10_0 is a good CCD-channel combination."
        truth = ("IR10", "0")
        self.assertEqual(hirise.get_ccdchannel(s), truth)

    def test_get_ccdchannel_bad(self):
        s = "There is a CCD, BG12, but no channel here."
        self.assertRaises(ValueError, hirise.get_ccdchannel, s)


class TestFromFile(unittest.TestCase):

    # HiRISE_img = Path('test-resources') / 'PSP_010502_2090_RED5_0.img'

    def setUp(self):
        named_cube = img.with_suffix(".cub")
        isis.hi2isis(img, to=named_cube)
        noname_cube = named_cube.with_name("noname.cub")
        shutil.copyfile(named_cube, noname_cube)

        self.paths = (
            (named_cube, "PSP_010502_2090", "PSP_010502_2090_RED5_0"),
            (noname_cube, "PSP_010502_2090", "PSP_010502_2090_RED5_0"),
            (
                Path("PSP_010502_2090_RED5.fake"),
                "PSP_010502_2090",
                "PSP_010502_2090_RED5",
            ),
        )

    def tearDown(self):
        with contextlib.suppress(FileNotFoundError):
            for p in self.paths:
                p[0].unlink()
            Path("print.prt").unlink()

    def test_get_ObsID(self):
        for t in self.paths:
            with self.subTest(path=t[0], ObsID=t[1]):
                oid = hirise.get_ObsID_fromfile(t[0])
                self.assertEqual(t[1], str(oid))

    def test_get_ObsID_fail(self):
        p = Path("Not-a-file")
        self.assertRaises(ValueError, hirise.get_ObsID_fromfile, p)
