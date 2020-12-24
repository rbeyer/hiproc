#!/usr/bin/env python
"""This module has tests for the img functions."""

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

import collections
import unittest

# from pathlib import Path
from unittest.mock import call, mock_open, patch

import numpy as np

# import kalasiris as isis

import hiproc.img as img

# from .utils import resource_check as rc

# # Hardcoding this, but I sure would like a better solution.
# HiRISE_img = Path('test-resources') / 'PSP_010502_2090_RED5_0.img'
# img = HiRISE_img
#
#
# class TestResources(unittest.TestCase):
#     '''Establishes that the test image exists.'''
#
#     def test_resources(self):
#         (truth, test) = rc(img)
#         self.assertEqual(truth, test)

example_lut = (
    (0, 1108),
    (1109, 1125),
    (1126, 1143),
    (1144, 1160),
    (1161, 1177),
    (1178, 1194),
    (1195, 1212),
    (1213, 1229),
    (1230, 1246),
    (1247, 1263),
    (1264, 1281),
    (1282, 1298),
    (1299, 1315),
    (1316, 1332),
    (1333, 1349),
    (1350, 1367),
    (1368, 1384),
    (1385, 1401),
    (1402, 1418),
    (1419, 1436),
    (1437, 1453),
    (1454, 1470),
    (1471, 1487),
    (1488, 1505),
    (1506, 1522),
    (1523, 1539),
    (1540, 1556),
    (1557, 1574),
    (1575, 1591),
    (1592, 1618),
    (1619, 1654),
    (1655, 1690),
    (1691, 1726),
    (1727, 1762),
    (1763, 1798),
    (1799, 1834),
    (1835, 1870),
    (1871, 1906),
    (1907, 1942),
    (1943, 1978),
    (1979, 2014),
    (2015, 2050),
    (2051, 2086),
    (2087, 2122),
    (2123, 2158),
    (2159, 2194),
    (2195, 2230),
    (2231, 2266),
    (2267, 2302),
    (2303, 2339),
    (2340, 2375),
    (2376, 2411),
    (2412, 2447),
    (2448, 2483),
    (2484, 2519),
    (2520, 2555),
    (2556, 2591),
    (2592, 2627),
    (2628, 2663),
    (2664, 2699),
    (2700, 2735),
    (2736, 2771),
    (2772, 2807),
    (2808, 2843),
    (2844, 2879),
    (2880, 2915),
    (2916, 2951),
    (2952, 2987),
    (2988, 3023),
    (3024, 3059),
    (3060, 3094),
    (3095, 3130),
    (3131, 3165),
    (3166, 3200),
    (3201, 3236),
    (3237, 3271),
    (3272, 3306),
    (3307, 3341),
    (3342, 3377),
    (3378, 3412),
    (3413, 3447),
    (3448, 3483),
    (3484, 3518),
    (3519, 3553),
    (3554, 3588),
    (3589, 3624),
    (3625, 3659),
    (3660, 3694),
    (3695, 3730),
    (3731, 3765),
    (3766, 3800),
    (3801, 3835),
    (3836, 3871),
    (3872, 3906),
    (3907, 3941),
    (3942, 3977),
    (3978, 4012),
    (4013, 4047),
    (4048, 4082),
    (4083, 4118),
    (4119, 4153),
    (4154, 4188),
    (4189, 4224),
    (4225, 4259),
    (4260, 4294),
    (4295, 4329),
    (4330, 4365),
    (4366, 4400),
    (4401, 4435),
    (4436, 4471),
    (4472, 4506),
    (4507, 4542),
    (4543, 4578),
    (4579, 4614),
    (4615, 4650),
    (4651, 4686),
    (4687, 4722),
    (4723, 4759),
    (4760, 4795),
    (4796, 4831),
    (4832, 4867),
    (4868, 4903),
    (4904, 4939),
    (4940, 4976),
    (4977, 5012),
    (5013, 5048),
    (5049, 5084),
    (5085, 5120),
    (5121, 5156),
    (5157, 5193),
    (5194, 5229),
    (5230, 5265),
    (5266, 5301),
    (5302, 5337),
    (5338, 5374),
    (5375, 5410),
    (5411, 5446),
    (5447, 5482),
    (5483, 5518),
    (5519, 5554),
    (5555, 5591),
    (5592, 5627),
    (5628, 5663),
    (5664, 5699),
    (5700, 5735),
    (5736, 5771),
    (5772, 5808),
    (5809, 5844),
    (5845, 5880),
    (5881, 5916),
    (5917, 5952),
    (5953, 5988),
    (5989, 6025),
    (6026, 6062),
    (6063, 6099),
    (6100, 6136),
    (6137, 6173),
    (6174, 6210),
    (6211, 6247),
    (6248, 6284),
    (6285, 6321),
    (6322, 6358),
    (6359, 6396),
    (6397, 6433),
    (6434, 6470),
    (6471, 6507),
    (6508, 6544),
    (6545, 6581),
    (6582, 6618),
    (6619, 6655),
    (6656, 6692),
    (6693, 6729),
    (6730, 6771),
    (6772, 6817),
    (6818, 6863),
    (6864, 6910),
    (6911, 6956),
    (6957, 7002),
    (7003, 7049),
    (7050, 7095),
    (7096, 7141),
    (7142, 7187),
    (7188, 7234),
    (7235, 7280),
    (7281, 7326),
    (7327, 7373),
    (7374, 7419),
    (7420, 7465),
    (7466, 7522),
    (7523, 7590),
    (7591, 7657),
    (7658, 7724),
    (7725, 7792),
    (7793, 7859),
    (7860, 7926),
    (7927, 7994),
    (7995, 8061),
    (8062, 8128),
    (8129, 8196),
    (8197, 8263),
    (8264, 8331),
    (8332, 8398),
    (8399, 8465),
    (8466, 8533),
    (8534, 8600),
    (8601, 8667),
    (8668, 8735),
    (8736, 8802),
    (8803, 8869),
    (8870, 8937),
    (8938, 9017),
    (9018, 9109),
    (9110, 9202),
    (9203, 9295),
    (9296, 9387),
    (9388, 9480),
    (9481, 9573),
    (9574, 9665),
    (9666, 9758),
    (9759, 9850),
    (9851, 9943),
    (9944, 10036),
    (10037, 10128),
    (10129, 10221),
    (10222, 10314),
    (10315, 10406),
    (10407, 10514),
    (10515, 10638),
    (10639, 10761),
    (10762, 10885),
    (10886, 11008),
    (11009, 11132),
    (11133, 11255),
    (11256, 11379),
    (11380, 11502),
    (11503, 11626),
    (11627, 11749),
    (11750, 11873),
    (11874, 12040),
    (12041, 12252),
    (12253, 12464),
    (12465, 12676),
    (12677, 12888),
    (12889, 13100),
    (13101, 13312),
    (13313, 13541),
    (13542, 13788),
    (13789, 14035),
    (14036, 14282),
    (14283, 14529),
    (14530, 14776),
    (14777, 15270),
    (15271, 16011),
    (16012, 16382),
    (16383, 16383),
    (-9998, -9998),
)


class TestLUT_Table(unittest.TestCase):
    def setUp(self):
        self.nolut = img.LUT_Table([[0, 0]])
        self.lut = img.LUT_Table(example_lut)

    def test_init(self):
        self.assertEqual(self.lut.table[0], 0)
        self.assertEqual(self.lut.table[1], 1117)
        self.assertEqual(self.lut.table[100], 4136)
        self.assertEqual(self.lut.table[253], 16197)
        self.assertEqual(self.lut.table[254], 16383)
        self.assertEqual(self.lut.table[255], -9998)

        self.assertEqual(self.nolut.table, False)
        self.assertEqual(str(self.nolut), "No applied LUT.")

        self.assertRaises(ValueError, img.LUT_Table, [(1, 10), (2, 20)])
        self.assertRaises(ValueError, img.LUT_Table, [None] * 257)

    def test_unlut(self):
        self.assertEqual(self.nolut.unlut(23), 23)
        self.assertRaises(IndexError, self.lut.unlut, 256)
        self.assertEqual(self.lut.unlut(0), 0)
        self.assertEqual(self.lut.unlut(1), 1117)

    def test_lookup(self):
        self.assertEqual(self.nolut.lookup(23), 23)
        self.assertEqual(self.lut.lookup(50), 0)
        self.assertEqual(self.lut.lookup(1145), 3)
        self.assertEqual(self.lut.lookup(1500), 23)
        self.assertEqual(self.lut.lookup(16000), 252)
        self.assertEqual(self.lut.lookup(-9998), 255)
        self.assertEqual(self.lut.lookup(20000), 255)


class TestFunctions(unittest.TestCase):
    def setUp(self):
        Units = collections.namedtuple("Units", ["value", "units"])
        self.label = {
            "^foo": Units(50013, "BYTES"),
            "foo": {
                "SAMPLE_BITS": 8,
                "SAMPLE_TYPE": "MSB_UNSIGNED_INTEGER",
                "LINES": 2,
                "LINE_PREFIX_BYTES": 1,
                "LINE_SAMPLES": 5,
                "LINE_SUFFIX_BYTES": 1,
            },
            "INSTRUMENT_SETTING_PARAMETERS": {
                "MRO:LOOKUP_CONVERSION_TABLE": example_lut
            },
        }

    def test_byte_info(self):
        self.assertRaises(KeyError, img.byte_info, "bar", self.label)
        self.assertEqual(img.byte_info("foo", self.label), (1, "big", False))

    def test_object_asarray(self):
        # Can't get this to work, possibly because there's a bug in 3.6
        # such that mock_open doesn't properly function when read(n) is
        # called on it, which is what object_asarray() does.
        # lutted = [1, 2, 0, 1, 1, 2, 1, 1, 1, 0]
        # b = bytes()
        # for i in lutted:
        #     b += i.to_bytes(1, byteorder='big', signed=False)
        b = (1).to_bytes(1, byteorder="big", signed=False)

        with patch("hiproc.img.open", mock_open(read_data=b)):
            with patch("hiproc.img.pvl.load", return_value=self.label):
                np.testing.assert_array_equal(
                    np.array([[1117, 0, 0, 0, 0], [0, 0, 0, 0, 0]]),
                    img.object_asarray("dummy.img", "foo"),
                )

    @patch("hiproc.img.open", new_callable=mock_open)
    def test_overwrite_object(self, m_open):
        with patch("hiproc.img.pvl.load", return_value=self.label):
            arr = [[1111, 0, 1120, 0, 1125], [1126, 1140, -9998, 0, 0]]
            img.overwrite_object("dummy.img", "foo", np.array(arr))

        handle = m_open()
        self.assertEqual(
            handle.seek.call_args_list,
            [call(50012), call(1, 1), call(1, 1), call(1, 1), call(1, 1)],
        )

        lutted = [1, 0, 1, 0, 1, 2, 2, 255, 0, 0]
        call_list = list()
        for i in lutted:
            call_list.append(
                call(i.to_bytes(1, byteorder="big", signed=False))
            )

        self.assertEqual(handle.write.call_args_list, call_list)
