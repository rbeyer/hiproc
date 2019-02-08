#!/usr/bin/env python
"""This module has tests for the HiRISE utility functions."""

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

import unittest, sys
sys.path.append('../')
from hirise import *

class TestParse(unittest.TestCase):
    def test_parseObsID_good_embed(self):
        truth = 'ABC','123456','1234'
        s = 'This is an ObservationID: '+'_'.join(truth)
        tokens = parseObsID( s )
        self.assertTupleEqual( truth, tokens )

    def test_parseObsID_no_ObsID(self):
        s = 'There is no Observation ID here.'
        self.assertRaises( ValueError, parseObsID, s )

    def test_getObsID(self):
        truth = 'ABC_123456_1234'
        s = 'This is an ObservationID: '+truth
        self.assertEqual( truth, getObsID(s) )


class TestInit(unittest.TestCase):
    def setUp(self):
        self.tuples = ( ('PSP','005632','1225'),
                        ('ESP', '57866','1670'), 
                        ('ESP','034783','1850'), 
                        ('R01','000000','0000') )
        self.strings = ('PSP_005632_1225', 'ESP_57866_1670', 'ESP_034783_1850')

    def test_init_tuple(self):
        for t in self.tuples:
            with self.subTest():
                oid = ObservationID( *t )
                self.assertTupleEqual(t,(oid.phase,oid.orbit_number,oid.latesque))

    def test_init_string(self):
        for s in self.strings:
            with self.subTest():
                oid = ObservationID( s )
                self.assertEqual( s, str(oid) )

    def test_init_bad_tuples(self):
        tuples = ( ('ABCD','123456', '1234'),
                   ('AB',  '123456', '1234'),
                   ('ABC', '1234',   '1234'), 
                   ('ABC', '1234567','1234'), 
                   ('ABC', '12345',  '12345'), 
                   ('ABC', '12345',  '123') )
        for t in tuples:
            with self.subTest(t):
                self.assertRaises( ValueError, ObservationID, *t )

    def test_init_bad_strings(self):
        strings = ('ABC-123456-1234', 
                   'AB_123456_1234',  'ABCD_123456_1234',
                   'ABC_1234_1234',   'ABC_1234567_1234',
                   'ABC_123456_12345', 'ABC_123456_123',
                   'foobar')
        for s in strings:
            with self.subTest(s):
                self.assertRaises( ValueError, ObservationID, s )

    def test_init_too_many_args(self):
        self.assertRaises( IndexError, ObservationID,'ABC','123456','1234','extra' )

    def test_init_zero_args(self):
        self.assertRaises( IndexError, ObservationID )


class TestGetCCDchannel(unittest.TestCase):
    def setUp(self):
        self.ccds = ('RED0','RED1','RED2',
                     'RED3','RED4','RED5',
                     'RED6','RED7','RED8','RED9',
                     'IR10','IR11','BG12','BG13')
        self.channels = ('0','1')

    def test_getccd_good(self):
        s = 'This is a good CCD name: RED4'
        self.assertIn( getccd(s), self.ccds )

    def test_getccd_bad(self):
        s = 'There is no CCD name in here.'
        self.assertRaises( ValueError, getccd, s )

    def test_getccdchannel_good(self):
        s = 'IR10_0 is a good CCD-channel combination.'
        tokens = getccdchannel( s )
        self.assertIn( tokens[0], self.ccds )
        self.assertIn( tokens[1], self.channels )

    def test_getccdchannel_bad(self):
        s = 'There is a CCD, BG12, but no channel here.'
        self.assertRaises( ValueError, getccdchannel, s )
