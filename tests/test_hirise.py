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

import itertools
import unittest
# from pathlib import Path

from PyRISE import hirise


class TestParse(unittest.TestCase):

    def test_parseObsID_good_embed(self):
        truth = 'ESP', '123456', '1235'
        s = 'This is an ObservationID: ' + '_'.join(truth)
        tokens = hirise.parseObsID(s)
        self.assertTupleEqual(truth, tokens)

    def test_parseObsID_no_ObsID(self):
        s = 'There is no Observation ID here.'
        self.assertRaises(ValueError, hirise.parseObsID, s)

    def test_getObsID(self):
        truth = 'ESP_123456_1235'
        s = 'This is an ObservationID: ' + truth
        self.assertEqual(truth, hirise.getObsID(s))

    @unittest.skip('Thinking about just abolishing these functions.')
    def test_parseProdID(self):
        truth = (('ESP', '123456', '1235'),
                 ('PSP', '1100', '0000', 'IR', '10'),
                 ('ESP', '123456', '1235', 'RED', '5', '0'))
        strings = (('ESP_123456_1235'),
                   ('PSP_1100_0000_IR10'),
                   ('ESP_123456_1235_RED5_0'))
        for t, s in zip(truth, strings):
            tokens = hirise.parseProdID(s)
            self.assertTupleEqual(tokens, t)


class TestPhase(unittest.TestCase):

    def setUp(self):
        self.tuples = (('0', 'AEB', True),
                       ('0', 'TRA', False),
                       ('500', 'TRA', True),
                       ('1100', 'TRA', False),
                       ('1100', 'PSP', True),
                       ('12000', 'PSP', False),
                       ('12000', 'ESP', True),
                       ('999999999999', 'ESP', True))

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
                self.assertEqual(hirise.orbit_in_phase(t[0], t[1]), t[2])

    def test_getphase(self):
        for t in filter(lambda x: int(x[0]) and x[-1], self.tuples):
            with self.subTest():
                self.assertEqual(t[1], hirise.getphase(t[0]))

    def test_getphase_bad(self):
        for t in filter(lambda y: int(y[0]),
                        itertools.filterfalse(lambda x: x[-1],
                                              self.tuples)):
            with self.subTest():
                self.assertNotEqual(t[1], hirise.getphase(t[0]))

    def test_getphase_Error(self):
        for t in filter(lambda x: int(x[0]) == 0, self.tuples):
            with self.subTest():
                self.assertRaises(IndexError, hirise.getphase, t[0])


class TestObsID(unittest.TestCase):

    def setUp(self):
        self.tuples = (('PSP', '005632', '1225'),
                       ('ESP', '57866', '1670'),
                       ('ESP', '034783', '1850'),
                       ('AEB', '000000', '0000'))
        self.strings = ('PSP_005632_1225', 'ESP_057866_1670', 'ESP_034783_1850')

    def test_init_tuple(self):
        for t in self.tuples:
            with self.subTest():
                oid = hirise.ObservationID(*t)
                orbit = '{:0>6}'.format(t[1])
                truth = (t[0], orbit, t[2])
                self.assertTupleEqual(truth, (oid.phase, oid.orbit_number, oid.latesque))

    def test_init_string(self):
        for s in self.strings:
            with self.subTest():
                oid = hirise.ObservationID(s)
                self.assertEqual(s, str(oid))

    def test_init_bad_tuples(self):
        tuples = (('ABCD', '123456', '1235'),
                  ('AB',  '123456',  '1230'),
                  ('ESP', '1234',    '1235'),
                  ('ESP', '1234567', '1230'),
                  ('ESP', '12345',   '12345'),
                  ('ESP', '12345',   '123'))
        for t in tuples:
            with self.subTest(t):
                self.assertRaises(ValueError, hirise.ObservationID, *t)

    def test_init_bad_strings(self):
        strings = ('ABC-123456-1234',
                   'AB_123456_1234',  'ABCD_123456_1234',
                   'ABC_1234_1234',   'ABC_1234567_1234',
                   'ABC_123456_12345', 'ABC_123456_123',
                   'foobar')
        for s in strings:
            with self.subTest(s):
                self.assertRaises(ValueError, hirise.ObservationID, s)

    def test_init_too_many_args(self):
        self.assertRaises(IndexError, hirise.ObservationID,
                          'ABC', '123456', '1234', 'extra')

    def test_init_zero_args(self):
        self.assertRaises(IndexError, hirise.ObservationID)


class TestProdID(unittest.TestCase):

    def setUp(self):
        self.tuples = (('PSP', '005632', '1225'),
                       ('ESP', '57866', '1670'),
                       ('ESP', '034783', '1850'),
                       ('AEB', '000000', '0000'))
        self.strings = ('PSP_005632_1225', 'ESP_057866_1670', 'ESP_034783_1850')

    def test_init_tuple(self):
        for t in self.tuples:
            with self.subTest():
                oid = hirise.ProductID(*t)
                orbit = '{:0>6}'.format(t[1])
                truth = (t[0], orbit, t[2])
                self.assertTupleEqual(truth, (oid.phase, oid.orbit_number, oid.latesque))

    def test_init_string(self):
        for s in self.strings:
            with self.subTest():
                oid = hirise.ProductID(s)
                self.assertEqual(s, str(oid))

    def test_init_bad_tuples(self):
        tuples = (('ABCD', '123456', '1235'),
                  ('AB',  '123456',  '1230'),
                  ('ESP', '1234',    '1235'),
                  ('ESP', '1234567', '1230'),
                  ('ESP', '12345',   '12345'),
                  ('ESP', '12345',   '123'))
        for t in tuples:
            with self.subTest(t):
                self.assertRaises(ValueError, hirise.ProductID, *t)

    def test_init_bad_strings(self):
        strings = ('ABC-123456-1234',
                   'AB_123456_1234',  'ABCD_123456_1234',
                   'ABC_1234_1234',   'ABC_1234567_1234',
                   'ABC_123456_12345', 'ABC_123456_123',
                   'foobar')
        for s in strings:
            with self.subTest(s):
                self.assertRaises(ValueError, hirise.ProductID, s)

    def test_init_too_many_args(self):
        self.assertRaises(ValueError, hirise.ProductID,
                          'ABC', '123456', '1234', 'extra')

    def test_init_zero_args(self):
        self.assertRaises(IndexError, hirise.ProductID)


class TestStrings(unittest.TestCase):

    def test_str_(self):
        test = 'This is an Obs Id: ESP_057866_1670'
        truth = 'ESP_057866_1670'
        self.assertEqual(truth, hirise.ObservationID(test).__str__())


class TestGetters(unittest.TestCase):

    def setUp(self):
        self.ccds = (('RED', '0'),
                     ('RED', '5'),
                     ('IR', '10'),
                     ('IR', '11'),
                     ('BG', '12'),
                     ('BG', '13'))
        self.channels = ('0', '1')

    def test_getccd_good(self):
        s = 'This is a good CCD name: RED4'
        truth = 'RED4'
        self.assertIn(hirise.getccd(s), truth)

    def test_getccd_bulk(self):
        for c in self.ccds:
            with self.subTest(c):
                s = ''.join(c)
                self.assertEquals(hirise.getccd(s), s)

    def test_getccd_bad(self):
        s = 'There is no CCD name in here.'
        self.assertRaises(ValueError, hirise.getccd, s)

    def test_getccdname(self):
        for c in self.ccds:
            with self.subTest(c):
                s = ''.join(c)
                self.assertEquals(hirise.getccdname(s), c[0])

    def test_getccdname_bad(self):
        s = 'There is no CCD name in here: ESP_057866_1670_YEL5_0'
        self.assertRaises(ValueError, hirise.getccd, s)

    def test_getccdnumber(self):
        for c in self.ccds:
            with self.subTest(c):
                s = ''.join(c)
                self.assertEquals(hirise.getccdnumber(s), c[1])

    def test_getccdnamenumber(self):
        for c in self.ccds:
            with self.subTest(c):
                s = ''.join(c)
                self.assertEquals(hirise.getccdnamenumber(s), c)

    def test_getccdchannel_good(self):
        s = 'IR10_0 is a good CCD-channel combination.'
        truth = ('IR10', '0')
        self.assertEquals(hirise.getccdchannel(s), truth)

    def test_getccdchannel_bad(self):
        s = 'There is a CCD, BG12, but no channel here.'
        self.assertRaises(ValueError, hirise.getccdchannel, s)
