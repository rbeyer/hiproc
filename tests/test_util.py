#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for the `util` module."""

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

import argparse
import logging
import unittest
from pathlib import Path

import PyRISE.util as util

from .utils import resource_check as rc

# Hardcoding these, but I sure would like a better solution.
HiRISE_img = Path('test-resources') / 'PSP_010502_2090_RED5_0.img'
img = HiRISE_img


class TestResources(unittest.TestCase):
    '''Establishes that the test image exists.'''

    def test_resources(self):
        (truth, test) = rc(img)
        self.assertEqual(truth, test)


class TestUtil(unittest.TestCase):

    def test_parent_parser(self):
        self.assertIsInstance(util.parent_parser(), argparse.ArgumentParser)

    def test_logging(self):
        util.set_logging('WARNING')
        logger = logging.getLogger()
        self.assertEqual(30, logger.getEffectiveLevel())

    def test_path_w_suffix(self):
        self.assertEqual('bar.foo',
                         str(util.path_w_suffix('.foo', 'bar.cub')))
        self.assertEqual('foo.foo',
                         str(util.path_w_suffix('foo.foo', 'bar.cub')))

    def test_pid_path_w_suffix(self):
        self.assertRaises(FileNotFoundError, util.pid_path_w_suffix,
                          'foo.foo', 'dummy')
        self.assertRaises(FileNotFoundError, util.pid_path_w_suffix, '.foo',
                          'PSP_010502_2090_RED5_0.dummy')
        self.assertEqual(img,
                         util.pid_path_w_suffix('.img',
                                                img.with_suffix('.foo.foo')))

    def test_get_path(self):
        self.assertRaises(ValueError, util.get_path, 'dummy')
        self.assertRaises(TypeError, util.get_path, 'dummy', 42)
        self.assertRaises(NotADirectoryError, util.get_path, 'dummy',
                          'not_a_dir')
        self.assertRaises(FileNotFoundError, util.get_path, 'dummy', 'PyRISE')
        self.assertRaises(FileNotFoundError, util.get_path, 'dummy',
                          Path('PyRISE'))
        self.assertRaises(FileNotFoundError, util.get_path, 'dummy',
                          ['PyRISE', 'tests'])
        self.assertEqual(Path('PyRISE/util.py'), util.get_path('util.py',
                                                               'PyRISE'))

    def test_conf_check_strings(self):
        self.assertIsNone(util.conf_check_strings('foo', ('YES', 'NO'), 'YES'))
        self.assertRaises(AssertionError,
                          util.conf_check_strings, 'foo', ('YES', 'NO'),
                          'MAYBE')

    def test_conf_check_count(self):
        self.assertIsNone(util.conf_check_count('foo', 3, 'what', ['one',
                                                                   'two',
                                                                   'three']))
        self.assertRaises(AssertionError,
                          util.conf_check_count, 'foo', 3, 'what', ['one, two'])

    def test_conf_check_bounds(self):
        self.assertIsNone(util.conf_check_bounds('foo', (0.1, 1), '0.5'))
        self.assertRaises(AssertionError,
                          util.conf_check_bounds, 'foo', (0.1, 1), '1.5')
