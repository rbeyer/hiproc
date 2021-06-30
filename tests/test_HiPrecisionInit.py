#!/usr/bin/env python
"""This module has tests for the HiRISE HiPrecisionInit functions."""

# Copyright 2020-2021, Ross A. Beyer (rbeyer@seti.org)
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
import unittest
from unittest.mock import patch

import hiproc.HiPrecisionInit as hpi


class TestMock(unittest.TestCase):
    def test_arg_parser(self):
        p = hpi.arg_parser()
        self.assertIsInstance(p, argparse.ArgumentParser)

    @patch("hiproc.HiPrecisionInit.pvl.load")
    @patch("hiproc.HiPrecisionInit.arg_parser")
    @patch("hiproc.HiPrecisionInit.util.set_logger")
    @patch("hiproc.HiPrecisionInit.util.main_exceptions")
    @patch("hiproc.HiPrecisionInit.print")
    @patch(
        "hiproc.HiPrecisionInit.check",
        return_value=(
            [True, True],
            [5, 5],
            1
        )
    )
    def test_main(self, m_check, m_print, m_except, m_logger, m_parse, m_load):
        with patch("hiproc.HiPrecisionInit.check", return_value=(
            [True, True],
            [5, 5],
            1
        )):
            self.assertIsNone(hpi.main())

        with patch("hiproc.HiPrecisionInit.check", return_value=(
            [False, False],
            [5, 5],
            1
        )):
            self.assertIsNone(hpi.main())

    @patch("hiproc.HiPrecisionInit.sstats.Polyfit", return_value=("_", 5, "_"))
    def test_needs_HiJACK(self, sstats):
        need, avediff = hpi.needs_HiJACK("dummy_path", 1)
        self.assertEqual(5, avediff)
        self.assertTrue(need)
        need, _ = hpi.needs_HiJACK("dummy_path", 10)
        self.assertFalse(need)

    @patch("hiproc.HiPrecisionInit.needs_HiJACK", return_value=(True, 5))
    def test_check(self, needs):
        truth_thresh = 1
        yes_HiJACK, avediffs, thresh = hpi.check(
            ["dummy1", "dummy2"],
            dict(HiPrecisionInit=dict(
                Mean_Jitter_Magnitude_Threshold=truth_thresh
            ))
        )
        self.assertEqual([True, True], yes_HiJACK)
        self.assertEqual([5, 5], avediffs)
        self.assertEqual(truth_thresh, thresh)
