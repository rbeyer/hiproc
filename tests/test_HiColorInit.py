#!/usr/bin/env python
"""This module has tests for the HiRISE HiccdStitch functions."""

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

# import contextlib
import unittest
from pathlib import Path
from unittest.mock import patch

import hiproc.HiColorInit as hci


def getkey(cube, group, key):
    values = {
        "ProductId": None,
        "Summing": 2,
        "Lines": 1024,
        "Samples": 1024,
        "TDI": 64,
    }
    return values[key]


class TestHiColorCube(unittest.TestCase):
    @patch("hiproc.HiColorInit.isis.getkey_k", side_effect=getkey)
    def test_init(self, mock_getkey):
        c = hci.HiColorCube("dummy/PSP_010502_2090_RED5_0")
        self.assertTrue(c.tdi, 64)


class TestMock(unittest.TestCase):
    @patch("hiproc.HiColorInit.isis.getkey_k", side_effect=getkey)
    def test_separate_ccds(self, mock_getkey):
        c04 = hci.HiColorCube("dummy/PSP_010502_2090_RED4")
        c05 = hci.HiColorCube("dummy/PSP_010502_2090_RED5")
        c10 = hci.HiColorCube("dummy/PSP_010502_2090_IR10")
        c11 = hci.HiColorCube("dummy/PSP_010502_2090_IR11")
        c12 = hci.HiColorCube("dummy/PSP_010502_2090_BG12")
        c13 = hci.HiColorCube("dummy/PSP_010502_2090_BG13")

        (red4, red5, ir10, ir11, bg12, bg13) = hci.separate_ccds(
            [c04, c05, c10, c11, c12, c13]
        )
        self.assertEqual(c04, red4)
        self.assertEqual(c05, red5)
        self.assertEqual(c10, ir10)
        self.assertEqual(c11, ir11)
        self.assertEqual(c12, bg12)
        self.assertEqual(c13, bg13)

    @patch("hiproc.HiColorInit.isis.reduce")
    @patch("hiproc.HiColorInit.isis.enlarge")
    @patch("hiproc.HiColorInit.isis.handmos")
    @patch("hiproc.HiColorInit.isis.editlab")
    @patch("hiproc.HiColorInit.Path.unlink")
    @patch("hiproc.HiColorInit.isis.getkey_k", side_effect=getkey)
    def test_HiColorInit_reduce1(
        self, m_getkey, m_unlink, m_editlab, m_handmos, m_enlarge, m_reduce
    ):
        red = Path("dummy/PSP_010502_2090_RED4")
        ir = Path("dummy/PSP_010502_2090_IR10")
        bg = Path("dummy/PSP_010502_2090_BG12")
        suffix = ".dummy.cub"
        hci.HiColorInit([red, ir, bg], suffix, keep=True)
        self.assertEqual(2, m_reduce.call_count)
        self.assertAlmostEqual(m_reduce.call_args[1]["sscale"], 1.0006)
        m_enlarge.assert_not_called()
        self.assertEqual(2, m_editlab.call_count)
        self.assertTrue(str(m_editlab.call_args[0][0]).endswith(suffix))
        self.assertAlmostEqual(m_handmos.call_args[1]["nlines"], 1024)
        m_unlink.assert_not_called()

    @patch("hiproc.HiColorInit.isis.reduce")
    @patch("hiproc.HiColorInit.isis.enlarge")
    @patch("hiproc.HiColorInit.isis.handmos")
    @patch("hiproc.HiColorInit.isis.editlab")
    @patch("hiproc.HiColorInit.Path.unlink")
    @patch("hiproc.HiColorInit.isis.getkey_k", side_effect=getkey)
    def test_HiColorInit_reduce2(
        self, m_getkey, m_unlink, m_editlab, m_handmos, m_enlarge, m_reduce
    ):
        red = Path("dummy/PSP_010502_2090_RED4")
        ir = Path("dummy/PSP_010502_2090_IR10")
        s = ".dummy.cub"
        hci.HiColorInit([red, ir], s)
        m_reduce.assert_called_once()
        self.assertAlmostEqual(1.0006, m_reduce.call_args[1]["sscale"])
        m_enlarge.assert_not_called()
        m_unlink.assert_called_once()
