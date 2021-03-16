#!/usr/bin/env python
"""This module has tests for the HiRISE HiJitReg functions."""

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
from unittest.mock import patch
from pathlib import Path

import hiproc.HiColorInit as hci
import hiproc.HiSlither as sli


def getkey(cube, group, key):
    values = {
        "ProductId": None,
        "Summing": 2,
        "Lines": 1024,
        "Samples": 1024,
        "TDI": 64,
        "StitchedProductIds": [
            "PSP_010502_2090_RED4_0",
            "PSP_010502_2090_RED4_1",
        ],
    }
    return values[key]


class TestHiSlither(unittest.TestCase):
    @patch("hiproc.HiColorInit.isis.getkey_k", side_effect=getkey)
    def test_cube_check(self, m_getkey):
        r = hci.HiColorCube("dummy/PSP_010502_2090_RED5.HiStitch.balance.cub")
        i = hci.HiColorCube(
            "dummy/PSP_010502_2090_IR11.HiStitch.balance.precolor.cub"
        )
        b = hci.HiColorCube(
            "dummy/PSP_010502_2090_BG13.HiStitch.balance.precolor.cub"
        )
        self.assertTrue(sli.cube_check(r, i, b))
        self.assertTrue(sli.cube_check(r, None, b))
        self.assertFalse(sli.cube_check(None, None, None))
        self.assertFalse(sli.cube_check(None, i, b))
        self.assertFalse(sli.cube_check(r, i, None))

    @patch("hiproc.HiColorInit.isis.getkey_k", side_effect=getkey)
    def test_get_slither_path(self, m_getkey):
        c = hci.HiColorCube(
            "dummy/PSP_010502_2090_IR10.HiStitch.balance.precolor.cub"
        )
        self.assertEqual(
            sli.get_slither_path(c),
            Path("dummy/PSP_010502_2090_IR10.slither.cub"),
        )

    @patch("hiproc.HiSlither.isis.editlab")
    @patch("hiproc.HiSlither.isis.mask")
    @patch("hiproc.HiColorInit.isis.getkey_k", side_effect=getkey)
    def test_make_dummy_IR(self, m_getkey, m_mask, m_editlab):
        r = hci.HiColorCube("dummy/PSP_010502_2090_RED4.HiStitch.balance.cub")
        b = hci.HiColorCube(
            "dummy/PSP_010502_2090_BG12.HiStitch.balance.precolor.cub"
        )
        self.assertEqual(
            sli.make_dummy_IR(r, b),
            Path("dummy/PSP_010502_2090_IR10.slither.cub"),
        )

    @patch("hiproc.HiSlither.isis.slither")
    @patch("hiproc.HiColorInit.isis.getkey_k", side_effect=getkey)
    def test_run_slither(self, m_getkey, m_slither):
        c = hci.HiColorCube(
            "dummy/PSP_010502_2090_IR10.HiStitch.balance.precolor.cub"
        )
        sli.run_slither(c)
        m_slither.called_once()

    @patch("hiproc.HiSlither.isis.trim")
    @patch("hiproc.HiSlither.isis.hicubeit")
    @patch("hiproc.HiSlither.run_slither")
    @patch("hiproc.HiColorInit.isis.getkey_k", side_effect=getkey)
    def test_process_set(self, m_getkey, m_sli, m_hicubeit, m_trim):
        r = hci.HiColorCube("dummy/PSP_010502_2090_RED5.HiStitch.balance.cub")
        i = hci.HiColorCube(
            "dummy/PSP_010502_2090_IR11.HiStitch.balance.precolor.cub"
        )
        b = hci.HiColorCube(
            "dummy/PSP_010502_2090_BG13.HiStitch.balance.precolor.cub"
        )
        self.assertEqual(
            sli.process_set(r, i, b, keep=True),
            Path("dummy/PSP_010502_2090_COLOR5.cub"),
        )
        _, trim_args = m_trim.call_args
        self.assertEqual(trim_args["left"], 3)
