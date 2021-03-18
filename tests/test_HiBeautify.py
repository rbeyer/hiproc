#!/usr/bin/env python
"""This module has tests for the HiRISE HiBeautify functions."""

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

import unittest
from unittest.mock import call, patch
from pathlib import Path

import hiproc.HiBeautify as hbeaut


def getkey(cube, group, key):
    values = {
        "ObservationId": "PSP_010502_2090",
        "ProductId": None,
        "Summing": 2,
        "Lines": 1024,
        "Samples": 1024,
        "Bands": 3,
        "TDI": 64,
        "Center": "(900, 700, 500) <NANOMETERS>",
        "SourceProductId": ("PSP_010502_2090_RED4_0, PSP_010502_2090_RED4_1"),
    }
    return values[key]


class TestHiBeautify(unittest.TestCase):
    # @patch("hiproc.HiColorNorm.ColorCube.get_binning", return_value=2)
    # @patch("hiproc.HiColorNorm.isis.getkey_k", side_effect=getkey)
    # def setUp(self, m_getkey, m_get_binning):
    #     c4 = hcn.ColorCube("dummy/PSP_010502_2090_COLOR4.HiColorNorm")
    #     c5 = hcn.ColorCube("dummy/PSP_010502_2090_COLOR5.HiColorNorm")
    #     self.cubes = [c4, c5]

    @patch("hiproc.HiColorNorm.ColorCube.get_binning", return_value=2)
    @patch("hiproc.HiColorNorm.isis.getkey_k", side_effect=getkey)
    @patch("hiproc.HiBeautify.isis.cubeit_k")
    @patch("hiproc.HiBeautify.isis.algebra")
    @patch("hiproc.HiBeautify.isis.handmos")
    @patch("hiproc.HiBeautify.isis.editlab")
    def test_HiBeautify(
        self,
        m_editlab,
        m_handmos,
        m_algebra,
        m_cubeit_k,
        m_getkey,
        m_get_binning
    ):
        cubes = [
            Path("dummy/PSP_010502_2090_COLOR4.HiColorNorm"),
            Path("dummy/PSP_010502_2090_COLOR5.HiColorNorm")
        ]
        conf = {
            "Beautify": {
                "Synthetic_A_Coefficient": 2,
                "Synthetic_B_Coefficient": 0.3,
            }
        }
        outirb = Path("outirb.cub")
        outrgb = Path("outrgb.cub")
        hbeaut.HiBeautify(cubes, conf, outirb, outrgb, keep=True)

        self.assertEqual(
            m_handmos.call_args_list,
            [
                call(
                    cubes[0],
                    create="Y",
                    mosaic=outirb,
                    nbands=3,
                    nlines=1024,
                    nsamp=2024.0,
                    outband=1,
                    outline=1,
                    outsample=1,
                ),
                call(
                    cubes[1],
                    mosaic=outirb,
                    outband=1,
                    outline=1,
                    outsample=1001,
                ),
            ],
        )

        args, kwargs = m_algebra.call_args
        red = f"{outirb}+2"
        bg = f"{outirb}+3"
        blue = m_algebra.call_args[1]["to"]
        self.assertEqual(m_algebra.call_args[0][0], bg)
        self.assertEqual(m_algebra.call_args[1]["A"], 2)
        self.assertEqual(m_algebra.call_args[1]["B"], 0.3)
        self.assertEqual(m_algebra.call_args[1]["from2"], red)
        self.assertEqual(m_algebra.call_args[1]["op"], "subtract")
        self.assertTrue(str(blue).endswith("_B.cub"))

        self.assertEqual(
            m_cubeit_k.call_args_list, [call([red, bg, blue], to=outrgb)]
        )
