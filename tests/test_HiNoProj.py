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
import pkg_resources
import unittest
from pathlib import Path
from unittest.mock import call
from unittest.mock import mock_open
from unittest.mock import patch

import pvl

# import hiproc.hirise as hirise
import hiproc.HiNoProj as hnp

conf = pvl.load(
    pkg_resources.resource_stream(
        "hiproc",
        "data/HiNoProj.conf"
    )
)


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
        "StitchedProductIds": (
            "PSP_010502_2090_RED4_0, PSP_010502_2090_RED4_1"
        ),
    }
    return values[key]


class TestHiNoProjCube(unittest.TestCase):
    @patch("hiproc.HiNoProj.isis.getkey_k", side_effect=getkey)
    def test_init(self, mock_getkey):
        c = hnp.Cube("dummy/PSP_010502_2090_RED5_0")
        self.assertIsNone(c.next_path)


class TestConf(unittest.TestCase):
    def test_conf_check(self):
        self.assertIsNone(hnp.conf_check(conf["HiNoProj"]))


class TestHiNoProj(unittest.TestCase):
    @patch("hiproc.HiNoProj.Path.with_suffix")
    @patch(
        "hiproc.HiNoProj.pvl.loads",
        return_value={
            "UniversalGroundRange": {
                "MinimumLatitude": -66,
                "MaximumLatitude": 67,
            }
        },
    )
    @patch("hiproc.HiNoProj.isis.camrange")
    @patch("hiproc.HiNoProj.isis.spiceinit")
    @patch("hiproc.HiNoProj.shutil.copyfile")
    @patch("hiproc.HiNoProj.isis.getkey_k", side_effect=getkey)
    def test_is_polar(
        self, mock_getkey, m_copy, m_spice, m_cam, m_pvl, m_path
    ):
        c1 = hnp.Cube("dummy/PSP_010502_2090_RED4.HiStitch.balance.cub")
        c2 = hnp.Cube("dummy/PSP_010502_2090_RED3.HiStitch.balance.cub")

        self.assertFalse(hnp.is_polar([c1, c2], 87, "tt"))
        self.assertTrue(hnp.is_polar([c1, c2], 60, "tt"))

    @patch("hiproc.HiNoProj.isis.handmos")
    @patch("hiproc.HiNoProj.isis.getkey_k", side_effect=getkey)
    def test_handmos_side(self, mock_getkey, m_hand):
        c2 = hnp.Cube("dummy/PSP_010502_2090_RED2.HiStitch.balance.cub")
        c3 = hnp.Cube("dummy/PSP_010502_2090_RED3.HiStitch.balance.cub")
        c4 = hnp.Cube("dummy/PSP_010502_2090_RED4.HiStitch.balance.cub")
        c5 = hnp.Cube("dummy/PSP_010502_2090_RED5.HiStitch.balance.cub")
        c6 = hnp.Cube("dummy/PSP_010502_2090_RED6.HiStitch.balance.cub")

        cubes = (c2, c3, c4, c5, c6)

        for i, c in enumerate(cubes):
            cubes[i].next_path = c.path
            cubes[i].line_offset = 10
            cubes[i].samp_offset = 10

        hnp.handmos_side(cubes, c4, "dummy/mosaic.cub", True)
        hnp.handmos_side(cubes, c4, "dummy/mosaic.cub", False)
        self.assertEqual(
            m_hand.call_args_list,
            [
                call(
                    c3.path,
                    mosaic="dummy/mosaic.cub",
                    outband=1,
                    outline=11,
                    outsample=11,
                    priority="ontop",
                ),
                call(
                    c2.path,
                    mosaic="dummy/mosaic.cub",
                    outband=1,
                    outline=21,
                    outsample=21,
                    priority="ontop",
                ),
                call(
                    c5.path,
                    mosaic="dummy/mosaic.cub",
                    outband=1,
                    outline=-9,
                    outsample=-9,
                    priority="beneath",
                ),
                call(
                    c6.path,
                    mosaic="dummy/mosaic.cub",
                    outband=1,
                    outline=-19,
                    outsample=-19,
                    priority="beneath",
                ),
            ],
        )

    # @patch('hiproc.HiNoProj.isis.getkey_k', side_effect=getkey)
    # def test_fix_kernel(self, mock_getkey):
    #     c5 = hnp.Cube('dummy/PSP_010502_2090_RED5.HiStitch.balance.cub')
    #     hnp.fix_kernel(c5)

    @patch("hiproc.HiNoProj.isis.editlab")
    @patch("hiproc.HiNoProj.handmos_side")
    @patch("hiproc.HiNoProj.shutil.copyfile")
    @patch(
        "hiproc.HiNoProj.open",
        mock_open(
            read_data="""
# Average Line Offset: 5
# Average Sample Offset: 5"""
        ),
    )
    @patch("hiproc.HiNoProj.isis.hijitreg")
    @patch("hiproc.HiNoProj.isis.noproj")
    @patch("hiproc.HiNoProj.isis.spicefit")
    @patch("hiproc.HiNoProj.isis.spiceinit")
    @patch("hiproc.HiNoProj.isis.getkey_k", side_effect=getkey)
    def test_HiNoProj(
        self,
        mock_getkey,
        m_spice,
        m_sfit,
        m_noproj,
        m_hijitreg,
        m_copy,
        m_handside,
        m_edit,
    ):
        c2 = Path("dummy/PSP_010502_2090_RED2.HiStitch.balance.cub")
        c3 = Path("dummy/PSP_010502_2090_RED3.HiStitch.balance.cub")
        c4 = Path("dummy/PSP_010502_2090_RED4.HiStitch.balance.cub")
        c5 = Path("dummy/PSP_010502_2090_RED5.HiStitch.balance.cub")
        c6 = Path("dummy/PSP_010502_2090_RED6.HiStitch.balance.cub")

        cubes = (c2, c3, c4, c5, c6)

        conf = {"HiNoProj": {
            "Shape": "SYSTEM"
        }}

        hnp.HiNoProj(cubes, conf, "dummy/mosaic.cub", 4, keep=True)
        self.assertEqual(len(cubes), m_spice.call_count)
        self.assertEqual(len(cubes), m_sfit.call_count)
        self.assertEqual(len(cubes), m_noproj.call_count)
        self.assertEqual(len(cubes) - 1, m_hijitreg.call_count)
        self.assertEqual(len(cubes) + 1, m_copy.call_count)
        self.assertEqual(2, m_handside.call_count)
        self.assertEqual(
            m_edit.call_args_list,
            [
                call(
                    Path("dummy/mosaic.cub"),
                    grpname="Instrument",
                    keyword="ImageJitterCorrected",
                    option="addkey",
                    value=0,
                ),
                call(
                    Path("dummy/mosaic.cub"),
                    grpname="Archive",
                    keyword="ProductId",
                    option="modkey",
                    value="PSP_010502_2090_RED",
                ),
                call(
                    Path("dummy/mosaic.cub"),
                    grpname="Instrument",
                    keyword="MatchedCube",
                    option="modkey",
                    value="PSP_010502_2090_RED4",
                ),
                call(
                    Path("dummy/mosaic.cub"),
                    grpname="Archive",
                    keyword="SourceProductId",
                    option="ADDKEY",
                    value="(PSP_010502_2090_RED4_0, "
                    "PSP_010502_2090_RED4_1, PSP_010502_2090_RED4_0, "
                    "PSP_010502_2090_RED4_1, PSP_010502_2090_RED4_0, "
                    "PSP_010502_2090_RED4_1, PSP_010502_2090_RED4_0, "
                    "PSP_010502_2090_RED4_1, PSP_010502_2090_RED4_0, "
                    "PSP_010502_2090_RED4_1)",
                ),
            ],
        )
