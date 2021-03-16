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
import subprocess
import unittest
from pathlib import Path
from unittest.mock import call
from unittest.mock import mock_open
from unittest.mock import patch

import pvl

# import hiproc.hirise as hirise
import hiproc.HiccdStitch as hcs

conf = pvl.load(
    pkg_resources.resource_stream(
        "hiproc",
        "data/HiccdStitch.conf"
    )
)


def getkey(cube, group, key):
    values = {
        "ProductId": None,
        "Summing": 2,
        "Lines": 1024,
        "Samples": 1024,
        "Special_Processing_Flag": "NOMINAL",
    }
    return values[key]


class TestHiccdStitchCube(unittest.TestCase):
    @patch("hiproc.HiStitch.isis.getkey_k", side_effect=getkey)
    def test_init(self, mock_getkey):
        c = hcs.HiccdStitchCube("dummy/PSP_010502_2090_RED5_0")
        self.assertTrue(c.ns, 1024)

    @patch("hiproc.HiStitch.isis.getkey_k", side_effect=getkey)
    def test_set_cubenorm_lines(self, mock_getkey):

        c = hcs.HiccdStitchCube("dummy/PSP_010502_2090_RED5_0")
        self.assertTrue(
            (1.5, 1023.0), c.set_cubenorm_lines(1, 1, None, None, 1)
        )
        self.assertRaises(
            ValueError, c.set_cubenorm_lines, 1, 1, 2000, None, 2
        )
        self.assertRaises(
            ValueError, c.set_cubenorm_lines, 1, 1, None, 2000, 2
        )
        self.assertRaises(
            RuntimeError, c.set_cubenorm_lines, 1, 1, 1000, 800, 2
        )

    @patch("hiproc.HiStitch.isis.getkey_k", side_effect=getkey)
    def test_set_balance(self, mock_getkey):
        c = hcs.HiccdStitchCube("dummy/PSP_010502_2090_RED5_0")
        c.set_balance(1, 1, {2: [1, 1]}, 1)
        self.assertEqual(2, c.ss_balance_left)
        self.assertEqual(1, c.ns_balance_left)
        self.assertEqual(1023, c.ss_balance_right)
        self.assertEqual(1, c.ns_balance_right)
        self.assertEqual(1, c.sl_balance)
        self.assertEqual(0, c.nl_balance)

    @patch("hiproc.HiStitch.isis.getkey_k", side_effect=getkey)
    def test_set_ls_path(self, mock_getkey):
        c = hcs.HiccdStitchCube("dummy/PSP_010502_2090_RED5_0")
        p = "dummy-path"
        c.set_ls_path(p)
        self.assertEqual(p, c.ls_path)
        self.assertEqual(p, c.lm_path)

    @patch("hiproc.HiStitch.isis.getkey_k", side_effect=getkey)
    def test_set_rs_path(self, mock_getkey):
        c = hcs.HiccdStitchCube("dummy/PSP_010502_2090_RED5_0")
        p = "dummy-path"
        c.set_rs_path(p)
        self.assertEqual(p, c.rs_path)
        self.assertEqual(p, c.rm_path)

    @patch("hiproc.HiStitch.isis.getkey_k", side_effect=getkey)
    def test_gather_from_db(self, mock_getkey):
        c = hcs.HiccdStitchCube("tmp/PSP_010502_2090_RED5_0")
        db = {"hical_status": "BadCal", "IMAGE_SIGNAL_TO_NOISE_RATIO": 5}
        with patch("hiproc.HiccdStitch.json.load", return_value=db):
            with patch(
                "hiproc.HiccdStitch.Path.glob", return_value=["d", "d"]
            ):
                with patch("hiproc.HiccdStitch.open", mock_open()):
                    c.gather_from_db()
                    self.assertEqual("BadCal", c.hical_status)
                    self.assertListEqual([5.0, 5.0], c.snr_list)


class TestConf(unittest.TestCase):
    def test_conf_check(self):
        self.assertIsNone(hcs.conf_check(conf))

    def test_make_area_dict(self):
        # print(hcs.make_area_dict(c['HiccdStitch']))
        truth = {
            1: [12, 36],
            2: [6, 18],
            3: [5, 15],
            4: [5, 9],
            8: [5, 5],
            16: [1, 2],
        }
        self.assertDictEqual(truth, hcs.make_area_dict(conf["HiccdStitch"]))


class TestBasic(unittest.TestCase):
    @patch("hiproc.HiStitch.isis.getkey_k", side_effect=getkey)
    def test_set_outpath(self, mock_getkey):
        c1 = hcs.HiccdStitchCube("dummy/PSP_010502_2090_RED5_0")
        c2 = hcs.HiccdStitchCube("dummy/PSP_010502_2090_RED3_1")
        diff = hcs.HiccdStitchCube("tmp/PSP_005632_1225_RED5_1")

        self.assertRaises(ValueError, hcs.set_outpath, "dummy", [c1, c2, diff])

        out = Path("foo")
        self.assertEqual(out, hcs.set_outpath(out, [c1, c2]))

        truth = Path("dummy/PSP_010502_2090_RED.foo")
        self.assertEqual(truth, hcs.set_outpath(".foo", [c1, c2]))

    @patch("hiproc.HiStitch.isis.enlarge")
    @patch("hiproc.HiStitch.isis.crop")
    @patch("hiproc.HiStitch.isis.getkey_k", side_effect=getkey)
    def test_crop_and_scale(self, mock_getkey, mock_crop, moc_enlarge):
        c1 = hcs.HiccdStitchCube("dummy/PSP_010502_2090_RED5_0")
        c2 = hcs.HiccdStitchCube("dummy/PSP_010502_2090_RED3_1")
        c1.smag = 1
        c1.lmag = 1
        c2.smag = 1
        c2.lmag = 1
        new_c, to_del = hcs.crop_and_scale([c1, c2])
        p = Path("dummy/PSP_010502_2090_RED5_0.left.crop.cub")
        self.assertEqual(p, new_c[0].ls_path)
        self.assertIn(p, to_del)

    @patch("hiproc.HiStitch.isis.getkey_k", side_effect=getkey)
    def test_get_group_i(self, mock_getkey):
        c3 = hcs.HiccdStitchCube("dummy/PSP_010502_2090_RED3")
        c5 = hcs.HiccdStitchCube("dummy/PSP_010502_2090_RED5")
        c6 = hcs.HiccdStitchCube("dummy/PSP_010502_2090_RED6")
        c3.rstats = 1
        c3.lstats = 1
        c5.rstats = 1
        c5.lstats = 1
        c6.lstats = 1
        c6.rstats = 1
        pairs = hcs.get_group_i([c5, c3, c6])
        self.assertListEqual([3], pairs[0][1])
        self.assertListEqual([5, 6], pairs[1][1])

    # Not testing hcs.get_stats(), because it just straight up calls ISIS
    # routines.

    @patch("hiproc.HiStitch.isis.getkey_k", side_effect=getkey)
    def test_get_correction(self, mock_getkey):
        c5 = hcs.HiccdStitchCube("dummy/PSP_010502_2090_RED5")
        c6 = hcs.HiccdStitchCube("dummy/PSP_010502_2090_RED6")
        c5.correction = 2
        c5.rstats = 2
        c6.lstats = 2
        self.assertEqual(2, hcs.get_correction(c6, c5, "MULTIPLY", 6))
        self.assertEqual(1, hcs.get_correction(c6, c5, "MULTIPLY", 0))
        self.assertEqual(0, hcs.get_correction(c6, c5, "ADD", 0))

    @patch("hiproc.HiStitch.isis.getkey_k", side_effect=getkey)
    def test_get_normalization(self, mock_getkey):
        c1 = hcs.HiccdStitchCube("dummy/PSP_010502_2090_RED1")
        c2 = hcs.HiccdStitchCube("dummy/PSP_010502_2090_RED2")
        c3 = hcs.HiccdStitchCube("dummy/PSP_010502_2090_RED3")
        c4 = hcs.HiccdStitchCube("dummy/PSP_010502_2090_RED4")
        c5 = hcs.HiccdStitchCube("dummy/PSP_010502_2090_RED5")
        c6 = hcs.HiccdStitchCube("dummy/PSP_010502_2090_RED6")
        c1.correction = 1
        c2.correction = 2
        c3.correction = 3
        c4.correction = 4
        c5.correction = 5
        c6.correction = 6
        self.assertEqual(
            5,
            hcs.get_normalization(
                [c1, c2, c3, c4, c5, c6], [1, 2, 3, 4, 5, 6], -1, [5, 11, 13]
            ),
        )
        self.assertEqual(
            3, hcs.get_normalization([c2, c3, c4], [2, 3, 4], -2, [5, 11, 13])
        )


class TestMock(unittest.TestCase):
    def setUp(self):
        self.cubenorm_data = """    Band  RowCol    ValidPoints        Average         Median         StdDev        Minimum        Maximum
       1       1           4000       0.154037       0.154375     0.00549645       0.132305       0.173675
       1       2           4000       0.151859       0.152177     0.00559696       0.127681       0.172393"""
        hiconf = {
            "HiccdStitch_Bin1_Skip_Top_Lines": 1,
            "HiccdStitch_Bin1_Skip_Bot_Lines": 1,
            "HiccdStitch_Bin01_Area": [1, 1],
            "HiccdStitch_Bin02_Area": [1, 1],
            "HiccdStitch_Bin03_Area": [1, 1],
            "HiccdStitch_Bin04_Area": [1, 1],
            "HiccdStitch_Bin08_Area": [1, 1],
            "HiccdStitch_Bin16_Area": [1, 1],
            "HiccdStitch_Balance": True,
            "HiccdStitch_Balance_Correction": "MULTIPLY",
            "HiccdStitch_Cubenorm_Method": "DIVIDE",
            "HiccdStitch_Interpolation": "BILINEAR",
            "HiccdStitch_Normalization_Minimum": 0.0,
            "HiccdStitch_Normalization_Maximum": 1.5,
            "HiccdStitch_Control_CCD": [5, 11, 13],
            "HiccdStitch_SNR_Threshold": 50,
            "HiccdStitch_Clean": "KEEP",
            "HiccdStitch_Version_Enable": False,
            "HiccdStitch_Balance_Method": "AVERAGE"
        }
        self.conf = {"HiccdStitch": hiconf}

    @patch("hiproc.HiStitch.isis.getkey_k", side_effect=getkey)
    def test_GetImageDims(self, mock_getkey):
        c1 = hcs.HiccdStitchCube("dummy/PSP_010502_2090_RED5_0")
        c2 = hcs.HiccdStitchCube("dummy/PSP_010502_2090_RED3_1")
        cubes = [c1, c2]
        new_c = hcs.GetImageDims(cubes, self.conf, None, None)
        self.assertEqual(1, new_c[0].smag)
        self.assertEqual(1, new_c[0].lmag)

    def test_AnalyzeStats(self):
        m = mock_open(read_data=self.cubenorm_data)
        m.return_value.__iter__ = lambda self: self
        m.return_value.__next__ = lambda self: next(iter(self.readline, ""))
        with patch("hiproc.HiccdStitch.open", m):
            self.assertAlmostEqual(
                2.371842000000029e-06,
                hcs.AnalyzeStats("dummy_cubenormout", "dummy_filtered"),
                places=10,
            )

    @patch("hiproc.HiStitch.isis.PathSet.unlink")
    @patch("hiproc.HiStitch.isis.getkey_k", side_effect=getkey)
    @patch("hiproc.HiStitch.isis.cubenorm")
    @patch("hiproc.HiStitch.isis.crop")
    def test_CubeNormStep(
        self, mock_crop, mock_cubenorm, mock_getkey, mock_PathSetunlink
    ):
        m = mock_open(read_data=self.cubenorm_data)
        m.return_value.__iter__ = lambda self: self
        m.return_value.__next__ = lambda self: next(iter(self.readline, ""))
        with patch("hiproc.HiccdStitch.open", m):
            c = hcs.HiccdStitchCube("dummy/PSP_010502_2090_RED5_0")
            new_c = hcs.CubeNormStep(c, self.conf["HiccdStitch"])
            self.assertEqual(
                Path("dummy/PSP_010502_2090_RED5_0.cubenorm.cub"),
                new_c.nextpath,
            )

    @patch("hiproc.HiStitch.isis.algebra")
    @patch("hiproc.HiStitch.isis.getkey_k", side_effect=getkey)
    def test_make_balance(self, mock_getkey, mock_algebra):
        c = hcs.HiccdStitchCube("dummy/PSP_010502_2090_RED5_0")
        c.correction = 1
        hcs.make_balance(
            c, self.conf["HiccdStitch"], Path("dummy_balance.cub")
        )
        mock_algebra.assert_has_calls(
            [
                call(
                    Path("dummy/PSP_010502_2090_RED5_0"),
                    a=1,
                    c=0,
                    operator="unary",
                    to="dummy_balance.cub+SignedWord+0.0:1.5",
                )
            ]
        )

        mock_algebra.reset_mock()
        new_conf = self.conf["HiccdStitch"]
        new_conf["HiccdStitch_Balance_Correction"] = "ADD"
        hcs.make_balance(c, new_conf, Path("dummy_balance.cub"))
        mock_algebra.assert_has_calls(
            [
                call(
                    Path("dummy/PSP_010502_2090_RED5_0"),
                    a=1,
                    c=1,
                    operator="unary",
                    to="dummy_balance.cub+SignedWord+0.0:1.5",
                )
            ]
        )

    @patch(
        "hiproc.HiccdStitch.pvl.loads",
        return_value={"Results": {"Average": 1}},
    )
    @patch("hiproc.HiStitch.isis.PathSet.unlink")
    @patch("hiproc.HiccdStitch.isis.stats")
    @patch("hiproc.HiccdStitch.isis.enlarge")
    @patch("hiproc.HiccdStitch.isis.crop")
    @patch("hiproc.HiccdStitch.isis.algebra")
    @patch("hiproc.HiccdStitch.isis.mask")
    @patch("hiproc.HiccdStitch.isis.getkey_k", side_effect=getkey)
    def test_BalanceStep(
        self,
        m_getkey,
        m_mask,
        m_algebra,
        m_crop,
        m_enlarge,
        m_stats,
        m_unlink,
        m_loads,
    ):
        c1 = hcs.HiccdStitchCube("dummy/PSP_010502_2090_RED1")
        c2 = hcs.HiccdStitchCube("dummy/PSP_010502_2090_RED2")
        c3 = hcs.HiccdStitchCube("dummy/PSP_010502_2090_RED3")
        c4 = hcs.HiccdStitchCube("dummy/PSP_010502_2090_RED4")
        c5 = hcs.HiccdStitchCube("dummy/PSP_010502_2090_RED5")
        c6 = hcs.HiccdStitchCube("dummy/PSP_010502_2090_RED6")

        cubes = hcs.BalanceStep(
            [c1, c2, c3, c4, c5, c6], self.conf["HiccdStitch"]
        )
        self.assertEqual(1, cubes[0].lstats)
        self.assertEqual(1, cubes[0].correction)
        self.assertEqual(
            Path("dummy/PSP_010502_2090_RED1.balance.cub"), cubes[0].nextpath
        )

    @patch("hiproc.HiccdStitch.isis.editlab")
    @patch("hiproc.HiccdStitch.isis.getkey_k", side_effect=getkey)
    def test_SpecialProcessingFlags(self, m_getkey, m_editlab):
        p = Path("dummy/PSP_010502_2090_RED1")
        c = hcs.HiccdStitchCube(p)
        hcs.SpecialProcessingFlags(c)
        m_editlab.assert_has_calls(
            [
                call(
                    p,
                    grpname="Instrument",
                    keyword="Special_Processing_Flag",
                    option="MODKEY",
                    value="NOMINAL",
                )
            ]
        )
        m_editlab.reset_mock()
        c.hical_status = "BadCal"
        hcs.SpecialProcessingFlags(c)
        m_editlab.assert_has_calls(
            [
                call(
                    p,
                    grpname="Instrument",
                    keyword="Special_Processing_Flag",
                    option="MODKEY",
                    value="BADCAL",
                )
            ]
        )
        m_editlab.reset_mock()
        c.cubenormstep = True
        hcs.SpecialProcessingFlags(c)
        m_editlab.assert_has_calls(
            [
                call(
                    p,
                    grpname="Instrument",
                    keyword="Special_Processing_Flag",
                    option="MODKEY",
                    value="CUBENORM",
                )
            ]
        )
        with patch(
            "hiproc.HiccdStitch.isis.getkey_k",
            side_effect=subprocess.CalledProcessError(cmd="foo", returncode=1),
        ):
            m_editlab.reset_mock()
            hcs.SpecialProcessingFlags(c)
            m_editlab.assert_has_calls(
                [
                    call(
                        p,
                        grpname="Instrument",
                        keyword="Special_Processing_Flag",
                        option="ADDKEY",
                        value="CUBENORM",
                    )
                ]
            )

    @patch("hiproc.HiccdStitch.isis.getkey_k", side_effect=getkey)
    def test_SNR_Check(self, m_getkey):
        c1 = hcs.HiccdStitchCube("dummy/PSP_010502_2090_RED1")
        c2 = hcs.HiccdStitchCube("dummy/PSP_010502_2090_RED2")
        self.assertIsNone(hcs.SNR_Check([c1, c2], 50))
        c1.bin = c2.bin = 1
        c1.snr_list = [5]
        c2.snr_list = [3, 4]
        with self.assertLogs(level="WARN"):
            self.assertIsNone(hcs.SNR_Check([c1, c2], 50))

    @patch("hiproc.HiccdStitch.SNR_Check")
    @patch("hiproc.HiccdStitch.isis.hiccdstitch")
    @patch("hiproc.HiccdStitch.SpecialProcessingFlags")
    @patch(
        "hiproc.HiccdStitch.BalanceStep",
        side_effect=lambda a, b, keep=False: a,
    )
    @patch("hiproc.HiccdStitch.CubeNormStep", side_effect=lambda a, b, c: a)
    @patch("hiproc.HiStitch.isis.PathSet.unlink")
    @patch("hiproc.HiccdStitch.isis.getkey_k", side_effect=getkey)
    def test_HiccdStitch(
        self,
        m_getkey,
        m_unlink,
        m_CubeNormStep,
        m_Balance,
        m_Special,
        m_hiccdstitch,
        m_SNR,
    ):
        c1 = hcs.HiccdStitchCube("dummy/PSP_010502_2090_RED1")
        c2 = hcs.HiccdStitchCube("dummy/PSP_010502_2090_RED2")
        p = "out-dummy/foo.cub"
        hcs.HiccdStitch([c1, c2], p, self.conf)
        m_CubeNormStep.assert_not_called()
        m_Balance.assert_has_calls(
            [call([c1, c2], self.conf["HiccdStitch"], keep=False)]
        )
        m_Special.assert_has_calls([call(c1), call(c2)])
        m_hiccdstitch.assert_called_once()
        m_SNR.assert_has_calls([call([c1, c2], 50)])
