#!/usr/bin/env python
"""This module has tests for the HiRISE HiColorNorm functions."""

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
from unittest.mock import call, patch, mock_open
from pathlib import Path

import pvl

import hiproc.hirise as hirise
import hiproc.HiColorNorm as hcn


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


class TestColorCube(unittest.TestCase):
    @patch("hiproc.HiColorNorm.ColorCube.get_binning", return_value=2)
    @patch("hiproc.HiColorNorm.isis.getkey_k", side_effect=getkey)
    def setUp(self, m_getkey, m_get_binning):
        self.c = hcn.ColorCube("dummy/PSP_010502_2090_COLOR5")

    def test_init(self):
        self.assertEqual(2, self.c.red_bin)

    def test_get_boxcar_size(self):
        self.assertEqual(1, self.c.get_boxcar_size(2))
        self.assertEqual(3, self.c.get_boxcar_size(3))
        self.assertEqual(5, self.c.get_boxcar_size(5))

    def test_set_crop_lines(self):
        conf = {
            "HiColorNorm": {
                "HiColorNorm_Crop_Top": 1000,
                "HiColorNorm_Crop_Bot": 2000,
            }
        }
        self.assertEqual((1, 1024), self.c.set_crop_lines(conf))

        self.c.lines = 50000
        self.assertEqual((1001, 47000), self.c.set_crop_lines(conf))


class TestHiColorNorm(unittest.TestCase):
    @patch("hiproc.HiColorNorm.ColorCube.get_binning", return_value=2)
    @patch("hiproc.HiColorNorm.isis.getkey_k", side_effect=getkey)
    def setUp(self, m_getkey, m_get_binning):
        c4 = hcn.ColorCube("dummy/PSP_010502_2090_COLOR4")
        c5 = hcn.ColorCube("dummy/PSP_010502_2090_COLOR5")
        self.cubes = [c4, c5]

    @unittest.skip("Tests on a real file.")
    def test_conf_check(self):
        conf_path = Path("data") / "HiColorNorm.conf"
        c = pvl.load(str(conf_path))
        self.assertIsNone(hcn.conf_check(c))

    def test_set_outpath(self):
        self.assertEqual(
            Path("unchanged"), hcn.set_outpath("unchanged", self.cubes)
        )
        self.assertEqual(
            Path("dummy/PSP_010502_2090_COLOR.cub"),
            hcn.set_outpath("_COLOR.cub", self.cubes),
        )
        t_getkey = getkey
        with patch("hiproc.HiColorNorm.ColorCube.get_binning", return_value=2):
            with patch(
                "hiproc.HiColorNorm.isis.getkey_k", side_effect=t_getkey
            ):
                with patch(
                    "hiproc.hirise.get_ObsID_fromfile",
                    return_value=hirise.ObservationID("ESP_060680_3180"),
                ):
                    oddball = hcn.ColorCube("dummy/ESP_060680_3180_COLOR4")
                    newlist = self.cubes + [oddball]
                    self.assertRaises(
                        ValueError, hcn.set_outpath, "unchanged", newlist
                    )

    @patch("hiproc.HiColorNorm.open", new_callable=mock_open)
    def test_FurrowCheck(self, m_open):
        self.assertFalse(hcn.FurrowCheck(self.cubes, None))

        with patch(
            "hiproc.HiColorNorm.json.load", return_value={"zapped": True}
        ):
            self.assertTrue(
                hcn.FurrowCheck(
                    self.cubes, ["dummy_db_path1", "dummy_db_path2"]
                )
            )
        with patch(
            "hiproc.HiColorNorm.json.load", return_value={"zapped": False}
        ):
            self.assertFalse(
                hcn.FurrowCheck(
                    self.cubes, ["dummy_db_path1", "dummy_db_path2"]
                )
            )

    @patch("hiproc.HiColorNorm.isis.handmos")
    def test_make_LR_mosaic(self, m_handmos):
        hcn.make_LR_mosaic("left.cub", "right.cub", 10, "mosaic.cub", 100, 20)
        expected = [
            call(
                "left.cub",
                create="YES",
                mosaic="mosaic.cub",
                nbands=1,
                nlines=100,
                nsamp=20,
                outband=1,
                outline=1,
                outsamp=1,
            ),
            call(
                "right.cub",
                mosaic="mosaic.cub",
                outband=1,
                outline=1,
                outsamp=11,
            ),
        ]
        self.assertListEqual(m_handmos.call_args_list, expected)

    @patch(
        "hiproc.HiColorNorm.csv.DictReader",
        return_value=[{"Average": 4}, {"Average": 5}],
    )
    @patch("hiproc.HiColorNorm.open", new_callable=mock_open)
    @patch("hiproc.HiColorNorm.isis.cubenorm")
    def test_cubenorm_stats(self, m_cubenorm, m_open, m_reader):
        self.assertAlmostEqual(
            0.70710678,
            hcn.cubenorm_stats(
                "crpmos.cub", "mos.cub", "mosnrm.cub", keep=True
            ),
        )
        expected = [
            call(
                "crpmos.cub",
                direction="column",
                format="table",
                norm="average",
                stats=Path("crpmos.cubenorm.txt"),
            ),
            call(
                "mos.cub",
                direction="column",
                format="table",
                fromstats=Path("crpmos.cubenorm.txt"),
                norm="average",
                preserve=True,
                statsource="table",
                to="mosnrm.cub",
            ),
        ]
        self.assertListEqual(m_cubenorm.call_args_list, expected)

    @patch("hiproc.HiColorNorm.isis.handmos")
    @patch("hiproc.HiColorNorm.isis.algebra")
    @patch("hiproc.HiColorNorm.shutil.copyfile")
    def test_make_unfiltered(self, m_copy, m_alg, m_hand):
        self.assertEqual(
            Path("in_UNFILTERED_COLOR4.cub"),
            hcn.make_unfiltered(
                "in_COLOR4.cub", "nrm.cub", "ttoken", "CC", 2, keep=True
            ),
        )
        self.assertEqual(
            m_alg.call_args_list,
            [
                call(
                    "nrm.cub",
                    from2="in_COLOR4.cub+2",
                    operator="MULTIPLY",
                    to=Path("in_COLOR4.ttoken_CC.algebra.cub"),
                )
            ],
        )
        self.assertEqual(
            m_hand.call_args_list,
            [
                call(
                    Path("in_COLOR4.ttoken_CC.algebra.cub"),
                    matchbandbin=False,
                    mosaic=Path("in_UNFILTERED_COLOR4.cub"),
                    outband=2,
                    outline=1,
                    outsample=1,
                )
            ],
        )

    @patch("hiproc.HiColorNorm.isis.lowpass")
    def test_lpfz_filtering(self, m_low):
        hcn.lpfz_filtering("from.cub", "to.cub", 3, 5)
        self.assertListEqual(
            m_low.call_args_list,
            [
                call(
                    "from.cub",
                    HIS=True,
                    HRS=True,
                    LIS=True,
                    LRS=True,
                    filter="OUTSIDE",
                    high=2.0,
                    lines=3,
                    low=0.0,
                    minimum=25,
                    minopt="PERCENTAGE",
                    null=True,
                    samples=5,
                    to="to.cub",
                )
            ],
        )

    @patch("hiproc.HiColorNorm.isis.lowpass")
    def test_lpfz_triplefilter(self, m_low):
        hcn.lpfz_triplefilter("from.cub", "to.cub", keep=True)
        self.assertListEqual(
            m_low.call_args_list,
            [
                call(
                    Path("from.cub"),
                    HIS=True,
                    HRS=True,
                    LIS=True,
                    LRS=True,
                    filter="OUTSIDE",
                    high=2.0,
                    lines=11,
                    low=0.0,
                    minimum=25,
                    minopt="PERCENTAGE",
                    null=True,
                    samples=5,
                    to=Path("from.z1.cub"),
                ),
                call(
                    Path("from.z1.cub"),
                    HIS=True,
                    HRS=True,
                    LIS=True,
                    LRS=True,
                    filter="OUTSIDE",
                    high=2.0,
                    lines=21,
                    low=0.0,
                    minimum=25,
                    minopt="PERCENTAGE",
                    null=True,
                    samples=9,
                    to=Path("from.z2.cub"),
                ),
                call(
                    Path("from.z2.cub"),
                    HIS=True,
                    HRS=True,
                    LIS=True,
                    LRS=True,
                    filter="OUTSIDE",
                    high=2.0,
                    lines=41,
                    low=0.0,
                    minimum=25,
                    minopt="PERCENTAGE",
                    null=True,
                    samples=11,
                    to="to.cub",
                ),
            ],
        )

    @patch("hiproc.HiColorNorm.isis.crop")
    @patch("hiproc.HiColorNorm.isis.mask")
    @patch("hiproc.HiColorNorm.isis.ratio")
    def test_per_color(self, m_rat, m_mask, m_crop):
        mask_p = Path("dummy/PSP_010502_2090_COLOR4.ttoken_IR.mask.cub")
        ratc_p = Path("dummy/PSP_010502_2090_COLOR4.ttoken_IR.ratcrop.cub")
        rati_p = Path("dummy/PSP_010502_2090_COLOR4.ttoken_IR.ratio.cub")

        self.assertTupleEqual(
            (mask_p, ratc_p),
            hcn.per_color(self.cubes[0], "ttoken", "IR", keep=True),
        )
        self.assertEqual(
            m_rat.call_args_list,
            [
                call(
                    den="dummy/PSP_010502_2090_COLOR4+2",
                    num="dummy/PSP_010502_2090_COLOR4+1",
                    to=rati_p,
                )
            ],
        )
        self.assertEqual(
            m_mask.call_args_list,
            [
                call(
                    rati_p,
                    mask=rati_p,
                    maximum=4.0,
                    minimum=0.0,
                    preserve="INSIDE",
                    to=mask_p,
                )
            ],
        )
        self.assertEqual(
            m_crop.call_args_list,
            [call(mask_p, line=None, nlines=None, to=ratc_p)],
        )

    @patch("hiproc.HiColorNorm.isis.handmos")
    @patch("hiproc.HiColorNorm.isis.algebra")
    @patch("hiproc.HiColorNorm.lpfz_triplefilter")
    @patch("hiproc.HiColorNorm.isis.lowpass")
    @patch("hiproc.HiColorNorm.make_unfiltered")
    @patch("hiproc.HiColorNorm.isis.crop")
    @patch("hiproc.HiColorNorm.cubenorm_stats", return_value=3)
    @patch("hiproc.HiColorNorm.make_LR_mosaic")
    def test_per_band(
        self,
        m_LR,
        m_cn,
        m_crop,
        m_unfilt,
        m_low,
        m_3filt,
        m_algebra,
        m_handmos,
    ):
        conf = {
            "HiColorNorm": {
                "HiColorNorm_Crop_Top": 1000,
                "HiColorNorm_Crop_Bot": 2000,
            }
        }
        for c in self.cubes:
            c.set_crop_lines(conf)
        self.assertEqual(
            3,
            hcn.per_band(
                self.cubes,
                Path("out_path.cub"),
                "ttoken",
                "IR",
                True,
                True,
                keep=True,
            ),
        )
        num_cubes = len(self.cubes)
        m_cn.assert_called_once()
        self.assertEqual(num_cubes, m_LR.call_count)
        self.assertEqual(num_cubes, m_crop.call_count)
        self.assertEqual(num_cubes, m_unfilt.call_count)
        self.assertEqual(num_cubes, m_low.call_count)
        self.assertEqual(num_cubes, m_3filt.call_count)
        self.assertEqual(num_cubes, m_algebra.call_count)
        self.assertEqual(num_cubes, m_handmos.call_count)

    @patch("hiproc.HiColorNorm.ColorCube.get_binning", return_value=2)
    @patch("hiproc.HiColorNorm.isis.hiccdstitch")
    @patch("hiproc.HiColorNorm.per_band", return_value=3)
    @patch(
        "hiproc.HiColorNorm.per_color",
        return_value=("CCmask.cub", "CCcrop.cub"),
    )
    @patch("hiproc.HiColorNorm.isis.cubeit")
    @patch("hiproc.HiColorNorm.Path.write_text")
    @patch("hiproc.HiColorNorm.isis.mask")
    @patch("hiproc.HiColorNorm.isis.getkey_k", side_effect=getkey)
    def test_HiColorNorm(
        self, m_getkey, m_mask, m_Pathwrite, m_cubeit, m_pc, m_pb, m_hcs, m_gb
    ):
        conf = {
            "HiColorNorm": {
                "HiColorNorm_Crop_Top": 1000,
                "HiColorNorm_Crop_Bot": 2000,
                "HiColorNorm_Make_Stitch": True,
                "HiColorNorm_NoiseFilter_IR10": True,
                "HiColorNorm_Normalization_Minimum": 0,
                "HiColorNorm_Normalization_Maximum": 16000,
            }
        }
        self.assertEqual(
            (3, 3),
            hcn.HiColorNorm(
                [
                    Path("dummy/PSP_010502_2090_COLOR4"),
                    Path("dummy/PSP_010502_2090_COLOR5")
                ],
                "out_path.cub",
                conf,
                make_unfiltered=True,
                keep=True,
            ),
        )
        self.assertEqual(6, m_mask.call_count)
        self.assertEqual(2, m_cubeit.call_count)
        self.assertEqual(4, m_pc.call_count)
        self.assertEqual(2, m_pb.call_count)
