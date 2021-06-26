#!/usr/bin/env python
"""This module has tests for the bitflips functions."""

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

import collections
import unittest
from pathlib import Path
from unittest.mock import mock_open, patch

import numpy as np
from scipy.signal import find_peaks

import kalasiris as isis

from hiproc import bitflips as bf

# from .utils import resource_check as rc

# # Hardcoding this, but I sure would like a better solution.
# HiRISE_img = Path('test-resources') / 'PSP_010502_2090_RED5_0.img'
# img = HiRISE_img
#
#
# class TestResources(unittest.TestCase):
#     '''Establishes that the test image exists.'''
#
#     def test_resources(self):
#         (truth, test) = rc(img)
#         self.assertEqual(truth, test)


class TestHist(unittest.TestCase):
    def setUp(self):
        HistRow = collections.namedtuple("HistRow", ["DN", "Pixels"])
        self.hist = [
            HistRow(1, 1),  # 0 *
            HistRow(2, 2),  # 1 **
            HistRow(5, 3),  # 2 ***
            HistRow(6, 4),  # 3 ****
            HistRow(7, 3),  # 4 ***
            HistRow(8, 2),  # 5 **
            HistRow(10, 1),  # 6 *
            HistRow(11, 2),  # 7 **
            HistRow(12, 2),  # 8 **
            # Biggest gap here between 12 and 20
            HistRow(20, 1),  # 9 *
        ]

    def test_find_select_idx_name(self):
        self.assertRaises(ValueError, bf.find_select_idx_name, 1, 1, True)
        self.assertEqual(bf.find_select_idx_name(1, 2, True), (-1, "max"))
        self.assertEqual(bf.find_select_idx_name(1, 2, False), (0, "max"))
        self.assertEqual(bf.find_select_idx_name(2, 1, True), (0, "min"))
        self.assertEqual(bf.find_select_idx_name(2, 1, False), (-1, "min"))

    def test_find_minima_index(self):
        pixel_counts = np.fromiter((int(x.Pixels) for x in self.hist), int)
        minima_i, _ = find_peaks(np.negative(pixel_counts))

        self.assertEqual(
            0, bf.find_minima_index(4, 0, minima_i, pixel_counts, default=0)
        )
        self.assertEqual(6, bf.find_minima_index(2, 8, minima_i, pixel_counts))

    def test_find_smart_window(self):

        hist_list = sorted(self.hist, key=lambda x: int(x.DN))
        pixel_counts = np.fromiter((int(x.Pixels) for x in hist_list), int)
        dn = np.fromiter((int(x.DN) for x in hist_list), int)

        self.assertEqual(
            (2, 12),
            bf.find_smart_window(dn, pixel_counts, 1, 12, 6, plot=False),
        )


class TestArrays(unittest.TestCase):
    def test_median_limit(self):
        data = np.array([[1, 2, 3], [3, 4, 5]])

        self.assertEqual(3, bf.median_limit(np.median(data), data))
        self.assertEqual(
            1.5, bf.median_limit(np.median(data), data, limit=2.9)
        )

    def test_medstd(self):
        data = [[1, 3, 3], [9, 5, 6]]

        nomask = np.ma.array(data, mask=np.ma.nomask)
        nm_vp = nomask.count(axis=0)
        nm_sd = np.std(nomask, axis=0)

        mask = np.ma.array(data, mask=[[1, 0, 0], [0, 0, 1]])
        ma_vp = mask.count(axis=0)
        ma_sd = np.std(mask, axis=0)

        self.assertEqual(1.5, bf.median_std(nm_vp, nm_sd))
        self.assertEqual(1, bf.median_std(ma_vp, ma_sd))

    def test_medstd_from_ma(self):
        data = [[1, 3, 3], [9, 5, 6]]

        nomask = np.ma.array(data, mask=np.ma.nomask)
        mask = np.ma.array(data, mask=[[1, 0, 0], [0, 0, 1]])
        self.assertEqual(1.5, bf.median_std_from_ma(nomask))
        self.assertEqual(1, bf.median_std_from_ma(mask))

    # Need to upgrade this test because of find_smart_index
    # def test_clean_array(self):
    #     data = [[1, 3, 3], [9, 5, 6]]

    #     arr = np.ma.array(data)
    #     np.testing.assert_array_equal(
    #         np.ma.masked_outside(arr, 2.9, 5.9), bf.clean_array(arr, width=1)
    #     )

    #     masked = np.ma.array(data, mask=[[1, 0, 0], [0, 0, 1]])
    #     np.testing.assert_array_equal(
    #         np.ma.masked_outside(arr, 2.9, 5.9),
    #         bf.clean_array(masked, width=1),
    #     )

    def test_min_max_ex(self):
        self.assertEqual(bf.min_max_ex(5, 2, 2, 400, 64), (1, 9, 16))
        self.assertEqual(bf.min_max_ex(5, 20, 2, 400, 64), (-35, 45, 20))
        self.assertEqual(bf.min_max_ex(100, 400, 2, 400, 64), (-700, 900, 400))
        self.assertEqual(bf.min_max_ex(100, 401, 2, 400, 64), (-28, 228, 16))

    # Need to upgrade: find_smart_window got complicated, this test is not.
    # def test_find_smart_window_from_ma(self):
    #     data = [[1, 3, 3, 4, 4], [4, 4, 9, 5, 6]]

    #     masked = np.ma.array(data, mask=[[1, 0, 0, 0, 0], [0, 0, 0, 0, 1]])
    #     self.assertEqual(
    #         (3, 9),
    #         bf.find_smart_window_from_ma(masked, width=1, axis=0, plot=False),
    #     )

    def test_mask_lists(self):
        data = [[1, 3, 3], [9, 5, 6]]
        arr = np.ma.array(data, mask=[[1, 0, 0], [0, 0, 1]])

        np.testing.assert_array_equal(
            [[0, 3, 3], [9, 5, 0]],
            bf.apply_special_pixels(arr, isis.specialpixels.UnsignedByte),
        )


class TestMock(unittest.TestCase):
    def setUp(self):
        self.cal = [
            [0, 1000, 1001, 991, 1000, 0],
            [1, 1000, 1002, 992, 1000, 1000],
            [2, 1000, 1003, 993, 1000, 1000],
            [3, 1000, 1004, 994, 1000, 1000],
            [4, 1000, 1005, 995, 1000, 1000],
            [5, 1000, 1011, 990, 1000, 1000],
            [6, 1000, 1010, 990, 1000, 1000],
            [7, 1000, 1010, 990, 1000, 1000],
            [8, 1000, 1010, 990, 1000, 1000],
            [9, 999, 1010, 990, 1001, 1001],
            [10, 1000, 1010, 990, 1000, 1000],
            [11, 1000, 1010, 990, 1000, 1000],
            [12, 1000, 1010, 990, 1000, 1000],
            [13, 1000, 1010, 990, 1000, 1000],
            [14, 1000, 1010, 990, 1000, 1000],
            [15, 1000, 1010, 990, 1000, 1000],
            [16, 1000, 1010, 990, 1000, 1000],
            [17, 1000, 1010, 990, 1000, 1000],
            [18, 1000, 1010, 990, 1000, 1000],
            [19, 1000, 1010, 990, 1000, 1000],
            [20, 1000, 1010, 990, 1000, 1000],
        ]
        self.cleaned = list()
        for row in self.cal:
            # Values less than 3 are less than the specialpixel Min
            # and the reverse-clocked area is only the first 20 lines.
            if 3 <= row[0] < 20:
                self.cleaned.append([0] + row[1:])
            else:
                self.cleaned.append(row)

    def test_clean_cal_tables(self):
        for notimpl in ("mask_area", "ramp_area"):
            self.assertRaises(
                NotImplementedError,
                bf.clean_cal_tables,
                "dummy-in.cub",
                "dummy-out.cub",
                **{notimpl: True}
            )

        masked = np.ma.masked_outside(np.array(self.cal), 3, 65522)

        np.testing.assert_array_equal(
            self.cleaned, bf.clean_cal_tables(masked, 1)
        )

    @patch("hiproc.bitflips.shutil.copy")
    @patch(
        "hiproc.bitflips.pvl.load",
        return_value={
            "IsisCube": {
                "Core": {"Pixels": {"Type": "UnsignedWord"}},
                "Instrument": {"Summing": 1},
            }
        },
    )
    @patch("hiproc.bitflips.open", new_callable=mock_open)
    @patch("hiproc.bitflips.isis.cube.overwrite_table")
    def test_clean_tables_from_cube(self, m_wtable, m_open, m_pvl, m_copy):
        for notimpl in ("mask_area", "ramp_area", "buffer_area", "dark_area"):
            self.assertRaises(
                NotImplementedError,
                bf.clean_tables_from_cube,
                "dummy-in.cub",
                "dummy-out.cub",
                **{notimpl: True}
            )

        with patch(
            "hiproc.bitflips.isis.cube.get_table",
            return_value={"Calibration": self.cal},
        ):
            bf.clean_tables_from_cube(
                Path("dummy-in.cub"),
                Path("dummy-out.cub"),
                width=5,
                rev_area=True,
            )

            # Need to build new test here.
            # self.assertEqual(
            #     self.cleaned, m_wtable.call_args[0][2]["Calibration"]
            # )

    @patch("hiproc.bitflips.clean_tables_from_cube")
    @patch(
        "hiproc.bitflips.pvl.load",
        return_value={
            "IsisCube": {"Core": {"Pixels": {"Type": "UnsignedWord"}}}
        },
    )
    @patch("hiproc.bitflips.isis.mask")
    def test_clean_cube(self, m_mask, m_load, m_clntbl):
        with patch(
            "hiproc.bitflips.gdal_array.LoadFile", return_value=self.cal
        ):
            bf.clean_cube(
                Path("dummy-in.cub"), Path("dummy-out.cub"), keep=True
            )

        self.assertEqual(m_mask.call_args[1]["minimum"], 990)
        self.assertEqual(m_mask.call_args[1]["maximum"], 1010)
