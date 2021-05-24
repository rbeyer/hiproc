#!/usr/bin/env python
"""This module has tests for the HiRISE EDR_Stats functions."""

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

import contextlib
import pkg_resources
import unittest
from pathlib import Path

import pvl

import kalasiris as isis
import hiproc.EDR_Stats as edr

from .utils import resource_check as rc

# Hardcoding these, but I sure would like a better solution.
HiRISE_img = Path("test-resources") / "PSP_010502_2090_RED5_0.img"
img = HiRISE_img
gains = pkg_resources.resource_filename(
    "hiproc", "data/EDR_Stats_gains_config.pvl"
)


class TestResources(unittest.TestCase):
    """Establishes that the test image exists."""

    def test_resources(self):
        files = (img, gains)
        for f in files:
            with self.subTest(filepath=f):
                (truth, test) = rc(f)
                self.assertEqual(truth, test)


class TestCheckLUT(unittest.TestCase):
    def test_check(self):
        self.assertEqual(edr.check_lut(img), 316)


class TestEDR_Stats(unittest.TestCase):
    def setUp(self):
        self.outfile = img.with_suffix(".TestEDR_Stats.cub")

    def tearDown(self):
        with contextlib.suppress(FileNotFoundError):
            self.outfile.unlink()
            Path("print.prt").unlink()

    def test_EDR_stats(self):
        h = edr.EDR_Stats(img, self.outfile, pvl.load(gains))
        self.assertIsInstance(h, dict)


class TestNeedHiCube(unittest.TestCase):
    def setUp(self):
        self.hicube = img.with_suffix(".TestNeedHiCube.cub")
        isis.hi2isis(img, to=self.hicube)

    def tearDown(self):
        with contextlib.suppress(FileNotFoundError):
            self.hicube.unlink()
            Path("print.prt").unlink()

    def test_parse_histat(self):
        h = edr.parse_histat(isis.histat(self.hicube).stdout)
        self.assertIsInstance(h, dict)

    def test_get_dncnt(self):
        c = edr.get_dncnt(self.hicube)
        self.assertEqual(c, 80)

    def test_calc_snr(self):
        histats = edr.parse_histat(isis.histat(self.hicube).stdout)
        histats["BINNING"] = isis.getkey_k(
            self.hicube, "Instrument", "Summing"
        )
        s = edr.calc_snr(self.hicube, pvl.load(gains), histats)
        self.assertAlmostEqual(s, 291.80442197)

    def test_tdi_bin_check(self):
        histats = edr.parse_histat(isis.histat(self.hicube).stdout)
        self.assertIsNone(edr.tdi_bin_check(self.hicube, histats))
        histats["PRODUCT_ID"] = "bogus id"
        self.assertIsNone(edr.tdi_bin_check(self.hicube, histats))

    def test_lut_check(self):
        histats = edr.parse_histat(isis.histat(self.hicube).stdout)
        self.assertIsNone(edr.lut_check(self.hicube, histats))
