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
from unittest.mock import call, patch

import numpy as np
from scipy.signal import find_peaks

import kalasiris as isis

from pyrise import bitflips as bf
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


class TestBasic(unittest.TestCase):

    def test_get_range(self):
        start = 2
        stop = 5
        truth = range(start, stop, 1)
        self.assertEqual(truth, bf.get_range(start, stop))

        start = 2.1
        stop = 4.8
        self.assertEqual(truth, bf.get_range(start, stop))

        start = 9.1
        truth = range(10, 4, -1)
        self.assertEqual(truth, bf.get_range(start, stop))


class TestHist(unittest.TestCase):

    def setUp(self):
        HistRow = collections.namedtuple('HistRow', ['DN', 'Pixels'])
        self.hist = [HistRow(1, 1),     # 0 *
                     HistRow(2, 2),     # 1 **
                     HistRow(5, 3),     # 2 ***
                     HistRow(6, 4),     # 3 ****
                     HistRow(7, 3),     # 4 ***
                     HistRow(8, 2),     # 5 **
                     HistRow(10, 1),    # 6 *
                     HistRow(11, 2),    # 7 **
                     HistRow(12, 2),    # 8 **
                     # Biggest gap here between 12 and 20
                     HistRow(20, 1)]    # 9 *

    def test_find_gap(self):
        self.assertEqual(13, bf.find_gap(self.hist, 1, 20))
        self.assertEqual(3, bf.find_gap(self.hist, 1, 20, findfirst=True))
        self.assertEqual(9, bf.find_gap(self.hist, 5, 12))

    def test_find_min_dn(self):
        self.assertEqual(1, bf.find_min_dn(self.hist, 1, 20))
        self.assertEqual(10, bf.find_min_dn(self.hist, 2, 20))
        self.assertEqual(20, bf.find_min_dn(self.hist, 20, 1))

    @patch('pyrise.bitflips.shutil.copyfile')
    @patch('pyrise.bitflips.isis.mask')
    @patch('pyrise.bitflips.isis.algebra')
    @patch('pyrise.bitflips.isis.handmos')
    def test_subtract_over_thresh(self, m_hand, m_alg, m_mask, m_copyfile):
        in_p = Path('in.dummy')
        out_p = Path('out.dummy')
        bf.subtract_over_thresh(in_p, out_p, thresh=1, delta=2, keep=True)

        m_copyfile.assert_called_once()

        m_mask.assert_called_once()
        t_path = in_p.with_suffix('.threshmask.cub')
        d_path = in_p.with_suffix('.delta.cub')
        mask_args = {'from': in_p, 'min': 1, 'to': t_path}
        self.assertEqual(m_mask.call_args_list,
                         [call(**mask_args)])

        m_alg.assert_called_once()
        self.assertEqual(m_alg.call_args_list,
                         [call(t_path, a=0, c=-2, from2=in_p, op='add',
                               to=d_path)])

        m_hand.assert_called_once()
        self.assertEqual(m_hand.call_args_list,
                         [call(d_path, mosaic=out_p)])

    def test_find_minima_index(self):
        pixel_counts = np.fromiter((int(x.Pixels) for x in self.hist), int)
        minima_i, _ = find_peaks(np.negative(pixel_counts))

        self.assertEqual(0,
                         bf.find_minima_index(4, 0, minima_i, pixel_counts))
        self.assertEqual(6,
                         bf.find_minima_index(2, 8, minima_i, pixel_counts))

    def test_find_smart_window(self):
        self.assertEqual((1, 10),
                         bf.find_smart_window(self.hist, 1, 12, 6, plot=False))


class TestMask(unittest.TestCase):

    def setUp(self):
        self.h = isis.Histogram('''Cube:           foo.cub
Band:           1
Average:        6490.68
Std Deviation:  166.739
Variance:       27801.9
Median:         1000
Mode:           6489
Skew:           0.0302975
Minimum:        3889
Maximum:        8230

Total Pixels:    2048000
Valid Pixels:    2048000
Null Pixels:     0
Lis Pixels:      0
Lrs Pixels:      0
His Pixels:      0
Hrs Pixels:      0


DN,Pixels,CumulativePixels,Percent,CumulativePercent
3889,1,1,4.88281e-05,4.88281e-05
3924,1,2,4.88281e-05,9.76563e-05
3960,2,4,9.76563e-05,0.000195313
3995,1,5,4.88281e-05,0.000244141
4030,4,9,0.000195313,0.000439453
4065,6,15,0.000292969,0.000732422
4101,12,27,0.000585937,0.00131836
4136,13,40,0.000634766,0.00195312
4171,11,51,0.000537109,0.00249023
4207,11,62,0.000537109,0.00302734
4242,10,72,0.000488281,0.00351562
4277,14,86,0.000683594,0.00419922
4312,4,90,0.000195313,0.00439453''')

    @patch('pyrise.HiCal.analyze_cubenorm_stats2', return_value=(1, 2))
    @patch('pyrise.bitflips.isis.cubenorm')
    def test_mask(self, m_cn, m_acns2):
        in_p = Path('in.dummy')
        out_p = Path('out.dummy')

        with patch('pyrise.bitflips.histogram', return_value=self.h):
            bf.mask(in_p, out_p)
