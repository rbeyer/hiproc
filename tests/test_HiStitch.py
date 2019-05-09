#!/usr/bin/env python
"""This module has tests for the HiRISE HiStitch functions."""

# Copyright 2019, Ross A. Beyer (rbeyer@seti.org)
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
import unittest
from pathlib import Path
from unittest.mock import call
from unittest.mock import mock_open
from unittest.mock import patch

import pvl

import PyRISE.hirise as hirise
import PyRISE.HiStitch as hs

conf_path = Path('resources') / 'HiStitch.conf'


class TestBasic(unittest.TestCase):

    def test_get_chids(self):
        cid_list = ('PSP_010502_2090_RED4_0',
                    'PSP_010502_2090_RED4_1',
                    'PSP_010502_2090_RED5_0')
        cube_list = list()
        for p in cid_list:
            cube_list.append(p + '.EDR_Stats.HiCal.cub')
        self.assertEquals(cid_list[0:2], tuple(map(str,
                                               hs.get_chids(cube_list[0:2]))))
        self.assertEquals(cid_list[-1], str(hs.get_chids((cube_list[-1],))[0]))
        self.assertRaises(ValueError, hs.get_chids, cube_list[1:3])
        with contextlib.suppress(FileNotFoundError):
            Path('print.prt').unlink()

    def test_set_outpath(self):
        d = Path('dummy-directory')
        output = Path('foo')
        chid = hirise.CCDID('PSP_010502_2090_RED5_0')
        self.assertEquals(output, hs.set_outpath(output, chid, d))
        suffix = '.foo'
        with_suffix = d / Path('PSP_010502_2090_RED5.foo')
        self.assertEquals(with_suffix, hs.set_outpath(suffix, chid, d))


class TestConf(unittest.TestCase):

    def test_conf_check(self):
        c = pvl.load(str(conf_path))
        self.assertIsNone(hs.conf_check(c))


class TestMock(unittest.TestCase):

    def test_sort_input_cubes(self):
        cub0 = 'PSP_010502_2090_RED4_0'
        cub1 = 'PSP_010502_2090_RED4_1'
        with patch('PyRISE.HiStitch.isis.getkey_k', side_effect=[cub0[-1],
                                                                 cub1[-1]]):
            self.assertEquals((cub0, cub1), hs.sort_input_cubes(cub0, cub1))

        with patch('PyRISE.HiStitch.isis.getkey_k', side_effect=[cub1[-1],
                                                                 cub0[-1]]):
            self.assertEquals((cub0, cub1), hs.sort_input_cubes(cub1, cub0))

        with patch('PyRISE.HiStitch.isis.getkey_k', return_value=0):
            self.assertRaises(RuntimeError, hs.sort_input_cubes, cub0, cub0)

        with patch('PyRISE.HiStitch.isis.getkey_k', side_effect=[0, 2]):
            self.assertRaises(RuntimeError, hs.sort_input_cubes, cub0, cub1)

    def test_sort_databases(self):
        chids = (hirise.ChannelID('PSP_010502_2090_RED4_0'),
                 hirise.ChannelID('PSP_010502_2090_RED4_1'))
        dbs = ({'PRODUCT_ID': str(chids[0])}, {'PRODUCT_ID': str(chids[1])})
        _ = 'dummy_path'
        with patch('PyRISE.HiStitch.open', mock_open(read_data='dummy')):
            with patch('PyRISE.HiStitch.json.load', side_effect=dbs):
                self.assertRaises(IndexError, hs.sort_databases, (_, _, _), chids)
            with patch('PyRISE.HiStitch.json.load', side_effect=dbs):
                self.assertRaises(LookupError, hs.sort_databases,
                                  (_, _), (chids[0], chids[0]))
            with patch('PyRISE.HiStitch.json.load', side_effect=dbs):
                self.assertEquals(dbs, hs.sort_databases((_, _), chids))
            with patch('PyRISE.HiStitch.json.load', side_effect=reversed(dbs)):
                self.assertEquals(dbs, hs.sort_databases((_, _), chids))

    @patch('PyRISE.HiStitch.isis.fx')
    @patch('PyRISE.HiStitch.isis.mask')
    @patch('PyRISE.HiStitch.isis.lowpass')
    @patch('PyRISE.HiStitch.isis.highpass')
    @patch('PyRISE.HiStitch.isis.algebra')
    @patch('PyRISE.HiStitch.isis.handmos')
    @patch('shutil.copyfile')
    def test_HiFurrow_Fix(self, mock_copyfile, mock_handmos, mock_algebra,
                          mock_highpass, mock_lowpass, mock_mask, mock_fx):
        with patch('PyRISE.HiStitch.isis.getkey_k', return_value=1):
            self.assertRaises(ValueError, hs.HiFurrow_Fix, 'dum_in', 'dum_out', 0)
        with patch('PyRISE.HiStitch.isis.getkey_k', side_effect=[2, 1, 1]):
            self.assertRaises(ValueError, hs.HiFurrow_Fix, 'dum_in', 'dum_out', 0)
        with patch('PyRISE.HiStitch.isis.getkey_k', side_effect=[2, 1024, 1024]):
            in_cube = 'dummy_in.cub'
            out_cube = 'dummy_out.cub'
            hs.HiFurrow_Fix(in_cube, out_cube, 1000, keep=True)

            mock_fx.assert_called_once()
            eqn = '(1*(sample<512)+ 1*(sample>513) + 0)'
            fx_path = mock_fx.call_args[1]['to']
            mock_fx.assert_called_once_with(equation=eqn,
                                            lines=1024, mode='OUTPUTONLY',
                                            samples=1024,
                                            to=fx_path)

            mask1_path = mock_mask.call_args_list[0][1]['to']
            mask2_path = mock_mask.call_args_list[1][1]['to']
            mock_calls = [call(Path(in_cube),
                               mask=fx_path,
                               max_=1, min_=1, preserve='INSIDE', spixels='NULL',
                               to=mask1_path),
                          call(Path(in_cube),
                               mask=fx_path,
                               max_=0, min_=0, preserve='INSIDE', spixels='NULL',
                               to=mask2_path)]
            mock_mask.assert_has_calls(mock_calls)

            lpf_path = mock_lowpass.call_args[1]['to']
            mock_lowpass.assert_called_once_with(mask1_path, his=False, hrs=False,
                                                 line=41, lis=False, null=True,
                                                 sample=5, to=lpf_path)

            hpf_path = mock_highpass.call_args[1]['to']
            mock_highpass.assert_called_once_with(mask2_path, line=41, sample=1,
                                                  to=hpf_path)

            alg_path = mock_algebra.call_args[1]['to']
            mock_algebra.assert_called_once_with(A=1.0, B=1.0, from2=hpf_path,
                                                 from_=lpf_path, operator='ADD',
                                                 to=alg_path)

            mock_handmos.assert_called_once_with(alg_path, create='NO', inband=1,
                                                 inline=1, insample=1,
                                                 mosaic=out_cube, outband=1,
                                                 outline=1, outsample=1)


class TestHiStitch(unittest.TestCase):

    def setUp(self):
        self.my_c = {'HiStitch_Furrow_Normalization': 4000,
                     'HiStitch_Balance': True,
                     'HiStitch_Equalize': False,
                     'HiStitch_Balance_Processing': (1, 1, 0, 1, 0, 0, 0,
                                                     0, 0, 1, 1, 0, 0, 0),
                     'HiStitch_LIS_Pixels': 10000,
                     'HiStitch_Balance_Bin_Mask_STD': (200.0, 200.0, 100.0,
                                                       100.0, 100.0),
                     'HiStitch_Balance_Bin_DarkPixel_STD': (70.0, 70.0, 70.0,
                                                            70.0, 70.0),
                     'HiStitch_Gap_Percent': 10.0,
                     'HiStitch_Bin01_Area': (7, 11),
                     'HiStitch_Bin02_Area': (7, 7),
                     'HiStitch_Bin04_Area': (7, 5),
                     'HiStitch_Control_Channel': (1, 0, 0, 0, 0, 0, 0,
                                                  0, 0, 0, 0, 0, 0, 0),
                     'HiStitch_Equalize_Width': 1501,
                     'HiStitch_Equalize_Correction': 'MULTIPLY'}
        self.my_dbs = list()
        self.my_dbs.append({'IMAGE_MEAN': 6491.34508964,
                            'LOW_SATURATED_PIXELS': 0,
                            'CAL_MASK_STANDARD_DEVIATION': 17.97,
                            'IMAGE_DARK_STANDARD_DEVIATION': 18.05,
                            'GAP_PIXELS_PERCENT': 0,
                            'BINNING': 2,
                            'zapped': False})
        self.my_dbs.append({'IMAGE_MEAN': 9.34508964,
                            'LOW_SATURATED_PIXELS': 0,
                            'CAL_MASK_STANDARD_DEVIATION': 17.97,
                            'IMAGE_DARK_STANDARD_DEVIATION': 18.05,
                            'GAP_PIXELS_PERCENT': 0,
                            'BINNING': 2,
                            'zapped': False})

    def test_set_flags(self):
        b = 1, 2, 4, 8, 16
        self.assertEquals((True, True, False), hs.set_flags(self.my_c,
                                                            self.my_dbs,
                                                            5, b.index(2),
                                                            6491.35))
        self.my_c['HiStitch_Equalize'] = True
        self.assertEquals((True, True, True), hs.set_flags(self.my_c,
                                                           self.my_dbs,
                                                           5, b.index(2),
                                                           6491.35))
        self.my_dbs[0]['zapped'] = True
        self.assertEquals((False, True, True), hs.set_flags(self.my_c,
                                                            self.my_dbs,
                                                            5, b.index(2),
                                                            6491.35))

    def test_HiStitchStep(self):
        cubes = ('dummy1.in', 'dummy2.in')
        out_c = 'dummy.out'

        with patch('PyRISE.HiStitch.isis.histitch') as mock:
            hs.HiStitchStep(cubes, out_c, 2, 5, self.my_c, False, False)
            mock.assert_called_with(from1=cubes[0], from2=cubes[1],
                                    to=out_c, balance=False)

            hs.HiStitchStep([cubes[0], ], out_c, 2, 5, self.my_c, False, False)
            mock.assert_called_with(from1=cubes[0],
                                    to=out_c, balance=False)

            hs.HiStitchStep(cubes, out_c, 2, 5, self.my_c, True, False)
            mock.assert_called_with(from1=cubes[0], from2=cubes[1],
                                    to=out_c, balance='TRUE',
                                    channel=0, seamsize=7, skip=7)

            hs.HiStitchStep(cubes, out_c, 4, 0, self.my_c, False, True)
            mock.assert_called_with(from1=cubes[0], from2=cubes[1],
                                    to=out_c, balance='EQUALIZE',
                                    channel=1, seamsize=5, skip=7,
                                    operator='MULTIPLY', width=1501)

    @patch('PyRISE.HiStitch.Path.rename')
    @patch('PyRISE.HiStitch.isis.fx')
    @patch('PyRISE.HiStitch.isis.mask')
    @patch('PyRISE.HiStitch.isis.lowpass')
    @patch('PyRISE.HiStitch.isis.highpass')
    @patch('PyRISE.HiStitch.isis.algebra')
    @patch('PyRISE.HiStitch.isis.handmos')
    @patch('PyRISE.HiStitch.isis.histitch')
    @patch('shutil.copyfile')
    def test_HiStitch(self, mock_histitch, mock_copyfile, mock_handmos,
                      mock_algebra, mock_highpass, mock_lowpass,
                      mock_mask, mock_fx, mock_rename):
        cubes = ('dummy1.in', 'dummy2.in')
        out_c = 'dummy.out'

        def getkey(cube, group, key):
            values = {'TruthChannel': 0,  # need to look at actual output
                      'BalanceRatio': 3,
                      'Summing': 2,
                      'Lines': 1024,
                      'Samples': 1024}
            return values[key]

        with patch('PyRISE.HiStitch.isis.getkey_k', side_effect=getkey):
            self.assertEquals((0, 3), hs.HiStitch(cubes, out_c,
                                                  self.my_c, self.my_dbs,
                                                  5, keep=True))

        with patch('PyRISE.HiStitch.isis.getkey_k', side_effect=getkey):
            self.assertEquals((None, None), hs.HiStitch([cubes, ], out_c,
                                                        self.my_c, self.my_dbs,
                                                        5, keep=True))
