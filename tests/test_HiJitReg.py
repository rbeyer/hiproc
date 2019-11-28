#!/usr/bin/env python
"""This module has tests for the HiRISE HiJitReg functions."""

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

# import contextlib
import unittest
from unittest.mock import call, mock_open, patch
from pathlib import Path

import pvl

import PyRISE.HiColorInit as hci
import PyRISE.HiJitReg as hjr


reg_def = {'AutoRegistration': {'PatternChip': {'Samples': 100, 'Lines': 60},
                                'SearchChip': {'Samples': 120, 'Lines': 100}}}

flat_tab_text = '''#  **** Registration Data ****
#   RegFile: tmp.regdef
#   OverlapSize:         1024    4000
#   Sample Spacing:   1024.0
#   Line Spacing:     100.0
#   Columns, Rows:    1 40
#   Corr. Algorithm:  MaximumCorrelation
#   Corr. Tolerance:  0.50
#   Total Registers:  38 of 40
#   Number Suspect:   0
#   Average Sample Offset: 0.4470  StdDev: 0.8409
#   Average Line Offset:   3.1888 StdDev: 0.6859

#  Column Headers and Data
            FromTime  FromSamp  FromLine           MatchTime MatchSamp MatchLine        RegSamp        RegLine   RegCorr      B0_Offset       B1_Slope   B_RCorr
  277977284.04784042       512       250  277977283.90413469       512       250       512.0000       247.0109  0.870027      -0.008024       0.887162  0.870215
  277977284.06594044       512       350  277977283.92223471       512       350       512.3985       347.7992  0.864923       0.005009       0.813267  0.882290
  277977284.08404046       512       450  277977283.94033474       512       450       512.0000       446.4992  0.892799       0.001585       0.832279  0.897579
  277977284.10214043       512       550  277977283.95843470       512       550       511.3900       546.8042  0.865275       0.010708       0.780005  0.867080
  277977284.12024045       512       650  277977283.97653472       512       650       511.4987       647.9978  0.892329       0.004874       0.812275  0.900263'''

cnet_pvl = pvl.loads('''Object = ControlNetwork
  NetworkId    = Null
  TargetName   = Mars
  UserName     = rbeyer
  Created      = 2019-08-19T15:15:59
  LastModified = 2019-08-19T15:15:59
  Description  = "Records s/c jitter between two adjacent HiRISE images"
  Version      = 5

  Object = ControlPoint
    PointType   = Free
    PointId     = "Row 0 Column 0"
    ChooserName = hijitreg
    DateTime    = 2019-08-19T15:16:02
    Ignore      = True

    Group = ControlMeasure
      SerialNumber = MRO/HIRISE/909172438:57141/IR10/2
      MeasureType  = Candidate
      ChooserName  = hijitreg
      DateTime     = 2019-08-19T15:16:02
      Sample       = 512.0
      Line         = 50.0
    End_Group

    Group = ControlMeasure
      SerialNumber = MRO/HIRISE/909172438:57141/RED4/2
      MeasureType  = Candidate
      ChooserName  = hijitreg
      DateTime     = 2019-08-19T15:16:02
      Sample       = 512.0
      Line         = 50.0
      Reference    = True
    End_Group
  End_Object

  Object = ControlPoint
    PointType   = Free
    PointId     = "Row 2 Column 0"
    ChooserName = hijitreg
    DateTime    = 2019-08-19T15:16:02

    Group = ControlMeasure
      SerialNumber   = MRO/HIRISE/909172438:57141/IR10/2
      MeasureType    = RegisteredPixel
      ChooserName    = hijitreg
      DateTime       = 2019-08-19T15:16:02
      Sample         = 512.0
      Line           = 247.01086232571
      SampleResidual = 0.0 <pixels>
      LineResidual   = 2.9891376742897 <pixels>
      GoodnessOfFit  = 0.87002737689241
      Reference      = True
    End_Group

    Group = ControlMeasure
      SerialNumber = MRO/HIRISE/909172438:57141/RED4/2
      MeasureType  = Candidate
      ChooserName  = hijitreg
      DateTime     = 2019-08-19T15:16:02
      Sample       = 512.0
      Line         = 250.0
      Reference    = True
    End_Group
  End_Object
End_Object
End
''')

jitter_params = {'COLS': 4, 'ROWS': 200, 'TOLERANCE': 0.5,
                 'PATTERN_SAMPLES': 100, 'PATTERN_LINES': 60,
                 'SEARCH_SAMPLES': 120, 'SEARCH_LINES': 100,
                 'SEARCHLONGER_SAMPLES': 130,
                 'SEARCHLONGER_LINES': 140,
                 'GROUP': 'HiJitReg'}


def getkey(cube, group, key):
            values = {'ProductId': None,
                      'Summing': 2,
                      'Lines': 1024,
                      'Samples': 1024,
                      'TDI': 64,
                      'StitchedProductIds': ['PSP_010502_2090_RED4_0',
                                             'PSP_010502_2090_RED4_1']}
            return values[key]


class TestJitterCube(unittest.TestCase):

    @patch('PyRISE.HiColorInit.isis.getkey_k', side_effect=getkey)
    def test_init(self, mock_getkey):
        c = hci.HiColorCube('dummy/PSP_010502_2090_IR10_0')
        j = hjr.JitterCube(c)
        self.assertEqual(j.cnet_path,
                         Path('dummy/PSP_010502_2090_RED4-IR10.control.pvl'))

    @patch('PyRISE.HiColorInit.isis.getkey_k', side_effect=getkey)
    def test_badinit(self, mock_getkey):
        c = hci.HiColorCube('dummy/PSP_010502_2090_RED4_0')
        self.assertRaises(KeyError, hjr.JitterCube, c)

    @patch('PyRISE.HiColorInit.isis.getkey_k', side_effect=getkey)
    def test_dict(self, mock_getkey):
        c = hjr.JitterCube('dummy/PSP_010502_2090_IR10_0')
        self.assertEqual(c['bin'], c.bin)
        self.assertEqual(c['ExcludeLimit'], 2)

    @patch('PyRISE.HiColorInit.isis.getkey_k', side_effect=getkey)
    def test_parseRegDefs(self, mock_getkey):
        c = hjr.JitterCube('dummy/PSP_010502_2090_IR10_0')
        with patch('PyRISE.HiJitReg.pvl.load', return_value=reg_def):
            c.parseRegDefs()
            self.assertEqual(c['PatternSamples'], 100)

    @patch('PyRISE.HiColorInit.isis.getkey_k', side_effect=getkey)
    def test_parseFlatTab(self, mock_getkey):
        c = hjr.JitterCube('dummy/PSP_010502_2090_IR10_0')
        with patch('PyRISE.HiJitReg.pvl.load', return_value=reg_def):
            c.parseRegDefs()
        with patch('PyRISE.HiJitReg.open', mock_open(read_data=flat_tab_text)):
            c.parseFlatTab()
        self.assertEqual(c['MatchedLineCount'], 5)
        self.assertEqual(c['AvgSampleOffset'], 0.4470)

    @patch('PyRISE.HiColorInit.isis.getkey_k', side_effect=getkey)
    def test_get_control_measures(self, m_getkey):
        self.assertEqual(len(hjr.JitterCube._get_control_measures(cnet_pvl)), 1)

    @patch('PyRISE.HiColorInit.isis.getkey_k', side_effect=getkey)
    def test_parseCNetPVL(self, mock_getkey):
        c = hjr.JitterCube('dummy/PSP_010502_2090_IR10_0')
        with patch('PyRISE.HiJitReg.pvl.load', return_value=reg_def):
            c.parseRegDefs()
        with patch('PyRISE.HiJitReg.open', mock_open(read_data=flat_tab_text)):
            c.parseFlatTab()
        with patch('PyRISE.HiJitReg.pvl.load', return_value=cnet_pvl):
            c.parseCNetPVL()
            self.assertFalse(c['CanSlither'])
            self.assertEqual(c['MatchedLineCount'], 1)

    @unittest.skip("Tests on a real file, time consuming.")
    def test_parseCNetPVL_file(self):
        c = hjr.JitterCube('tmp/PSP_010502_2090_IR10.HiStitch.balance.precolor.cub')
        c.parseRegDefs()
        c.parseFlatTab()
        c.parseCNetPVL('tmp/test.pvl')

    @patch('PyRISE.HiColorInit.isis.getkey_k', side_effect=getkey)
    def test_filterCNetPVL(self, mock_getkey):
        c = hjr.JitterCube('dummy/PSP_010502_2090_IR10_0')
        with patch('PyRISE.HiJitReg.pvl.load', return_value=reg_def):
            c.parseRegDefs()
        with patch('PyRISE.HiJitReg.open', mock_open(read_data=flat_tab_text)):
            c.parseFlatTab()
        c.IgnoredPoints.add('Row 2 Column 0')
        with patch('PyRISE.HiJitReg.pvl.load', return_value=cnet_pvl):
            with patch('PyRISE.HiJitReg.open', mock_open()) as m:
                c.filterCNetPVL()
                ignorecall = call().write(b'Ignore')
                self.assertEqual(len(tuple(filter(lambda x: x == ignorecall,
                                                  m.mock_calls))),
                                 2)


class TestHiJitReg(unittest.TestCase):

    @patch('PyRISE.HiColorInit.isis.getkey_k', side_effect=getkey)
    def setUp(self, mock_getkey):
        self.r = hci.HiColorCube('dummy/PSP_010502_2090_RED4_0')
        self.c = hci.HiColorCube('dummy/PSP_010502_2090_IR10_0')
        self.j = hjr.JitterCube(self.c)

        self.conf = {'AutoRegistration': {'ControlNet': {'Control_Cols': 4,
                                                         'Control_Lines': 200},
                                          'Algorithm': {'Tolerance': 0.5,
                                                        'Increment': 0.1},
                                          'PatternChip': {'Samples': 100,
                                                          'Lines': 60},
                                          'SearchChip': {'Samples': 120,
                                                         'Lines': 100},
                                          'SearchLongerChip': {'Samples': 130,
                                                               'Lines': 140}},
                     'Smoothing': {'Exclude_Limit': 2,
                                   'Badness_Limit': 1,
                                   'Boxcar_Length': 10}}

    def test_Analyze_Flat(self):
        with patch('PyRISE.HiJitReg.open',
                   mock_open(read_data=flat_tab_text)):
            with patch('PyRISE.HiJitReg.pvl.load',
                       side_effect=[reg_def, cnet_pvl, cnet_pvl]):
                self.assertEqual(hjr.Analyze_Flat(self.j, 1, 1), 0)

            with patch('PyRISE.HiJitReg.pvl.load',
                       side_effect=[reg_def, cnet_pvl, cnet_pvl]):
                self.j['CanSlither'] = True
                self.assertEqual(hjr.Analyze_Flat(self.j, 1, 1), 1)

    @patch('PyRISE.HiJitReg.Path')
    def test_write_regdef(self, m_Path):
        hjr.write_regdef('dummy/out.pvl', jitter_params)
        self.assertEqual(len(m_Path.mock_calls), 2)

    @patch('PyRISE.HiJitReg.isis.cnetbin2pvl')
    @patch('PyRISE.HiJitReg.isis.hijitreg')
    @patch('PyRISE.HiJitReg.Path')
    def test_run_HiJitReg(self, m_Path, m_hijitreg, m_cnetbin2pvl):
        hjr.run_HiJitReg(self.r.path, self.j, jitter_params, 'foo', keep=True)
        m_hijitreg.assert_called_once()
        m_cnetbin2pvl.assert_called_once()

    @patch('PyRISE.HiJitReg.isis.cnetbin2pvl')
    @patch('PyRISE.HiJitReg.isis.hijitreg')
    @patch('PyRISE.HiJitReg.Path.write_text')
    @patch('PyRISE.HiColorInit.isis.getkey_k', side_effect=getkey)
    def test_jitter_iter(self, m_getkey, m_Path, m_hijitreg, m_cnetbin2pvl):
        with patch('PyRISE.HiJitReg.open',
                   mock_open(read_data=flat_tab_text)):
            with patch('PyRISE.HiJitReg.pvl.load',
                       side_effect=[reg_def, cnet_pvl, reg_def, cnet_pvl]):
                self.assertFalse(hjr.jitter_iter(self.r, self.c, self.conf,
                                                 keep=True))

            with patch('PyRISE.HiJitReg.pvl.load',
                       side_effect=[reg_def, cnet_pvl, reg_def, cnet_pvl]):
                with patch('PyRISE.HiJitReg.Analyze_Flat', return_value=1):
                    self.assertTrue(hjr.jitter_iter(self.r, self.c, self.conf,
                                                    keep=True))

    @patch('PyRISE.HiJitReg.jitter_iter', return_value=True)
    def test_HiJitReg(self, m_jit_it):
        c04 = c05 = c10 = c11 = c12 = c13 = None
        conf = dict(dummy='yes')
        self.assertRaises(RuntimeError, hjr.HiJitReg, c04, c05, c10, c11, c12, c13, conf)
        c04 = 'RED4'
        c12 = 'BG12'
        c13 = 'BG13'
        self.assertListEqual(['BG12'], hjr.HiJitReg(c04, c05, c10, c11, c12, c13, conf))
        c05 = 'RED5'
        c10 = 'IR10'
        c11 = 'IR11'
        self.assertListEqual(['IR10', 'BG12', 'IR11', 'BG13'],
                             hjr.HiJitReg(c04, c05, c10, c11, c12, c13, conf))
