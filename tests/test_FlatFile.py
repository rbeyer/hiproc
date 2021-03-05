#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `FlatFile` class."""

# Copyright 2020, Ross A. Beyer (rbeyer@seti.org)
#
# Reuse is permitted under the terms of the license.
# The AUTHORS file and the LICENSE file are at the
# top level of this library.

import unittest

import hiproc.FlatFile as ff


class TestFlatFile(unittest.TestCase):
    def setUp(self):
        self.flat = ff.FlatFile(
            """#          Hijitreg ISIS Application Results
#    Coordinates are (Sample, Line) unless indicated
#           RunDate:  2012-04-11T20:02:49
#
#    ****  Image Input Information ****
#  FROM:  /HiRISE/Users/audrie/pipe_dev_sandbox/HiRISE/Data/ResolveJitter/ESP/ORB_016200_016299/ESP_016213_2315/ESP_016213_2315_RED4.prehijack.cub
#    Lines:       70000
#    Samples:     2048
#    FPSamp0:     0
#    SampOffset:  0
#    LineOffset:  0
#    CPMMNumber:  5
#    Summing:     1
#    TdiMode:     128
#    Channel:     2
#    LineRate:    0.00009700 <seconds>
#    TopLeft:           1       1
#    LowerRight:     2048   70000
#    StartTime:   2010-01-10T20:07:22.068 <UTC>
#    SCStartTime: 947621272:49930 <SCLK>
#    StartTime:   316426108.23449421 <seconds>

#  MATCH: /HiRISE/Users/audrie/pipe_dev_sandbox/HiRISE/Data/ResolveJitter/ESP/ORB_016200_016299/ESP_016213_2315/ESP_016213_2315_BG12.prehijack.cub
#    Lines:       70000
#    Samples:     2048
#    FPSamp0:     0
#    SampOffset:  0
#    LineOffset:  0
#    CPMMNumber:  4
#    Summing:     1
#    TdiMode:     128
#    Channel:     2
#    LineRate:    0.00009700 <seconds>
#    TopLeft:           1       1
#    LowerRight:     2048   70000
#    StartTime:   2010-01-10T20:07:21.972 <UTC>
#    SCStartTime: 947621272:43630 <SCLK>
#    StartTime:   316426108.13836384 <seconds>


#  **** Registration Data ****
#   RegFile: /HiRISE/Users/audrie/pipe_dev_sandbox/HiRISE/Data/ResolveJitter/ESP/ORB_016200_016299/ESP_016213_2315/ESP_016213_2315_BG12-RED4.regdef.pvl
#   OverlapSize:         2048   70000
#   Sample Spacing:   682.0
#   Line Spacing:     20.0
#   Columns, Rows:    3 3500
#   Corr. Algorithm:  MaximumCorrelation
#   Corr. Tolerance:  0.70
#   Total Registers:  10145 of 10500
#   Number Suspect:   0
#   Average Sample Offset: 0.3586  StdDev: 0.3605
#   Average Line Offset:   3.9526 StdDev: 1.5403

#  Column Headers and Data
            FromTime  FromSamp  FromLine           MatchTime MatchSamp MatchLine        RegSamp        RegLine   RegCorr      B0_Offset       B1_Slope   B_RCorr
  316426108.25476718       341       210  316426108.15863681       341       210       340.6611       204.2240  0.883177       0.024032       0.273207  0.894237
  316426108.25476718      1023       210  316426108.15863681      1023       210      1022.5484       204.4079  0.805408       0.020917       0.286734  0.839953
  316426108.25476718      1705       210  316426108.15863681      1705       210      1704.6887       204.6171  0.823995       0.026994       0.258996  0.853913
  316426108.25670719       341       230  316426108.16057682       341       230       340.7258       224.1064  0.868361       0.025117       0.266694  0.879799
  316426108.25670719      1023       230  316426108.16057682      1023       230      1022.5470       224.5099  0.807298       0.023772       0.269487  0.842616
  316426108.25670719      1705       230  316426108.16057682      1705       230      1704.6864       224.5820  0.819560       0.028520       0.250281  0.848002
  316426108.25864720       341       250  316426108.16251683       341       250       340.6485       244.1028  0.848406       0.026669       0.257074  0.869427
  316426108.25864720      1023       250  316426108.16251683      1023       250      1022.4197       244.4699  0.798717       0.024407       0.265783  0.836001
  316426108.25864720      1705       250  316426108.16251683      1705       250      1704.5953       244.6173  0.865055       0.022810       0.282049  0.885654
  316426108.26058722       341       270  316426108.16445684       341       270       340.6213       264.1180  0.798498       0.029366       0.241046  0.828585
  316426108.26058722      1023       270  316426108.16445684      1023       270      1022.3535       264.4123  0.790646       0.025977       0.257441  0.823464
  316426108.26058722      1705       270  316426108.16445684      1705       270      1704.6039       264.4040  0.860945       0.024428       0.272569  0.882096
  316426108.26252723       341       290  316426108.16639686       341       290       340.5806       284.2182  0.754269       0.032774       0.221183  0.794989
  316426108.26252723      1023       290  316426108.16639686      1023       290      1022.3008       284.0952  0.763710       0.028660       0.242425  0.791072
  316426108.26252723      1705       290  316426108.16639686      1705       290      1704.6190       284.3197  0.862750       0.024277       0.272852  0.885356
  316426108.26446718       341       310  316426108.16833681       341       310       340.6429       304.1541  0.734999       0.034928       0.209090  0.767223
  316426108.26446718      1705       310  316426108.16833681      1705       310      1704.6546       304.1873  0.817908       0.030489       0.237423  0.849629
  316426108.26640719       341       330  316426108.17027682       341       330       340.6055       324.1532  0.713359       0.039534       0.183031  0.748323
  316426108.26640719      1023       330  316426108.17027682      1023       330      1022.5223       324.4402  0.696912       0.034920       0.206285  0.745516
  316426108.26640719      1705       330  316426108.17027682      1705       330      1704.6477       324.1382  0.781064       0.034834       0.212454  0.820036
  316426108.26834720       341       350  316426108.17221683       341       350       340.6203       344.0387  0.681599       0.041267       0.173113  0.716820
  316426108.26834720      1023       350  316426108.17221683      1023       350      1022.4567       344.2591  0.694442       0.034860       0.206966  0.742064
  316426108.26834720      1705       350  316426108.17221683      1705       350      1704.6779       344.1707  0.754361       0.035685       0.207219  0.795146
  316426108.27028722      1023       370  316426108.17415684      1023       370      1022.5797       364.2458  0.693458       0.035562       0.202958  0.743528
  316426108.27028722      1705       370  316426108.17415684      1705       370      1704.7215       364.0940  0.783275       0.033000       0.221816  0.817121
  316426108.27222723       341       390  316426108.17609686       341       390       340.6190       383.9682  0.688795       0.039935       0.180578  0.735695
"""
        )

    def test_init_str(self):
        self.assertIsInstance(self.flat, ff.FlatFile)

    def test_dictlike(self):
        self.assertEqual(self.flat["RunDate"], "2012-04-11T20:02:49")
        self.assertEqual(self.flat["FROM"]["Lines"], "70000")
        # print(self.flat)

    def test_listlike(self):
        self.assertEqual(12, len(self.flat[0]))

    def test_contains(self):
        self.assertTrue("OverlapSize" in self.flat)

    def test_len(self):
        self.assertEqual(26, len(self.flat))
