#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `FlatFile` class."""

# Copyright 2020, Ross A. Beyer (rbeyer@seti.org)
#
# Reuse is permitted under the terms of the license.
# The AUTHORS file and the LICENSE file are at the
# top level of this library.

import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
import numpy.testing as npt

from hiproc.FlatFile import FlatFile
import hiproc.resolve_jitter as rj


class TestResolveJitter(unittest.TestCase):

    def setUp(self):
        self.flat_text = """#          Hijitreg ISIS Application Results
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
  316426108.27222723      1705       390  316426108.17609686      1705       390      1704.6824       384.0053  0.800740       0.031919       0.227514  0.832411
  316426108.27416718       341       410  316426108.17803681       341       410       340.5761       403.9824  0.711003       0.039963       0.180518  0.759016
  316426108.27416718      1705       410  316426108.17803681      1705       410      1704.7217       403.9904  0.794964       0.032618       0.223333  0.824529
  316426108.27610719       341       430  316426108.17997682       341       430       340.5468       424.0275  0.714079       0.036748       0.198811  0.757016
  316426108.27610719      1705       430  316426108.17997682      1705       430      1704.6834       424.0754  0.792730       0.031869       0.227960  0.820162
  316426108.27804720       341       450  316426108.18191683       341       450       340.5045       444.0639  0.738056       0.033717       0.216061  0.772310
  316426108.27804720      1705       450  316426108.18191683      1705       450      1704.6754       444.0586  0.808393       0.026047       0.261489  0.831419
"""

    def test_upper_power_of_two(self):
        self.assertEqual(rj.upper_power_of_two(5), 8)

    def test_set_file_path(self):
        loc = Path("foo")
        path1 = loc / "bar"
        path2 = Path("bah")
        self.assertEqual(rj.set_file_path(loc, path1), path1)
        self.assertEqual(rj.set_file_path(loc, path2), loc / path2)

    def test_create_matrices(self):
        flat_no_dupes = FlatFile("""#          Hijitreg ISIS Application Results
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
  316426108.25670719       341       230  316426108.16057682       341       230       340.7258       224.1064  0.868361       0.025117       0.266694  0.879799
  316426108.25864720       341       250  316426108.16251683       341       250       340.6485       244.1028  0.848406       0.026669       0.257074  0.869427
  316426108.26058722       341       270  316426108.16445684       341       270       340.6213       264.1180  0.798498       0.029366       0.241046  0.828585
  316426108.26252723       341       290  316426108.16639686       341       290       340.5806       284.2182  0.754269       0.032774       0.221183  0.794989
  316426108.26446718       341       310  316426108.16833681       341       310       340.6429       304.1541  0.734999       0.034928       0.209090  0.767223
  316426108.26640719       341       330  316426108.17027682       341       330       340.6055       324.1532  0.713359       0.039534       0.183031  0.748323
  316426108.26834720       341       350  316426108.17221683       341       350       340.6203       344.0387  0.681599       0.041267       0.173113  0.716820
  316426108.27222723       341       390  316426108.17609686       341       390       340.6190       383.9682  0.688795       0.039935       0.180578  0.735695
  316426108.27416718       341       410  316426108.17803681       341       410       340.5761       403.9824  0.711003       0.039963       0.180518  0.759016
  316426108.27610719       341       430  316426108.17997682       341       430       340.5468       424.0275  0.714079       0.036748       0.198811  0.757016
  316426108.27804720       341       450  316426108.18191683       341       450       340.5045       444.0639  0.738056       0.033717       0.216061  0.772310
""")
        nfft = rj.upper_power_of_two(12)
        dt = -0.09613037
        time_list = list()
        offx_list = list()
        offy_list = list()
        for row in flat_no_dupes:
            time_list.append(float(row["FromTime"]))
            offx_list.append(float(row["RegSamp"]) - float(row["FromSamp"]))
            offy_list.append(float(row["RegLine"]) - float(row["FromLine"]))
        tt = np.linspace(0, nfft - 1, nfft) / nfft

        t_arr = np.array(time_list)
        t0 = t_arr[0]
        nfftime = np.linspace(t0, t_arr[-1], nfft)
        duration = t_arr[-1] - t0

        xinterp, yinterp, x, y, ddt, overxx, overyy = rj.create_matrices(
            t_arr, np.array(offx_list), np.array(offy_list),
            dt, duration, t0, nfft, nfftime, tt
        )
        npt.assert_almost_equal(
            np.array([
                -0.3389, -0.2765867, -0.31849644, -0.3637385, -0.38710571,
                -0.4194, -0.36357846, -0.38133553, -0.38929053, -0.37971056,
                -0.37999117, -0.38073106, -0.40392277, -0.43590395,
                -0.46046941, -0.4955
            ]),
            xinterp
        )
        npt.assert_almost_equal(
            np.array([
                -5.776, -5.88606176, -5.89660349, -5.89438368, -5.86820246,
                -5.7818, -5.83900716, -5.84639643, -5.88186211, -5.97255114,
                -6.01086814, -6.03085679, -6.02570878, -6.00248216,
                -5.96458296, -5.9361
            ]),
            yinterp
        )
        npt.assert_almost_equal(
            np.array([
                -7.71832599e-01 + 0.j, 1.91855307e-03 - 0.030364j,
                1.37450069e-02 - 0.04498451j, 8.84287023e-03 - 0.03538853j,
                4.14559077e-04 - 0.01371299j, 1.00119877e-02 - 0.0115939j,
                1.96447927e-03 - 0.01359445j, 4.42185446e-03 + 0.00183916j,
                1.13939776e-02 + 0.j, 4.42185446e-03 - 0.00183916j,
                1.96447927e-03 + 0.01359445j, 1.00119877e-02 + 0.0115939j,
                4.14559077e-04 + 0.01371299j, 8.84287023e-03 + 0.03538853j,
                1.37450069e-02 + 0.04498451j, 1.91855307e-03 + 0.030364j
            ]),
            x
        )
        npt.assert_almost_equal(
            np.array([
                -1.18266834e+01 + 0.j, 7.83362427e-03 - 0.09056171j,
                3.55530799e-02 + 0.03217207j, 2.38036628e-02 - 0.00477377j,
                1.99110490e-02 - 0.00810523j, 4.66141255e-03 - 0.00175215j,
                2.34592026e-02 + 0.00620169j, 1.66323568e-02 - 0.00878693j,
                1.09746079e-02 + 0.j, 1.66323568e-02 + 0.00878693j,
                2.34592026e-02 - 0.00620169j, 4.66141255e-03 + 0.00175215j,
                1.99110490e-02 + 0.00810523j, 2.38036628e-02 + 0.00477377j,
                3.55530799e-02 - 0.03217207j, 7.83362427e-03 + 0.09056171j
            ]),
            y
        )
        npt.assert_almost_equal(
            np.array([
                0., 5.47072499, 4.65826468, 3.84580436, 3.03334404, 2.22088373,
                1.40842341, 0.5959631
            ]),
            ddt
        )
        npt.assert_almost_equal(
            np.array([
                [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
                 np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
                 np.nan, np.nan],
                [3.43348564e-02, 2.50579620e-02, 1.19662200e-02,
                 -2.94727045e-03, -1.74120657e-02, -2.92260318e-02,
                 -3.65905995e-02, -3.83845802e-02, -3.43348564e-02,
                 -2.50579620e-02, -1.19662200e-02, 2.94727045e-03,
                 1.74120657e-02, 2.92260318e-02, 3.65905995e-02,
                 3.83845802e-02],
                [1.44341709e-02, -1.03013686e-02, -2.90025061e-02,
                 -3.07143688e-02, -1.44341709e-02, 1.03013686e-02,
                 2.90025061e-02, 3.07143688e-02, 1.44341709e-02,
                 -1.03013686e-02, -2.90025061e-02, -3.07143688e-02,
                 -1.44341709e-02, 1.03013686e-02, 2.90025061e-02,
                 3.07143688e-02],
                [2.07973398e-03, -1.70523385e-02, -1.51310288e-02,
                 5.47155040e-03, 1.93187722e-02, 9.31439771e-03,
                 -1.21898408e-02, -1.86440980e-02, -2.07973398e-03,
                 1.70523385e-02, 1.51310288e-02, -5.47155040e-03,
                 -1.93187722e-02, -9.31439771e-03, 1.21898408e-02,
                 1.86440980e-02],
                [-5.78745261e-04, -6.84526323e-03, 5.78745261e-04,
                 6.84526323e-03, -5.78745261e-04, -6.84526323e-03,
                 5.78745261e-04, 6.84526323e-03, -5.78745261e-04,
                 -6.84526323e-03, 5.78745261e-04, 6.84526323e-03,
                 -5.78745261e-04, -6.84526323e-03, 5.78745261e-04,
                 6.84526323e-03],
                [-7.88064431e-03, -4.64290642e-05, 7.91617957e-03,
                 -6.01235248e-03, -3.31452421e-03, 8.54917948e-03,
                 -3.22873448e-03, -6.07801309e-03, 7.88064431e-03,
                 4.64290642e-05, -7.91617957e-03, 6.01235248e-03,
                 3.31452421e-03, -8.54917948e-03, 3.22873448e-03,
                 6.07801309e-03],
                [-8.98355133e-03, 2.36355000e-03, 5.64098686e-03,
                 -1.03411101e-02, 8.98355133e-03, -2.36355000e-03,
                 -5.64098686e-03, 1.03411101e-02, -8.98355133e-03,
                 2.36355000e-03, 5.64098686e-03, -1.03411101e-02,
                 8.98355133e-03, -2.36355000e-03, -5.64098686e-03,
                 1.03411101e-02],
                [7.83214697e-04, 2.38315858e-03, -5.18671757e-03,
                 7.20064583e-03, -8.11834103e-03, 7.80009240e-03,
                 -6.29435042e-03, 3.83035064e-03, -7.83214697e-04,
                 -2.38315858e-03, 5.18671757e-03, -7.20064583e-03,
                 8.11834103e-03, -7.80009240e-03, 6.29435042e-03,
                 -3.83035064e-03]
            ]),
            overxx
        )
        npt.assert_almost_equal(
            np.array([
                [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
                 np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
                 np.nan, np.nan],
                [1.01349184e-01, 7.28216597e-02, 3.32076982e-02,
                 -1.14618344e-02, -5.43864066e-02, -8.90311414e-02,
                 -1.10121692e-01, -1.14447213e-01, -1.01349184e-01,
                 -7.28216597e-02, -3.32076982e-02, 1.14618344e-02,
                 5.43864066e-02, 8.90311414e-02, 1.10121692e-01,
                 1.14447213e-01],
                [-3.30146687e-02, -2.38776946e-02, -7.53490856e-04,
                 2.28120976e-02, 3.30146687e-02, 2.38776946e-02,
                 7.53490856e-04, -2.28120976e-02, -3.30146687e-02,
                 -2.38776946e-02, -7.53490856e-04, 2.28120976e-02,
                 3.30146687e-02, 2.38776946e-02, 7.53490856e-04,
                 -2.28120976e-02],
                [-1.10248501e-02, -1.04642838e-02, 3.01583405e-03,
                 1.27725033e-02, 6.75981671e-03, -7.59876353e-03,
                 -1.25756585e-02, -2.02622881e-03, 1.10248501e-02,
                 1.04642838e-02, -3.01583405e-03, -1.27725033e-02,
                 -6.75981671e-03, 7.59876353e-03, 1.25756585e-02,
                 2.02622881e-03],
                [-1.01750840e-02, -3.51325297e-03, 1.01750840e-02,
                 3.51325297e-03, -1.01750840e-02, -3.51325297e-03,
                 1.01750840e-02, 3.51325297e-03, -1.01750840e-02,
                 -3.51325297e-03, 1.01750840e-02, 3.51325297e-03,
                 -1.01750840e-02, -3.51325297e-03, 1.01750840e-02,
                 3.51325297e-03],
                [-2.76514220e-03, 1.31658513e-03, 1.75747157e-03,
                 -2.66169563e-03, 2.79702072e-04, 2.44762094e-03,
                 -2.15303003e-03, -7.99763089e-04, 2.76514220e-03,
                 -1.31658513e-03, -1.75747157e-03, 2.66169563e-03,
                 -2.79702072e-04, -2.44762094e-03, 2.15303003e-03,
                 7.99763089e-04],
                [-8.07945869e-03, 1.76689965e-02, -1.69082759e-02,
                 6.24291649e-03, 8.07945869e-03, -1.76689965e-02,
                 1.69082759e-02, -6.24291649e-03, -8.07945869e-03,
                 1.76689965e-02, -1.69082759e-02, 6.24291649e-03,
                 8.07945869e-03, -1.76689965e-02, 1.69082759e-02,
                 -6.24291649e-03],
                [-2.26212648e-02, 2.95800949e-02, -3.20356237e-02,
                 2.96140192e-02, -2.26839487e-02, 1.23004527e-02,
                 -4.43242297e-05, -1.22185522e-02, 2.26212648e-02,
                 -2.95800949e-02, 3.20356237e-02, -2.96140192e-02,
                 2.26839487e-02, -1.23004527e-02, 4.43242297e-05,
                 1.22185522e-02]
            ]),
            overyy
        )

    def test_filter_data(self):
        data = np.array([0.6611, 0.7258, 0.6485, 0.6213, 0.5806, 0.6429, 0.6055,
                         0.6203, 0.6190, 0.5761, 0.5468, 0.5045])
        nfft = rj.upper_power_of_two(len(data))

        filtered = rj.filter_data(nfft, 2 / nfft, data)
        npt.assert_almost_equal(
            np.array([
                0.67617421, 0.68494727, 0.65859713, 0.62266901, 0.60921676,
                0.61720618, 0.6176552, 0.61536194, 0.60570956, 0.57824275,
                0.54485826, 0.51752379
            ]),
            filtered)

    def test_parse_file(self):
        (
            t_arr, offx, offy, lines, dt, tdi, line_rate
        ) = rj.parse_file(self.flat_text, 11, 2, True)
        self.assertEqual(128, tdi)
        self.assertEqual(70000, lines)
        self.assertEqual(0.000097, line_rate)
        # Need to test arrays, too

    def test_mask_frequencies(self):
        (
            t_arr, offx, offy, lines, dt, tdi, line_rate
        ) = rj.parse_file(self.flat_text, 11, 2)
        nfft = rj.upper_power_of_two(lines / 20)
        tt = np.linspace(0, 1, nfft, endpoint=False)
        nfftime = np.linspace(t_arr[0], t_arr[-1], nfft)
        xinterp, yinterp, x, y, ddt, overxx, overyy = rj.create_matrices(
            t_arr, offx, offy,
            dt, t_arr[-1] - t_arr[0], t_arr[0], nfft, nfftime, tt
        )
        rj.mask_frequencies(0.01, ddt, overxx, overyy)
        # Not sure how to test these large arrays :(

    @unittest.skip("Slow to run.")
    @patch('hiproc.resolve_jitter.Path')
    def test_pixel_smear(self, m_path):
        # This test is now slow because it basically engages the whole
        # process.  Need to figure out how to mock this data better.
        (
            nfftime,
            sample,
            line,
            linerate,
            tdi,
            t0,
            t1,
            t2,
            t3,
            offx_filtered,
            offy_filtered,
            xinterp,
            yinterp,
            min_avg_error,
            min_k,
            jittercheckx,
            jitterchecky,
            rh0,
        ) = rj.resolve_jitter(
            self.flat_text,
            True,
            self.flat_text,
            True,
            self.flat_text,
            True,
            20,
        )

        (
            max_smear_s,
            max_smear_l,
            max_smear_mag,
            dysdx,
            dyldx,
            xi
        ) = rj.pixel_smear(
            nfftime, sample, line, linerate, tdi
        )

        self.assertAlmostEqual(-0.8419796875685424, max_smear_s)
        self.assertAlmostEqual(-9.052314552452234, max_smear_l)
        self.assertAlmostEqual(9.09026020017119, max_smear_mag)

    # def test_write_data_for_plotting(self):
    #     m = mock_open()
    #     col1 = [1.234, 2.345, 3.456]
    #     col2 = [0.12345, 0.23456, 0.34567, 0.45678]
    #     with patch('hiproc.resolve_jitter.open', m):
    #         self.assertRaises(
    #             IndexError,
    #             rj.write_data_for_plotting, "dummy-file",
    #             ['one', 'two', 'three'], col1, col2
    #         )

    #         rj.write_data_for_plotting(
    #         "dummy-file", ['one', 'two'], col1, col2)
    #         print(m.mock_calls)
