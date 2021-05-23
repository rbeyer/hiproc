#!/usr/bin/env python
"""This module has tests for the HiRISE HiCal functions."""

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
import csv
import pkg_resources
import unittest
from pathlib import Path
from unittest.mock import patch

import pvl

import kalasiris as isis
import kalasiris.version as isisversion
import hiproc.hirise as hirise
import hiproc.EDR_Stats as edr
import hiproc.HiCal as hc

from .utils import resource_check as rc

# Hardcoding these, but I sure would like a better solution.
HiRISE_imgs = (
    "PSP_010502_2090_RED3_0.img",
    "PSP_010502_2090_RED3_1.img",
    "PSP_010502_2090_RED4_0.img",
    "PSP_010502_2090_RED4_1.img",
    "PSP_010502_2090_RED5_0.img",
    "PSP_010502_2090_RED5_1.img",
    "PSP_010502_2090_IR10_0.img",
    "PSP_010502_2090_IR10_1.img",
    "PSP_010502_2090_IR11_0.img",
    "PSP_010502_2090_IR11_1.img",
    "PSP_010502_2090_BG12_0.img",
    "PSP_010502_2090_BG12_1.img",
    "PSP_010502_2090_BG13_0.img",
    "PSP_010502_2090_BG13_1.img",
)
test_resources = Path("test-resources")
imgs = list(map(test_resources.joinpath, HiRISE_imgs))

gains = pvl.load(
    pkg_resources.resource_stream(
        "hiproc",
        "data/EDR_Stats_gains_config.pvl"
    )
)
conf_path = pkg_resources.resource_filename("hiproc", "data/HiCal.conf")
conf = pvl.load(conf_path)
# hgf_conf = Path("data") / "HiGainFx.conf"
nf_conf = pvl.load(
    pkg_resources.resource_stream("hiproc", "data/NoiseFilter.conf")
)


class TestResources(unittest.TestCase):
    """Establishes that the test image exists."""

    def test_resources(self):
        for f in imgs:
            with self.subTest(filepath=f):
                (truth, test) = rc(f)
                self.assertEqual(truth, test)


class TestBasic(unittest.TestCase):
    def test_furrow_setup(self):
        self.assertEqual(hc.furrow_setup("RED0", 4)[0], 8000)
        self.assertRaises(KeyError, hc.furrow_setup, "RED5", 1)

    def test_samp_setup(self):
        chansamp = hc.chan_samp_setup(1, 2)
        self.assertEqual(chansamp.samp[1], 511)
        self.assertEqual(chansamp.ssamp, 1)
        self.assertEqual(chansamp.nsamp, 502)

    def test_cut_size(self):
        self.assertEqual((6, 6), hc.cut_size(0, 10))
        self.assertEqual((50, 6), hc.cut_size(0, 255))
        self.assertEqual((6, 40), hc.cut_size(1, 511))

    def test_set_lines(self):
        self.assertEqual((1, 500), hc.set_lines(0, 0, 1, 500))
        self.assertEqual((2, 3998), hc.set_lines(2, 2, 2, 4000))
        self.assertRaises(ZeroDivisionError, hc.set_lines, 0, 0, 0, 0)

    def test_process_this(self):
        flags = (
            1,
            1,
            1,
            1,
            0,
            0,
            1,
            1,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            1,
            1,
            1,
            1,
            0,
            0,
            0,
            0,
            0,
            0,
        )
        self.assertEqual(0, hc.process_this(("RED5", 0), flags))
        self.assertEqual(1, hc.process_this(("RED0", 1), flags))
        self.assertEqual(0, hc.process_this(("BG13", 1), flags))

    def test_set_flags(self):
        my_c = {
            "HiCal_Noise_Processing": (
                1,
                1,
                1,
                1,
                0,
                0,
                1,
                1,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                1,
                1,
                1,
                1,
                0,
                0,
                0,
                0,
                0,
                0,
            ),
            "HiCal_Noise_Bin_DarkPixel_STD": (20.0, 30.0, 70.0, 100.0, 100.0),
            "HiCal_Noise_Bin_Mask_STD": (20.0, 30.0, 70.0, 100.0, 100.0),
            "HiCal_Noise_LIS_Count": 1,
            "HiCal_HPF_Cubenorm": "DIVIDE",
        }
        my_db = {
            "IMAGE_DARK_STANDARD_DEVIATION": 18.04943853369,
            "CAL_MASK_STANDARD_DEVIATION": 17.974952301437,
            "LOW_SATURATED_PIXELS": 0,
        }
        b = 1, 2, 4, 8, 16
        self.assertEqual(
            (False, False, True),
            hc.set_flags(my_c, my_db, ("RED5", 0), b.index(2)),
        )
        my_db["LOW_SATURATED_PIXELS"] = 1
        my_c["HiCal_HPF_Cubenorm"] = "SUBTRACT"
        self.assertEqual(
            (True, True, False),
            hc.set_flags(my_c, my_db, ("RED1", 1), b.index(2)),
        )

    def test_getHistVal(self):

        if tuple(isisversion.version_info()[:3]) < (4, 3, 0):
            col = "DN"
        else:
            col = "MaxExclusive"

        histogram = isis.Histogram(
            f"""Total Pixels:    2048000
Null Pixels:     0
Lis Pixels:      0

{col},Pixels,CumulativePixels,Percent,CumulativePercent
3889,1,1,4.88281e-05,4.88281e-05
3924,1,2,4.88281e-05,9.76563e-05
3960,2,4,9.76563e-05,0.000195313
3995,1,5,4.88281e-05,0.000244141
4030,4,9,0.000195313,0.000439453
6841,17215,2020810,0.840576,98.6724
6887,9753,2030563,0.476221,99.1486
7258,460,2045676,0.0224609,99.8865
7304,364,2046040,0.0177734,99.9043
7350,279,2046319,0.013623,99.9179
7759,108,2047734,0.00527344,99.987
7826,104,2047838,0.00507813,99.9921
7893,78,2047916,0.00380859,99.9959
7961,37,2047953,0.00180664,99.9977
8028,21,2047974,0.00102539,99.9987
8095,12,2047986,0.000585937,99.9993
8163,12,2047998,0.000585937,99.9999
8230,2,2048000,9.76563e-05,100"""
        )
        conf = dict(
            NoiseFilter_HighEnd_Percent=99.999,
            NoiseFilter_Hard_Tolmax=1.5,
            NoiseFilter_Hard_HighEnd_Percent=99.9,
        )
        self.assertEqual((0, 8095), hc.getHistVal(histogram, conf))

        bad_conf = dict(conf)
        bad_conf["NoiseFilter_HighEnd_Percent"] = 150
        self.assertRaises(ValueError, hc.getHistVal, histogram, bad_conf)

        histogram.dictionary["Lis Pixels"] = int(histogram["Total Pixels"]) / 2
        self.assertEqual((50, 7304), hc.getHistVal(histogram, conf))

    def test_FurrowCheck(self):
        vpnts = [1000, 1000, 1000]
        self.assertFalse(hc.FurrowCheck(vpnts, 0))

        vpnts0 = [1] + vpnts
        self.assertTrue(hc.FurrowCheck(vpnts0, 0))
        self.assertFalse(hc.FurrowCheck(vpnts0, 1))

        vpnts1 = vpnts + [1]
        self.assertTrue(hc.FurrowCheck(vpnts1, 1))
        self.assertFalse(hc.FurrowCheck(vpnts1, 0))

    def test_Cubenorm_Filter_filter_boxfilter(self):
        mylist = [4000, 4000, 0, 4000, 4000, 3800, 4000, 4000]

        truthlist = [
            3974.1882476041565,
            3974.0124192049366,
            3973.604595524187,
            3973.0384085658607,
            3972.4160735379446,
            3971.849930095303,
            3971.442163554957,
            3971.266366671061,
        ]

        self.assertEqual(
            truthlist, hc.Cubenorm_Filter_filter_boxfilter(mylist, 2, 50)
        )


class TestMock(unittest.TestCase):
    def test_get_bins_fromfiles(self):
        cubes = list()
        for p in imgs:
            cub = p.with_suffix(".PathOnly.cub")
            cubes.append(cub)
        with patch("hiproc.HiCal.isis.getkey_k", return_value="2"):
            with patch("pathlib.Path.glob", return_value=cubes):
                d = hc.get_bins_fromfiles(cubes[0])
                self.assertEqual(len(cubes) / 2, len(d))

    def test_check_destripe(self):
        bins = dict(RED3=2, RED4=2, RED5=2, IR10=4, IR11=4, BG12=4, BG13=4)
        powered = (
            "Off, Off, Off, On, On, On, On, On, On, On, Off, Off, Off, Off"
        )

        self.assertEqual(False, hc.check_destripe("dummy", 0, True, True))
        self.assertEqual(True, hc.check_destripe("dummy", 2, True, True))
        self.assertEqual(False, hc.check_destripe("dummy", 0, False, False))
        with patch("hiproc.HiCal.get_bins_fromfiles", return_value=bins):
            with patch("hiproc.HiCal.isis.getkey_k", return_value=powered):
                self.assertEqual(
                    False, hc.check_destripe("dummy", 2, None, None)
                )


class TestConf(unittest.TestCase):
    def test_conf_check(self):
        self.assertIsNone(hc.conf_check(conf))


class TestNeedCubenormStatsFile(unittest.TestCase):
    def setUp(self):
        cube = imgs[0].with_suffix(".TestHiCal_TestNeedCNstatsfile.cub")
        self.statsfile = cube.with_suffix(".stats")
        self.output = Path("dummy.tab")

        isis.hi2isis(imgs[0], to=cube)
        isis.cubenorm(cube, stats=self.statsfile)
        cube.unlink()

    def tearDown(self):
        with contextlib.suppress(FileNotFoundError):
            self.statsfile.unlink()
            self.output.unlink()
            Path("print.prt").unlink()

    def test_analyze_cubenorm_stats(self):
        self.assertEqual(
            (2901.0, 10491.599999999999),
            hc.analyze_cubenorm_stats(self.statsfile, 2)
        )

    @patch("hiproc.HiCal.isis.cubenormfile.DictWriter")
    @patch("hiproc.HiCal.NoiseFilter_cubenorm_writer")
    def test_NoiseFilter_cubenorm_edit(self, m_writer, m_DictWriter):
        conf = dict(
            NoiseFilter_Zap_Fraction=0.4, NoiseFilter_Nonvalid_Fraction=0.90
        )

        hc.NoiseFilter_cubenorm_edit(
            self.statsfile, self.output, 0, 2, conf, True
        )
        # print(m_writer.call_args_list)

    def test_Cubenorm_Filter_filter(self):
        vpnts = list()
        averages = list()
        with open(self.statsfile) as csvfile:
            reader = csv.DictReader(csvfile, dialect=isis.cubenormfile.Dialect)
            for row in reader:
                vpnts.append(int(row.pop("ValidPoints")))
                averages.append(float(row.pop("Average")))

        self.assertAlmostEqual(
            1.0036515,
            hc.Cubenorm_Filter_filter(averages, 5, 50, 0, True, vpnts, True)[
                -1
            ],
            6,
        )

    def test_Cubenorm_Filter(self):
        with patch("hiproc.HiCal.csv.DictWriter"):
            t = hc.Cubenorm_Filter(
                self.statsfile, self.output, False, 5, False, 0
            )
            self.assertAlmostEqual(15.4995326, t[0], 6)
            self.assertFalse(t[1])


class TestNeedISISCube(unittest.TestCase):
    def setUp(self):
        self.cube = imgs[0].with_suffix(".TestHiCal_TestNeedISISCube.cub")
        isis.hi2isis(imgs[0], to=self.cube)
        self.pid = hirise.get_ChannelID_fromfile(self.cube)
        self.binning = int(isis.getkey_k(self.cube, "Instrument", "Summing"))

    def tearDown(self):
        with contextlib.suppress(FileNotFoundError):
            self.cube.unlink()
            Path("print.prt").unlink()

    def test_furrow_nulling(self):
        ccdchan = (self.pid.get_ccd(), int(self.pid.channel))
        outcube = Path("test_furrow_nulling-out.cub")
        self.assertFalse(
            hc.furrow_nulling(self.cube, outcube, self.binning, ccdchan, False)
        )
        outcube.unlink()

    def test_mask(self):
        outcube = Path("test_mask-out.cub")
        self.assertIsNone(
            hc.mask(self.cube, outcube, 1200, 16383, self.binning, False)
        )
        outcube.unlink()

    def test_run_hical(self):
        myconf = dict(HiCal=None, NoiseFilter=None)
        myconf["HiCal"] = dict(
            HiCal_Normalization_Minimum="0.0",
            HiCal_Normalization_Maximum="1.5",
            HiCal_ISIS_Conf="hical.pipelines.conf",
            HiCal_ISIS_Conf_Noise="hical.noise.pipelines.conf",
        )
        myconf["NoiseFilter"] = dict(
            NoiseFilter_Raw_Min="1200", NoiseFilter_Raw_Max="16383"
        )
        outcube = Path("test_run_hical-out.cub")

        self.assertEqual(
            "Standard",
            hc.run_hical(
                self.cube,
                outcube,
                myconf,
                conf_path,
                3,
                3,
                self.binning,
                True,
                keep=False,
            ),
        )
        outcube.unlink()

    # def test_HiGainFx(self):
    #     outcube = Path('test_run_HiGainFx-out.cub')
    #     self.assertIsNone(hc.HiGainFx(self.cube, outcube, Path('data'),
    #                                   '0001', keep=True))
    #     # outcube.unlink()

    def test_highlow_destripe(self):
        myconf = dict(
            NoiseFilter_LPF_Line=501,
            NoiseFilter_LPF_Samp=9,
            NoiseFilter_LPF_Minper=5,
            NoiseFilter_HPF_Line=501,
            NoiseFilter_HPF_Samp=1,
            NoiseFilter_HPF_Minper=5,
        )
        outcube = Path("test_highlow_destripe-out.cub")
        self.assertIsNone(
            hc.highlow_destripe(
                self.cube,
                self.cube,
                outcube,
                myconf,
                isisnorm="",
                lnull=True,
                lhrs=True,
                lhis=True,
                llrs=True,
                llis=True,
                keep=False,
            )
        )
        outcube.unlink()

    def test_NoiseFilter_noisefilter(self):
        myconf = dict(
            NoiseFilter_Minimum_Value=0,
            NoiseFilter_Noise_Samp=7,
            NoiseFilter_Noise_Line=7,
        )
        outcube = Path("test_NoiseFilter_noisefilter-out.cub")
        self.assertIsNone(
            hc.NoiseFilter_noisefilter(
                self.cube,
                outcube,
                flattol=1,
                conf=myconf,
                maxval=4000,
                tolmin=3.5,
                tolmax=3.5,
            )
        )
        outcube.unlink()

    def test_NoiseFilter(self):
        myconf = dict(
            NoiseFilter_HighEnd_Percent=99.999,
            NoiseFilter_Hard_Tolmax=1.5,
            NoiseFilter_Hard_Tolmin=1.5,
            NoiseFilter_Hard_HighEnd_Percent=99.9,
            NoiseFilter_Zap_Fraction=0.4,
            NoiseFilter_Nonvalid_Fraction=0.90,
            NoiseFilter_LPF_Line=501,
            NoiseFilter_LPF_Samp=9,
            NoiseFilter_LPF_Minper=5,
            NoiseFilter_HPF_Line=501,
            NoiseFilter_HPF_Samp=1,
            NoiseFilter_HPF_Minper=5,
            NoiseFilter_Tolmin=3.5,
            NoiseFilter_Tolmax=3.5,
            NoiseFilter_Hard_Filtering=5,
            NoiseFilter_Flattol=1,
            NoiseFilter_Minimum_Value=0,
            NoiseFilter_Noise_Samp=7,
            NoiseFilter_Noise_Line=7,
            NoiseFilter_LPFZ_Line=5,
            NoiseFilter_LPFZ_Samp=5,
        )
        outcube = Path("test_NoiseFilter-out.cub")
        self.assertIsNone(
            hc.NoiseFilter(
                self.cube,
                outcube,
                myconf,
                minimum=None,
                maximum=None,
                zapc=False,
                keep=False,
            )
        )
        outcube.unlink()

    def test_Hidestripe(self):
        to_del = isis.PathSet()
        calcube = to_del.add(Path("test_Hidestripe-out.hical.cub"))
        isis.hical(self.cube, to=calcube)
        to_del.add(Path(str(self.pid)).with_suffix(".hical.log"))
        outcube = to_del.add(Path("test_Hidestripe-out.cub"))
        samps = int(isis.getkey_k(calcube, "Dimensions", "Samples"))

        self.assertRaises(
            KeyError,
            hc.Hidestripe,
            self.cube,
            outcube,
            self.binning,
            minimum=0.0,
            maximum=1.5,
            hidcorr="ADD",
            line_samples=samps,
            keep=False,
        )

        self.assertAlmostEqual(
            0.000101402295171637,
            hc.Hidestripe(
                calcube,
                outcube,
                self.binning,
                minimum=0.0,
                maximum=1.5,
                hidcorr="ADD",
                line_samples=samps,
                keep=False,
            ),
        )
        to_del.unlink()


class TestHiCal(unittest.TestCase):
    def setUp(self):
        self.cube = imgs[0].with_suffix(".TestHiCal.cub")
        self.pid = hirise.get_ChannelID_fromfile(self.cube)
        self.db = edr.EDR_Stats(imgs[0], self.cube, gains)
        self.binning = int(isis.getkey_k(self.cube, "Instrument", "Summing"))
        self.conf = conf
        # self.conf['HiGainFx'] = pvl.load(str(hgf_conf))['HiGainFx']
        self.conf["NoiseFilter"] = nf_conf["NoiseFilter"]

    def tearDown(self):
        with contextlib.suppress(FileNotFoundError):
            self.cube.unlink()
            Path("print.prt").unlink()

    def test_HiCal(self):
        outcube = Path("test_HiCal-out.cub")
        # ccdchan = (self.pid.get_ccd(), self.pid.channel)
        db = hc.HiCal(
            self.cube,
            outcube,
            self.db,
            self.conf,
            conf_path,
            bin2=False,
            bin4=False,
            keep=False,
        )
        # print(db)
        self.assertAlmostEqual(
            db["HIGH_PASS_FILTER_CORRECTION_STANDARD_DEVIATION"],
            0.00204738,
        )
        self.assertIsNone(db["DESTRIPED_DIFFERENCE_STANDARD_DEVIATION"])
        self.assertFalse(db["zapped"])
        self.assertEqual(db["hical_status"], "Standard")
        outcube.unlink()
