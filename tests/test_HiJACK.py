#!/usr/bin/env python
"""This module has tests for the HiRISE HiJACK functions."""

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
from pathlib import Path
from unittest.mock import call

# from unittest.mock import mock_open
from unittest.mock import patch
from unittest.mock import Mock

# import hiproc.hirise as hirise
import hiproc.HiJitReg as hjr
import hiproc.HiNoProj as hnp
import hiproc.HiJACK as hjk

hijackconf_path = Path("data") / "HiJACK.conf"
resjitconf_path = Path("data") / "ResolveJitter.conf"


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
        "StitchedProductIds": (
            "PSP_010502_2090_RED4_0, PSP_010502_2090_RED4_1"
        ),
    }
    return values[key]


class TestHiJACK(unittest.TestCase):
    @patch("hiproc.HiNoProj.isis.getkey_k", side_effect=getkey)
    def setUp(self, getkey_k):
        self.r4 = hnp.Cube("dummy/PSP_010502_2090_RED4.HiStitch.balance.cub")
        self.r3 = hnp.Cube("dummy/PSP_010502_2090_RED3.HiStitch.balance.cub")
        self.r5 = hnp.Cube("dummy/PSP_010502_2090_RED5.HiStitch.balance.cub")
        self.i1 = hnp.Cube("dummy/PSP_010502_2090_IR10.HiStitch.balance.cub")

        self.sets = dict(
            CCDs="RED3,RED4,RED5,BG12",
            Alt_1_CCDs="RED4,RED5,RED6,IR11",
            Alt_2_CCDs="RED3,RED4,RED5,IR10",
            Alt_3_CCDs="IR10,RED4,RED5,BG12",
        )

    @patch("hiproc.HiJACK.shutil.move")
    @patch("hiproc.HiJACK.isis.handmos")
    @patch("hiproc.HiJACK.isis.hijitreg")
    @patch("hiproc.HiJACK.shutil.copy")
    @patch("hiproc.HiJACK.isis.reduce")
    @patch("hiproc.HiJACK.isis.enlarge")
    def test_match_red(
        self, m_enlarge, m_reduce, m_copy, m_hjr, m_hand, m_move
    ):

        cubes = [self.r4, self.r3, self.r5, self.i1]
        for c in cubes:
            c.next_path = c.path.with_suffix(".next.cub")
        cubes.sort()

        hjk.match_red(cubes, self.r5, "dummy/flat")
        m_enlarge.assert_not_called()
        m_reduce.assert_called_once()
        self.assertEqual(3, m_copy.call_count)
        self.assertEqual(
            m_hjr.call_args_list,
            [
                call(
                    self.r4.next_path,
                    flat="dummy/flat",
                    match=self.r5.next_path,
                    row=20,
                )
            ],
        )

    def test_find_common(self):
        s1 = set(("RED3", "RED4", "RED5", "BG12"))
        s2 = set(("RED4", "RED5", "RED6", "IR11"))
        s3 = set(("RED3", "RED4", "RED5", "IR10"))
        s4 = set(("IR10", "RED4", "RED5", "BG12"))

        self.assertEqual("RED4", hjk.find_common(s1))
        self.assertEqual("RED5", hjk.find_common(s2))
        self.assertEqual("RED4", hjk.find_common(s3))
        self.assertEqual("RED4", hjk.find_common(s4))
        self.assertRaises(ValueError, hjk.find_common, set(("RED1", "RED2")))

    def test_ccd_set_generator(self):
        g = hjk.ccd_set_generator(self.sets)
        self.assertEqual(next(g), {"BG12", "RED3", "RED5", "RED4"})
        self.assertEqual(next(g), {"IR11", "RED6", "RED5", "RED4"})
        self.assertEqual(next(g), {"IR10", "RED3", "RED5", "RED4"})
        self.assertEqual(next(g), {"BG12", "IR10", "RED5", "RED4"})
        self.assertRaises(StopIteration, next, g)

    def test_determine_resjit_set(self):
        cubes = [self.r4, self.r3, self.r5, self.i1]
        cubes.sort()

        (rcubes, common) = hjk.determine_resjit_set(
            cubes, hjk.ccd_set_generator(self.sets)
        )
        self.assertEqual(self.r4, common)
        rcubes.sort()
        self.assertEqual(cubes, rcubes)

    @patch("hiproc.HiJitReg.run_HiJitReg")
    @patch("hiproc.HiNoProj.isis.getkey_k", side_effect=getkey)
    def test_make_flats(self, m_get, m_rhjr):
        cubes = [self.r4, self.r3, self.r5, self.i1]
        for c in cubes:
            c.next_path = c.path.with_suffix(".next.cub")
        cubes.sort()

        conf = {
            "AutoRegistration": {
                "ControlNet": {
                    "Control_Cols_Red": 1,
                    "Control_Cols_Color": 3,
                    "Control_Lines": 20,
                },
                "Algorithm": {"Tolerance": 0.7, "Increment": 0.1, "Steps": 2},
                "PatternChipRed": {"Samples": 17, "Lines": 45},
                "PatternChipColor": {"Samples": 100, "Lines": 60},
                "SearchChipRed": {"Samples": 30, "Lines": 70},
                "SearchChipColor": {"Samples": 140, "Lines": 100},
                "SearchLongerChip": {"Samples": 130, "Lines": 140},
                "AnaylyzeFlat": {"Minimum_Good": 0.20},
            },
            "Smoothing": {
                "Exclude_Limit": 2,
                "Badness_Limit": 1,
                "Boxcar_Length": 10,
            },
        }

        with patch("hiproc.HiJitReg.JitterCube") as m_jit:
            m_path = Mock(spec_set=Path)
            m_jit().flattab_path = m_path
            s = hjk.make_flats(cubes, self.r5, conf, "tt", keep=True)
            self.assertEqual(s, [m_path] * (len(cubes) - 1))

        with patch("hiproc.HiJitReg.Analyze_Flat", return_value=1):
            s = hjk.make_flats(cubes, self.r5, conf, "tt", keep=True)
            flat_list = list()
            for x in (self.r3, self.r4, self.i1):
                p = hjr.JitterCube.get_pair_name(x, self.r5.get_ccd())
                flat_list.append(self.r5.path.parent / (p + ".flat.tab"))
            self.assertEqual(flat_list, s)

    @unittest.skip("ResolveJitter needs more work")
    @patch("hiproc.HiJACK.subprocess.run")
    def test_ResolveJitter(self, m_run):
        cubes = [self.r4, self.r3, self.r5, self.i1]
        for c in cubes:
            c.next_path = c.path.with_suffix(".next.cub")
        cubes.sort()

        flat_list = list()
        for x in (self.r3, self.r4, self.i1):
            p = hjr.JitterCube.get_pair_name(x, self.r5.get_ccd())
            flat_list.append(self.r5.path.parent / (p + ".flat.tab"))

        with patch(
            "hiproc.HiJACK.determine_resjit_set", return_value=(cubes, self.r5)
        ):
            with patch("hiproc.HiJACK.make_flats", return_value=flat_list):
                m_path = Mock(spec_set=Path)
                m_path.parent = self.r5.path.parent
                conf = dict(
                    ResolveJitter=self.sets,
                    AutoRegistration={"ControlNet": {"Control_Lines": 10}},
                )
                hjk.ResolveJitter(cubes, conf, m_path, "tt", keep=True)
                rjp = "/Users/rbeyer/software/HiPrecision_Oleg/resolveJitter"
                self.assertEqual(
                    m_run.call_args_list,
                    [
                        call(
                            [
                                rjp,
                                str(m_path.parent),
                                str(self.r5.get_obsid()),
                                "10",
                                flat_list[0].relative_to(m_path.parent),
                                "-1",
                                flat_list[1].relative_to(m_path.parent),
                                "-1",
                                flat_list[2].relative_to(m_path.parent),
                                "-1",
                            ],
                            check=True,
                        )
                    ],
                )

    @patch("hiproc.HiNoProj.handmos_side")
    @patch("hiproc.HiNoProj.fix_labels")
    def test_mosaic_dejittered(self, m_fix, m_hand):
        out_p = "output_path"
        prodid = "product_id"
        self.assertRaises(
            ValueError,
            hjk.mosaic_dejittered,
            [self.r4, self.r3, self.r5],
            out_p,
            prodid,
        )
        hjk.mosaic_dejittered([self.r4, self.r5], out_p, prodid)
        self.assertEqual(
            m_hand.call_args_list,
            [call([self.r4, self.r5], self.r5, out_p, left=True)],
        )
        self.assertEqual(
            m_fix.call_args_list,
            [call([self.r4, self.r5], out_p, str(self.r5), prodid)],
        )

    @patch("hiproc.HiNoProj.isis.getkey_k", side_effect=getkey)
    @patch("hiproc.HiJACK.Path.mkdir")
    @patch("hiproc.HiJACK.ResolveJitter")
    @patch("hiproc.HiNoProj.fix_labels")
    @patch("hiproc.HiJACK.isis.PathSet")
    @patch("hiproc.HiJACK.plot_flats")
    @patch("hiproc.HiJACK.isis.fromlist.temp")
    @patch("hiproc.HiJACK.isis.cubeit")
    @patch("hiproc.HiJACK.isis.handmos")
    @patch("hiproc.HiJACK.mosaic_dejittered")
    @patch("hiproc.HiJACK.shutil.copyfile")
    @patch("hiproc.HiJACK.isis.hijitter")
    @patch("hiproc.HiNoProj.copy_and_spice")
    @patch("hiproc.HiNoProj.conf_check")
    @patch("hiproc.HiJACK.pvl.load")
    @patch("hiproc.HiJACK.match_red")
    def test_HiJACK(
        self,
        m_match,
        m_pvl,
        m_ccheck,
        m_cns,
        m_hijit,
        m_copy,
        m_mosdejit,
        m_handmos,
        m_cubeit,
        m_fromlist,
        m_plot,
        m_PathSet,
        m_fixlab,
        m_ResJit,
        m_mkdir,
        m_getkey
    ):
        confd = "confdir"
        outd = "outdir"
        with patch(
            "hiproc.HiNoProj.add_offsets",
            side_effect=[
                (
                    [self.r4, self.r5],
                    [Path("RED5-RED4flat1"), Path("RED3-RED4flat2")],
                ),
                ([self.i1, self.r3], ["irflat1", "irflat2"]),
                ([self.i1, self.r3], ["bgflat1", "bgflat2"]),
            ],
        ):
            hjk.HiJACK(
                [
                    Path("dummy/PSP_010502_2090_RED3.HiStitch.balance.cub"),
                    Path("dummy/PSP_010502_2090_RED4.HiStitch.balance.cub"),
                    Path("dummy/PSP_010502_2090_RED5.HiStitch.balance.cub"),
                    Path("dummy/PSP_010502_2090_IR10.HiStitch.balance.cub")
                ],
                confd,
                outd,
                base_ccd_number=4,
                plot=False,
                keep=True,
            )

            m_match.assert_called_once()
            self.assertEqual(
                m_pvl.call_args_list,
                [
                    call(confd + "/ResolveJitter.conf"),
                    call(confd + "/HiJACK.conf"),
                ],
            )
            m_ccheck.assert_called_once()
            self.assertEqual(m_cns.call_count, 4)
            m_hijit.assert_called_once()

            oid = str(self.r5.get_obsid())
            r5_p = Path(
                f"{outd}/"
                + self.r5.path.with_suffix("")
                .with_suffix(".dejittered.cub")
                .name
            )
            r3_p = Path(
                f"{outd}/"
                + self.r3.path.with_suffix("")
                .with_suffix(".dejittered.cub")
                .name
            )
            rnp = Path(f"{outd}/{oid}_RED.NOPROJ.cub")
            r45np = Path(f"{outd}/{oid}_RED4-5.NOPROJ.cub")
            irnp = Path(f"{outd}/{oid}_IR.NOPROJ.cub")
            bgnp = Path(f"{outd}/{oid}_BG.NOPROJ.cub")
            self.assertEqual(
                m_copy.call_args_list,
                [
                    call(r5_p, rnp),
                    call(None, r45np),
                    call(None, irnp),
                    call(None, bgnp),
                ],
            )

            self.assertEqual(
                m_mosdejit.call_args_list,
                [
                    call([self.r4, self.r5], r45np, f"{oid}_RED4-5"),
                    call([self.i1, self.r3], irnp, f"{oid}_IR"),
                    call([self.i1, self.r3], bgnp, f"{oid}_BG"),
                ],
            )

            r4_p = Path(
                f"{outd}/"
                + self.r4.path.with_suffix("")
                .with_suffix(".dejittered.cub")
                .name
            )
            self.assertEqual(
                m_handmos.call_args_list,
                [
                    call(
                        None,
                        mosaic=rnp,
                        outband=1,
                        outline=1,
                        outsample=1,
                        priority="ontop",
                    )
                ],
            )

            irbnp = Path(f"{outd}/{oid}_IRB.NOPROJ.cub")
            self.assertEqual(m_cubeit.call_args_list[0][1]["proplab"], r45np)
            self.assertEqual(m_cubeit.call_args_list[0][1]["to"], irbnp)

            i1_p = Path(
                f"{outd}/"
                + self.i1.path.with_suffix("")
                .with_suffix(".dejittered.cub")
                .name
            )

            self.assertEqual(
                m_fromlist.call_args_list,
                [
                    call(
                        [
                            self.r3.path,
                            self.r4.path,
                            self.r5.path,
                            self.i1.path,
                        ]
                    ),
                    call([r3_p, r4_p, r5_p, i1_p]),
                    call([irnp, r45np, bgnp]),
                ],
            )

            m_plot.called_once()
            m_PathSet.called_once()

            self.assertEqual(
                m_fixlab.call_args_list,
                [call([self.r4, self.r5], rnp, self.r5, f"{oid}_RED")],
            )

            self.assertEqual(
                m_ResJit.call_args_list[0][0][0],
                [self.r3, self.r4, self.r5, self.i1],
            )
            self.assertEqual(
                m_ResJit.call_args_list[0][0][2],
                # Path(f"{outd}/{oid}_jitter_cpp.txt"),
                Path(f"{outd}/{oid}_jitter_py.txt"),
            )
