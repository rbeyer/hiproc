#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for the `util` module."""

# Copyright 2019-2020, Ross A. Beyer (rbeyer@seti.org)
#
# Reuse is permitted under the terms of the license.
# The LICENSE file is at the top level of this library.

import argparse
import logging
import unittest
from pathlib import Path

import hiproc.util as util

from .utils import resource_check as rc

# Hardcoding these, but I sure would like a better solution.
HiRISE_img = Path("test-resources") / "PSP_010502_2090_RED5_0.img"
img = HiRISE_img


class TestResources(unittest.TestCase):
    """Establishes that the test image exists."""

    def test_resources(self):
        (truth, test) = rc(img)
        self.assertEqual(truth, test)


class TestUtil(unittest.TestCase):
    def test_parent_parser(self):
        self.assertIsInstance(util.parent_parser(), argparse.ArgumentParser)

    def test_logging(self):
        util.set_logger(loglvl="WARNING")
        logger = logging.getLogger()
        self.assertEqual(30, logger.getEffectiveLevel())

    def test_path_w_suffix(self):
        self.assertEqual("bar.foo", str(util.path_w_suffix(".foo", "bar.cub")))
        self.assertEqual(
            "foo.foo", str(util.path_w_suffix("foo.foo", "bar.cub"))
        )

    def test_pid_path_w_suffix(self):
        self.assertRaises(
            FileNotFoundError, util.pid_path_w_suffix, "foo.foo", "dummy"
        )
        self.assertRaises(
            FileNotFoundError,
            util.pid_path_w_suffix,
            ".foo",
            "PSP_010502_2090_RED5_0.dummy",
        )
        self.assertEqual(
            img, util.pid_path_w_suffix(".img", img.with_suffix(".foo.foo"))
        )

    def test_get_path(self):
        self.assertRaises(ValueError, util.get_path, "dummy")
        self.assertRaises(TypeError, util.get_path, "dummy", 42)
        self.assertRaises(
            NotADirectoryError, util.get_path, "dummy", "not_a_dir"
        )
        self.assertRaises(FileNotFoundError, util.get_path, "dummy", "hiproc")
        self.assertRaises(
            FileNotFoundError, util.get_path, "dummy", Path("hiproc")
        )
        self.assertRaises(
            FileNotFoundError, util.get_path, "dummy", ["hiproc", "tests"]
        )
        self.assertEqual(
            Path("hiproc/util.py"), util.get_path("util.py", "hiproc")
        )

    def test_conf_check_strings(self):
        self.assertIsNone(util.conf_check_strings("foo", ("YES", "NO"), "YES"))
        self.assertRaises(
            AssertionError,
            util.conf_check_strings,
            "foo",
            ("YES", "NO"),
            "MAYBE",
        )

    def test_conf_check_count(self):
        self.assertIsNone(
            util.conf_check_count("foo", 3, "what", ["one", "two", "three"])
        )
        self.assertRaises(
            AssertionError,
            util.conf_check_count,
            "foo",
            3,
            "what",
            ["one, two"],
        )

    def test_conf_check_bounds(self):
        self.assertIsNone(util.conf_check_bounds("foo", (0.1, 1), "0.5"))
        self.assertRaises(
            AssertionError, util.conf_check_bounds, "foo", (0.1, 1), "1.5"
        )

    def test_pause_slicer(self):
        self.assertEqual(slice(4, 7), util.pause_slicer(5, 3))
        self.assertEqual(slice(2, 5), util.pause_slicer(5, -3))
