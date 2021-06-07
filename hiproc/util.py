#!/usr/bin/env python
"""This module contains PyRISE utility functions."""

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

import argparse
import collections
import logging
import os
import subprocess
import sys
import textwrap
from pathlib import Path

import hiproc
import hiproc.hirise as hirise


# These are primarily for use by HiCal and lisfix, primarily meant to be
# used by pause_slicer()
# 1st pixel = index 1 for these
ch_pause = (
    (252, 515, 778),  # Channel 0 pause point sample locations
    (247, 510, 773),  # Channel 1 pause point sample locations
)
ch_width = (
    (17, 17, 17),  # Number of pixels to cut from pause point
    (-17, -17, -17),
)


class main_exceptions:
    """This is a context manager that does the "standard handling"
    of exceptions at the top level in the main() functions in all of
    the programs in this library.  This is just a way to keep that code
    in one place rather than having to repeat the identical try-except blocks
    in each program.

    Using this context manager will cause information about exceptions to
    be printed to sys.stderr (or an alternate file-like of your choosing)
    and may cause sys.exit() to be called.

    Use it like this::

     with util.main_exceptions(args.verbose):
         isis.cubeit(fromlist=f, to='stacked.cub')

    """

    def __init__(self, verbose=0, file=sys.stderr, logger=None):
        self.verbose = verbose
        self.file = file
        self.logger = logger

    def __enter__(self):
        return None

    def __exit__(self, exc_type, exc_val, exc_tb):
        if exc_type is not None:
            if issubclass(exc_type, subprocess.CalledProcessError):
                print(
                    textwrap.dedent(
                        f"""\
                        Had a subprocess error:
                        {' '.join(exc_val.cmd)}
                        {exc_val.stdout}
                        {exc_val.stderr}
                        """
                    ),
                    file=self.file
                )
                if self.verbose >= 2:
                    return False
                else:
                    sys.exit(exc_val.returncode)
            if issubclass(exc_type, Warning):
                if self.logger is None:
                    print(exc_val, file=self.file)
                else:
                    self.logger.warning(exc_val)
                return True
            else:
                print("Error:", file=self.file)
                print(exc_val, file=self.file)
                if self.verbose >= 2:
                    return False
                else:
                    sys.exit(1)

        return


def parent_parser() -> argparse.ArgumentParser:
    """Returns a parent parser with common arguments for PyRISE programs."""
    parent = argparse.ArgumentParser(add_help=False)
    parent.add_argument(
        "-v",
        "--verbose",
        action="count",
        default=0,
        help="Displays additional information as processing progresses."
    )
    parent.add_argument(
        "-l",
        "--log",
        required=False,
        help=argparse.SUPPRESS,
        # help="The log level to show for this program, can "
        # "be a named log level or a numerical level.",
    )
    parent.add_argument(
        "--logfile",
        required=False,
        help="The log file to write log messages to in addition "
        "to the terminal.",
    )
    parent.add_argument(
        "-k",
        "--keep",
        required=False,
        default=False,
        action="store_true",
        help="Normally, the program will clean up any "
        "intermediary files, but if this option is given, it "
        "won't.",
    )
    parent.add_argument(
        '--version',
        action='version',
        version=f"hiproc version {hiproc.__version__}",
        help="Show hiproc version number."
    )
    return parent


def set_logger(verblvl=None, filename=None, loglvl=0) -> None:
    """Sets the log level and configuration for applications."""
    logger = logging.getLogger(__name__.split(".")[0])
    if loglvl is not None:
        if isinstance(loglvl, int):
            lvl = loglvl
        else:
            lvl = getattr(logging, loglvl.upper(), logging.WARNING)
    else:
        lvl_dict = {0: logging.WARNING, 1: logging.INFO, 2: logging.DEBUG}
        if verblvl in lvl_dict:
            lvl = lvl_dict[verblvl]
        else:
            lvl = lvl_dict[max(lvl_dict.keys())]

    logger.setLevel(lvl)

    ch = logging.StreamHandler()
    ch.setLevel(lvl)

    if lvl < 20:  # less than INFO
        formatter = logging.Formatter("%(name)s - %(levelname)s: %(message)s")
        # kalasiris logger
        k_logger = logging.getLogger("kalasiris")
        k_logger.setLevel(lvl)
        k_logger.addHandler(ch)
    else:
        formatter = logging.Formatter("%(levelname)s: %(message)s")

    ch.setFormatter(formatter)
    logger.addHandler(ch)

    if filename is not None:
        fh = logging.FileHandler(filename)
        fh.setLevel(lvl)
        fh.setFormatter(formatter)
        logger.addHandler(fh)
        if lvl < 20:
            k_logger.addHandler(fh)

    return


def path_w_suffix(in_path: str, template_path: os.PathLike) -> Path:
    """If the input starts with a '.' assume it is a suffix and return the
    template with the suffix replaced, otherwise return the input."""
    if in_path.startswith("."):
        return Path(template_path).with_suffix(in_path)
    else:
        return Path(in_path)


def pid_path_w_suffix(in_path: str, template_path: os.PathLike) -> Path:
    """A little extra twist to look for the db file."""
    p = path_w_suffix(in_path, template_path)
    if p.exists():
        return p
    elif in_path.startswith("."):
        pid = hirise.get_ChannelID_fromfile(template_path)
        t_path = Path(template_path)
        if t_path.is_dir():
            d = t_path
        else:
            d = t_path.parent
        pid_path = d / Path(str(pid)).with_suffix(in_path)
        if pid_path.exists():
            return pid_path
        else:
            raise FileNotFoundError(f"Could not find {pid_path}")
    else:
        raise FileNotFoundError(f"Could not find {p}")


def get_path(in_path: os.PathLike, search=None) -> Path:
    """Returns a path that can resovled in the search path
    (or list of paths)."""
    in_p = Path(in_path)
    if in_p.exists():
        return in_p

    search_paths = list()
    if isinstance(search, str):
        search_paths.append(Path(search))
    elif isinstance(search, Path):
        search_paths.append(search)
    elif isinstance(search, collections.abc.Sequence):
        for s in search:
            search_paths.append(Path(s))
    elif search is None:
        raise ValueError("You must provide a path or list of paths to search.")
    else:
        raise TypeError(
            f"Unfortunately, {search} isn't an os.PathLike or a list of them."
        )

    for p in search_paths:
        if p.is_dir():
            out_p = p / in_p
            if out_p.exists():
                return out_p
        else:
            raise NotADirectoryError(f"{p} is not a directory.")
    else:
        raise FileNotFoundError(f"Could not find {in_p} in {search_paths}.")


def conf_check_strings(
    conf_name: str, choices: tuple, conf_value: str
) -> None:
    assert (
        conf_value in choices
    ), f"The {conf_name} parameter can be {choices}, but was {conf_value}"


def conf_check_bool(conf_name: str, conf_value: bool) -> None:
    assert isinstance(conf_value, bool), (
        f"The {conf_name} parameter must be boolean but is {type(conf_value)} "
        f"with value {conf_value}."
    )


def conf_check_count(
    conf_name: str, count: int, what: str, conf_value: list
) -> None:
    assert len(conf_value) == count, (
        f"The {conf_name} parameter must have {count} entries, one for each "
        f"{what}, but it was {conf_value}"
    )


def conf_check_bounds(conf_name: str, bounds: tuple, conf_value: str) -> None:
    assert bounds[0] <= float(conf_value) <= bounds[1], (
        "The {} parameter must be between {} and {} inclusive, "
        "but was {}".format(conf_name, *bounds, conf_value)
    )


def pause_slicer(samp: int, width: int) -> slice:
    """Returns a slice object which satisfies the range of indexes for a pause
    point.

    The incoming numbers for samp are 1-based pixel numbers, so must
    subtract 1 to get a list index.
    The width values are the number of pixels to affect, including the
    pause point pixel.  If positive they start with the pause point pixel
    and count 'up.'  If negative, they start with the pause point pixel
    and count 'down.'
    """
    # We don't need to protect for indices less than zero or greater than the
    # length of the list, because slice objects can take values that would not
    # be valid for item access.
    if width > 0:
        s_start = samp - 1
        s_stop = s_start + width
    else:
        s_start = samp + width
        s_stop = samp
    return slice(s_start, s_stop)
