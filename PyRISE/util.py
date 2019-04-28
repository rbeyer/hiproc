#!/usr/bin/env python
"""This module contains PyRISE utility functions."""

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

import argparse
import collections
import logging
import os
from pathlib import Path

import PyRISE.hirise as hirise


def parent_parser() -> argparse.ArgumentParser:
    '''Returns a parent parser with common arguments for PyRISE programs.'''
    parent = argparse.ArgumentParser(add_help=False)
    parent.add_argument('--db',         required=False, default='.HiCat.json',
                        help="The .json file to use.  Optionally, if it starts "
                        "with a '.' it is considered an extension and will be "
                        "swapped with the input file's extension to find the "
                        ".json file to use.")
    parent.add_argument('-l', '--log',  required=False, default='WARNING',
                        help="The log level to show for this program, can be a "
                        "named log level or a numerical level.")
    parent.add_argument('-k', '--keep', required=False, default=False,
                        help="Normally, the program will clean up any "
                        "intermediary files, but if this option is given, it "
                        "won't.")
    return parent


def set_logging(i) -> None:
    '''Sets the log level and basic configuration.'''

    log_level = 0
    if isinstance(i, int):
        log_level = i
    else:
        log_level = getattr(logging, i.upper(), logging.WARNING)

    logging.basicConfig(format='%(levelname)s: %(message)s', level=log_level)
    return


def path_w_suffix(in_path: str, template_path: os.PathLike) -> Path:
    """If the input starts with a '.' assume it is a suffix and return the
       template with the suffix replaced, otherwise return the input."""
    if in_path.startswith('.'):
        return Path(template_path).with_suffix(in_path)
    else:
        return Path(in_path)


def pid_path_w_suffix(in_path: str, template_path: os.PathLike) -> Path:
    """A little extra twist to look for the db file."""
    p = path_w_suffix(in_path, template_path)
    if p.exists():
        return p
    elif(in_path.startswith('.')):
        pid = hirise.get_ProdID_fromfile(template_path)
        t_path = Path(template_path)
        d = None
        if t_path.is_dir():
            d = t_path
        else:
            d = t_path.parent
        pid_path = d / Path(str(pid)).with_suffix(in_path)
        if pid_path.exists():
            return pid_path
        else:
            raise FileNotFoundError(f'Could not find {pid_path}')
    else:
        raise FileNotFoundError(f'Could not find {p}')


def get_path(in_path: os.PathLike, search=None) -> Path:
    '''Returns a path that can resovled in the search path (or list of paths).'''
    inPath = Path(in_path)
    if inPath.exists():
        return inPath

    search_paths = list()
    if isinstance(search, str):
        search_paths.append(Path(search))
    elif isinstance(search, Path):
        search_paths.append(search)
    elif isinstance(search, collections.abc.Sequence):
        for s in search:
            search_paths.append(Path(s))
    elif search is None:
        raise ValueError('You must provide a path or list of paths to search.')
    else:
        raise TypeError(f"Unfortunately, {search} isn't an os.PathLike "
                        " or a list of them.")

    for p in search_paths:
        if p.is_dir():
            outPath = p / inPath
            if outPath.exists():
                return outPath
        else:
            raise NotADirectoryError(f'{p} is not a directory.')
    else:
        raise FileNotFoundError(f'Could not find {inPath} in {search_paths}.')


def conf_check_strings(conf_name: str, choices: tuple, conf_value: str) -> None:
    assert conf_value in choices, ('The {} parameter can be {}, '
                                   'but was {}'.format(conf_name,
                                                       str(choices),
                                                       conf_value))


def conf_check_bool(conf_name: str, conf_value: bool) -> None:
    assert isinstance(conf_value, bool), ('The {} parameter must be '
                                          'boolean but is {} with value '
                                          '{}'.format(conf_name,
                                                      type(conf_value),
                                                      conf_value))


def conf_check_count(conf_name: str, count: int, what: str, conf_value: list) -> None:
    assert len(conf_value) == count, ('The {} parameter must have {} '
                                      'entries, one for each {}, but '
                                      'it was {}'.format(conf_name,
                                                         str(count),
                                                         what,
                                                         conf_value))


def conf_check_bounds(conf_name: str, bounds: tuple, conf_value: str) -> None:
    assert (float(conf_value) >= bounds[0] and
            float(conf_value) <= bounds[1]), ('The {} parameter must be '
                                              'between {} and {} inclusive, '
                                              'but was {}'.format(conf_name,
                                                                  *bounds,
                                                                  conf_value))
