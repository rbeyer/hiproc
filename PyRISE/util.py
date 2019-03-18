#!/usr/bin/env python
"""This module contains utility functions."""

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

import collections
from pathlib import Path


class PathSet(collections.abc.MutableSet):
    """A class for containing a set of Path objects."""

    def __init__(self, iterable, temp_token=None):
        for value in iterable:
            if not isinstance(value, Path):
                raise TypeError('only accepts Path or pathlike objects')
        super().__init__(iterable)

    def add(self, elem):
        '''This variation on add() returns the element.'''
        if not isinstance(elem, Path):
            raise TypeError('only accepts Path or pathlike objects')
        super().add(elem)
        return elem

    def unlink(self):
        '''Just runs Path.unlink() on all members.'''
        for p in self:
            p.unlink()
