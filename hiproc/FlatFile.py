#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2020, Ross A. Beyer (rbeyer@seti.org)
#
# Reuse is permitted under the terms of the license.
# The AUTHORS file and the LICENSE file are at the
# top level of this library.

import collections
import csv
import itertools

from pathlib import Path


class FlatFile(collections.abc.Sequence):
    """Reads the output from the flat file created by hijitreg and provides
    it as a sequence.

    The resulting FlatFile object primarily behaves like a list, where
    that list represents the rows of the ISIS hijitreg output.  Each of those
    rows is a :func:`collections.namedtuple` which contains the elements
    of each row, referenced by their column names.

    The FlatFile object also has some dictionary-like capabilities, in
    order to get at the values listed in the ISIS hijitreg output in the
    section before the numerical output.
    """

    def __init__(self, flatinfo):
        self.flatinfo = flatinfo

        try:  # Assume it is a string?
            (self.dictionary, self.headers, self.the_list) = self.parse(
                flatinfo
            )
        except ValueError:
            # Maybe it is a file path?
            (self.dictionary, self.headers, self.the_list) = self.parse(
                Path(flatinfo).read_text()
            )

    def __str__(self):
        return str(self.dictionary)

    def __repr__(self):
        return f"{self.__class__.__name__}('{self.flatinfo}')"

    def __len__(self):
        return len(self.the_list)

    def __getitem__(self, key):
        try:
            return self.dictionary[key]
        except KeyError:
            return self.the_list[key]

    def __iter__(self):
        return self.the_list.__iter__()

    def __contains__(self, item):
        if item in self.dictionary:
            return True
        else:
            return item in self.the_list

    def keys(self):
        """Gets the keys from the initial portion of the flat output file.

           These will be items like 'FROM', 'MATCH', 'RegFile', etc.
        """
        return self.dictionary.keys()

    def values(self):
        """Gets the values from the initial portion of the flat output file."""
        return self.dictionary.values()

    @staticmethod
    def parse(flatinfo: str) -> tuple:
        """Takes a string (expecting the output of ISIS ``hijitreg``), and
        parses the output.

        A three-element namedtuple is returned: the first element
        is a *dictionary* of the name:value information at the
        top of the file, the second element is a *list* of the
        of the fields that decorate the top of the data
        rows, and the third element is a *list* of FlatRow
        :func:`collections.namedtuple` objects that represent each
        row of the flat output.

        The contents of the file that results from ISIS ``hijitreg``
        look like this::

            #          Hijitreg ISIS Application Results
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
            [... other lines like this ...]


            #  MATCH: /HiRISE/Users/audrie/pipe_dev_sandbox/HiRISE/Data/ResolveJitter/ESP/ORB_016200_016299/ESP_016213_2315/ESP_016213_2315_BG12.prehijack.cub
            #    Lines:       70000
            #    Samples:     2048
            #    FPSamp0:     0
            #    SampOffset:  0
            [... other lines like this ...]

            #   Total Registers:  10145 of 10500
            #   Number Suspect:   0
            #   Average Sample Offset: 0.3586  StdDev: 0.3605
            #   Average Line Offset:   3.9526 StdDev: 1.5403

            #  Column Headers and Data
                        FromTime  FromSamp  FromLine           MatchTime MatchSamp MatchLine        RegSamp        RegLine   RegCorr      B0_Offset       B1_Slope   B_RCorr
              316426108.25476718       341       210  316426108.15863681       341       210       340.6611       204.2240  0.883177       0.024032       0.273207  0.894237
              316426108.25476718      1023       210  316426108.15863681      1023       210      1022.5484       204.4079  0.805408       0.020917       0.286734  0.839953
              316426108.25476718      1705       210  316426108.15863681      1705       210      1704.6887       204.6171  0.823995       0.026994       0.258996  0.853913
            [... more comma-separated lines like above ...]

        But this function sees it like this::

            #          Hijitreg ISIS Application Results
            #    Coordinates are (Sample, Line) unless indicated
            #           k:  v
            #
            #    ****  Image Input Information ****
            #  k: v_FROM
            #  [here the key in the parsed dict will be 'FROM']
            #  [the value of the 'FROM' key will be another dictionary]
            #  [that dictionary will contain a 'path' key with the filename]
            #  [the other key-value pairs will come from these items]
            #    Lines:       70000
            #    Samples:     2048
            #    FPSamp0:     0
            #    SampOffset:  0
            #    LineOffset:  0
            #    CPMMNumber:  5
            #    Summing:     1
            [... other lines like this ...]

            #  k: v_MATCH
            # [again the MATCH key will have a similar dictionary of values.

            [... other lines like this ...]

            #   k: v
            #   k: v
            #   k: v k_StdDev: v
            #   k: v k_StdDev: v

            #  Column Headers and Data
            h h h h h h h h h h h h
            n n n n n n n n n n n n
            n n n n n n n n n n n n
            n n n n n n n n n n n n
            [... more whitespace-separated lines like above ...]

        Where each of the letters above is a string value that the parser
        reads.

        First, it takes all of the ``k`` and ``v`` elements and saves them
        as the keys and values in the returned dictionary.
        For the special FROM and MATCH keys, their values will each be a
        dictionary with their printed elements and a special key named
        path which will contain the string after "FROM:" or "MATCH:".

        If there are "StdDev" elements on a line that follows a key, they
        will be saved with a key value of the previous key plus " StdDev"
        such that a line like this::

            Average Sample Offset: 0.3586  StdDev: 0.3605

        Will be turned into two entries in the dictionary, one with a key
        value of "Average Sample Offset" and one with
        "Average Sample Offset StdDev".

        Second, the ``h`` elements are returned as the list of
        fieldnames.

        Third, it reads the lines with ``n`` and stores them as
        ``namedtuples`` in the returned list.
        """
        d = dict()
        flat_rows = []

        # This would be so much easier if they just wrote out JSON or something.
        # Gong to have to get a little weird to deal with the FROM and MATCH
        # dicts.
        indent_level = 0
        prev_dict = None
        for line in filter(
            lambda x: ":" in x, str(flatinfo).splitlines()
        ):
            (k, colon, v) = line.lstrip("#").partition(":")
            leading_spaces = len(k) - len(k.lstrip(' '))
            the_key = k.strip()
            the_value = v.strip()
            if leading_spaces < indent_level:
                prev_dict = None

            # Does this have a weird trailing standard deviation?
            val, stdsep, std = the_value.partition("StdDev:")
            if stdsep != "":
                the_value = val.strip()

            if prev_dict is not None:
                d[prev_dict].setdefault(the_key, the_value)

            else:
                if the_key == "FROM" or the_key == "MATCH":
                    prev_dict = the_key
                    the_value = {"path": the_value}

                d.setdefault(the_key, the_value)

            # Since this had a standard deviation, add it as a key:
            if stdsep != "":
                d.setdefault(the_key + " " + stdsep.strip(":"), std.strip())

            indent_level = leading_spaces

        if len(d) == 0:
            raise ValueError(
                "Could not extract any header information from the flat file."
            )

        # Read the data rows.
        dialect = csv.Dialect
        dialect.delimiter = ' '
        dialect.skipinitialspace = True
        dialect.quoting = csv.QUOTE_NONE
        dialect.lineterminator = '\n'

        reader = csv.DictReader(itertools.filterfalse(lambda x:
                                                      x.startswith('#') or
                                                      x.isspace() or
                                                      len(x) == 0,
                                                      flatinfo.splitlines()),
                                dialect=dialect)

        for row in reader:
            flat_rows.append(row)

        FlatParsed = collections.namedtuple(
            "FlatParsed", ["info", "fieldnames", "data"]
        )
        return FlatParsed(d, reader.fieldnames, flat_rows)
