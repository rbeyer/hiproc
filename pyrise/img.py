#!/usr/bin/env python
"""Functions to work with HiRISE .img files.
"""

# Copyright 2020, Ross A. Beyer (rbeyer@seti.org)
#
# Reuse is permitted under the terms of the license.
# The LICENSE file is at the top level of this library.

import os
import statistics

import numpy as np
import pvl
import kalasiris as isis


class LUT_Table:
    """HiRISE Look-Up Table class.

    This is very specific to the LUT tables in HiRISE EDRs.
    """

    def __init__(self, lut):
        if len(lut) != 256 and len(lut) != 1:
            raise ValueError(
                "HiRISE Lookup Tables must have 1 or 256 entries."
            )

        if len(lut) == 1 and lut == [[0, 0]]:
            # A single-item (0, 0) LUT indicates that
            # no LUT has been applied to these data.
            self.table = False
        else:
            # The final two table entries must have identical values:
            for i in (-1, -2):
                if lut[i][0] != lut[i][1]:
                    text = {-1: "final", -2: "second to last"}
                    raise ValueError(
                        f"The {text[i]} element of the table must have "
                        f"identical values, and this does not: {lut[i]}"
                    )

            # Pre-compute the lookups:
            self.table = list()
            for i, pair in enumerate(lut):
                if pair[0] == pair[1] or i + 1 == len(lut):
                    self.table.append(pair[0])
                else:
                    self.table.append(
                        int(statistics.mean((lut[i][0], lut[i + 1][0])))
                    )

            # Ensure the first lut element is zero
            if lut[0][0] == 0:
                self.table[0] = 0
            else:
                raise ValueError(
                    f"The first LUT pair doesn't have a zero: {lut[0]}"
                )

        self.lut = lut

    def __str__(self):
        if self.table:
            return str(self.table)
        else:
            return "No applied LUT."

    def unlut(self, pixel: int) -> int:
        if self.table:
            return self.table[pixel]
        else:
            # That means there's a single pair "(0, 0)" which indicates
            # that no lut has been applied, so this function should no-op.
            return pixel

    def lookup(self, dn: int) -> int:
        if self.table:
            for i, pair in enumerate(self.lut):
                if pair[0] <= dn <= pair[1]:
                    return i
            else:
                return 255

        else:
            # See comment in unlut() above.
            return dn

    def specialpixels(self):
        """A kalasiris.specialpixels.SpecialPixels named tuple is
        returned which describes the values in this LUT.

        These are slightly different than the ISIS special pixel values,
        and ISIS special pixels values can't be used to deal with values
        that are read from these EDRs, however, those values are provided
        for your use by this function.
        """
        # Based on the FixDns8() and FixDns16() functions in
        # https://github.com/USGS-Astrogeology/ISIS3/blob/dev/isis/src/mro/apps/hi2isis/main.cpp
        if self.table:
            return isis.specialpixels.SpecialPixels(
                Min=1,
                Null=self.unlut(255),
                Lrs=0,
                Lis=0,
                His=self.unlut(254),
                Hrs=self.unlut(255),
                Max=self.lut[-3][1],
            )
        else:
            return isis.specialpixels.SpecialPixels(
                Min=1,
                Null=65535,
                Lrs=0,
                Lis=0,
                His=16383,
                Hrs=65535,
                Max=16382,
            )


def byte_info(name: str, label: dict) -> tuple:
    if "SAMPLE_BITS" in label[name]:
        if label[name]["SAMPLE_BITS"] == 16:
            byte_width = 2
        elif label[name]["SAMPLE_BITS"] == 8:
            byte_width = 1
        else:
            raise ValueError(
                "Don't know how to handle '{}'".format(
                    label[name]["SAMPLE_BITS"]
                )
                + f"as a value of SAMPLE_BITS in {name}, should be 8 or 16."
            )
    else:
        raise ValueError("Can't determine how to set the width.")

    if label[name]["SAMPLE_TYPE"].startswith("MSB"):
        byte_order = "big"
    else:
        byte_order = "little"

    if "UNSIGNED" in label[name]["SAMPLE_TYPE"]:
        byte_signed = False
    else:
        byte_signed = True

    return byte_width, byte_order, byte_signed


def object_asarray(file_path: os.PathLike, name: str) -> np.ndarray:
    """Return the data object of *name* from *file_path* as a numpy array.

    Data objects are found by having both a a pointer parameter (like ^*name*)
    and an object in the label with *name*.

    Doesn't work for objects that have columns.  Since such objects don't
    have SAMPLE_BITS and SAMPLE_TYPE fields, this will result in a KeyError.
    """
    label = pvl.load(str(file_path))
    width, b_order, b_signed = byte_info(name, label)

    lut = LUT_Table(
        label["INSTRUMENT_SETTING_PARAMETERS"]["MRO:LOOKUP_CONVERSION_TABLE"]
    )

    table = list()
    with open(file_path, mode="rb") as f:
        # Since The first byte is location 1 (not 0) in PDS-counting.
        f.seek(int(label[f"^{name}"].value) - 1)
        for line in range(label[name]["LINES"]):
            row = list()
            # Skip over the prefix bytes.
            f.seek(label[name]["LINE_PREFIX_BYTES"], os.SEEK_CUR)
            for samp in range(label[name]["LINE_SAMPLES"]):
                dn = lut.unlut(
                    int.from_bytes(
                        f.read(width), byteorder=b_order, signed=b_signed
                    )
                )
                row.append(dn)
            table.append(row)
            # Skip over the suffix bytes.
            f.seek(label[name]["LINE_SUFFIX_BYTES"], os.SEEK_CUR)

    return np.array(table)


def overwrite_object(file_path: os.PathLike, name: str, arr: np.ndarray):
    """The file at *file_path* will have its *name* object overwritten
    with the data in *array*.

    Doesn't work for objects that have columns.  Since such objects don't
    have SAMPLE_BITS and SAMPLE_TYPE fields, this will result in a KeyError.
    """
    label = pvl.load(str(file_path))
    width, b_order, b_signed = byte_info(name, label)

    if arr.shape != (label[name]["LINES"], label[name]["LINE_SAMPLES"]):
        raise ValueError(
            f"The object {name} has different dimensions "
            f"({label[name]['LINES']}, "
            f"{label[name]['LINE_SAMPLES']})than "
            f"the array ({arr.shape})."
        )

    lut = LUT_Table(
        label["INSTRUMENT_SETTING_PARAMETERS"]["MRO:LOOKUP_CONVERSION_TABLE"]
    )

    with open(file_path, mode="r+b") as f:
        # Since The first byte is location 1 (not 0) in PDS-counting.
        f.seek(int(label[f"^{name}"].value) - 1)
        for line in range(arr.shape[0]):
            # Skip over the prefix bytes.
            f.seek(label[name]["LINE_PREFIX_BYTES"], os.SEEK_CUR)
            for samp in range(arr.shape[1]):
                f.write(
                    lut.lookup(arr[line][samp]).to_bytes(
                        width, byteorder=b_order, signed=b_signed
                    )
                )
            # Skip over the suffix bytes.
            f.seek(label[name]["LINE_SUFFIX_BYTES"], os.SEEK_CUR)
    return
