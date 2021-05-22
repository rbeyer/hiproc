#!/usr/bin/env python
"""Normalizes HiRISE Color products.

This program will normalize the BG and IR color across the left and right
halves of each set of color CCDs. It uses the HiSlither COLOR4 and COLOR5 cubes.


Data Flow
---------
Input Products:

- ``COLOR4`` and ``COLOR5`` files which are the result of HiSlither.

Output Products:

- creates UNFILTERED versions of the COLOR4 and COLOR5 files.
- creates HiColorNorm versions of the COLOR4 and COLOR5 files.

"""

# Copyright 2006-2020, Arizona Board of Regents on behalf of the Lunar and
# Planetary Laboratory at the University of Arizona.
#   - Orignal Perl program.
#
# Copyright 2020, Ross A. Beyer (rbeyer@seti.org)
#   - Elements of this Python program are are based on the original Perl
#     but the logic here is rewritten from scratch to emulate functionality.
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
#
# This program is based on HiColorNorm version 1.6.2 (2020/04/28),
# and on the Perl HiColorNorm program ($Revision: 1.28 $
#                                      $Date: 2020/04/28 16:17:31 $)
# by Eric Eliason as an employee of the University of Arizona.

import argparse
import collections
import csv
import json
import logging
import os
import pkg_resources
import re
import shutil
import statistics
from datetime import datetime
from itertools import repeat
from pathlib import Path

import pvl

import kalasiris as isis
import hiproc.hirise as hirise
import hiproc.util as util

logger = logging.getLogger(__name__)


class ColorCube(hirise.ObservationID):
    """A class for HiRISE multiband COLOR cubes."""

    def __init__(self, pathlike, dbs=None):
        self.path = Path(pathlike)
        match = re.search(r"COLOR(\d)", str(pathlike))
        if match:
            self.ccdnumber = match.group(1)
        else:
            raise ValueError(f"Could not extract a COLOR from {pathlike}")

        super().__init__(hirise.get_ObsID_fromfile(self.path))
        self.lines = int(isis.getkey_k(self.path, "Dimensions", "Lines"))
        self.samps = int(isis.getkey_k(self.path, "Dimensions", "Samples"))
        self.bands = int(isis.getkey_k(self.path, "Dimensions", "Bands"))

        if self.bands != 3:
            raise Exception("{} must have 3 bands.".format(self.path))

        # centers = isis.getkey_k(self.path, 'BandBin', 'Center')
        # match_c = re.search(r'900,\s+700,\s+500', centers)
        if not match:
            raise Exception(
                f"{self.path} must have BandBin/Center of 900, 700, 500."
            )
        self.band = {"IR": 1, "RED": 2, "BG": 3}

        self.pids = isis.getkey_k(
            self.path, "Mosaic", "SourceProductId"
        ).split(", ")
        self.ir_bin = self.get_binning("IR", dbs)
        self.red_bin = self.get_binning("RED", dbs)
        # self.bg_bin = get_binning('BG')

        self.mask_path = {"IR": None, "BG": None}
        self.crop_path = {"IR": None, "BG": None}
        self.nrm_path = {"IR": None, "BG": None}

        self.crop_sline = None
        self.crop_lines = None

        self.final_path = None

    def __str__(self):
        return "_".join([super().__str__(), self.get_ccd()])

    def __repr__(self):
        return f"{self.__class__.__name__}('{self.path}')"

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return super().__eq__(other) and self.ccdnumber == other.ccdnumber
        return False

    def __lt__(self, other):
        if isinstance(other, self.__class__):
            if super().__eq__(other):
                return int(self.ccdnumber) < int(other.ccdnumber)
            else:
                return super().__lt__(other)
        else:
            return NotImplemented

    def get_ccd(self) -> str:
        return f"COLOR{self.ccdnumber}"

    def get_obsid(self) -> hirise.ObservationID:
        return hirise.ObservationID(
            self.phase, self.orbit_number, self.target
        )

    def set_crop_lines(self, conf):
        sline = conf["HiColorNorm"]["HiColorNorm_Crop_Top"] + 1
        nline = (
            self.lines
            - conf["HiColorNorm"]["HiColorNorm_Crop_Top"]
            - conf["HiColorNorm"]["HiColorNorm_Crop_Bot"]
        )
        if nline < 5000:
            sline = 1
            nline = self.lines

        self.crop_sline = sline
        self.crop_lines = nline
        return (sline, nline)

    def get_binning_from(self, color_code, pids=None, db_paths=None):
        dbs = list()
        if pids is None:
            pids = self.pids

        parent = self.path.parent
        if db_paths is None:
            db_paths = list()
            for pid in pids:
                db_paths.extend(list(parent.glob(f"{pid}*.json")))

        for pid in pids:
            if color_code in pid:
                if len(db_paths) > 0:
                    for p in db_paths:
                        with open(p, "r") as f:
                            dbs.append(json.load(f))
                else:
                    for p in parent.glob(f"{pid}*.cub"):
                        temp_d = dict()
                        temp_d["PRODUCT_ID"] = isis.getkey_k(
                            p, "Archive", "ProductId"
                        )
                        temp_d["BINNING"] = int(
                            isis.getkey_k(p, "Instrument", "Summing")
                        )
                        dbs.append(temp_d)
        return dbs

    def get_binning(self, color_code, dbs=None, pids=None, db_paths=None):
        # Get the binning values for the component red, ir, and bg
        # observations.
        if pids is None:
            pids = self.pids

        if dbs is None:
            dbs = self.get_binning_from(color_code, pids, db_paths)

        bins = collections.Counter()
        for db in dbs:
            if (
                "PRODUCT_ID" in db
                and db["PRODUCT_ID"] in pids
                and color_code in db["PRODUCT_ID"]
            ):
                bins[db["BINNING"]] += 1

        if len(bins) != 1:
            raise Exception(
                "Gathered {} different ".format(len(bins))
                + f"binning values for the {color_code} channel, "
                "must be an error."
            )
        return int([*bins][0])

    def get_boxcar_size(self, color_bin: int) -> int:
        if color_bin > self.red_bin * 2:
            return 5
        elif color_bin > self.red_bin:
            return 3
        else:
            return 1


def arg_parser():
    parser = argparse.ArgumentParser(
        description=__doc__,
        parents=[util.parent_parser()],
        conflict_handler="resolve",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-o", "--output",
        required=False,
        default="_COLOR.cub",
        help="The filename to be used for the output color cube.  If it "
             "begins with an underscore ('_') it will be assumed to be a "
             "suffix that will be appended to a name derived from the "
             "observation. Default: %(default)s"
    )
    parser.add_argument(
        "-c",
        "--conf",
        required=False,
        type=argparse.FileType('r'),
        default=pkg_resources.resource_stream(
            __name__,
            'data/HiColorNorm.conf'
        ),
        help="Path to the HiColorNorm config file.  Defaults to "
             "HiColorNorm.conf distributed with the library."
    )
    parser.add_argument(
        "-n",
        "--nounfiltered",
        action="store_false",
        dest="Make_Unfiltered",
        help="Stops creation of an unfiltered cube.",
    )
    parser.add_argument(
        "cubes",
        type=Path,
        nargs=2,
        metavar="COLOR-cube",
        help="The COLOR4.cub and COLOR5.cub files"
    )
    return parser


def main():
    args = arg_parser().parse_args()

    util.set_logger(args.verbose, args.logfile, args.log)

    with util.main_exceptions(args.verbose):
        (ir_ratio, bg_ratio) = HiColorNorm(
            args.cubes,
            args.output,
            pvl.load(args.conf),
            make_unfiltered=args.Make_Unfiltered,
            keep=args.keep,
        )

    # Original Perl connects to HiCat to query the binning for CCDs
    #   and also whether there were furrows in the IR or BG components
    #   Finally, there's an insert into the Color_Processing_Statistics
    #       table that puts in the
    #       HICOLORNORM_RATIO_CORRECTION_STANDARD_DEVIATION for each FILTER
    # At the moment, we're not doing that here.


def conf_check(conf: dict) -> None:
    """Various checks on parameters in the configuration."""

    t = "HiColorNorm_NoiseFilter_IR10"
    util.conf_check_bool(t, conf["HiColorNorm"][t])

    t = "HiColorNorm_Make_Stitch"
    util.conf_check_bool(t, conf["HiColorNorm"][t])

    # This is not set via conf files, but via command line parameter.
    # t = 'HiColorNorm_Make_Unfiltered_Cube'
    # util.conf_check_bool(t, conf['HiColorNorm'][t])

    t = "HiColorNorm_Normalization_Minimum"
    util.conf_check_bounds(t, (-16384.0, 16364.0), conf["HiColorNorm"][t])

    t = "HiColorNorm_Normalization_Maximum"
    util.conf_check_bounds(t, (-16384.0, 16364.0), conf["HiColorNorm"][t])

    return


def set_outpath(out: os.PathLike, cubes) -> Path:
    # Check that they all have the same Obs ID:
    if len(set(c.get_obsid() for c in cubes)) != 1:
        raise ValueError(
            "These cube files don't all have the same "
            f"Observation ID: {cubes}"
        )

    if str(out).startswith("_"):
        return cubes[0].path.parent / Path(str(cubes[0].get_obsid()) + out)
    else:
        return Path(out)


def HiColorNorm(
    cubes: list,
    output,
    conf: dict,
    make_unfiltered=True,
    db_list=None,
    keep=False
):
    logger.info(f"HiColorNorm start: {cubes}")

    # GetConfigurationParameters()
    conf_check(conf)

    cubes = list(map(ColorCube, cubes, repeat(db_list)))
    cubes.sort()

    outcub_path = set_outpath(output, cubes)

    temp_token = datetime.now().strftime("HiColorNorm-%y%m%d%H%M%S")
    out_p = Path(outcub_path)

    furrow_flag = FurrowCheck(cubes)

    to_del = isis.PathSet()

    for i, _ in enumerate(cubes):
        cubes[i].set_crop_lines(conf)

    for i, c in enumerate(cubes):
        # Protect the processing from erroneous/spurious pixels
        mask_list = list()
        for b in (1, 2, 3):
            tmp_p = to_del.add(
                c.path.with_suffix(f".{temp_token}.temp{b}.cub")
            )
            isis.mask(
                f"{c.path}+{b}",
                mask=f"{c.path}+{b}",
                to=tmp_p,
                minimum=0.0,
                maximum=2.0,
                preserve="INSIDE",
            )
            mask_list.append(tmp_p)

        c.final_path = c.path.with_suffix(".HiColorNorm.cub")
        isis.cubeit_k(mask_list, to=c.final_path)

        (cubes[i].mask_path["IR"], cubes[i].crop_path["IR"]) = per_color(
            c, temp_token, "IR", keep=keep
        )
        (cubes[i].mask_path["BG"], cubes[i].crop_path["BG"]) = per_color(
            c, temp_token, "BG", keep=keep
        )

    ir_ratio_stddev = per_band(
        cubes, out_p, temp_token, "IR", furrow_flag, make_unfiltered, keep=keep
    )
    bg_ratio_stddev = per_band(
        cubes, out_p, temp_token, "BG", furrow_flag, make_unfiltered, keep=keep
    )

    if conf["HiColorNorm"]["HiColorNorm_Make_Stitch"]:
        # listpath = to_del.add(c.path.with_suffix(f'.{temp_token}.list.txt'))
        # listpath.write_text(
        #   '\n'.join(str(c.final_path) for c in cubes) + '\n')

        with isis.fromlist.temp([str(c.final_path) for c in cubes]) as f:
            isis.hiccdstitch(fromlist=f, to=out_p)

        for c in cubes:
            to_del.add(c.final_path)

    if not keep:
        to_del.unlink()
        for c in cubes:
            for cc in ("IR", "BG"):
                c.mask_path[cc].unlink()
                c.crop_path[cc].unlink()
                c.nrm_path[cc].unlink()

    logger.info("HiColorNorm done.")
    return ir_ratio_stddev, bg_ratio_stddev


def FurrowCheck(cubes, db_path=None) -> bool:
    # If any of these cubes have been flagged as having furrows
    if db_path is None:
        db_path = list()
        for c in cubes:
            parent = c.path.parent
            for pid in c.pids:
                if "RED" in pid:
                    continue
                for p in parent.glob(f"{pid}*.json"):
                    db_path.append(p)
    dbs = list()
    for p in db_path:
        with open(p) as f:
            dbs.append(json.load(f))

    furrow = any(bool(x.get("zapped", False)) for x in dbs)
    return furrow


def per_color(cube, temp_token, color_code, keep=False):

    tt = f"{temp_token}_{color_code}"

    numerator = "{}+{}".format(cube.path, cube.band[color_code])
    denominator = "{}+{}".format(cube.path, cube.band["RED"])

    # Generate the IR/RED and BG/RED ratios for each of the COLOR products
    rat_p = cube.path.with_suffix(f".{tt}.ratio.cub")
    isis.ratio(num=numerator, den=denominator, to=rat_p)

    # mask out invalid pixels
    # Possible future update: Make mosaic of ratio files then run cubnorm
    # correction on these for the unfiltered products, to avoid any null
    # pixels created during the ratio process.
    mask_p = cube.path.with_suffix(f".{tt}.mask.cub")
    isis.mask(
        rat_p,
        mask=rat_p,
        to=mask_p,
        minimum=0.0,
        maximum=4.0,
        preserve="INSIDE",
    )

    # Generate the crop files
    crop_p = cube.path.with_suffix(f".{tt}.ratcrop.cub")
    isis.crop(mask_p, to=crop_p, line=cube.crop_sline, nlines=cube.crop_lines)

    if not keep:
        rat_p.unlink()

    return mask_p, crop_p


def make_LR_mosaic(
    left_path, right_path, left_samps, mosaic_path, lines, samps
):
    isis.handmos(
        left_path,
        mosaic=mosaic_path,
        create="YES",
        nlines=lines,
        nsamp=samps,
        nbands=1,
        outsamp=1,
        outline=1,
        outband=1,
    )
    isis.handmos(
        right_path,
        mosaic=mosaic_path,
        outsamp=(left_samps + 1),
        outline=1,
        outband=1,
    )
    return


def cubenorm_stats(crpmos, mos, mosnrm, keep=False):
    # get cubenorm statistics on the crop mosaic files
    stats_p = Path(crpmos).with_suffix(".cubenorm.txt")
    isis.cubenorm(
        crpmos,
        stats=stats_p,
        format="table",
        direction="column",
        norm="average",
    )
    averages = list()
    with open(stats_p) as csvfile:
        reader = csv.DictReader(csvfile, dialect=isis.cubenormfile.Dialect)
        for row in reader:
            averages.append(float(row["Average"]))
    stddev = statistics.stdev(averages)

    # Now apply the statistics to the non-crop image
    isis.cubenorm(
        mos,
        fromstats=stats_p,
        to=mosnrm,
        preserve=True,
        statsource="table",
        format="table",
        direction="column",
        norm="average",
    )

    if not keep:
        stats_p.unlink()

    return stddev


def make_unfiltered(
    in_path, nrm_path, temp_token, color_code, color_band, keep=False
):
    tt = f"{temp_token}_{color_code}"
    in_p = Path(in_path)
    unfiltered_p = in_p.with_name(
        in_p.name.replace("COLOR", "UNFILTERED_COLOR")
    )
    shutil.copyfile(in_p, unfiltered_p)
    # run ISIS algebra program to create unfiltered
    # (normalized-IR/RED)*RED files
    alg_unfiltered_path = in_p.with_suffix(f".{tt}.algebra.cub")
    isis.algebra(
        nrm_path,
        from2="{}+2".format(in_p),
        to=alg_unfiltered_path,
        operator="MULTIPLY",
    )
    # Update the output file with the non-filtered normalized IR and BG bands
    isis.handmos(
        alg_unfiltered_path,
        mosaic=unfiltered_p,
        outsample=1,
        outline=1,
        matchbandbin=False,
        outband=color_band,
    )
    if not keep:
        alg_unfiltered_path.unlink()
    return unfiltered_p


def lpfz_triplefilter(
    from_path: os.PathLike, to_path: os.PathLike, keep=False
) -> None:
    to_del = isis.PathSet()
    from_p = Path(from_path)
    z1 = to_del.add(from_p.with_suffix(".z1.cub"))
    lpfz_filtering(from_p, z1, 11, 5)

    z2 = to_del.add(from_p.with_suffix(".z2.cub"))
    lpfz_filtering(z1, z2, 21, 9)

    lpfz_filtering(z2, to_path, 41, 11)

    if not keep:
        to_del.unlink()

    return


def lpfz_filtering(
    from_path: os.PathLike, to_path: os.PathLike, boxl: int, boxs: int
) -> None:
    isis.lowpass(
        from_path,
        to=to_path,
        lines=boxl,
        samples=boxs,
        filter="OUTSIDE",
        minopt="PERCENTAGE",
        minimum=25,
        low=0.00,
        high=2.0,
        null=True,
        HRS=True,
        HIS=True,
        LIS=True,
        LRS=True,
    )


def per_band(
    cubes, out_p, temp_token, color_code, furrow_flag, unfiltered, keep=False
) -> float:
    to_del = isis.PathSet()

    # Generate the handmos of both the ratio and the croped ratio files
    if len(cubes) == 2:
        maxlines = max(c.lines for c in cubes)
        maxsamps = sum(c.samps for c in cubes)

        # Generate the mosaic of the ratio image
        mos_p = to_del.add(
            out_p.with_suffix(f".{temp_token}.{color_code}.mos.cub")
        )
        make_LR_mosaic(
            cubes[0].mask_path[color_code],
            cubes[1].mask_path[color_code],
            cubes[0].samps,
            mos_p,
            maxlines,
            maxsamps,
        )

        # Generate the mosaic of the crop ratio image
        maxlines = max(c.crop_lines for c in cubes)
        crpmos_p = to_del.add(
            out_p.with_suffix(f".{temp_token}.{color_code}.moscrop.cub")
        )
        make_LR_mosaic(
            cubes[0].crop_path[color_code],
            cubes[1].crop_path[color_code],
            cubes[0].samps,
            crpmos_p,
            maxlines,
            maxsamps,
        )
    else:
        mos_p = cubes[0].mask_path[color_code]
        crpmos_p = cubes[0].crop_path[color_code]

    mosnrm_p = to_del.add(
        out_p.with_suffix(f".{temp_token}.{color_code}.mosnorm.cub")
    )
    ratio_stddev = cubenorm_stats(crpmos_p, mos_p, mosnrm_p, keep=keep)

    # from the handmos cubes, pull out the individual CCDs with crop
    if len(cubes) == 2:
        samp = 1
        for i, c in enumerate(cubes):
            cubes[i].nrm_path[color_code] = c.path.with_suffix(
                f".{temp_token}.{color_code}.norm.cub"
            )
            isis.crop(
                mosnrm_p,
                to=cubes[i].nrm_path[color_code],
                sample=samp,
                nsamples=c.samps,
            )
            samp += c.samps
    else:
        # If there was only one file then all we need to do is rename files
        cubes[0].nrm_path[color_code] = mosnrm_p

    if unfiltered:
        # Create the unfiltered color cubes, if needed
        for c in cubes:
            make_unfiltered(
                c.path,
                c.nrm_path[color_code],
                temp_token,
                color_code,
                c.band[color_code],
                keep=keep,
            )

    # low pass filter the files
    for c in cubes:
        # If $lpfz=1 then we're only interpolating null pixels caused by
        # noise or furrow removal

        lowpass_args = dict(
            minopt="PERCENT",
            minimum=25,
            low=0.00,
            high=2.0,
            null=True,
            hrs=True,
            lis=True,
            lrs=True,
        )

        # As the comment below states, I *think* that the OP had a copy
        # and paste error here, because even though it kept track of both
        # the IR and BG boxcar dimensions, they were always identical.
        logger.warning(
            "The original Perl calculates the boxcar size for "
            "both the IR and BG based only on the ratio of the "
            "IR to the RED."
        )
        if c.get_boxcar_size(c.ir_bin) == 1:
            lowpass_args["lines"] = 3
            lowpass_args["samples"] = 3
            lowpass_args["filter"] = "OUTSIDE"
        else:
            lowpass_args["lines"] = 1
            lowpass_args["samples"] = 1

        lpf_path = to_del.add(
            c.path.with_suffix(f".{temp_token}.{color_code}.lpf.cub")
        )
        isis.lowpass(c.nrm_path[color_code], to=lpf_path, **lowpass_args)

        # Perform lpfz filters to interpolate null pixels due to furrows or
        #   bad pixels
        if furrow_flag:
            lpfz_path = c.path.with_suffix(
                f".{temp_token}.{color_code}.lpfz.cub"
            )
            lpfz_triplefilter(lpf_path, lpfz_path, temp_token, keep=keep)
        else:
            lpfz_path = lpf_path

        # run ISIS algebra program to created (normalized-IR/RED)*RED files
        alg_path = to_del.add(
            c.path.with_suffix(f".{temp_token}.{color_code}.algebra.cub")
        )
        isis.algebra(
            lpfz_path, from2=f"{c.path}+2", to=alg_path, operator="MULTIPLY"
        )

        # Update the output file with the normalized IR and BG bands
        isis.handmos(
            alg_path,
            mosaic=c.final_path,
            outsample=1,
            outline=1,
            outband=c.band[color_code],
            matchbandbin=False,
        )

    if not keep:
        to_del.unlink()

    return ratio_stddev


if __name__ == "__main__":
    main()
