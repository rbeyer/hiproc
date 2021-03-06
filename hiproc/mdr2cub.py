#!/usr/bin/env python
"""Convert one of Alan's MDR files (in DN or I/F) to a cube file (in I/F)."""

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

# This program is based on equations in make_IOF.pro, by Alan Delamere,
# 2019-05-22, and on equations in the public domain ISIS program hical.

import argparse
import logging
import pkg_resources
import shutil
import statistics
import sys
from pathlib import Path
from osgeo import gdal

import pvl
import spiceypy

import hiproc.hirise as hirise
import hiproc.util as util
import kalasiris as isis

logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-o", "--output", required=False, default=".mdr.iof.cub"
    )
    parser.add_argument("-e", "--edr", required=True)
    parser.add_argument(
        "-c",
        "--conf",
        required=False,
        type=argparse.FileType('r'),
        default=pkg_resources.resource_stream(
            __name__,
            'data/hical.pipelines.conf'
        ),
    )
    parser.add_argument("mdr", metavar="MDR_file")
    parser.add_argument(
        "-l",
        "--log",
        required=False,
        default="WARNING",
        help="The log level to show for this program, can "
        "be a named log level or a numerical level.",
    )
    parser.add_argument(
        "-k",
        "--keep",
        required=False,
        default=False,
        action="store_true",
        help="Normally, the program will clean up any "
        "intermediary files, but if this option is given, it "
        "won't.",
    )

    args = parser.parse_args()

    util.set_logger(args.verbose, args.logfile, args.log)

    edr_path = Path(args.edr)
    mdr_path = Path(args.mdr)

    to_del = isis.PathSet()

    h2i_path = to_del.add(edr_path.with_suffix(".hi2isis.cub"))

    out_path = util.path_w_suffix(args.output, edr_path)

    # The first thing Alan's program did was to crop the image down to only the
    # 'imaging' parts.  We're not doing that so the resultant file has a
    # geometry similar to what comes out of ISIS hical.

    # Convert the EDR to a cube file
    isis.hi2isis(edr_path, to=h2i_path)

    # Convert Alan's MDR to a cube file
    mdr_cub_path = to_del.add(mdr_path.with_suffix(".alan.cub"))
    logger.info(f"Running gdal_translate {mdr_path} -of ISIS3 {mdr_cub_path}")
    gdal.Translate(str(mdr_cub_path), str(mdr_path), format="ISIS3")

    h2i_s = int(isis.getkey_k(h2i_path, "Dimensions", "Samples"))
    h2i_l = int(isis.getkey_k(h2i_path, "Dimensions", "Lines"))
    mdr_s = int(isis.getkey_k(mdr_cub_path, "Dimensions", "Samples"))
    mdr_l = int(isis.getkey_k(mdr_cub_path, "Dimensions", "Lines"))

    if h2i_s != mdr_s:
        label = pvl.load(str(h2i_path))
        hirise_cal_info = get_one(
            label, "Table", "HiRISE Calibration Ancillary"
        )

        buffer_pixels = get_one(hirise_cal_info, "Field", "BufferPixels")[
            "Size"
        ]
        dark_pixels = get_one(hirise_cal_info, "Field", "DarkPixels")["Size"]
        rev_mask_tdi_lines = hirise_cal_info["Records"]

        if h2i_s + buffer_pixels + dark_pixels == mdr_s:
            logger.info(
                f"The file {mdr_cub_path} has "
                f"{buffer_pixels + dark_pixels} more sample pixels "
                f"than {h2i_path}, assuming those are dark and "
                "buffer pixels and will crop accordingly."
            )
            if h2i_l + rev_mask_tdi_lines != mdr_l:
                logger.critical(
                    'Even assuming this is a "full" channel '
                    "image, this has the wrong number of lines. "
                    f"{mdr_cub_path} should have "
                    f"{h2i_l + rev_mask_tdi_lines}, but "
                    f"has {mdr_l} lines. Exiting"
                )
                sys.exit()
            else:
                crop_path = to_del.add(mdr_cub_path.with_suffix(".crop.cub"))
                # We want to start with the next pixel (+1) after the cal
                # pixels.
                isis.crop(
                    mdr_cub_path,
                    to=crop_path,
                    sample=buffer_pixels + 1,
                    nsamples=h2i_s,
                    line=rev_mask_tdi_lines + 1,
                )
                mdr_cub_path = crop_path
                mdr_l = int(isis.getkey_k(mdr_cub_path, "Dimensions", "Lines"))

        else:
            logger.critical(
                f"The number of samples in {h2i_path} ({h2i_s}) "
                f"and {mdr_cub_path} ({mdr_s}) are different. "
                "Exiting."
            )
            sys.exit()

    if h2i_l != mdr_l:
        logger.critical(
            f"The number of lines in {h2i_path} ({h2i_l}) "
            f"and {mdr_cub_path} ({mdr_l}) are different. "
            "Exiting."
        )
        sys.exit()

    # Convert the EDR to the right bit type for post-HiCal Pipeline:
    h2i_16b_p = to_del.add(h2i_path.with_suffix(".16bit.cub"))
    isis.bit2bit(
        h2i_path,
        to=h2i_16b_p,
        bit="16bit",
        clip="minmax",
        minval=0,
        maxval=1.5,
    )
    shutil.copyfile(h2i_16b_p, out_path)

    # If it is a channel 1 file, Alan mirrored it so that he could process
    # the two channels in an identical way (which we also took advantage
    # of above if the buffer and dark pixels were included), so we need to
    # mirror it back.
    cid = hirise.get_ChannelID_fromfile(h2i_16b_p)
    if cid.channel == "1":
        mirror_path = to_del.add(mdr_cub_path.with_suffix(".mirror.cub"))
        isis.mirror(mdr_cub_path, to=mirror_path)
        mdr_cub_path = mirror_path

    # Is the MDR in DN or I/F?
    maximum_pxl = float(
        pvl.loads(isis.stats(mdr_cub_path).stdout)["Results"]["Maximum"]
    )
    if maximum_pxl < 1.5:
        logger.info("MDR is already in I/F units.")
        mdr_16b_p = to_del.add(mdr_cub_path.with_suffix(".16bit.cub"))
        isis.bit2bit(
            mdr_cub_path,
            to=mdr_16b_p,
            bit="16bit",
            clip="minmax",
            minval=0,
            maxval=1.5,
        )
        isis.handmos(mdr_16b_p, mosaic=out_path)
    else:
        logger.info("MDR is in DN units and will be converted to I/F.")

        fpa_t = statistics.mean(
            [
                float(
                    isis.getkey_k(
                        h2i_16b_p, "Instrument", "FpaPositiveYTemperature"
                    )
                ),
                float(
                    isis.getkey_k(
                        h2i_16b_p, "Instrument", "FpaNegativeYTemperature"
                    )
                ),
            ]
        )
        print(f"fpa_t {fpa_t}")

        conf = pvl.load(args.conf)

        tdg = t_dep_gain(get_one(conf["Hical"], "Profile", cid.ccdname), fpa_t)
        suncorr = solar_correction()
        sclk = isis.getkey_k(
            h2i_16b_p, "Instrument", "SpacecraftClockStartCount"
        )
        target = isis.getkey_k(h2i_16b_p, "Instrument", "TargetName")
        suncorr = solar_correction(sunDistanceAU(sclk, target))
        sed = float(
            isis.getkey_k(h2i_16b_p, "Instrument", "LineExposureDuration")
        )
        zbin = get_one(conf["Hical"], "Profile", "GainUnitConversion")[
            "GainUnitConversionBinFactor"
        ]

        # The 'ziof' name is from the ISIS HiCal/GainUnitConversion.h, it is a
        # divisor in the calibration equation.
        print(f"zbin {zbin}")
        print(f"tdg {tdg}")
        print(f"sed {sed}")
        print(f"suncorr {suncorr}")
        ziof = zbin * tdg * sed * 1e-6 * suncorr
        eqn = f"\(F1 / {ziof})"  # noqa W605

        mdriof_p = to_del.add(mdr_cub_path.with_suffix(".iof.cub"))
        to_s = "{}+SignedWord+{}:{}".format(mdriof_p, 0, 1.5)
        isis.fx(f1=mdr_cub_path, to=to_s, equ=eqn)

        isis.handmos(mdriof_p, mosaic=out_path)

    if not args.keep:
        to_del.unlink()


def solar_correction(au=1.498) -> float:
    # au = 1.498  # Placeholder from Alan
    s = 1.5 / au  # Not sure about this, comes from the ISIS code
    return s * s


def sunDistanceAU(time: str, target: str) -> float:
    """Returns distance in AU between Sun and observed body from MRO."""

    base_kernel_path = Path(isis.environ["ISIS3DATA"]) / "base" / "kernels"
    lsk = sorted(Path(base_kernel_path / "lsk").glob("naif*.tls"))[-1]
    pck = sorted(Path(base_kernel_path / "spk").glob("de*.bsp"))[-1]
    sat = sorted(Path(base_kernel_path / "spk").glob("mar*.bsp"))[-1]

    sclk = sorted(
        Path(
            Path(isis.environ["ISIS3DATA"]) / "mro" / "kernels" / "sclk"
        ).glob("MRO_SCLKSCET.*.65536.tsc")
    )[-1]

    spiceypy.furnsh([str(lsk), str(pck), str(sat), str(sclk)])

    et = spiceypy.scs2e(-74999, time)

    targ = target.lower()
    if targ == "sky" or targ == "cal" or targ == "phobos" or targ == "deimos":
        targ = "mars"

    (sunv, lt) = spiceypy.spkpos(targ, et, "J2000", "LT+S", "sun")

    sunkm = spiceypy.vnorm(sunv)
    # Return in AU units
    return sunkm / 1.49597870691e8


def t_dep_gain(profile: dict, t: float) -> float:
    """Given the profile, and the FPA temperature in C,
    calculate the temperature dependent gain."""
    # Equivalent to getTempDepGain() in ISIS HiCal/GainUnitConversion.h
    # These equations are really g * (1 + t - baseT) * Q * absgainTDI
    # where these variables are from the hical.pipelines.conf file.
    #   g = FilterGainCorrection (DN/s) - Alan's numbers don't match
    #   baseT = IoverFbasetemperature (deg C)  - Alan's numbers are the same
    #   Q = QEpercentincreaseperC (1/deg C) - Alan's numbers are the same
    #   absgainTDI = AbsGain_TDI128 (? units) - Alan's numbers don't match
    #
    # These are Alan's numbers from the clean11 era.  Need to hunt down
    # authoritative values and/or verify hical.pipelines.conf which had a
    # syntax error.
    # zgain = dict(RED=157709797, IR=56467454, BG=115121269)
    # baseT = dict(RED=18.9, IR=18.9, BG=18.9)
    # QEpcntC = dict(RED=0.0005704, IR=0.002696, BG=0.00002295)
    # absgainTDI = dict(RED=6.37688, IR=6.99017, BG=7.00042)
    zgain = profile["FilterGainCorrection"]
    baseT = profile["IoverFbasetemperature"]
    QEpcntC = profile["QEpercentincreaseperC"]
    absgainTDI = profile["AbsGain_TDI128"]

    return zgain * (1 + (t - baseT) * QEpcntC * absgainTDI)


def get_one(conf, thing: str, name: str) -> dict:
    """Return the item for the named thing (of which there are multiple)."""
    for it in conf.getlist(thing):
        if it["Name"] == name:
            return it
