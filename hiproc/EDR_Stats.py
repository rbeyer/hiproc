#!/usr/bin/env python
"""Create an ISIS cube from a HiRISE EDR .img file and record some
statistics.

These are the functionalities that the EDR_Stats pipeline performs that are
reproduced here:

* Convert a HiRISE EDR Product to an ISIS cube for subsequent pipeline
  processing using the ISIS system for much of the work.
* Gather image statistics about the observation and write them out
  to a .json file for other programs to use.
* Calculate the Signal-to-Noise Ratio (SNR).

In general, this program is the first, simple step of the HiRISE processing
chain, and prepares data for the next step: HiCal.


Data Flow
---------
Input Products:

- Any HiRISE PDS EDR, typically ending in ``.img``.

Output Products:

- A ``.cub`` file for each ``.img`` file processed.
- A ``.json`` file for each ``.img`` file processed containing summary
    information about the image.

"""

# Copyright 2004-2020, Arizona Board of Regents on behalf of the Lunar and
# Planetary Laboratory at the University of Arizona.
#   - Orignal Perl program.
#
# Copyright 2020-2021, Ross A. Beyer (rbeyer@seti.org)
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
# This Python program is based on EDR_Stats version 2.18.2 (2021/04/09)
# and on the Perl EDR_Stats program ($Revision: 1.44 $
#                                    $Date: 2021/04/09 22:56:18 $
# and on the check_lut.pm program (2.4, 2020/04/28)
# by Eric Eliason, Audrie Fennema, R. King, and Richard Leis as employees of
# the University of Arizona.

import argparse
import concurrent.futures
import hashlib
import json
import logging
import math
import os
import pkg_resources
import sys
from itertools import repeat
from pathlib import Path

import pvl.new as pvl  # Need the new PVL objects, which are pickleable.
import kalasiris as isis

import hiproc.hirise as hirise
import hiproc.util as util

logger = logging.getLogger(__name__)


def arg_parser():
    parser = argparse.ArgumentParser(
        description=__doc__,
        parents=[util.parent_parser()],
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-o", "--output",
        required=False,
        default=".EDR_Stats.cub",
        help="Output filename.  Optionally, if it starts with a '.' it is "
             "considered a suffix and will be swapped with the input"
             "file's suffix to determine the file to write. "
             "Default: %(default)s"
    )
    parser.add_argument(
        "--db",
        required=False,
        default=".HiCat.json",
        help="The .json file to output.  Optionally, if it "
        "starts with a '.' it is considered a suffix and "
        "will be swapped with the input file's suffix to "
        "determine the file to write. Default: %(default)s",
    )
    parser.add_argument(
        "--histmin",
        required=False,
        default=0.01,
        help="The minimum percentage above which DN are counted.  This only "
             "affects the value of STD_DN_LEVELS in the output .json file. "
             "If it is lower, more DN will be counted. Default: %(default)s"
    )
    parser.add_argument(
        "--histmax",
        required=False,
        default=99.99,
        help="The minimum percentage below which DN are counted. See --histmin "
             "for more information. Default: %(default)s"
    )
    parser.add_argument(
        "-g",
        "--gains",
        required=False,
        type=argparse.FileType('r'),
        default=pkg_resources.resource_stream(
            __name__,
            'data/EDR_Stats_gains_config.pvl'
        ),
        help="Path to the gains config PVL file.  Defaults to "
             "EDR_Stats_gains_config.pvl distributed with the library.",
    )
    parser.add_argument(
        "--max_workers",
        default=None,
        type=int,
        help="If more than one image is provided to work on, this program "
             "will engage multiprocessing to parallelize the work.  This "
             "multiprocessing will default to the number of processors on the "
             "machine.  If you want to throttle this to use less resources on "
             "your machin, indicate the number of processors you want to use."
    )
    parser.add_argument(
        "img",
        metavar="some.img-file",
        nargs="+",
        help="More than one can be listed here.",
    )
    return parser


def main():
    args = arg_parser().parse_args()

    util.set_logger(args.verbose, args.logfile, args.log)

    if len(args.img) > 1:
        if not args.output.startswith("."):
            logger.critical(
                "With more than one input IMG file, the --output must start "
                f"with a period, and it does not: {args.output}"
            )
            sys.exit()

        if not args.db.startswith("."):
            logger.critical(
                "With more than one input IMG file, the --db must start with "
                f"a period, and it does not: {args.db}"
            )
            sys.exit()

    gainsinfo = pvl.load(args.gains)

    with util.main_exceptions(args.verbose):
        if len(args.img) == 1:
            # No need to fire up the multiprocessing for just one image.
            db_path = write_json(
                EDR_Stats(
                    args.img[0],
                    util.path_w_suffix(args.output, args.img[0]),
                    gainsinfo,
                    args.histmin,
                    args.histmax,
                    keep=args.keep
                ),
                args.db,
                args.img[0],
            )
            logger.info(f"Wrote {db_path}")
        else:
            with concurrent.futures.ProcessPoolExecutor(
                max_workers=args.max_workers
            ) as executor:
                for img, histats, in zip(args.img, executor.map(
                    EDR_Stats,
                    args.img,
                    map(util.path_w_suffix, repeat(args.output), args.img),
                    repeat(gainsinfo),
                    repeat(args.histmin),
                    repeat(args.histmax),
                    repeat(args.keep),
                )):
                    db_path = write_json(histats, args.db, img)
                    logger.info(f"Wrote {db_path}")
    return


def write_json(d: dict, outpath: str, template_path: os.PathLike) -> Path:
    """
    Writes out a Python dict as JSON to *outpath*.

    :param d: Python dictionary to serialize to a JSON file.
    :param outpath: If it starts with a "." assume it is a suffix that
        should be swapped with the suffix on *img* to get the output filename,
        otherwise use as an outpath.
    :param template_path: Pathlike that could have its suffix replaced.
    """
    json_path = util.path_w_suffix(outpath, template_path)

    with open(json_path, "w") as f:
        json.dump(d, f, indent=0, sort_keys=True)

    return json_path


def EDR_Stats(
    img: os.PathLike,
    out_path: os.PathLike,
    gainsinfo: dict,
    histmin=0.01,
    histmax=99.99,
    keep=False,
) -> dict:
    cid = hirise.get_ChannelID_fromfile(img)

    logger.info(f"{cid}: EDR_Stats start: {img}")
    try:
        logger.info(
            f"{cid}: The LUT for this file is: " + str(check_lut(img))
        )
    except KeyError as err:
        logger.error(
            f"{cid}: The LUT header area is either corrupted or has a gap."
        )
        raise err

    # Convert to .cub
    isis.hi2isis(img, to=out_path)

    histat_complete = isis.histat(
        out_path,
        useoffsets=True,
        leftimage=0,
        rightimage=1,
        leftcalbuffer=3,
        rightcalbuffer=1,
        leftcaldark=3,
        rightcaldark=1,
        leftbuffer=3,
        rightbuffer=1,
        leftdark=3,
        rightdark=1,
    )
    histats = parse_histat(histat_complete.stdout)

    # Get some info from the new cube:
    histats["PRODUCT_ID"] = isis.getkey_k(out_path, "Archive", "ProductId")
    histats["IMAGE_LINES"] = int(
        isis.getkey_k(out_path, "Dimensions", "Lines")
    )
    histats["LINE_SAMPLES"] = int(
        isis.getkey_k(out_path, "Dimensions", "Samples")
    )
    histats["BINNING"] = int(isis.getkey_k(out_path, "Instrument", "Summing"))

    histats["STD_DN_LEVELS"] = get_dncnt(out_path, histmin, histmax, keep=keep)
    histats["IMAGE_SIGNAL_TO_NOISE_RATIO"] = calc_snr(
        out_path, gainsinfo, histats, cid=cid
    )
    histats["GAP_PIXELS_PERCENT"] = (
        histats["GAP_PIXELS"]
        / (int(histats["IMAGE_LINES"]) * int(histats["LINE_SAMPLES"]))
    ) * 100.0

    tdi_bin_check(out_path, histats, cid=cid)
    lut_check(out_path, histats)

    logger.info(f"{cid}: EDR_Stats done: {out_path}")
    return histats


def parse_histat(pvltext: str) -> dict:
    """Parse the output of histat into a dictionary"""

    p = pvl.loads(pvltext)
    d = dict()

    # Image Area Statistics
    d["IMAGE_MEAN"] = p["IMAGE"]["Average"]
    d["IMAGE_STANDARD_DEVIATION"] = p["IMAGE"]["StandardDeviation"]
    d["IMAGE_MINIMUM"] = p["IMAGE"]["Minimum"]
    d["IMAGE_MAXIMUM"] = p["IMAGE"]["Maximum"]
    d["GAP_PIXELS"] = p["IMAGE"]["NullPixels"]
    d["LOW_SATURATED_PIXELS"] = p["IMAGE"]["LisPixels"]
    d["HIGH_SATURATED_PIXELS"] = p["IMAGE"]["HisPixels"]

    # Calibration Reverse-readout Statistics
    d["CAL_REVERSE_MEAN"] = p["CAL_REVERSE"]["Average"]
    d["CAL_REVERSE_STANDARD_DEVIATION"] = p["CAL_REVERSE"]["StandardDeviation"]
    d["CAL_REVERSE_MINIMUM"] = p["CAL_REVERSE"]["Minimum"]
    d["CAL_REVERSE_MAXIMUM"] = p["CAL_REVERSE"]["Maximum"]
    d["CAL_REVERSE_LIS"] = p["CAL_REVERSE"]["LisPixels"]

    # Calibration Mask Statistics
    d["CAL_MASK_MEAN"] = p["CAL_MASK"]["Average"]
    d["CAL_MASK_STANDARD_DEVIATION"] = p["CAL_MASK"]["StandardDeviation"]
    d["CAL_MASK_MINIMUM"] = p["CAL_MASK"]["Minimum"]
    d["CAL_MASK_MAXIMUM"] = p["CAL_MASK"]["Maximum"]

    # Calibration Ramp Statistics
    d["CAL_RAMP_MEAN"] = p["CAL_RAMP"]["Average"]
    d["CAL_RAMP_STANDARD_DEVIATION"] = p["CAL_RAMP"]["StandardDeviation"]
    d["CAL_RAMP_MINIMUM"] = p["CAL_RAMP"]["Minimum"]
    d["CAL_RAMP_MAXIMUM"] = p["CAL_RAMP"]["Maximum"]

    # Image Dark Reference Statistics
    d["IMAGE_DARK_MEAN"] = p["IMAGE_DARK"]["Average"]
    d["IMAGE_DARK_STANDARD_DEVIATION"] = p["IMAGE_DARK"]["StandardDeviation"]
    d["IMAGE_DARK_MINIMUM"] = p["IMAGE_DARK"]["Minimum"]
    d["IMAGE_DARK_MAXIMUM"] = p["IMAGE_DARK"]["Maximum"]

    # Image Buffer Area
    d["IMAGE_BUFFER_MEAN"] = p["IMAGE_BUFFER"]["Average"]
    d["IMAGE_BUFFER_STANDARD_DEVIATION"] = p["IMAGE_BUFFER"][
        "StandardDeviation"
    ]
    d["IMAGE_BUFFER_MINIMUM"] = p["IMAGE_BUFFER"]["Minimum"]
    d["IMAGE_BUFFER_MAXIMUM"] = p["IMAGE_BUFFER"]["Maximum"]
    d["IMAGE_BUFFER_LIS"] = p["IMAGE_BUFFER"]["LisPixels"]

    # Calibration Image Dark Reference
    d["CAL_DARK_MEAN"] = p["CAL_DARK"]["Average"]
    d["CAL_DARK_STANDARD_DEVIATION"] = p["CAL_DARK"]["StandardDeviation"]
    d["CAL_DARK_MINIMUM"] = p["CAL_DARK"]["Minimum"]
    d["CAL_DARK_MAXIMUM"] = p["CAL_DARK"]["Maximum"]

    # Calibration Image Buffer Area
    d["CAL_BUFFER_MEAN"] = p["CAL_BUFFER"]["Average"]
    d["CAL_BUFFER_STANDARD_DEVIATION"] = p["CAL_BUFFER"]["StandardDeviation"]
    d["CAL_BUFFER_MINIMUM"] = p["CAL_BUFFER"]["Minimum"]
    d["CAL_BUFFER_MAXIMUM"] = p["CAL_BUFFER"]["Maximum"]
    d["CAL_BUFFER_LIS"] = p["CAL_BUFFER"]["LisPixels"]

    # Calibration Dark Ramp Area
    d["CAL_DARK_RAMP_MEAN"] = p["CAL_DARK_RAMP"]["Average"]
    d["CAL_DARK_RAMP_STANDARD_DEVIATION"] = p["CAL_DARK_RAMP"][
        "StandardDeviation"
    ]
    d["CAL_DARK_RAMP_MINIMUM"] = p["CAL_DARK_RAMP"]["Minimum"]
    d["CAL_DARK_RAMP_MAXIMUM"] = p["CAL_DARK_RAMP"]["Maximum"]

    # Image Post Ramp Area
    d["IMAGE_POST_RAMP_MEAN"] = p["IMAGE_POSTRAMP"]["Average"]
    d["IMAGE_POST_RAMP_STANDARD_DEVIATION"] = p["IMAGE_POSTRAMP"][
        "StandardDeviation"
    ]
    d["IMAGE_POST_RAMP_MINIMUM"] = p["IMAGE_POSTRAMP"]["Minimum"]
    d["IMAGE_POST_RAMP_MAXIMUM"] = p["IMAGE_POSTRAMP"]["Maximum"]

    return d


def get_dncnt(cub: os.PathLike, hmin=0.01, hmax=99.99, keep=False) -> int:
    """Extract DN count from the histogram of a cub file"""
    # I'm not sure about this method.
    # The code below is the exact logic that the original Perl program has,
    # but this is just counting the number of histogram bins
    # that are within the boundaries, not the number of DN.
    # And the # of bins is automatically computed by isis.hist,
    # so could be different for each cube.
    # logger.info(get_dncnt.__doc__)

    histfile = Path(cub).with_suffix(".hist")
    if not histfile.is_file():
        isis.hist(cub, to=histfile)

    h = isis.Histogram(histfile)

    count: int = 0
    for row in h:
        if hmin <= float(row.CumulativePercent) <= hmax:
            count += 1

    if not keep:
        histfile.unlink()
    return count


def calc_snr(
    cub: os.PathLike, gainsinfo: dict, histats: dict, cid=None
) -> float:
    """Calculate the signal to noise ratio."""
    if cid is None:
        cid = hirise.get_ChannelID_fromfile(cub)
    # logger.info(f"{cid}: " + calc_snr.__doc__)

    ccdchan = f"{cid.get_ccd()}_{cid.channel}"

    # gainspvl = pvl.load(str(gainsfile))
    gain = float(gainsinfo["Gains"][ccdchan]["Bin" + str(histats["BINNING"])])

    img_mean = float(histats["IMAGE_MEAN"])
    lis_pixels = float(histats["LOW_SATURATED_PIXELS"])
    buf_mean = float(histats["IMAGE_BUFFER_MEAN"])

    snr = -9999
    r = 90
    # Note from original file about r:
    # 150 e *Changed value to 90 e- 1/31/2012 to bring closer
    # to HIPHOP value for read noise. SM

    if 0 == lis_pixels and img_mean > 0.0 and buf_mean > 0.0:
        s = (img_mean - buf_mean) * gain
        snr = s / math.sqrt(s + r * r)
        logger.info(f"{cid}: Calculation of Signal/Noise Ratio:")
        logger.info(f"{cid}: \tIMAGE_MEAN:        {img_mean}")
        logger.info(f"{cid}: \tIMAGE_BUFFER_MEAN: {buf_mean}")
        logger.info(f"{cid}: \tR (electrons/DN):  {r}")
        logger.info(f"{cid}: \tGain:              {gain}")
        logger.info(f"{cid}: Signal/Noise ratio: {snr}")

    return snr


def check_lut(img: os.PathLike):
    """Checks whether a stored look up table (LUT) matches a known LUT."""
    # Original author of this function was Robert King, in December 2006.
    # logger.info(check_lut.__doc__)

    img_pvl = pvl.load(str(img))
    lut_type = img_pvl["INSTRUMENT_SETTING_PARAMETERS"][
        "MRO:LOOKUP_TABLE_TYPE"
    ]
    if "STORED" in lut_type:
        lut = dict()
        # dec.28.06 zero-filled LUT (ie, the image was not LUT'ed)
        lut["897256b6709e1a4da9daba92b6bde39ccfccd8c1"] = None

        # dec.19.06 build
        lut["35c8318042da3c30949f3ecb9f3876de524f158f"] = 300
        lut["874f6c26fe21577bc88f8986ec5643ff43747908"] = 301
        lut["02ae8c3834160614a65190cfee95fd62a2227447"] = 302
        lut["919d8f16363d1161b1843783c713564110cef5cc"] = 303
        lut["086adb2152c5f47741ce16e128ad3533fff1fa31"] = 304
        lut["020b40e1aab312ddf9fce61d5eb7cacd33bbeaec"] = 305
        lut["b45de94586e5ab539cc1c0b72d0e0267765e3b69"] = 306
        lut["8e4cd2ba1646d117afb9af658b9f762fcb0bcc48"] = 307
        lut["6227fbb77374afe54ba5de26a5af2cfb954f7783"] = 308
        lut["31d88a70dc4a6c9cf4991eb3f7a07c3b5df624a9"] = 309
        lut["166d23b2c7e71acebfaf07798a697670bff55df9"] = 310
        lut["d0f29aae50c21e443d7d188e7937b249fd8ba447"] = 311
        lut["a44541a295899f7fe55bf038dee9d5c130190c08"] = 312
        lut["58914a9cf58f4abbf7adc60357d8daa0d8104794"] = 313
        lut["bd416e8935556cd22e6caae1e9541ea42fc32fe7"] = 314
        lut["332a862b5963d55716bc00185906ab0195372e03"] = 315
        lut["b7d5169700fa7925adfc9ab4463c4c0c64cb9c31"] = 316
        lut["de35fb5719fabf3da400de015cf2c3155bbfe1ac"] = 317
        lut["55b78c9ca2822a6f87b3c3884d95034d4abbb725"] = 318
        lut["2aaf567e760eb313eca50cff9b3ed74f035af0d5"] = 319
        lut["43a31577737b140eecc327d09f0900c182dddbe3"] = 320
        lut["79d29aa8e13dd1991f810894821ebc518bcf3ff3"] = 321
        lut["dca0da083569fc25b910611f7d43a72738923d59"] = 322
        lut["809c4a81505eced9203e4648d7905a3de9c71f7e"] = 323
        lut["349bf694a6b11ee3e47956189aebeab234be4e9e"] = 324
        lut["063a2d3ed43bb1b47e73a38fbe49091ac4183c3c"] = 325
        lut["9ea3b154b83925f65c36888d0981b84eed35baa6"] = 326
        lut["4ca66adf72eefe25bc53edc801c1c71e99cf4ab9"] = 327

        # may.16.05 build
        lut["23f9889760466065c78b070aa83097229c391ba1"] = 200
        lut["d46c82527435d5b8a3b4a329f4a80ec80ed3a6e5"] = 201
        lut["a814e4a5a679d1ca710c92b745f79c6d2a1c2f43"] = 202
        lut["b67f6fb3adb08c25955878e5f685846c3ede6194"] = 203
        lut["7c03fd37caac821554e6c9ebe593ca60200cee4b"] = 204
        lut["6c21207a2e78517a8e069f99efd8f3c549ace95d"] = 205
        lut["b21a8cb5dfee2be1eaf5e7c8ce2b7b30c7f385ef"] = 206
        lut["801e83161358bf9a4a00a76c8970454521c7aff8"] = 207
        lut["7223df1f14d33a7a57b96d2d9d39b6ae526a1f8b"] = 208
        lut["90b01c7c6055ea580293771601e6f31cc966671c"] = 209
        lut["a23833f7421ae17c066010b1150feb696b7dbb2c"] = 210
        lut["9b20c521c050c37faa959d84e6f1031f4ad42ac6"] = 211
        lut["80b6d126c6ed96d1635bb0755ed701cefb001c72"] = 212
        lut["598df42422e6167304db3d9e2613912a62f80931"] = 213
        lut["d02cda160d162001fbd515a609837320058572c6"] = 214
        lut["876a7779d8418e1a0cb9c9a4312cb65f68fd6d0e"] = 215
        lut["fbc4af8337bd8e68659ffa9bd016ce3d9aa1892b"] = 216
        lut["2948ffacdc179dba578fcf162234a61f3ef5e19c"] = 217
        lut["010344d1afc116b92a3e5a5c3884a5c3d4f1997d"] = 218
        lut["7904b0a9c14194b6209f02161be15fe6d1416926"] = 219
        lut["a57657074cb1059b4529854e2b64db750b78ec93"] = 220
        lut["e55530adfd7b4f7bb9715bff81870f32c528c559"] = 221
        lut["e8ed81aa0c81a4d54b43d8a94fcd75133585a32c"] = 222
        lut["901e6113448d3de51fff452c67d4b9f93f700a9c"] = 223
        lut["258f9f498804b3bcfb4fea679df5aa528b71021b"] = 224
        lut["71305d279ae2f2b2915fe33ceff59bd1038f6df5"] = 225
        lut["53d91396bee5b7f8c4026343e4145b1b54269a47"] = 226
        lut["07c0e1a8453a53232ccdb09d183612dd33e28ecd"] = 227

        # sep.04.04 build
        lut["3b2797fbf2588e6815f5cef7fa33ee4d3d9ef7db"] = 100
        lut["37e4bf9213d992ce974caa936211099fbf63befc"] = 101
        lut["da37672307cb95fe51b381fad59c68c537f7dff4"] = 102
        lut["23f795926899e5bf1f32f0126418be4e91ac0f94"] = 103
        lut["3f1835f9a832cd5c31ed4eb8276a964dcbf7ce12"] = 104
        lut["95f2043d5eb69c2c5d987bccf8ffd6714cbfc333"] = 105
        lut["2cd93e0418b89329ddf94daf44bfe605884d2e7f"] = 106
        lut["51817ef40212588e80155d03e19fcdc586ef7781"] = 107
        lut["8d2ec79837e6cb684a54068218f225262cc392e3"] = 108
        lut["9a9802f68a3828a0003a302b553364fd46cc1797"] = 109
        lut["f2a728884a478c0826a8c52700a6fdc53cf645e7"] = 110
        lut["247b510cfc32ec50c2029fef7a9e2dc1e906f807"] = 111
        lut["9d37b56fec19fecc73e90c183c74a17de7e29bd1"] = 112
        lut["99e4a80b7acc14d0c903f45a379cd43d14ab471e"] = 113
        lut["c9cc0c2dce9f053a1fb00fac37540c27bbe287c3"] = 114
        lut["d1f7ec6a8a8ba2a249c04b8a05571bc352fbefd9"] = 115
        lut["9544bb6b8c664413b76d8045f0291caf74d30654"] = 116
        lut["d81ed931335a16a84c5b5335d533ea54018bcc75"] = 117
        lut["e5dc08252003d7e9b4efa7ee6491387b4e0912c0"] = 118
        lut["c4bfc2c546e588e3036333d485ede76d3619d407"] = 119
        lut["ac3e23a5bd22d6300e46e7d30671bac3593858d5"] = 120
        lut["687d908ad9821d90a5efb6a9661bce6620a21630"] = 121
        lut["16df8d98efe62be7e18a1a83e29e15dcd70390d9"] = 122
        lut["802b52209dda8c0d20a0889cb8dd4d1432bc08f2"] = 123
        lut["b054f3a6bd8ddf33c7361b73aef1357fbf0dc894"] = 124
        lut["bfdb8f8a90e3f7e1d8495a878a42cfb790bb7e88"] = 125
        lut["49584c404d0eb35f96d947b7498d90f8deb9138c"] = 126
        lut["386ef62256956ef1d42ab53fe5c08abdb6faadd2"] = 127

        with open(img, mode="rb") as f:
            f.seek(32768 + 800)
            buf = f.read(16384)

        return lut[hashlib.sha1(buf).hexdigest()]

    return None


def tdi_bin_check(cube: os.PathLike, histats: dict, cid=None):
    """This function only logs warnings and returns nothing."""

    if cid is None:
        try:
            cid = histats['PRODUCT_ID']
        except KeyError:
            cid = hirise.get_ChannelID_fromfile(cube)

    # TDI and binning check
    if float(histats["IMAGE_MEAN"]) >= 8000:
        logger.warning(
            f"{cid}: "
            f"Channel mean greater than 8000 (TDI or binning too high)."
        )
    elif float(histats["IMAGE_MEAN"]) < 2500:
        tdi = isis.getkey_k(cube, "Instrument", "Tdi")
        if tdi == "32" or tdi == "64":
            logger.warning(f"{cid}: TDI too low.")
    return


def lut_check(cube: os.PathLike, histats: dict, cid=None):
    # LUT check
    lut = int(isis.getkey_k(cube, "Instrument", "LookupTableNumber"))
    orbit_number = int(isis.getkey_k(cube, "Archive", "OrbitNumber"))
    threshhold = dict()
    if lut != -9998:
        if orbit_number > 65881:
            # After orbit 65881, RED1 RED2 RED3 moved to to DN 900-1000 offset

            threshhold["RED0"] = (
                (6814, 22),
                (5341, 23),
                (3869, 24),
                (3133, 25),
                (2397, 26),
                (1200, 27),
            )
            threshhold["RED1"] = (
                (6619, 1),
                (5116, 2),
                (3614, 3),
                (2863, 4),
                (2112, 5),
                (900, 6),
            )
            threshhold["RED2"] = threshhold["RED1"]
            threshhold["RED3"] = threshhold["RED1"]
            threshhold["RED4"] = (
                (6684, 8),
                (5191, 9),
                (3699, 10),
                (2953, 11),
                (2207, 12),
                (1000, 13)
            )
            threshhold["RED5"] = (
                (6749, 15),
                (5266, 16),
                (3784, 17),
                (3043, 18),
                (2302, 19),
                (1100, 20),
            )
            threshhold["RED6"] = threshhold["RED4"]
            threshhold["RED7"] = threshhold["RED4"]
            threshhold["RED8"] = threshhold["RED0"]
            threshhold["RED9"] = threshhold["RED4"]
            threshhold["IR10"] = threshhold["RED1"]
            threshhold["IR11"] = threshhold["RED4"]
            threshhold["BG12"] = threshhold["RED1"]
            threshhold["BG13"] = threshhold["RED1"]
        elif orbit_number > 13057:
            # After orbit 13057, IR10 moved to DN 900-1000 offset
            threshhold["RED0"] = (
                (6814, 22),
                (5341, 23),
                (3869, 24),
                (3133, 25),
                (2397, 26),
                (1200, 27),
            )
            threshhold["RED1"] = (
                (6684, 8),
                (5191, 9),
                (3699, 10),
                (2953, 11),
                (2207, 12),
                (1000, 13),
            )
            threshhold["RED2"] = threshhold["RED1"]
            threshhold["RED3"] = threshhold["RED1"]
            threshhold["RED4"] = threshhold["RED1"]
            threshhold["RED5"] = (
                (6749, 15),
                (5266, 16),
                (3784, 17),
                (3043, 18),
                (2302, 19),
                (1100, 20),
            )
            threshhold["RED6"] = threshhold["RED1"]
            threshhold["RED7"] = threshhold["RED1"]
            threshhold["RED8"] = threshhold["RED0"]
            threshhold["RED9"] = threshhold["RED1"]
            threshhold["IR10"] = (
                (6619, 1),
                (5116, 2),
                (3614, 3),
                (2863, 4),
                (2112, 5),
                (900, 6),
            )
            threshhold["IR11"] = threshhold["RED1"]
            threshhold["BG12"] = threshhold["IR10"]
            threshhold["BG13"] = threshhold["IR10"]
        elif orbit_number > 11710:
            # After orbit 11710, RED6 moved to DN 1000-1100 offset
            threshhold["RED0"] = (
                (6814, 22),
                (5341, 23),
                (3869, 24),
                (3133, 25),
                (2397, 26),
                (1200, 27),
            )
            threshhold["RED1"] = (
                (6684, 8),
                (5191, 9),
                (3699, 10),
                (2953, 11),
                (2207, 12),
                (1000, 13),
            )
            threshhold["RED2"] = threshhold["RED1"]
            threshhold["RED3"] = threshhold["RED1"]
            threshhold["RED4"] = threshhold["RED1"]
            threshhold["RED5"] = (
                (6749, 15),
                (5266, 16),
                (3784, 17),
                (3043, 18),
                (2302, 19),
                (1100, 20),
            )
            threshhold["RED6"] = threshhold["RED1"]
            threshhold["RED7"] = threshhold["RED1"]
            threshhold["RED8"] = threshhold["RED0"]
            threshhold["RED9"] = threshhold["RED1"]
            threshhold["IR10"] = threshhold["RED1"]
            threshhold["IR11"] = threshhold["RED1"]
            threshhold["BG12"] = (
                (6619, 1),
                (5116, 2),
                (3614, 3),
                (2863, 4),
                (2112, 5),
                (900, 6),
            )
            threshhold["BG13"] = threshhold["BG12"]
        elif orbit_number > 2660:
            # New LUT table after orbit 2660
            threshhold["RED0"] = (
                (6814, 22),
                (5341, 23),
                (3869, 24),
                (3133, 25),
                (2397, 26),
                (1200, 27),
            )
            threshhold["RED1"] = (
                (6684, 8),
                (5191, 9),
                (3699, 10),
                (2953, 11),
                (2207, 12),
                (1000, 13),
            )
            threshhold["RED2"] = threshhold["RED1"]
            threshhold["RED3"] = threshhold["RED1"]
            threshhold["RED4"] = threshhold["RED1"]
            threshhold["RED5"] = (
                (6749, 15),
                (5266, 16),
                (3784, 17),
                (3043, 18),
                (2302, 19),
                (1100, 20),
            )
            threshhold["RED6"] = (
                (6619, 1),
                (5116, 2),
                (3614, 3),
                (2863, 4),
                (2112, 5),
                (900, 6),
            )
            threshhold["RED7"] = threshhold["RED1"]
            threshhold["RED8"] = threshhold["RED0"]
            threshhold["RED9"] = threshhold["RED1"]
            threshhold["IR10"] = threshhold["RED1"]
            threshhold["IR11"] = threshhold["RED1"]
            threshhold["BG12"] = threshhold["RED6"]
            threshhold["BG13"] = threshhold["RED6"]
        else:
            # Original LUTs prior to 2661
            threshhold["RED0"] = (
                (14057, 21),
                (11676, 22),
                (9295, 23),
                (7152, 24),
                (5248, 25),
                (3343, 26),
                (1200, 27),
            )
            threshhold["RED1"] = (
                (13857 + 1, 14),
                (11476 + 1, 15),
                (9095 + 1, 16),
                (6952 + 1, 17),
                (5048 + 1, 18),
                (3143 + 1, 19),
                (1000, 20),
            )
            threshhold["RED2"] = threshhold["RED1"]
            threshhold["RED3"] = threshhold["RED1"]
            threshhold["RED4"] = threshhold["RED1"]
            threshhold["RED5"] = threshhold["RED1"]
            threshhold["RED6"] = (
                (13657 + 1, 7),
                (11276 + 1, 8),
                (8895 + 1, 9),
                (6752 + 1, 10),
                (4848 + 1, 11),
                (2943 + 1, 12),
                (800, 13),
            )
            threshhold["RED7"] = threshhold["RED1"]
            threshhold["RED8"] = threshhold["RED0"]
            threshhold["RED9"] = threshhold["RED1"]
            threshhold["IR10"] = threshhold["RED1"]
            threshhold["IR11"] = threshhold["RED1"]
            threshhold["BG12"] = threshhold["RED6"]
            threshhold["BG13"] = threshhold["RED6"]

        ccd = hirise.ChannelID(
            isis.getkey_k(cube, "Archive", "ProductId")
        ).get_ccd()
        if cid is None:
            cid = ccd

        for (th, ex) in threshhold[ccd]:
            if float(histats["IMAGE_MEAN"]) >= th:
                lut_diff = lut - ex
                if lut_diff >= 1 or lut_diff <= -1:
                    if lut_diff > 0:
                        direction = "to the right"
                    else:
                        direction = "to the left"
                    logger.warning(
                        f"{cid}: LUT is {lut_diff} column(s) ({direction}) "
                        f"from ideal settings - image overcompressed."
                    )
                break
        else:
            logger.warning(
                f"{cid}: DN value, {histats['IMAGE_MEAN']}, lower than lowest "
                f"DN value with defined LUT for this channel."
            )
    return


if __name__ == "__main__":
    main()
