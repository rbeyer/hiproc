#!/usr/bin/env python
"""JitPlot converts the control network that HiJitReg creates and outputs
a simple text data file that can be plotted.
"""

# Copyright 2004-2020, Arizona Board of Regents on behalf of the Lunar and
# Planetary Laboratory at the University of Arizona.
#   - Orignal Java programs.
#
# Copyright 2020, Ross A. Beyer (rbeyer@seti.org)
#   - Elements of this Python program are are based on the original Java
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
# This program is based on HiColor version 5.4.2 (2020/02/14)
# and on the Java programs:
# - JitPlot.java ($Id: JitPlot.java,v 1.9 2020/02/14 22:46:48 $)
# - JitParser.java ($Id: JitParser.java,v 1.12 2020/02/14 22:46:48 $)
# by Guy McArthur and Michael Wendell as employees of the University of
# Arizona.

import argparse
import logging
import statistics
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import pvl
import hiproc.util as util
import hiproc.HiColorInit as hicolor
import hiproc.HiJitReg as hjr

logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, parents=[util.parent_parser()]
    )
    parser.add_argument(
        "-c",
        "--conf",
        required=False,
        default=Path(__file__).resolve().parent.parent
        / "data"
        / "HiJitReg.conf",
    )
    parser.add_argument(
        "cubes", metavar="balance.precolor.cub files", nargs="+"
    )

    args = parser.parse_args()

    util.set_logger(args.verbose, args.logfile, args.log)

    conf = pvl.load(str(args.conf))

    cubes = list(map(hicolor.HiColorCube, args.cubes))
    (red4, red5, ir10, ir11, bg12, bg13) = hicolor.separate_ccds(cubes)

    plt.ioff()
    fig, axes = plt.subplots(1, 4, sharex=True, sharey=True)

    axes[0].set_ylabel("Line")
    fig.text(0.5, 0.04, "Error Magnitude (pixels)", ha="center", va="center")

    for c, ax in zip([ir10, ir11, bg12, bg13], axes):
        if c is None:
            continue

        j = hjr.JitterCube(c, conf)
        j.reset()

        (accepted, rejected, smoothed) = filter_and_smooth(j)

        ax.set_title(j.get_ccd())
        size = 5
        if len(accepted) > 0:
            acceptA = np.array(accepted)
            ax.scatter(
                acceptA[:, 1],
                acceptA[:, 0],
                s=size,
                label="Accepted",
                facecolors="none",
                edgecolors="black",
            )
        if len(rejected) > 0:
            rejectedA = np.array(rejected)
            ax.scatter(
                rejectedA[:, 1],
                rejectedA[:, 0],
                c="red",
                s=size,
                label="Rejected",
            )
        if len(smoothed) > 0:
            smoothedA = np.array(smoothed)
            ax.scatter(
                smoothedA[:, 1],
                smoothedA[:, 0],
                c="blue",
                s=size,
                label="Smoothed",
            )
    fig.suptitle(str(j.get_obsid()))
    bottom, top = plt.ylim()
    plt.ylim(top, bottom)
    plt.xlim(0, 10)
    plt.legend(loc="upper right", bbox_to_anchor=(1, -0.05), ncol=3)
    plt.show()
    return


def filter_and_smooth(j: hjr.JitterCube) -> tuple:
    accepted = list()
    rejected = list()
    smoothed = list()

    for i, cm in enumerate(j.control_measures):
        start = 0
        end = int(i + j["BoxcarLength"] / 2)

        if i > j["BoxcarLength"] / 2:
            start = int(i - j["BoxcarLength"] / 2)

        error_median = statistics.median(
            map(lambda x: x["ErrorMagnitude"], j.control_measures[start:end])
        )

        point = (cm["Row"] * j["LineSpacing"], cm["ErrorMagnitude"])
        smoothed.append((point[0], error_median))

        if (
            abs(cm["ErrorMagnitude"] - error_median) > j["ExcludeLimit"]
            or cm["PointId"] in j.IgnoredPoints
        ):
            rejected.append(point)
            logger.info("Rejecting {}, {}".format(cm["Row"], cm["Column"]))
        else:
            accepted.append(point)

    return accepted, rejected, smoothed
