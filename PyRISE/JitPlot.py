#!/usr/bin/env python
"""JitPlot converts the control network that HiJitReg creates and outputs
a simple text data file that can be plotted.
"""

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


# This program is based on HiColor version 1.99 2017/10/10
# and on the Java programs
#   JitPlot: $Id: JitPlot.java,v 1.7 2013/01/02 18:31:38 guym Exp $
# and
#   JitParser: $Id: JitParser.java,v 1.9 2013/04/29 20:55:55 guym Exp $
# by Guy McArthur
# which are Copyright(C) 2004 Arizona Board of Regents, under the GNU GPL.
#
# Since that suite of software is under the GPL, none of it can be directly
# incorporated in this program, since I wish to distribute this software
# under the Apache 2 license.  Elements of this software (written in an
# entirely different language) are based on that software but rewritten
# from scratch to emulate functionality.

import argparse
import logging
import statistics
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import pvl
import pyrise.util as util
import pyrise.HiColorInit as hicolor
import pyrise.HiJitReg as hjr


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     parents=[util.parent_parser()])
    parser.add_argument('-c', '--conf',    required=False,
                        default=Path(__file__).resolve().parent.parent /
                        'data' / 'HiJitReg.conf')
    parser.add_argument('cubes', metavar="balance.precolor.cub files",
                        nargs='+')

    args = parser.parse_args()

    util.set_logging(args.log)

    conf = pvl.load(str(args.conf))

    cubes = list(map(hicolor.HiColorCube, args.cubes))
    (red4, red5, ir10, ir11, bg12, bg13) = hicolor.separate_ccds(cubes)

    plt.ioff()
    fig, axes = plt.subplots(1, 4, sharex=True, sharey=True)

    axes[0].set_ylabel('Line')
    fig.text(0.5, 0.04, 'Error Magnitude (pixels)', ha='center', va='center')

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
            ax.scatter(acceptA[:, 1], acceptA[:, 0], s=size,
                       label='Accepted', facecolors='none', edgecolors='black')
        if len(rejected) > 0:
            rejectedA = np.array(rejected)
            ax.scatter(rejectedA[:, 1], rejectedA[:, 0], c='red', s=size,
                       label='Rejected')
        if len(smoothed) > 0:
            smoothedA = np.array(smoothed)
            ax.scatter(smoothedA[:, 1], smoothedA[:, 0], c='blue', s=size,
                       label='Smoothed')
    fig.suptitle(str(j.get_obsid()))
    bottom, top = plt.ylim()
    plt.ylim(top, bottom)
    plt.xlim(0, 10)
    plt.legend(loc='upper right', bbox_to_anchor=(1, -0.05), ncol=3)
    plt.show()
    return


def filter_and_smooth(j: hjr.JitterCube) -> tuple:
    accepted = list()
    rejected = list()
    smoothed = list()

    for i, cm in enumerate(j.control_measures):
        start = 0
        end = int(i + j['BoxcarLength'] / 2)

        if i > j['BoxcarLength'] / 2:
            start = int(i - j['BoxcarLength'] / 2)

        error_median = statistics.median(map(lambda x: x['ErrorMagnitude'],
                                             j.control_measures[start:end]))

        point = (cm['Row'] * j['LineSpacing'], cm['ErrorMagnitude'])
        smoothed.append((point[0], error_median))

        if(abs(cm['ErrorMagnitude'] - error_median) > j['ExcludeLimit'] or
           cm['PointId'] in j.IgnoredPoints):
            rejected.append(point)
            logging.info('Rejecting {}, {}'.format(cm['Row'], cm['Column']))
        else:
            accepted.append(point)

    return(accepted, rejected, smoothed)
