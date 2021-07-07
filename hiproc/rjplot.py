#!/usr/bin/env python
"""Creates interactive matplotlib plots from the output of
resolve_jitter, namely the <ObsID>_jitter_plot_[cpp|py].txt files and smear.txt
files.  It can produce interactive matplotlib plots of the jitter or smear,
and it can difference two sets of "jitter_plot.txt" files.
"""

import argparse
import re

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import hiproc.hirise as hirise


def read_data(path: Path) -> dict:
    lines = path.read_text().splitlines()
    labels = lines[0].lstrip('# ').split()
    lists = [[] for i in range(len(labels))]
    for line in lines[1:]:
        tokens = line.split()
        for i, token in enumerate(tokens):
            if token != "nan":
                lists[i].append(float(token))

    return dict(zip(labels, lists))


def read_smear_data(path: Path):
    s = list()
    el = list()
    e = list()
    with open(path) as f:
        for line in f:
            if line.startswith("#") or line.isspace():
                continue
            else:
                tokens = line.strip().split()
                s.append(float(tokens[0]))
                el.append(float(tokens[1]))
                e.append(float(tokens[2]))
    return s, el, e


def get_titles(path: Path) -> list:
    plt_p = path.with_suffix(".plt")

    text = plt_p.read_text()

    t = list()
    for i in (1, 2, 3):
        # m = re.search(fr"filePath{i}\s+=\s+'.*/?(\S+)\.flat\.tab'", text)
        m = re.search(fr"filePath{i}\s+=\s+'(\S+)'", text)
        if m is None:
            raise ValueError(f"filePath{i} could not be found in {plt_p}")
        else:
            # t.append(m.group(1))
            p = Path(m.group(1))
            # The expectation is that this will end in .flat.tab, so x2
            t.append(p.with_suffix("").with_suffix("").stem)

    return t


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("file", type=Path, nargs="+")

    args = parser.parse_args()

    plt.ioff()
    fig = plt.figure(constrained_layout=True)

    oid = hirise.get_ObsID_fromfile(args.file[0])

    if "smear" in args.file[0].name:
        data = [read_smear_data(args.file[0])]

        if len(args.file) == 2:
            data.append(read_smear_data(args.file[1]))
            gs = fig.add_gridspec(3, 1)
        else:
            gs = fig.add_gridspec(2, 1)

        axtop = fig.add_subplot(gs[0, 0])
        axtop.set_title(f"{oid} Smear")
        axtop.set_xlabel("ET")
        axtop.set_ylabel("Sample")

        axbot = fig.add_subplot(gs[1, 0])
        axbot.set_xlabel("ET")
        axbot.set_ylabel("Line")

        for d, color, m_size in zip(data, ("red", "blue"), (None, 2.5)):
            samp, line, et = d
            axtop.plot(et, samp, "o", c=color, ms=m_size)
            axbot.plot(et, line, "o", c=color, ms=m_size)

        if len(args.file) == 2:
            axdiff = fig.add_subplot(gs[2, 0])
            axdiff.set_ylabel("Diff")

            samp1, line1, et1 = data[0]
            samp2, line2, et2 = data[1]

            l2interp = np.interp(et1, et2, line2)
            s2interp = np.interp(et1, et2, samp2)

            linediff = line1 - l2interp
            sampdiff = samp1 - s2interp

            ldmax = np.max(np.abs(linediff))
            sdmax = np.max(np.abs(sampdiff))

            axdiff.plot(
                et1, linediff, c="purple", label=f"Line Diff Max: {ldmax}"
            )
            axdiff.plot(
                et1, sampdiff, c="orange", label=f"Sample Diff Max: {sdmax}"
            )
            axdiff.legend()

    else:  # regular jitter files

        titles = get_titles(args.file[0])
        data = [read_data(args.file[0])]

        if len(args.file) == 2:
            t2 = get_titles(args.file[1])
            if titles != t2:
                raise ValueError(
                    f"The titles in file1 ({titles}) do not match the titles in "
                    f"file 2 ({t2})"
                )
            data.append(read_data(args.file[1]))

        gs = fig.add_gridspec(3, 6)

        ax00 = fig.add_subplot(gs[0, 0:2])
        ax00.set_title(titles[0])

        ax01 = fig.add_subplot(gs[0, 2:4])
        ax01.set_title(titles[1])

        ax02 = fig.add_subplot(gs[0, 4:])
        ax02.set_title(titles[2])

        ax10 = fig.add_subplot(gs[1, 0:2])
        ax11 = fig.add_subplot(gs[1, 2:4])
        ax12 = fig.add_subplot(gs[1, 4:])

        ax20 = fig.add_subplot(gs[2, 0:3])
        ax20.set_title("Cross-Track Jitter")

        ax21 = fig.add_subplot(gs[2, 3:])
        ax21.set_title("Down-Track Jitter")

        ax00.set_ylabel('X, Cross-Track')

        for ax, time, off, interp, jittercheck in zip(
            (ax00, ax01, ax02, ax11, ax12),
            ("t1_shift", "t2_shift", "t3_shift", "t2_shift", "t3_shift"),
            ("offx1", "offx2", "offx3", "offy2", "offy3"),
            ("xinterp1", "xinterp2", "xinterp3", "yinterp2", "yinterp3"),
            (
                "jittercheckx1_shift",
                "jittercheckx2_shift",
                "jittercheckx3_shift",
                "jitterchecky2_shift",
                "jitterchecky3_shift"
            )
        ):
            for d, color, m_size in zip(data, ("red", "pink"), (None, 2.5)):
                ax.plot(d[time], d[off], "o", c=color, ms=m_size)

            for d, color, l_style, l_width in zip(
                data, ("green", "lime"), (None, "--"), (None, 0.5)
            ):
                ax.plot(d["ET_shift"], d[interp], c=color, ls=l_style, lw=l_width)

            for d, color in zip(data, ("yellow", "gold")):
                ax.plot(d["ET_shift"], d[jittercheck], c=color)

        ax10.set_ylabel('Y, Down-Track')
        for d, color, m_size, lab in zip(
            data, ("red", "pink"), (None, 2.5), ("Offset", "Offset 2")
        ):
            ax10.plot(d["t1_shift"], d["offy1"], "o", c=color, ms=m_size, label=lab)

        for d, color, l_style, l_width, lab in zip(
            data,
            ("green", "lime"),
            (None, "--"),
            (None, 0.5),
            ("Interp", "Interp 2")
        ):
            ax10.plot(
                d["ET_shift"],
                d["yinterp1"],
                c=color,
                ls=l_style,
                lw=l_width,
                label=lab
            )

        for d, color, lab in zip(
            data, ("yellow", "gold"), ("jittercheck", "jittercheck 2")
        ):
            ax10.plot(d["ET_shift"], d["jitterchecky1_shift"], c=color, label=lab)
        ax10.legend()

        for ax, ls in zip((ax20, ax21), ("Sample", "Line")):
            for d, color, l_style in zip(data, ("blue", "cyan"), (None, "--")):
                ax.plot(d["ET_shift"], d["Sample"], c=color, ls=l_style)

    plt.show()
