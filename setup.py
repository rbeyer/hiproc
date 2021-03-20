#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

requirements = [
    "gdal>=3.0.1",
    "kalasiris>=1.8.0",
    "matplotlib>=3.2.1",
    "numpy>=1.18.1",
    "pvl>=1.0.1",
    "scipy>=1.4.1"
]

setup(
    entry_points={
        'console_scripts': [
            'EDR_Stats=hiproc.EDR_Stats:main',
            'HiBeautify=hiproc.HiBeautify:main',
            'HiCal=hiproc.HiCal:main',
            'HiColorInit=hiproc.HiColorInit:main',
            'HiColorNorm=hiproc.HiColorNorm:main',
            'HiJACK=hiproc.HiJACK:main',
            'HiJitReg=hiproc.HiJitReg:main',
            'HiNoProj=hiproc.HiNoProj:main',
            'HiPrecisionInit=hiproc.HiPrecisionInit:main',
            'HiSlither=hiproc.HiSlither:main',
            'HiStitch=hiproc.HiStitch:main',
            'HiccdStitch=hiproc.HiccdStitch:main',
            'JitPlot=hiproc.JitPlot:main',
            'SlitherStats=hiproc.SlitherStats:main',
            'bitflips=hiproc.bitflips:main',
            'hiproc=hiproc.hiproc:main',
            'lisfix=hiproc.lisfix:main',
            'mdr2cub=hiproc.mdr2cub:main',
            'resolve_jitter=hiproc.resolve_jitter:main'
        ],
    },
    install_requires=requirements,
    include_package_data=True,
    package_data={
        "hiproc": ["data/*"],
    },
    keywords='HiRISE',
    packages=find_packages(include=['hiproc', 'hiproc.*']),
    test_suite='tests',
    zip_safe=False,
)
