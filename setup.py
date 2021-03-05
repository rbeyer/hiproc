#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('CHANGELOG.rst') as change_file:
    changelog = change_file.read()

requirements = []

setup(
    author="Ross A. Beyer",
    author_email='rbeyer@seti.org',
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: Apache Software License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    description="A library to help process HiRISE EDRs with ISIS.",
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
            'hiproc=hiproc.hirise_proc:main',
            'mdr2cub=hiproc.mdr2cub:main',
        ],
    },
    install_requires=requirements,
    long_description=readme + '\n\n' + changelog,
    include_package_data=True,
    keywords='hiproc',
    name='hiproc',
    packages=find_packages(include=['hiproc', 'hiproc.*']),
    test_suite='tests',
    url='https://github.com/rbeyer/hiproc',
    version='0.5.0',
    zip_safe=False,
)
