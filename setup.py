#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = []

setup(
    author="Ross A. Beyer",
    author_email='rbeyer@seti.org',
    python_requires='>=3.5',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: BSD License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    description="A library to help process HiRISE EDRs with ISIS.",
    entry_points={
        'console_scripts': [
            'EDR_Stats=pyrise.EDR_Stats:main',
            'HiBeautify=pyrise.HiBeautify:main',
            'HiCal=pyrise.HiCal:main',
            'HiColorInit=pyrise.HiColorInit:main',
            'HiColorNorm=pyrise.HiColorNorm:main',
            'HiJACK=pyrise.HiJACK:main',
            'HiJitReg=pyrise.HiJitReg:main',
            'HiNoProj=pyrise.HiNoProj:main',
            'HiPrecisionInit=pyrise.HiPrecisionInit:main',
            'HiSlither=pyrise.HiSlither:main',
            'HiStitch=pyrise.HiStitch:main',
            'HiccdStitch=pyrise.HiccdStitch:main',
            'JitPlot=pyrise.JitPlot:main',
            'SlitherStats=pyrise.SlitherStats:main',
            'bitflips=pyrise.bitflips:main',
            'hirise_proc=pyrise.hirise_proc:main',
            'mdr2cub=pyrise.mdr2cub:main',
        ],
    },
    install_requires=requirements,
    license="BSD license",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='pyrise',
    name='pyrise',
    packages=find_packages(include=['pyrise', 'pyrise.*']),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/rbeyer/pyrise',
    version='0.1.0',
    zip_safe=False,
)
