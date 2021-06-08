======
hiproc
======

.. image:: https://readthedocs.org/projects/hiproc/badge/?version=latest
        :target: https://hiproc.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status

.. image:: https://img.shields.io/pypi/v/hiproc.svg
        :target: https://pypi.python.org/pypi/hiproc
        :alt: PyPI version

.. image:: https://anaconda.org/conda-forge/hiproc/badges/version.svg
        :target: https://anaconda.org/conda-forge/hiproc
        :alt: Conda version


A library to help process HiRISE EDRs with ISIS.


* Free software: Apache 2 License

.. * Documentation: https://hiproc.readthedocs.io.
.. * `PlanetaryPy`_ Affiliate Package (someday).


Features
--------

* TODO: Complete testing against Perl Pipelines.


External Dependencies
---------------------
These programs use as much 'vanilla' Python 3 as possible.

However, it does depend on the following:

- pvl library (https://pvl.readthedocs.io)
- kalasiris library (https://kalasiris.readthedocs.io)
- gdal
- numpy
- scipy
- matplotlib


Warning !
---------

The algorithms based on the HiRISE Processing Pipelines were emulated
and tested locally, but the results of each pipeline **have not**
been tested directly against the results of the HiRISE Processing
Pipelines, and this warning will remain until I have done so.  As
a result, I would not particularly 'trust' anything produced by
these programs at this time, and consider these algorithms a
work-in-progress.

These programs have been tested against their upstream Perl counterparts:

- EDR_Stats: Verified!
    Really just runs ``hi2isis`` so no surprise here.

- HiCal: Verified.
    Upstream is undergoing change, needs to be re-verified once upstream
    settles down.

- HiStitch: not verified
- HiccdStitch: not verified
- HiColorInit: not verified
- HiJitReg: not verified
- HiSlither: not verified
- HiColorNorm: not verified
- HiBeautify: not verified
- HiPrecisionInit: not verified
- HiNoProj: not verified
- HiJACK: not verified

Documentation
-------------
Full `documentation for hiproc is available <https://hiproc.readthedocs.io/en/latest/>`_,
including information on the processing flow of the various available programs, and
each program is self-documenting via their `-h` argument.

Due to the interaction with ISIS and GDAL, please read the installation instructions
carefully.

Contributing
------------

Feedback, issues, and contributions are always gratefully welcomed. See the
contributing guide for details on how to help and setup a development
environment.


Naming
------

The ISIS software has a number of processing or "proc" programs
(`mocproc`, `thmproc`, etc.) that are meant to be run to process
raw images to higher-level, more usable versions.  Naming this
library `hiproc` is an echo to that. There is a `hiproc` program
that is available after installation that provides a streamlined
one-stop-program, but this package provides a great deal more.


.. _PlanetaryPy: https://github.com/planetarypy
