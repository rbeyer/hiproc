======
hiproc
======


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

The HiJACK program also requires the ``resolveJitter`` program that
is still buggy, and is not working reliably.

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

- HiCal: Verified. (if ``HiGainFx()`` enabled)
    Upstream is undergoing change.  HiGainFx really shouldn't be
    applied, so it is commented out here.  We're also working
    to integrate the bitflip cleaning into the upstream Perl,
    so this is in flux.


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
