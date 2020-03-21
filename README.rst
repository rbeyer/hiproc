======
pyrise
======


..  image:: https://img.shields.io/pypi/v/pyrise.svg
        :target: https://pypi.python.org/pypi/pyrise

..  image:: https://img.shields.io/travis/rbeyer/pyrise.svg
        :target: https://travis-ci.org/rbeyer/pyrise

..  image:: https://readthedocs.org/projects/pyrise/badge/?version=latest
        :target: https://pyrise.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status




A library to help process HiRISE EDRs with ISIS.


* Free software: Apache 2 License
.. * Documentation: https://pyrise.readthedocs.io.
.. * `PlanetaryPy`_ Affiliate Package (someday).


Features
--------

* TODO


External Dependencies
---------------------
These programs use as much 'vanilla' Python 3 as possible.

However, it does depend on the following:

- pvl library (https://pvl.readthedocs.io)
- kalasiris library (https://kalasiris.readthedocs.io)
- numpy
- scipy
- matplotlib

The HiJACK program also requires an external ``resolveJitter``
program that has not been publicly released, but isn't that far
off.  There is a MATLAB version that has the appropriate licensing,
and there is a C++ version.  The C++ version could also be made to
have the appropriate licensing, it just hasn't gone through a release
process.  Maybe I'll write it in Python, too, and distribute it here.

Warning !
---------

The algorithms based on the HiRISE Processing Pipelines were emulated
and tested locally, but the results of each pipeline **have not**
been tested directly against the results of the HiRISE Processing
Pipelines, and this warning will remain until I have done so.  As
a result, I would not particularly 'trust' anything produced by
these programs at this time, and consider these algorithms a
work-in-progress.


Details
-------
The image processing pipelines that the HiRISE team operates
internally to produce higher order products are more than just the
'simple' programs available in ISIS.  Those processes that run in
the HiRISE Operations Center (HiROC) are a complicated dance of
primarily Perl and ISIS run by a custom job management system, all
of which interacts with the HiRISE catalog (HiCat) database.

This makes the HiROC system excellent for processing the Gigabytes
of new data that arrive daily from Mars, and allows the team to
perform massive reprocessing of the entire data set, as needed, and
to produce on a large scale a variety of derived data products.

However, that same complexity makes it difficult to reproduce exactly
what that system is doing on a small scale.

The programs here are meant to replicate the HiRISE processing chain
on a local scale, so that individual algorithms and processes can
be investigated, without needing a massive data processing system and
lots of infrastructure.

The programs here have similar names to HiRISE pipelines (hence the
perhaps strange intercapped naming conventions), but only focus on
the data processing.  The HiRISE pipeline programs do a lot of other
tasks relevant to being part of a massive ground data system, and
clearly, those functionalities aren't replicated here.

The HiROC system begins by watching the MRO project's raw data server for
new products with the ``FEI_Watchdog`` program, and then the HiDog pipeline
fetches those products down to HiROC and the the EDRgen pipeline converts
them into ``.img`` EDR products.

Since that is the most basic form of the data available from the PDS, we
will start there, and assume that you have downloaded a set of EDR ``.img``
files from the PDS.

As a final note, this library currently uses ``.json`` files to manage
passing information between programs, instead of a relational database system.


Contributing
------------

Feedback, issues, and contributions are always gratefully welcomed. See the
contributing guide for details on how to help and setup a development
environment.


.. _PlanetaryPy: https://github.com/planetarypy
