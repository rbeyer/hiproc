=========
Changelog
=========

All notable changes to this project will be documented in this file.

The format is based on `Keep a Changelog <https://keepachangelog.com/en/1.0.0/>`_,
and this project adheres to `Semantic Versioning <https://semver.org/spec/v2.0.0.html>`_.

When updating this file, please add an entry for your change under
Unreleased_ and one of the following headings:

- Added - for new features.
- Changed - for changes in existing functionality.
- Deprecated - for soon-to-be removed features.
- Removed - for now removed features.
- Fixed - for any bug fixes.
- Security - in case of vulnerabilities.

If the heading does not yet exist under Unreleased_, then add it
as a 3rd level heading, underlined with pluses (see examples below).

When preparing for a public release add a new 2nd level heading,
underlined with dashes under Unreleased_ with the version number
and the release date, in year-month-day format (see examples below).


Unreleased
----------

0.8.1 (2021-07-06)
------------------

Added
+++++
- Additional logging of the HiPrecisionInit output during a hiproc run, which was previously
  obscured.

Fixed
+++++
- In some legitimate cases, some elements returned from ISIS histat would be missing (e.g. in the
  case of no valid pixels in the CAL_BUFFER area, there is no returned value of "Average"), but
  were not handled properly (Issue #3).  More robust handling was added.


0.8.0 (2021-06-30)
------------------

Changed
+++++++
- bitflips.py: Removed the span check when considering end values, which may affect how
  the algorithm deals with "mostly good" DN histograms, such that it should now correctly
  handle mesa-shaped DN histograms and not chop them in the middle.

0.7.1 (2021-06-21)
------------------

Fixed
+++++
- The new-in-0.7.0 multiprocessing requires that the multidict package be installed
  (for pvl to use), and it wasn't included in the requirements.  If you have 0.7.0,
  a simple "conda install multidict" should get you working.


0.7.0 (2021-06-15)
------------------

Added
+++++
- Improved documentation
- environment.yml for conda development
- Added Python-based multiprocessing capabilities
  to EDR_Stats, HiCal, and hiproc.

Changed
+++++++
- Incorporated changes from upstream to HiCal to make bitflip and lisfix settings
  on by default.

Fixed
+++++
- resolve_jitter output header lines properly commented, now doesn't bomb HiJACK's run
  of hijitter.
- now, additionally supporting ISIS 4.1.1 through 5.0.0
- Upstream addressed a bug in analyze_cubenorm_stats() that we had noticed,
  so now fixed here, too.
- MANIFEST.in did not properly include the "data" directory, so that was a problem.

0.6.1 (2021-03-23)
------------------

Added
+++++
- lisfix will return a non-zero exit code if it chooses not to fix the input cube.

Changed
+++++++
- Updated documentation in various places.


0.6.0 (2021-03-18)
------------------

Added
+++++
- More complete documentation for the programs and their parameters.

Changed
+++++++
- Implemented better handling for configuration files so that they
  will get distributed with the package.


0.5.0 (2021-03-05)
------------------

Added
+++++
* lisfix: Added the lisfix module.
* bitflips: Sometimes the very end of the histogram (although not a formal minima) is the
  appropriate choice, and is now considered.
* bitflips: Added capability to ignore minor maxima at the ends of the histogram to "roll down"
  to a better solution.

Changed
+++++++
* name change of project from pyrise to hiproc.
* bitflips: Changed the default medstd_limit from 300 to 400 DN.

Fixed
+++++
* bitflips: There were a variety of edge cases that resulted in errors.  The appropriate guardrails, handlers,
  and recovery logic has now been added.
* HiStitch: The equalize and balance parameters cannot both be true.

0.4.0 (2020-09-22)
------------------
* Tremendous amount of re-working in bitflips to improve
  performance.
* Format cleanup

0.3.0 (2020-05-16)
------------------
* Confirmed that EDR_Stats and HiCal produce identical output cubes.

0.2.0 (2020-05-06)
------------------
* Updated with bit-flip correction.

0.1.0 2020-03-21
----------------
* First shared on GitHub
