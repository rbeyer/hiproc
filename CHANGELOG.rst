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


0.5.0 (2021-03-05)
------------------

Added
+++++
* lisfix: Added the lisfix module.
* bitflips: Sometimes the very end of the histogram (although not a formal minima) is the appropriate choice, and
  is now considered.
* bitflips: Added capability to ignore minor maxima at the ends of the histogram to "roll down" to a better solution.

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
