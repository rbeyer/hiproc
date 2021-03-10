==================
Processing Details
==================

The image processing pipelines that the HiRISE team operates
internally to produce higher order products are more than just the
'simple' programs available in ISIS.  Those processes that run in
the HiRISE Operations Center (HiROC) are a complicated dance of
primarily Perl and ISIS run by a custom job management system, all
of which interacts with the HiRISE database (HiCat).

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
perhaps strange CamelCase naming conventions), but only focus on
the data processing.  The HiRISE pipeline programs of the same name
do a lot of other tasks relevant to being part of a massive ground
data system, and clearly, those functionalities aren't replicated
here.

The HiROC system begins by watching the MRO project's raw data server for
new products with the ``FEI_Watchdog`` program, and then the HiDog pipeline
fetches those products down to HiROC and the the EDRgen pipeline converts
them into ``.img`` EDR products.

Since that is the most basic form of the data available from the PDS, we
will start there, and assume that you have downloaded a set of EDR ``.img``
files from the PDS.

As a final note, this library currently uses ``.json`` files to manage
passing information between programs, instead of a relational database system.

-------------
Pipeline Flow
-------------

Details of the individual programs in this library are detailed on their own pages,
this is an overview of the processing flow.

There are a variety of configuration files that are included with this library that
were designed to be used in the HiROC pipelines themselves, or are meant to be used
by ISIS programs that the pipelines run.  As such, they contain more information than
is needed by the programs in this library, but rather than redesign them, the programs
in this library are simply built to use them, and extract what they need.  This allows
these programs to use the *exact same* configuration files that are used by the HiROC
pipelines.

In general, the initial flow is always the same:

1. ``EDR_Stats``
2. ``HiCal``
    It should be noted that ``bitflips`` and ``lisfix`` are separate programs, and can
    be run separately, but their functionality is engaged by ``HiCal``, when appropriate.
3. ``HiStitch``
4. ``HiccdStitch``

Depending on what your goals are for processing, you could stop here.  However, if you
want to create color mosaics, or engage the "high precision" processing, a few more
preparatory steps are needed:

5. ``HiColorInit``
6. ``HiJitReg``
    If you want to visualize the output of HiJitReg, you can use the JitPlot program
    to do so.
