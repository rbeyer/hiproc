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
CamelCase naming conventions that aren't exactly Pythonic), but
only focus on the data processing.  The HiRISE pipeline programs
of the same name do a lot of other tasks relevant to being part of
a massive ground data system, and clearly, those functionalities
aren't replicated here.

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

At this point, you will have "stitched" and "balance" files for
each CCD.  Depending on what your goals are for processing, you
could stop here.  You can also get to this point via just using
``hiproc`` which runs the above four steps.

However, if you want to create color mosaics, or
engage the "high precision" processing, a few more steps are needed:

5. ``HiColorInit``
6. ``HiJitReg``
    If you want to visualize the output of HiJitReg, you can use the JitPlot program
    to do so.
7. ``HiSlither``
    If you want a summary of the HiSlither results, or to visualize them, you can use
    the SlitherStats program.

If you want to create composite color mosaics from HiRISE data, then you would use
these steps

* ``HiColorNorm``
* ``HiBeautify``

After running these steps, you will have IRB and RGB mosaics of the central color
HiRISE CCDs.  Alternately, you can just run ``hiproc -c`` which runs all of the above
steps from scratch, or just the additional steps needed, if you've already run plain
``hiproc``.

For "precision" processing, do the following:

* ``HiPrecisionInit`` to determine if you need to run HiJACK or just HiNoProj
* ``HiNoProj`` or ``HiJACK``

Alternately, you can just run ``hiproc -p`` (or ``hiproc -j`` if you want to force
HiJACK processing), which does the above from scratch, or just the additional steps
needed, if you've already run hiproc to produce basic or color products.
