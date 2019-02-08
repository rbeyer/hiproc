The image processing pipelines that the HiRISE team operates
internally to produce higher order products are more than just the
'simple' programs available in ISIS3.  Those processes that run in
the HiRISE Operations Center (HiROC) are a complicated dance of
primarily Perl and ISIS3 run by a custom job management system, all
of which is based on the HiRISE catalog (HiCat) database.

This makes the HiROC system excellent for processing the Gigabytes of new
data that arrive daily from Mars, and allow the team to perform massive
reprocessing of the entire data set, as needed, and to produce on a large
scale a variety of derived data products.

However, that same complexity makes it difficult to reproduce exactly
what that system is doing on a small scale.

The programs in this repo are meant to replicate the HiRISE processing
chain on a small scale, so that individual algorithms and processes can
be investigated, without needed a massive data processing system.

The programs here have similar names to HiRISE pipelines, but only focus
on the data processing.  The HiRISE pipeline programs do a lot of other 
tasks relevant to being part of a massive ground data system, and clearly,
those functionalities aren't replicated here.

The HiROC system begins by watching the MRO project's raw data server for
new products with the FEI_Watchdog program, and then the HiDog pipeline 
fetches those products down to HiROC and the the EDRgen pipeline converts 
them into .img EDR products.

Since that is the most basic form of the data available from the PDS, we 
will start there, and assume that you have downloaded a set of EDR .img
files from the PDS.

We will use simple sqlite files instead of a full-up relational database.

# EDR_Stats

* Convert a HiRISE EDR Product to an ISIS cube for subsequent pipeline processing using the ISIS system for much of the work.
* Gather image statistics about the observation and update HiCat's EDR_Products table with those statistics
* Create an image histogram and store it as an ASCII file for use by other applications.
* Continue the pipeline processing by adding table entries to HiCat's HiCal and HiPdsVal sources tables.
* Calculate the Signal-to-Noise Ratio (SNR).
