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

# External Dependencies

These programs use as much 'vanilla' Python 3 as possible.

However, the use of the pvl library (https://pvl.readthedocs.io).


# EDR_Stats

* Convert a HiRISE EDR Product to an ISIS cube for subsequent pipeline processing using the ISIS system for much of the work.
* Gather image statistics about the observation and update HiCat's EDR_Products table with those statistics
* Create an image histogram and store it as an ASCII file for use by other applications.
* Calculate the Signal-to-Noise Ratio (SNR).

# HiCal

* Perform a radiometric calibration correction and conversion to "I/F" units. The calibration corrects for the variable gain and dark current of each CCD detector column.
* As part of the radiometric calibration, the current pipeline implementation applies a separate gain line-drift correction. This original line-drift correction is used in lieu of the correction found in the current ISIS "hical" program as the original correction algorithm has been demonstrated to work better.
* Perform furrow correction. The image columns at the channel join are checked for furrows. Pixels in the furrow region (at the channel join) whose DN values have gone above a threshold are set to the null pixel value as these pixels can not be calibrated. If a RED-filter image experiences furrowing, a rare event usually caused by an improperly commanded observation, then the furrowed pixels will be permanently set to null pixels for all of the HiRISE standard and extras products. For BG and IR-filter images that make up the color products, the furrowed pixels will be interpolated in the HiColorNorm pipeline step.
* If an observation is determined to have furrows then an entry is made in HiCat's Tags table to indicate the level of furrowing that has occurred. The comment field in the Tag table entry contains a number indicating the percent of the first furrow column that had pixel values above the threshold furrow value.
* Due to HiRISE instrument instability problems, the HiCal pipeline performs a noise reduction procedure to reduce the number of bad pixels in an image observation. The noise correction is applied when the standard deviation of the dark pixel or mask regions exceed a threshold, or if the number of LIS pixels exceeds a threshold.
* A high-pass filter "cubenorm" step is applied to the calibrated image. Due to camera instabilities, residual vertical stripping often exists in the imaging that is corrected by this empirical method. The average standard deviation of this change is calculated and stored in HiCat.
* The ISIS program "hidestripe" is applied to suppress horizontal stripping seen whenever an observation is acquired using mixed-binning commanding. The standard deviation of this change is calculated and started in HiCat.
* Browse and thumbnail images of the EDR channel products are created.
* Entries for the browse and thumbnail EDR products are added to HiCat's Extras_Products table.
* An entry is made to the HiStitch pipeline sources table when two channel products of the same CCD have successfully gone through the HiCal pipeline step.
