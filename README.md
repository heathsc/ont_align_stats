# ont_align_stats

Utility for calculating mapping statistics for (mostly) ONT reads from indexed BAM files.

## Changes

 - 0.3.1 Correct double counting of reads that span multiple regions.  Prevent reads from non-mappable regions being 
omitted from the total counts
 - 0.3.0 Collection of coverage stats now working
 - 0.2.0 Add collection of basic stats, % mapped, split mapping stats
 - 0.1.0 First commit
