# ont_align_stats

Utility for calculating mapping statistics for (mostly) ONT reads from indexed BAM files.

## Changes

 - 0.3.2 Include unmapped reads in statistics.  Move to r_htslib 0.9.3. which allows recovering unmapped reads.
 - 0.3.1 Correct double counting of reads that span multiple regions.  Prevent reads from non-mappable regions being 
omitted from the total counts
 - 0.3.0 Collection of coverage stats now working
 - 0.2.0 Add collection of basic stats, % mapped, split mapping stats
 - 0.1.0 First commit
