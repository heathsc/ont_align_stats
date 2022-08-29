# ont_align_stats

Utility for calculating mapping statistics for (mostly) ONT reads from indexed BAM files.

- [Introduction](#intro)
- [Installation](#install)
- [General usage](#usage)
    - [Multithreading options](#thread)
    - [Command line options](#cli)
- [Output](#output)
- [Changes](#changes)

## <a name="intro"></a>Introduction

## <a name="install"></a>Installation

To compile you will need an up-to-date copy of rust.  This can be
installed locally following the instructions [here](https://www.rust-lang.org/learn/get-started).  
Note that if you have rust already installed you should update it
using ``rustup update`` before trying to compile baldur.

You will also need to have [htslib](https://github.com/samtools/htslib) installed in a place
where the rust compiler can find it.  Note that *baldur* has been tested with versions of htslib from 1.12 to 1.15, but other recent versions should also work.

Clone the repository and then from the baldur directory
use cargo to compile the application:
```
    git clone https://github.com/heathsc/ont_align_stats.git
    cd ont_align_stats
    cargo build --release
```
If you have an error during linkage about not being able to find libhts then you will need to specify the installation location of libhts
during the build process:

    RUSTFLAGS="-L /opt/share/htslib/1.14/lib/"  cargo build --release

After successful the executable will be found in target/release/.  It
should be copied somewhere where it can be found by the shell. 
If libhts is not installed in a standard place then it may not be found 
when ont_align_stats is run; in this case the location of libhts should be added to the users LD_LIBRARY_PATH shell variable.

Once installed, basic help can be found by invoking ont_align_stats with
the -h flag.

## <a name="usage"></a>General usage

ont_align_stats is invoked with a SAM/BAM/CRAM input file containing aligned sequence reads, and will generate a JSON file
with alignment statistics.  While intended for non-paired long reads it should work on all 
types of sequence data, however no account is taken of paired reads (and in particular any overlap
between the two reads in a pair will be counted twice for the coverage statistics).

The input file should normally be sorted by genomic coordinates and indexed; ont_align_stats will work
without an index file but will generally require much more memory (~12 GB for a human genome) while memory
usage for an indexed input file will depend on the threading options chosen but will generally be <1 GB.  Note in
particular that if the input is sorted on genomic coordinates but an index is not available then
the --n-tasks options will not be effective.

By default, the entire genome is processed; this can be altered using the --region, --regions and --mappability options.
The --region argument allows a genome region to be specified on the command line in the standard way
used with [samtools](https://github.com/samtools/htslib) i.e., chr2 will select all of chromosome 2,
chr6:40000000 or chr6:40,000,000 will select everything on chromosome 6 from 40 MBases until the end 
of the chromosome (note the optional commas), and chr14:135,000,000-150,000,000 selects
the region on chromosome 14 between 135 and 150 Mbases.  Multiple --region options can be given to
specify multiple genomic regions, however if a large set of regions is required then the --regions option
allows the specification of a BED file containing the desired regions.  Note that the supplied BED file can
be compressed and will be decompressed transparently as long as a suitable tool is present in the users PATH.

The --mappability option also specifies a BED file with regions to consider, but there are
subtle differences between the --mappability and --regions options.  The --mappability regions list
regions that have previously been annotated as being uniquely mappable; this allows the coverage
statistics to be calculated just on the mappable fraction of the genome.  The mappability BED file can often
contain many regions (potentially millions), some of which are very small. The --regions option is
optimized for the case of large regions (on the order of megabases) whereas the --mappability option will handle
a large number of very small regions efficiently.  Supplying the mappability BED file 
as an argument to the --regions option will generally lead to ont_align_stats running 
very slowly.

If both the --region (or --regions) option and the --mappability option are provided then the 
intersect of the two sets of regions is considered for the coverage statistics.

### <a name="thread"></a>Multithreading options

Increased performance can be achieved using the multithreading options.  There are two such options:
--n-tasks (or -n) and --hts-threads (or -t).  
The --hts-threads option determines the number of additional threads for each hts file reading 
task; libhts is used to read SAM/BAM/CRAM files, and can use multiple threads to increase the
performance of the file reading/decoding.  The behaviour of the --n-tasks option depends on 
whether the input file is indexed or not.  The --hts-threads and --n-tasks options are set to 
0 and 1 respectively by default, so to obtain any
benefit from multi-threading these options should be used.

If the file is indexed and the requested number of tasks is >1
then the input file will be opened multiple times and read from in parallel.  The optimal 
combination of --n-tasks and --hts-threads will depend on the computer system and should be found
by doing some preliminary testing.  On a 8 (physical) core test laptop with SSDs, saturation of the CPUs (so 
16 virtual cores running at almost 100%) can be achieved with n-tasks set to the number of physical cores and 
hts-threads set to the number of virtual cores, but your mileage might vary.  Note that if n-tasks is set
to >1 the hts threads are shared across tasks.

If the file is not indexed the input file will only be opened once irrespective of 
the --n-tasks setting; however the genomic regions to be considered will be split between the number
of specified tasks so some parallelism can be achieved.  This is not as effective as with an indexed 
input file, however setting n-tasks to the number of physical cores and threads-per-reader to at least the same value works well with the test laptop 
described above.  Note that if the file is sorted in genomic order (but not indexed) then increasing the number of tasks will have little effect.
It would be much better in this case to generate the index before running ont_align_stats.

### <a name="cli"></a>Command line options

ont_align_stats has several command line options for controlling its operation.

| Short | Long           | Description                                                | Default           |
|-------|----------------|------------------------------------------------------------|-------------------|
 | q     | min-maxq       | Minimum MAXQ for calculating mapping statistics            | 10                |
| b     | min-qual       | Minimum base quality for coverage statistics               | 10                |
| r     | region         | Genomic region                                             |                   |
| R     | regions        | BED file containing genomic regions                        |                   |
|||||
| d     | dir            | Set output directory                                       | Current directory |
| p     | prefix         | Set prefix for output files                                | ont_align_stats   |
| z     | compress-gzip  | Compress output file with gzip                             |                   |
| j     | compress-bzip2 | Compress output file with bzip2                            |                   |
| x     | compress-xz    | Compress output file with xz                               |                   |
|||||
| m     | mappability    | Mappability BED file                                       | 1                 |
| T     | reference      | Reference FASTA file (for CRAM files)                      | 1                 |
|||||
| n     | n-tasks        | Number of parallel read tasks                              | 1                 |
| t     | hts-threads    | Number of additional threads for hts readers               | 1                 |
|||||
| l     | log-level      | Set log level [none, error, warn, info, debug, trace]      | info              |
|       | timestamp      | Prepend log entries with timestamp [none, sec, ms, us, ns] | none              |
|       | quiet          | Suppress all logging                                       |                   |
| h     | help           | Print help information                                     |                   |
| V     | version        | Print version information                                  |                   |

## <a name="changes"></a>Output
The main output file from ont_align_stats is a JSON file containing at the root level a map with (currently) two elements, metadata and stats.
The metadata element is a map with a number of key:value pairs giving information about the command, the machine used, option selected etc. 
The stats element is also a map with the following elements: counts, mapped_pctg, primary_mapq, read_len, n_splits, and coverage.
These are described below:

 - counts: set of raw counts of either mappings, reads or bases.
 - mapped_pctg: distribution of the number of reads where the given % of the read is mapped (i.e., ignoring clipped regions)
 - primary_mapq: distribution of the number of primary mappings with the given MAPQ value.
 - read_len: distribution of the number of reads with a given read length.
 - n_splits: distribution of the number of reads with a given number of split mappings (MAPQ > 0).
 - coverage: map with the number of mappable bases with a given coverage as well as the % of mappable bases with at least the given coverage.

## <a name="changes"></a>Changes

 - 0.9.0 Switch prefix short option to 'p' from 'P'.  Add output compression options
 - 0.8.0 Switch to using thread pools for hts readers for indexed reading when the file is opened multiple times.  In this way the threads are shared across tasks so the total number of threads is reduced.
 - 0.7.0 Change output so there is a metadata section and the stats section
 - 0.6.1 Optimizations, particularly for the non-indexed case
 - 0.6.0 Change coverage output in JSON file so that we report the proportion of genome covered with at least the given coverage level as well as the raw number of bases
 - 0.5.3 Resolve some problems with read selection in rare cases
 - 0.5.2 Fix region selection bug
 - 0.5.1 Clean up code.  Correct errors in CLI interface.
 - 0.5.0 Output in JSON format.  File name generated using the dir and prefix options
 - 0.4.0 Non-indexed operation works (but is slower and uses much more memory than indexed operation)
 - 0.3.5 Refactor in preparation for implementing non-indexed operation
 - 0.3.4 Use Coverage struct for temporary coverage statistics using almost 4 times less memory than the previous Option< usize > solution
 - 0.3.3 Update to using Rust 1.63.  Switch from crossbeam to std::thread for scoped threads
 - 0.3.2 Include unmapped reads in statistics.  Move to r_htslib 0.9.3. which allows recovering unmapped reads.
 - 0.3.1 Correct double counting of reads that span multiple regions.  Prevent reads from non-mappable regions being 
omitted from the total counts
 - 0.3.0 Collection of coverage stats now working
 - 0.2.0 Add collection of basic stats, % mapped, split mapping stats
 - 0.1.0 First commit
