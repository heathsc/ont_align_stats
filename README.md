# ont_align_stats

Utility for calculating mapping statistics for (mostly) ONT reads from indexed BAM files.

- [Introduction](#intro)
- [Installation](#install)
- [General usage](#usage)
    - [Command line options](#cli)
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

Once installed, basic help can be found by invoking ont_align_stats with
the -h flag.

## <a name="usage"></a>General usage

### <a name="cli"></a>Command line options

ont_align_stats has several command line options for controlling its operation.

| Short | Long               | Description                                                | Default           |
|-------|--------------------|------------------------------------------------------------|-------------------|
 | q     | min-maxq           | Minimum MAXQ for calculating mapping statistics            | 10                |
| b     | min-qual           | Minimum base quality for coverage statistics               | 10                |
| r     | region             | Genomic region                                             |                   |
| R     | regions            | BED file containing genomic regions                        |                   |
|||||
| d     | dir                | Set output directory                                       | Current directory |
| P     | prefix             | Set prefix for output files                                | ont_align_stats   |
|||||
| m     | mappability        | Mappability BED file                                       | 1                 |
| T     | reference          | Reference FASTA file (for CRAM files)                      | 1                 |
|||||
| n     | n-tasks            | Number of parallel read tasks                              | 1                 |
| t     | threads-per-reader | Number of thread per SAM/BAM/CRAM reader                   | 1                 |
|||||
| l     | log-level          | Set log level [none, error, warn, info, debug, trace]      | info              |
|       | timestamp          | Prepend log entries with timestamp [none, sec, ms, us, ns] | none              |
|       | quiet              | Suppress all logging                                       |                   |
| h     | help               | Print help information                                     |                   |
| V     | version            | Print version information                                  |                   |

## <a name="changes"></a>Changes

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
