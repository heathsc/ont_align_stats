use std::path::PathBuf;

use clap::{Arg, ArgAction, ArgGroup, Command, command, value_parser};

use super::LogLevel;

pub(super) fn cli_model() -> Command {
    command!()
        .arg(
            Arg::new("input")
                .value_parser(value_parser!(PathBuf))
                .value_name("INPUT")
                .required(true)
                .help("Input FASTQ/SAM/BAM/CRAM file(s)"),
        )
        .next_help_heading("Processing")
        .arg(
            Arg::new("n_tasks")
                .short('n')
                .long("n-tasks")
                .value_parser(value_parser!(u64).range(1..))
                .value_name("INT")
                .default_value("1")
                .help("No. parallel read tasks"),
        )
        .arg(
            Arg::new("hts_threads")
                .short('t')
                .long("hts-threads")
                .value_parser(value_parser!(u64).range(1..))
                .value_name("INT")
                .help("No. threads for hts reading [default: No. physical cores]"),
        )
        .arg(
            Arg::new("max_block_size")
                .long("max-block-size")
                .value_parser(value_parser!(u64).range(1000..))
                .value_name("INT")
                .default_value("25000000")
                .hide(true)
                .help("Maximum block size"),
        )
        .arg(
            Arg::new("bam_rec_thread_buffer")
                .long("bam-rec-thread-buffer")
                .value_parser(value_parser!(u64).range(128..))
                .value_name("INT")
                .default_value("1024")
                .hide(true)
                .help("Size of bam rec buffer"),
        )
        .next_help_heading("Genome")
        .arg(
            Arg::new("mappability")
                .long("mappability")
                .short('m')
                .value_parser(value_parser!(PathBuf))
                .value_name("FILE")
                .help("Mappability bed file"),
        )
        .arg(
            Arg::new("reference")
                .long("reference")
                .short('T')
                .value_parser(value_parser!(PathBuf))
                .value_name("FILE")
                .help("Reference FASTA file"),
        )
        .next_help_heading("Filtering")
        .arg(
            Arg::new("region")
                .long("region")
                .short('r')
                .value_parser(value_parser!(String))
                .value_name("REGION")
                .action(clap::ArgAction::Append)
                .help("Genomic region"),
        )        
        .arg(
            Arg::new("regions")
                .long("regions")
                .short('R')
                .value_parser(value_parser!(PathBuf))
                .value_name("FILE")
                .help("BED file with genomic regions"),
        )
        .arg(
            Arg::new("min_mapq")
                .long("min-mapq")
                .short('q')
                .value_parser(value_parser!(u8))
                .value_name("INT")
                .default_value("10")
                .help("Minimum mapq for mapping statistics"),
        )
        .arg(
            Arg::new("min_qual")
                .long("min-qual")
                .short('b')
                .value_parser(value_parser!(u8))
                .value_name("INT")
                .default_value("10")
                .help("Minimum base quality for coverage statistics"),
        )
        .arg(
            Arg::new("min_read_len")
                .long("min-read-len")
                .short('L')
                .value_parser(value_parser!(u64))
                .value_name("INT")
                .help("Minimum read length for coverage statistics"),
        )
        .arg(
            Arg::new("max_read_len")
                .long("max-read-len")
                .short('R')
                .value_parser(value_parser!(u64))
                .value_name("INT")
                .help("Maximum read length for coverage statistics"),
        )
        .next_help_heading("Output")
        .arg(
            Arg::new("prefix")
                .long("prefix")
                .short('p')
                .value_parser(value_parser!(String))
                .value_name("STRING")
                .default_value("ont_align_stats")
                .help("Prefix for output files"),
        ) 
        .arg(
            Arg::new("compress_gzip")
                .action(ArgAction::SetTrue)
                .long("compress-gzip")
                .short('z')
                .help("Compress output file with gzip"),
        )
        .arg(
            Arg::new("compress_bzip2")
                .action(ArgAction::SetTrue)
                .long("compress-bzip2")
                .short('j')
                .help("Compress output file with bzip2"),
        )
        .arg(
            Arg::new("compress_zstd")
                .action(ArgAction::SetTrue)
                .long("compress-zstd")
                .short('Z')
                .help("Compress output file with zstd"),
        )
        .arg(
            Arg::new("compress_xz")
                .action(ArgAction::SetTrue)
                .long("compress-xz")
                .short('x')
                .help("Compress output file with xz"),
        )
        .group(ArgGroup::new("compress")
            .args(["compress_gzip", "compress_bzip2", "compress_zstd", "compress_xz"])
            .multiple(false))
        .arg(
            Arg::new("dir")
                .long("dir")
                .short('d')
                .value_parser(value_parser!(PathBuf))
                .value_name("PATH")
                .help("Output directory [default: current directory]"),
        )
        .next_help_heading("Misc")
        .arg(
            Arg::new("timestamp")
                .short('X')
                .long("timestamp")
                .value_parser(value_parser!(stderrlog::Timestamp))
                .value_name("GRANULARITY")
                .default_value("none")
                .help("Prepend log entries with a timestamp"),
        )
        .arg(
            Arg::new("loglevel")
                .short('l')
                .long("loglevel")
                .value_name("LOGLEVEL")
                .value_parser(value_parser!(LogLevel))
                .ignore_case(true)
                .default_value("info")
                .help("Set log level"),
        )
        .arg(
            Arg::new("quiet")
                .action(ArgAction::SetTrue)
                .long("quiet")
                .conflicts_with("loglevel")
                .help("Silence all output"),
        )
}
