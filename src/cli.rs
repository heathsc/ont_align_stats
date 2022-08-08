use std::{
    collections::HashMap,
    fmt::Display,
    num::NonZeroUsize,
    path::{Path, PathBuf},
    str::FromStr,
};

use crate::config::Config;
use crate::regions::Regions;
use anyhow::{Context, Error};
use clap::{crate_version, Arg, ArgMatches, Command};

use crate::input::open_input;
use crate::mappability::Mappability;

/// Set up stderrlog using options from clap::ArgMatches
///
/// Relevant options:
///
///   loglevel: can be error, warn, info, trace or none.  Default is info
///   quiet: turns off all loggin irrespective of loglevel
///   timestamp: sets timestamp options for logger
///
///   loglevel none is a synonym for the quiet option
///
fn init_log(m: &ArgMatches) {
    let verbose = m
        .value_of("loglevel")
        .map(|s| match s.to_lowercase().as_str() {
            "error" => Some(0),
            "warn" => Some(1),
            "info" => Some(2),
            "debug" => Some(3),
            "trace" => Some(4),
            "none" => None,
            _ => Some(0),
        })
        .unwrap_or(Some(2));

    let quiet = verbose.is_none() || m.is_present("quiet");
    let ts = m
        .value_of("timestamp")
        .map(|v| stderrlog::Timestamp::from_str(v).expect("Invalid value for timestamp"))
        .unwrap_or(stderrlog::Timestamp::Off);

    stderrlog::new()
        .quiet(quiet)
        .verbosity(verbose.unwrap_or(0))
        .timestamp(ts)
        .init()
        .unwrap();
}

fn parse<T>(s: &str) -> anyhow::Result<T>
where
    T: FromStr,
    T::Err: Display + Send + Sync + 'static,
{
    s.parse::<T>()
        .map_err(|e| Error::msg(format!("Parse error {}", e)))
}

fn cli_model() -> ArgMatches {
    Command::new("ont_align_stats").version(crate_version!()).author("Simon Heath")
        .about("ont_align_stats is a utility for extracting mapping stats for long ONT reads from SAM/BAM/CRAM")
        .arg(
            Arg::new("quiet")
                .long("quiet")
                .help("Silence all output")
                .conflicts_with("loglevel")
        )
        .arg(
            Arg::new("timestamp")
                .long("timestamp")
                .takes_value(true).value_name("GRANULARITY")
                .possible_values(&["none", "error", "warn", "info", "debug", "trace"])
                .ignore_case(true).default_value("none")
                .help("Prepend log entries with a timestamp")
        )
        .arg(
            Arg::new("loglevel")
                .short('l').long("loglevel")
                .takes_value(true).value_name("LOGLEVEL")
                .possible_values(&["none", "error", "warn", "info", "debug", "trace"])
                .ignore_case(true).default_value("info")
                .help("Set log level")
        )
        .arg(
            Arg::new("mappability")
                .short('m').long("mappability")
                .takes_value(true).value_name("PATH")
                .help("Path to mappability bed file")            
        )
        .arg(
            Arg::new("reference")
                .short('T').long("reference")
                .takes_value(true).value_name("PATH")
                .help("Path to reference fasta file")
        )
        .arg(
            Arg::new("region")
                .short('r').long("region")
                .takes_value(true).value_name("REGION").multiple_occurrences(true)
                .help("Genomic region")
        )
        .arg(
            Arg::new("regions")
                .short('R').long("regions")
                .takes_value(true).value_name("FILE")
                .conflicts_with("region")
                .help("BED file with genomic regions")
        )
        .arg(
            Arg::new("min_maxq")
                .short('q').long("min-pcnt")
                .takes_value(true).value_name("NUMBER")
                .default_value("10")
                .help("Minimum % samples passing filter for entry to be output")
        )
        .arg(
            Arg::new("prefix")
                .short('P').long("prefix")
                .takes_value(true).value_name("STRING")
                .default_value("bedmethyl")
                .help("Set prefix for output files")
        )
        .arg(
            Arg::new("dir")
                .short('d').long("dir")
                .takes_value(true).value_name("DIR")
                .help("Set output directory [default: current directory]")
        )
        .arg(
            Arg::new("n_tasks")
                .short('n').long("n-tasks")
                .takes_value(true).value_name("INT")
                .default_value("1")
                .help("No. of parallel read tasks")
        )
        .arg(
            Arg::new("threads_per_task")
                .short('t').long("threads-per-task")
                .takes_value(true).value_name("INT")
                .default_value("1")
                .help("No. of threads per read tasks")
        )
        .arg(
            Arg::new("max_block_size")
                .long("max-block-size").hidden(true)
                .takes_value(true).value_name("INT")
                .default_value("5000000")
        )
        .arg(
            Arg::new("max_gap")
                .long("max-gap").hidden(true)
                .takes_value(true).value_name("INT")
                .default_value("100000")
        )
        .arg(
            Arg::new("no_index")
                .long("no-index").hidden(true)
        )
        .arg(
            Arg::new("input")
                .takes_value(true).value_name("FILE")
                .required(true)
                .help("Input file (SAM/BAM/CRAM)")
        )
        .get_matches()
}

pub fn handle_cli() -> anyhow::Result<(Config, PathBuf, Regions)> {
    let m = cli_model();
    init_log(&m);

    debug!("Checking input file");

    // Handle Input file
    let input = m.value_of("input").unwrap();
    let no_index = m.is_present("no_index");
    let reference = m.value_of("reference").map(Path::new);

    let hts = open_input(input, no_index, reference, 1)?;
    let indexed = hts.has_index();

    // Get list of sequence names and lengths from file header
    let seq = hts.seq_names();
    let lengths = hts.seq_lengths();
    assert_eq!(seq.len(), lengths.len());

    debug!("Handling command line inputs");

    // Parse region options
    let mut regions = if let Some(f) = m.value_of("regions") {
        Some(Regions::from_file(f))
    } else {
        m.values_of("region").map(Regions::from_str_vec)
    }
    .unwrap_or_else(|| Regions::from_str_vec(&seq))?;

    // Threads argument
    let n_tasks = usize::from(
        parse::<NonZeroUsize>(
            m.value_of("n_tasks")
                .expect("Missing default n_tasks value"),
        )
        .with_context(|| "Error parsing threads option")?,
    );

    // If multithreading and file is indexed, split up regions into chunks of at most max_block_size
    if indexed {
        let max = parse::<usize>(
            m.value_of("max_block_size")
                .expect("Missing default max-block-size value"),
        )
        .with_context(|| "Error parsing max_block_size option")?;
        let max_block_size = if max >= 1000 {
            max
        } else {
            warn!("Requested max-block-size too low: setting to 1000");
            1000
        };
        let len_hash: HashMap<&str, usize> =
            seq.iter().copied().zip(lengths.iter().copied()).collect();
        regions.split_regions(max_block_size, len_hash);
    }

    // If mappability option set, add mappability data to regions
    let map = if let Some(s) = m.value_of("mappability") {
        let mp = Mappability::from_file(s)?;
        let max_gap = parse::<usize>(
            m.value_of("max_gap")
                .expect("Missing default max-gap value"),
        )
        .with_context(|| "Error parsing max_gap option")?;
        regions.add_mappability(&mp, max_gap);
        Some(mp)
    } else {
        None
    };

    let mut cfg = Config::default();
    cfg.set_indexed(indexed);
    cfg.set_n_tasks(n_tasks);
    cfg.set_mappability(map);
    cfg.set_min_mapq(m.value_of_t("min_maxq").unwrap());
    cfg.set_threads_per_task(m.value_of_t("threads_per_task").unwrap());
    if let Some(s) = reference {
        cfg.set_reference(s)
    }

    if let Some(s) = m.value_of("prefix") {
        cfg.set_prefix(s)
    }
    if let Some(s) = m.value_of("dir") {
        cfg.set_dir(s)
    }
    Ok((cfg, PathBuf::from(input), regions))
}