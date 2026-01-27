use clap::ArgMatches;
use super::LogLevel;

/// Initialize logging from command line arguments
pub fn init_log(m: &ArgMatches) {
    let verbose = m
        .get_one::<LogLevel>("loglevel")
        .copied()
        .expect("Missing default log level");
    let quiet = verbose.is_none() || m.get_flag("quiet");
    let ts = m
        .get_one::<stderrlog::Timestamp>("timestamp")
        .copied()
        .unwrap_or(stderrlog::Timestamp::Off);

    stderrlog::new()
        .quiet(quiet)
        .verbosity(verbose.get_level())
        .timestamp(ts)
        .init()
        .unwrap();
}