use std::{path::PathBuf, time::Instant};

use m_htslib::region::RegionList;

mod getters;
mod make_config;

#[derive(Copy, Clone, Debug, Default)]
pub enum CompressOpt {
    #[default]
    None,
    Gzip,
    Bzip2,
    Xz,
    Zstd,
}

pub struct Config {
    input: PathBuf,
    region_list: RegionList,

    // Output Options
    prefix: String,
    dir: Option<PathBuf>,
    compress: CompressOpt,

    // Read Select options
    // Minimum MAPQ for a mapping to be considered for the coverage stats
    min_mapq: u8,
    // Minimum base quality taken into account for the coverage stats
    min_qual: u8,
    // Minimum read length taken into account for the coverage stats
    min_read_len: Option<usize>,
    // Maximum read length taken into account for the coverage stats
    max_read_len: Option<usize>,

    // Operation options
    n_tasks: usize,
    hts_threads: usize,
    bam_rec_thread_buffer: usize,

    // Reference path
    reference: Option<PathBuf>,

    // Starting instant
    start_time: Instant,
}
