use std::{
    path::{Path, PathBuf},
    time::{Duration, Instant},
};

#[derive(Copy, Clone, Debug)]
pub enum CompressOpt {
    None,
    Gzip,
    Bzip2,
    Xz,
    Zstd,
}

impl Default for CompressOpt {
    fn default() -> Self {
        Self::None
    }
}

#[derive(Default)]
pub struct Config {
    // Output Options
    prefix: Option<String>,
    dir: Option<PathBuf>,
    compress: CompressOpt,

    // Read Select options
    // Minimum MAPQ for a mapping to be considered for the coverage stats
    min_mapq: usize,
    // Minimum base quality taken into account for the coverage stats
    min_qual: usize,
    // Minimum read length taken into account for the coverage stats
    min_read_len: Option<usize>,
    // Maximum read length taken into account for the coverage stats
    max_read_len: Option<usize>,

    // Operation options
    n_tasks: usize,
    hts_threads: usize,
    bam_rec_thread_buffer: usize,
    non_index_buffer_size: usize,

    // Is the file indexed?
    indexed: bool,

    // Bisulfite mode
    bisulfite: bool,

    // Reference path
    reference: Option<PathBuf>,

    // Starting instant
    start_time: Option<Instant>,
}

impl Config {
    pub fn set_prefix(&mut self, s: &str) {
        self.prefix = Some(s.to_owned())
    }
    pub fn prefix(&self) -> Option<&str> {
        self.prefix.as_deref()
    }
    pub fn set_dir<P: AsRef<Path>>(&mut self, s: P) {
        self.dir = Some(s.as_ref().to_owned())
    }
    pub fn dir(&self) -> Option<&Path> {
        self.dir.as_deref()
    }
    pub fn set_reference<P: AsRef<Path>>(&mut self, s: P) {
        self.reference = Some(s.as_ref().to_owned())
    }
    pub fn reference(&self) -> Option<&Path> {
        self.reference.as_deref()
    }
    pub fn set_min_mapq(&mut self, x: usize) {
        self.min_mapq = x
    }
    pub fn set_min_qual(&mut self, x: usize) {
        self.min_qual = x
    }
    pub fn set_min_read_len(&mut self, x: Option<usize>) {
        self.min_read_len = x
    }
    pub fn set_max_read_len(&mut self, x: Option<usize>) {
        self.max_read_len = x
    }
    pub fn min_mapq(&self) -> usize {
        self.min_mapq
    }
    pub fn min_qual(&self) -> usize {
        self.min_qual
    }
    pub fn min_read_len(&self) -> Option<usize> {
        self.min_read_len
    }
    pub fn max_read_len(&self) -> Option<usize> {
        self.max_read_len
    }
    pub fn set_hts_threads(&mut self, x: usize) {
        self.hts_threads = x
    }
    pub fn set_start_time(&mut self, ins: Instant) {
        self.start_time = Some(ins)
    }
    pub fn elapsed(&self) -> Option<Duration> {
        self.start_time.as_ref().map(|ins| ins.elapsed())
    }
    pub fn hts_threads(&self) -> usize {
        self.hts_threads
    }
    pub fn set_n_tasks(&mut self, x: usize) {
        self.n_tasks = x
    }
    pub fn set_bam_rec_thread_buffer(&mut self, x: usize) {
        self.bam_rec_thread_buffer = x
    }
    pub fn set_non_index_buffer_size(&mut self, x: usize) {
        self.non_index_buffer_size = x
    }
    pub fn set_compress(&mut self, x: CompressOpt) {
        self.compress = x
    }
    pub fn compress(&self) -> CompressOpt {
        self.compress
    }
    pub fn n_tasks(&self) -> usize {
        self.n_tasks
    }
    pub fn bam_rec_thread_buffer(&self) -> usize {
        self.bam_rec_thread_buffer
    }
    pub fn non_index_buffer_size(&self) -> usize {
        self.non_index_buffer_size
    }
    pub fn set_indexed(&mut self, x: bool) {
        self.indexed = x
    }
    pub fn indexed(&self) -> bool {
        self.indexed
    }
    pub fn set_bisulfite(&mut self, x: bool) {
        self.bisulfite = x
    }
    pub fn bisulfite(&self) -> bool {
        self.bisulfite
    }
}
