use std::path::{Path, PathBuf};

#[derive(Default)]
pub struct Config {
    // Output Options
    prefix: Option<String>,
    dir: Option<PathBuf>,

    // Read Select options
    // Minimum MAPQ for a mapping to be considered for the coverage stats
    min_mapq: usize,
    // Minimum base quality taken into account for the coverage stats
    min_qual: usize,

    // Operation options
    n_tasks: usize,
    threads_per_task: usize,
    bam_rec_thread_buffer: usize,

    // Is the file indexed?
    indexed: bool,

    // Reference path
    reference: Option<PathBuf>,
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
    pub fn min_mapq(&self) -> usize {
        self.min_mapq
    }
    pub fn min_qual(&self) -> usize {
        self.min_qual
    }
    pub fn set_threads_per_task(&mut self, x: usize) {
        self.threads_per_task = x
    }
    pub fn threads_per_task(&self) -> usize {
        self.threads_per_task
    }
    pub fn set_n_tasks(&mut self, x: usize) {
        self.n_tasks = x
    }
    pub fn set_bam_rec_thread_buffer(&mut self, x: usize) {
        self.bam_rec_thread_buffer = x
    }
    pub fn n_tasks(&self) -> usize {
        self.n_tasks
    }
    pub fn bam_rec_thread_buffer(&self) -> usize {
        self.bam_rec_thread_buffer
    }
    pub fn set_indexed(&mut self, x: bool) {
        self.indexed = x
    }
    pub fn indexed(&self) -> bool {
        self.indexed
    }
}
