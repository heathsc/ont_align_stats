use crate::mappability::Mappability;
use std::path::{Path, PathBuf};

#[derive(Default)]
pub struct Config {
    // Output Options
    prefix: Option<String>,
    dir: Option<PathBuf>,

    // Read Select options
    min_mapq: usize,

    // Operation options
    n_tasks: usize,
    threads_per_task: usize,

    // Is the file indexed?
    indexed: bool,

    // Mappability regions
    mappability: Option<Mappability>,

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
    pub fn min_maxq(&self) -> usize {
        self.min_mapq
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
    pub fn n_tasks(&self) -> usize {
        self.n_tasks
    }
    pub fn set_mappability(&mut self, map: Option<Mappability>) {
        self.mappability = map
    }
    pub fn mappability(&self) -> Option<&Mappability> {
        self.mappability.as_ref()
    }
    pub fn set_indexed(&mut self, x: bool) {
        self.indexed = x
    }
    pub fn indexed(&self) -> bool {
        self.indexed
    }
}
