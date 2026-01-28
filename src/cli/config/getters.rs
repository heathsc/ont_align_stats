use std::{time::Duration, path::Path};

use m_htslib::region::RegionList;

use super::{Config, CompressOpt};

impl Config {
    pub fn input(&self) -> &Path {
        self.input.as_ref()
    }
    
    pub fn region_list(&self) -> &RegionList {
        &self.region_list
    }
    
    pub fn prefix(&self) -> &str {
        self.prefix.as_str()
    }

    pub fn dir(&self) -> Option<&Path> {
        self.dir.as_deref()
    }

    pub fn reference(&self) -> Option<&Path> {
        self.reference.as_deref()
    }

    pub fn min_mapq(&self) -> u8 {
        self.min_mapq
    }

    pub fn min_qual(&self) -> u8 {
        self.min_qual
    }

    pub fn min_read_len(&self) -> Option<usize> {
        self.min_read_len
    }

    pub fn max_read_len(&self) -> Option<usize> {
        self.max_read_len
    }

    pub fn elapsed(&self) -> Duration {
        self.start_time.elapsed()
    }

    pub fn hts_threads(&self) -> usize {
        self.hts_threads
    }

    pub fn compress(&self) -> CompressOpt {
        self.compress
    }

    pub fn n_tasks(&self) -> usize {
        self.n_tasks
    }
}
