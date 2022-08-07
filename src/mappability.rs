use std::{collections::HashMap, io::BufRead, path::Path};

use anyhow::Context;
use compress_io::compress;

#[derive(Default)]
pub struct Mappability {
    ctg_regions: HashMap<Box<str>, Vec<[usize; 2]>>,
}

impl Mappability {
    fn add_region<S: AsRef<str>>(&mut self, ctg: S, start: usize, stop: usize) {
        assert!(stop >= start, "Invalid mappability region");
        let ctg = ctg.as_ref();
        if !self.ctg_regions.contains_key(ctg) {
            self.ctg_regions
                .insert(ctg.to_owned().into_boxed_str(), Vec::new());
        }
        let v = self.ctg_regions.get_mut(ctg).unwrap();
        v.push([start, stop]);
    }

    fn sort(mut self) -> anyhow::Result<Self> {
        for v in self.ctg_regions.values_mut() {
            v.sort_unstable();
            for w in v.windows(2) {
                if w[0][1] >= w[1][0] {
                    return Err(anyhow!("Mappability regions overlap"));
                }
            }
        }
        Ok(self)
    }

    pub fn from_file<P: AsRef<Path>>(name: P) -> anyhow::Result<Self> {
        let fname = name.as_ref();
        debug!(
            "Reading in mappabiity regions from file {}",
            fname.display()
        );
        let mut buf = String::new();
        let mut nreg = 0;
        let mut f = compress::CompressIo::new()
            .path(fname)
            .bufreader()
            .with_context(|| "Could not open mappability input file")?;
        let mut mreg = Self::default();
        loop {
            buf.clear();

            if f.read_line(&mut buf)? == 0 {
                break;
            }
            let fd: Vec<_> = buf.trim().split('\t').collect();

            if fd.len() >= 3 {
                match (fd[1].parse::<usize>(), fd[2].parse::<usize>()) {
                    (Ok(start), Ok(stop)) if start <= stop => {
                        mreg.add_region(fd[0], start, stop);
                        nreg += 1;
                    }
                    _ => return Err(anyhow!("Invalid line from mappability file: {}", buf)),
                }
            }
        }
        debug!(
            "Mappability file read in: {} contigs and {} regions",
            mreg.ctg_regions.len(),
            nreg
        );
        mreg.sort()
    }

    pub fn ctg_regions<S: AsRef<str>>(&self, ctg: S) -> Option<&Vec<[usize; 2]>> {
        self.ctg_regions.get(ctg.as_ref())
    }
}
