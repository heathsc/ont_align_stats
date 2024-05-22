use std::collections::btree_map::{Iter, Values};
use std::{
    cmp::{Ord, Ordering},
    collections::{BTreeMap, HashMap},
    fmt,
    io::BufRead,
    rc::Rc,
};

use crate::mappability::Mappability;
use anyhow::Context;
use compress_io::compress;
use lazy_static::lazy_static;
use r_htslib::Sequence;
use regex::{Regex, RegexSet};

lazy_static! {
    // Matches when the contig is disambiguated using brackets i.e.., {chr2}:20000-50000
    // The Regex for the contig name comes from the VCF4.3 spec
    static ref RE_REGION1: Regex = Regex::new(r#"^[{]([0-9A-Za-z!#$%&+./:;?@^_|~-][0-9A-Za-z!#$%&*+./:;=?@^_|~-]*)[}]:?([0-9,]+)?-?([0-9,]+)?"#).unwrap();

    // Matches when the contig is present without brackets i.e., chr2:20000-30000
    static ref RE_REGION2: Regex = Regex::new(r#"^([0-9A-Za-z!#$%&+./:;?@^_|~-][0-9A-Za-z!#$%&*+./;=?@^_|~-]*):?([0-9,]+)?-?([0-9,]+)?"#).unwrap();

    // RegexSet to check for valid contig
    static ref RE_SET_CONTIG: RegexSet = RegexSet::new(&[
        r"^[{]([0-9A-Za-z!#$%&+./:;?@^_|~-][0-9A-Za-z!#$%&*+./:;=?@^_|~-]*)[}]$",
        r"^([0-9A-Za-z!#$%&+./:;?@^_|~-][0-9A-Za-z!#$%&*+./;=?@^_|~-]*)$"
    ]).unwrap();
}

#[derive(Debug)]
enum ChkRegOverlap {
    Yes(Region),
    No(Region, Region),
}

fn parse_usize_with_commas(s: &str) -> Option<usize> {
    s.replace(',', "").parse::<usize>().ok()
}

fn parse_range(s1: &str, s2: &str) -> Option<(usize, usize)> {
    match (parse_usize_with_commas(s1), parse_usize_with_commas(s2)) {
        (Some(x), Some(y)) => Some((x, y)),
        _ => None,
    }
}

fn parse_region_str(s: &str) -> Option<(&str, usize, Option<usize>)> {
    if let Some(cap) = RE_REGION1.captures(s).or_else(|| RE_REGION2.captures(s)) {
        match (cap.get(1), cap.get(2), cap.get(3)) {
            (Some(c), None, None) => Some((c.as_str(), 1, None)),
            (Some(c), Some(p), None) => {
                parse_usize_with_commas(p.as_str()).map(|x| (c.as_str(), x.max(1), None))
            }
            (Some(c), None, Some(q)) => {
                parse_usize_with_commas(q.as_str()).map(|x| (c.as_str(), 1, Some(x)))
            }
            (Some(c), Some(p), Some(q)) => parse_range(p.as_str(), q.as_str())
                .map(|(x, y)| (c.as_str(), x.max(1), Some(y.max(1)))),
            _ => None,
        }
    } else {
        None
    }
}

#[derive(Debug, Clone, Eq, PartialEq)]
pub struct Region {
    start: usize,
    end: Option<usize>,
    mappability: Option<Vec<[usize; 2]>>,
}

impl Ord for Region {
    fn cmp(&self, other: &Self) -> Ordering {
        match self.start.cmp(&other.start) {
            Ordering::Equal => match (self.end.as_ref(), other.end.as_ref()) {
                (Some(x), Some(y)) => x.cmp(y),
                (Some(_), None) => Ordering::Less,
                (None, Some(_)) => Ordering::Greater,
                _ => Ordering::Equal,
            },
            x => x,
        }
    }
}

impl PartialOrd for Region {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl fmt::Display for Region {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match (self.start, self.end) {
            (0, Some(y)) => {
                write!(f, ":-{}", y)
            }
            (x, None) => {
                write!(f, ":{}-", x)
            }
            (x, Some(y)) => {
                write!(f, ":{}-{}", x, y)
            }
        }
    }
}

impl Region {
    // Assume that self <= other
    fn check_overlap(mut self, other: Self) -> ChkRegOverlap {
        if self > other {
            return other.check_overlap(self);
        }
        match (self.end.as_ref(), other.end.as_ref()) {
            (Some(x), Some(y)) => {
                if other.start <= x + 1 {
                    self.end = Some(*x.max(y));
                    ChkRegOverlap::Yes(self)
                } else {
                    ChkRegOverlap::No(self, other)
                }
            }
            (Some(x), None) => {
                if other.start <= x + 1 {
                    self.end = None;
                    ChkRegOverlap::Yes(self)
                } else {
                    ChkRegOverlap::No(self, other)
                }
            }
            (None, Some(_)) => ChkRegOverlap::Yes(self),
            (None, None) => ChkRegOverlap::Yes(self),
        }
    }

    pub fn mappability(&self) -> Option<&Vec<[usize; 2]>> {
        self.mappability.as_ref()
    }

    pub fn start(&self) -> usize {
        self.start
    }

    pub fn end(&self) -> Option<usize> {
        self.end
    }

    pub fn length(&self) -> Option<usize> {
        self.end.map(|x| x + 1 - self.start)
    }
}

#[derive(Default)]
pub struct Regions {
    ctg_reg: BTreeMap<Rc<str>, Vec<Region>>,
    ctg_vec: Vec<Option<Rc<str>>>,
}

impl Regions {
    fn add_contig(&mut self, ctg: &str) {
        if !self.ctg_reg.contains_key(ctg) {
            let ctg = Rc::from(ctg);
            self.ctg_reg.insert(ctg, Vec::new());
        }
    }

    fn add_region_from_str(&mut self, region_str: &str) -> anyhow::Result<()> {
        if let Some((ctg, start, end)) = parse_region_str(region_str) {
            self.add_region_from_parts(ctg, start, end)
        } else {
            Err(anyhow!("Could not parse regions string {}", region_str))
        }
    }

    fn add_region_from_parts(
        &mut self,
        ctg: &str,
        start: usize,
        end: Option<usize>,
    ) -> anyhow::Result<()> {
        if let Some(x) = end {
            if x < start {
                return Err(anyhow!(
                    "Invalid region {}:{}-{} Start point after end point",
                    ctg,
                    start,
                    x
                ));
            }
        }
        self.add_contig(ctg);
        let region = Region {
            start,
            end,
            mappability: None,
        };
        self.ctg_reg.get_mut(ctg).unwrap().push(region);
        Ok(())
    }

    pub fn from_str_vec<S: AsRef<str>, I: IntoIterator<Item = S>>(
        reg_str: I,
    ) -> anyhow::Result<Self> {
        let mut regions = Self::default();

        for reg_str in reg_str.into_iter() {
            regions.add_region_from_str(reg_str.as_ref())?
        }
        Ok(regions.merge())
    }

    pub fn from_file<P: AsRef<str>>(file: P) -> anyhow::Result<Self> {
        let fname = file.as_ref();
        let mut buf = String::new();
        let mut f = compress::CompressIo::new()
            .path(fname)
            .bufreader()
            .with_context(|| "Could not open regions file for input")?;
        let mut regions = Self::default();

        loop {
            buf.clear();

            if f.read_line(&mut buf)? == 0 {
                break;
            }
            let fd: Vec<_> = buf.trim().split('\t').collect();

            if fd.len() > 1 {
                if !RE_SET_CONTIG.is_match(fd[0]) {
                    return Err(anyhow!("Illegal contig name {}", fd[0]));
                }
                let start = fd.get(1).and_then(|s| s.parse::<usize>().ok()).unwrap_or(0);
                let stop = fd.get(2).and_then(|s| s.parse::<usize>().ok());
                regions.add_region_from_parts(fd[0], start, stop)?;
            }
        }
        Ok(regions.merge())
    }

    fn merge(mut self) -> Self {
        for v in self.ctg_reg.values_mut() {
            v.sort_unstable();
            let mut nv = Vec::new();
            let mut prev: Option<Region> = None;
            for r in v.drain(..) {
                prev = Some(if let Some(p) = prev.take() {
                    match p.check_overlap(r) {
                        ChkRegOverlap::Yes(a) => a,
                        ChkRegOverlap::No(a, b) => {
                            nv.push(a);
                            b
                        }
                    }
                } else {
                    r
                });
            }
            if let Some(p) = prev.take() {
                nv.push(p)
            }
            *v = nv;
        }
        self
    }

    pub fn split_regions(&mut self, max_block_size: usize, len_hash: HashMap<&str, usize>) {
        debug!(
            "Splitting regions into block size of maximum {}",
            max_block_size
        );

        for (ctg, v) in self.ctg_reg.iter_mut() {
            let l = *len_hash.get(ctg.as_ref()).expect("Missing sequence length");
            let mut nvec = Vec::new();
            trace!("split_regions - {}", ctg);
            for reg in v.drain(..) {
                let end = reg.end.unwrap_or(l).min(l);
                let start = reg.start.max(1);
                assert!(reg.start <= end, "Invalid region");
                let size = end - start + 1;
                if size <= max_block_size {
                    trace!("Region {} (original)", reg);
                    nvec.push(reg)
                } else {
                    trace!("Splitting Region {}:{}-{}", ctg, start, end);
                    let mut ns = 1 + size / max_block_size;
                    let mut x = start;
                    while ns > 0 {
                        assert!(x <= end);
                        let remainder = end - x;
                        let s = remainder / ns;
                        let reg1 = Region {
                            start: x,
                            end: Some(x + s),
                            mappability: None,
                        };
                        trace!("Region {} (split)", reg1);
                        nvec.push(reg1);
                        x += s + 1;
                        ns -= 1;
                    }
                    assert_eq!(x, end + 1);
                }
            }
            v.append(&mut nvec);
        }
    }

    pub fn add_mappability(&mut self, map: &Mappability, max_gap: usize) {
        debug!("Adding mappability to regions");
        let mut ncreg = BTreeMap::new();
        for (ctg, v) in self.ctg_reg.iter_mut() {
            if let Some(mapv) = map.ctg_regions(ctg) {
                let mut nv = Vec::new();
                let mut ix = 0;
                for mut reg in v.drain(..) {
                    if ix < mapv.len() {
                        trace!("Looking for mappable regions in {}", reg);
                        for mreg in mapv[ix..].iter() {
                            if mreg[1] >= reg.start || reg.end.map(|x| mreg[0] > x).unwrap_or(false)
                            {
                                break;
                            }
                            ix += 1;
                        }
                    }
                    if ix >= mapv.len() {
                        reg.mappability = Some(Vec::new());
                        nv.push(reg);
                        continue;
                    }
                    let mut map_reg: Vec<[usize; 2]> = Vec::new();
                    for mreg in mapv[ix..].iter() {
                        let y = reg.end.unwrap_or(mreg[1]);
                        if mreg[1] >= reg.start && mreg[0] <= y {
                            let i = map_reg.len();
                            if i > 1 && mreg[0] - map_reg[i - 1][1] > max_gap {
                                let reg1 = Region {
                                    start: reg.start,
                                    end: Some(map_reg[i - 1][1]),
                                    mappability: Some(map_reg),
                                };
                                trace!("(A) Adding region {}.  Subregions: {}", reg1, i);
                                reg.start = reg1.end.unwrap() + 1;
                                nv.push(reg1);
                                map_reg = Vec::new();
                            }
                            map_reg.push([mreg[0].max(reg.start), mreg[1].min(y)]);
                        } else {
                            break;
                        }
                    }
                    if reg.end.map(|x| reg.start < x).unwrap_or(true) {
                        trace!("(B) Adding region {}.  Subregions: {}", reg, map_reg.len());
                        reg.mappability = Some(map_reg);
                        nv.push(reg);
                    }
                }
                if !nv.is_empty() {
                    trace!("Inserting regions ({}) for {}", nv.len(), ctg);
                    ncreg.insert(Rc::clone(ctg), nv);
                }
            } else {
                let nv: Vec<_> = v
                    .drain(..)
                    .map(|mut r| {
                        r.mappability = Some(Vec::new());
                        r
                    })
                    .collect();
                trace!("Inserting regions ({}) for {}", nv.len(), ctg);
                ncreg.insert(Rc::clone(ctg), nv);
            }
        }
        self.ctg_reg = ncreg;
    }

    pub fn add_tid_info(&mut self, seq_names: &[&str]) {
        let v: Vec<_> = seq_names
            .iter()
            .enumerate()
            .map(|(_, ctg)| {
                self.ctg_reg
                    .get_key_value(*ctg)
                    .map(|(ctg, _)| Rc::clone(ctg))
            })
            .collect();
        self.ctg_vec = v;
    }

    pub fn iter(&self) -> Iter<'_, Rc<str>, Vec<Region>> {
        self.ctg_reg.iter()
    }

    pub fn values(&self) -> Values<'_, Rc<str>, Vec<Region>> {
        self.ctg_reg.values()
    }

    pub fn len(&self) -> usize {
        self.ctg_reg.values().map(|v| v.len()).sum()
    }

    pub fn ctg_regions(&self, ctg: &str) -> Option<&Vec<Region>> {
        self.ctg_reg.get(ctg)
    }

    pub fn tid2ctg(&self, tid: usize) -> Option<&str> {
        self.ctg_vec[tid].as_ref().map(|s| s.as_ref())
    }

    pub fn fix_open_intervals(&mut self, seq: &[&str], len: &[usize]) {
        let shash: HashMap<_, _> = seq.iter().zip(len.iter()).map(|(s, x)| (*s, *x)).collect();
        for (ctg, reg_vec) in self.ctg_reg.iter_mut().map(|(c, v)| (c.as_ref(), v)) {
            for reg in reg_vec.iter_mut() {
                if reg.end.is_none() {
                    reg.end = Some(*shash.get(ctg).unwrap())
                }
            }
        }
    }

    pub fn ctg_regions_length(&self, ctg: &str) -> Option<usize> {
        self.ctg_reg.get(ctg).map(|v| {
            let mut l = 0;
            for reg in v.iter() {
                let rl = reg.length().expect("Open region found");
                l += rl;
            }
            l
        })
    }

    pub fn total_regions_Length(&self) -> usize {
        self.ctg_reg
            .keys()
            .map(|ctg| self.ctg_regions_length(ctg).expect("Missing contig length"))
            .sum::<usize>()
    }
}

pub fn find_overlapping_regions(rvec: &[Region], x: usize, y: usize) -> Vec<usize> {
    assert!(y >= x);
    let mut v = Vec::new();

    // Find first region overlapping x-y
    if let Some(i) = match rvec.binary_search_by_key(&x, |r| r.start()) {
        Ok(i) => Some(i),
        Err(i) => {
            if i > 0 {
                // Check if the read overlaps previous region
                if rvec[i - 1].end().map(|end| end >= x).unwrap_or(false) {
                    Some(i - 1)
                } else if i < rvec.len() && rvec[i].start() <= y {
                    Some(i)
                } else {
                    None
                }
            } else if rvec[i].start() <= y {
                Some(i)
            } else {
                None
            }
        }
    } {
        v.push(i);

        // Check if the following regions also overlap.  If so add to v
        for (ix, r) in rvec[i + 1..].iter().enumerate() {
            if r.start() <= y {
                v.push(ix + i + 1)
            } else {
                break;
            }
        }
    }

    v
}

// Copy of region information for a task
pub struct TaskRegion<'a> {
    ctg: &'a str,
    regions: Vec<(&'a Region, usize, Option<Sequence>)>, // usize parameter gives index into list of regions for this contig
}

impl<'a> TaskRegion<'a> {
    pub fn new(ctg: &'a str, regions: Vec<(&'a Region, usize, Option<Sequence>)>) -> Self {
        Self { ctg, regions }
    }

    pub fn ctg(&self) -> &str {
        self.ctg
    }

    pub fn regions(&self) -> &[(&Region, usize, Option<Sequence>)] {
        &self.regions
    }
}
