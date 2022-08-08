use std::{
    cmp::{Ord, Ordering},
    collections::{btree_map, BTreeMap, HashMap},
    fmt,
    io::BufRead,
    rc::Rc,
    slice,
};

use crate::mappability::Mappability;
use anyhow::Context;
use compress_io::compress;
use lazy_static::lazy_static;
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
            (Some(c), None, None) => Some((c.as_str(), 0, None)),
            (Some(c), Some(p), None) => {
                parse_usize_with_commas(p.as_str()).map(|x| (c.as_str(), x, None))
            }
            (Some(c), None, Some(q)) => {
                parse_usize_with_commas(q.as_str()).map(|x| (c.as_str(), 0, Some(x)))
            }
            (Some(c), Some(p), Some(q)) => {
                parse_range(p.as_str(), q.as_str()).map(|(x, y)| (c.as_str(), x, Some(y)))
            }
            _ => None,
        }
    } else {
        None
    }
}

#[derive(Debug, Clone, Eq, PartialEq)]
pub struct Region {
    contig: Rc<str>,
    start: usize,
    end: Option<usize>,
    mappability: Option<Vec<[usize; 2]>>,
}

impl Ord for Region {
    fn cmp(&self, other: &Self) -> Ordering {
        match self.contig.cmp(&other.contig) {
            Ordering::Equal => match self.start.cmp(&other.start) {
                Ordering::Equal => match (self.end.as_ref(), other.end.as_ref()) {
                    (Some(x), Some(y)) => x.cmp(y),
                    (Some(_), None) => Ordering::Less,
                    (None, Some(_)) => Ordering::Greater,
                    _ => Ordering::Equal,
                },
                x => x,
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

fn write_ctg(f: &mut fmt::Formatter, s: &str) -> fmt::Result {
    if f.alternate() && s.contains(':') {
        write!(f, "{{{}}}", s)
    } else {
        write!(f, "{}", s)
    }
}

impl fmt::Display for Region {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match (self.start, self.end) {
            (0, None) => write_ctg(f, self.contig.as_ref()),
            (0, Some(y)) => {
                write_ctg(f, self.contig.as_ref())?;
                write!(f, ":-{}", y)
            }
            (x, None) => {
                write_ctg(f, self.contig.as_ref())?;
                write!(f, ":{}-", x)
            }
            (x, Some(y)) => {
                write_ctg(f, self.contig.as_ref())?;
                write!(f, ":{}-{}", x, y)
            }
        }
    }
}

impl Region {
    pub fn contig(&self) -> &str {
        self.contig.as_ref()
    }

    // Assume that self <= other
    fn check_overlap(mut self, other: Self) -> ChkRegOverlap {
        if self > other {
            return other.check_overlap(self);
        }
        if self.contig != other.contig {
            ChkRegOverlap::No(self, other)
        } else {
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
    }

    pub fn mappability(&self) -> Option<&Vec<[usize; 2]>> {
        self.mappability.as_ref()
    }
}

#[derive(Default)]
pub struct Regions {
    ctg_reg: BTreeMap<Rc<str>, Vec<Region>>,
}

impl Regions {
    fn add_contig(&mut self, ctg: &str) -> Rc<str> {
        if let Some((k, _)) = self.ctg_reg.get_key_value(ctg) {
            Rc::clone(k)
        } else {
            let ctg = Rc::from(ctg);
            let ctg1 = Rc::clone(&ctg);
            self.ctg_reg.insert(ctg1, Vec::new());
            ctg
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
        let contig = self.add_contig(ctg);
        let region = Region {
            contig,
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
                            contig: Rc::clone(ctg),
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
                    trace!("Looking for mappable regions in {}", reg);
                    for mreg in mapv[ix..].iter() {
                        if mreg[1] >= reg.start || reg.end.map(|x| mreg[0] > x).unwrap_or(false) {
                            break;
                        }
                        ix += 1;
                    }
                    if ix >= mapv.len() {
                        break;
                    }
                    let mut map_reg: Vec<[usize; 2]> = Vec::new();
                    for mreg in mapv[ix..].iter() {
                        let y = reg.end.unwrap_or(mreg[1]);
                        if mreg[1] >= reg.start && mreg[0] <= y {
                            let i = map_reg.len();
                            if i > 1 && mreg[0] - map_reg[i - 1][1] > max_gap {
                                let reg1 = Region {
                                    contig: Rc::clone(ctg),
                                    start: reg.start.max(map_reg[0][0]),
                                    end: Some(map_reg[i - 1][1]),
                                    mappability: Some(map_reg),
                                };
                                trace!("(A) Adding region {}.  Subregions: {}", reg1, i);
                                nv.push(reg1);
                                map_reg = Vec::new();
                            }
                            map_reg.push([mreg[0].max(reg.start), mreg[1].min(y)]);
                        } else {
                            break;
                        }
                    }
                    let i = map_reg.len();
                    if i > 0 {
                        reg.start = reg.start.max(map_reg[0][0]);
                        reg.end = reg.end.map(|x| x.min(map_reg[i - 1][1]));
                        trace!("(B) Adding region {}.  Subregions: {}", reg, map_reg.len());
                        reg.mappability = Some(map_reg);
                        nv.push(reg);
                    }
                }
                if !nv.is_empty() {
                    trace!("Inserting regions ({}) for {}", nv.len(), ctg);
                    ncreg.insert(Rc::clone(ctg), nv);
                }
            }
        }
        self.ctg_reg = ncreg;
    }

    pub fn iter(&self) -> RegionIter {
        let mut reg_vec = self.ctg_reg.values();
        let regs = reg_vec.next().map(|v| v.iter());
        RegionIter { reg_vec, regs }
    }

    pub fn len(&self) -> usize {
        self.ctg_reg.values().map(|v| v.len()).sum()
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
}

pub struct RegionIter<'a> {
    reg_vec: btree_map::Values<'a, Rc<str>, Vec<Region>>,
    regs: Option<slice::Iter<'a, Region>>,
}

impl<'a> Iterator for RegionIter<'a> {
    type Item = &'a Region;

    fn next(&mut self) -> Option<Self::Item> {
        if self.regs.is_none() {
            None
        } else {
            match self.regs.as_mut().and_then(|v| v.next()) {
                Some(x) => Some(x),
                None => {
                    self.regs = self.reg_vec.next().map(|v| v.iter());
                    self.next()
                }
            }
        }
    }
}
