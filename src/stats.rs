use std::{collections::BTreeMap, ops::AddAssign};

use serde::{
    self,
    ser::{SerializeMap, Serializer},
    Serialize,
};

use crate::read::{bisulfite::BSStrand, coverage::Matches};

pub enum StatType {
    Mappings = 0,
    Reads,
    Unmapped,
    Reversed,
    Secondary,
    QcFail,
    Duplicate,
    TotalBases,
    MappedBases,
    TotalPairs,
    CorrectPairs,
    MateUnmapped,
    DifferentContigs,
    BadTemplateLength,
    ConvergentPair,
    DivergentPair,
    OrientationFF,
    OrientationFR,
    OrientationRF,
    OrientationRR,
    IllegalOrientation,
    BisulfiteC2T,
    BisulfiteG2A,
    OverlapBases,
}

const N_COUNTS: usize = (StatType::OverlapBases as usize) + 1;

const STAT_NAMES: [&str; N_COUNTS] = [
    "Mappings",
    "Reads",
    "Unmapped",
    "Reversed",
    "Secondary",
    "QcFail",
    "Duplicate",
    "TotalBases",
    "MappedBases",
    "TotalPairs",
    "CorrectPairs",
    "MateUnmapped",
    "DifferentContigs",
    "BadTemplateLength",
    "ConvergentPair",
    "DivergentPair",
    "OrientationFF",
    "OrientationFR",
    "OrientationRF",
    "OrientationRR",
    "IllegalOrientation",
    "BisulfiteC2T",
    "BisulfiteG2A",
    "OverlapBases",
];

fn serialize_counts<S>(ct: &[usize], serializer: S) -> Result<S::Ok, S::Error>
where
    S: Serializer,
{
    let l: usize = ct.iter().fold(0, |s, v| if *v > 0 { s + 1 } else { s });
    let mut map = serializer.serialize_map(Some(l))?;
    for (k, v) in STAT_NAMES.iter().zip(ct.iter()) {
        if *v > 0 {
            map.serialize_entry(*k, v)?;
        }
    }
    map.end()
}

fn serialize_vec<S, T>(ct: &[T], serializer: S) -> Result<S::Ok, S::Error>
where
    S: Serializer,
    T: Serialize,
{
    let mut map = serializer.serialize_map(Some(ct.len()))?;
    for (ix, val) in ct.iter().enumerate() {
        map.serialize_entry(&ix, &val)?;
    }
    map.end()
}

const BASE: [char; 4] = ['A', 'C', 'G', 'T'];

#[derive(Serialize)]
struct CountProp {
    count: u64,
    proportion: f64,
}

fn serialize_matches<S>(m: &Matches, serializer: S) -> Result<S::Ok, S::Error>
where
    S: Serializer,
{
    let mut map = serializer.serialize_map(Some(13))?;

    let mut tot_mm = 0;
    let mut tot = 0;
    for (i, v) in m.iter().enumerate() {
        tot += v[..4].iter().map(|x| *x as u64).sum::<u64>();
        for j in (0..4).filter(|j| *j != i) {
            let t = m.iter().map(|x| x[j] as u64).sum::<u64>() as f64;
            let count = v[j] as u64;
            let proportion = (count as f64) / t;
            let cp = CountProp { count, proportion };
            tot_mm += count;
            map.serialize_entry(format!("{}->{}", BASE[j], BASE[i]).as_str(), &cp)?
        }
    }
    let cp = CountProp {
        count: tot_mm,
        proportion: (tot_mm as f64) / (tot as f64),
    };
    map.serialize_entry("All mismatches", &cp)?;
    map.end()
}

fn serialize_mm_vec<S>(ct: &[usize], serializer: S) -> Result<S::Ok, S::Error>
where
    S: Serializer,
{
    let l: usize = ct.iter().fold(0, |s, v| if *v > 0 { s + 1 } else { s });
    let mut map = serializer.serialize_map(Some(l))?;
    for (ix, val) in ct.iter().enumerate() {
        if *val > 0 {
            let key = format!("{:.1}", (ix as f64) / 10.0);
            map.serialize_entry(&key, &val)?;
        }
    }
    map.end()
}

fn serialize_cov<S>(ct: &BTreeMap<usize, usize>, serializer: S) -> Result<S::Ok, S::Error>
where
    S: Serializer,
{
    let mut map = serializer.serialize_map(Some(ct.len()))?;
    let tot: usize = ct.values().copied().sum();
    let mut rest = tot;
    let ftot = tot as f64;

    for (ix, val) in ct.iter() {
        let z = (rest as f64) / ftot;
        rest -= *val;
        map.serialize_entry(&ix, &(*val, z))?;
    }
    map.end()
}

#[derive(Default, Debug, Copy, Clone)]
pub struct BaseCounts {
    counts: [u64; 4],
}

impl AddAssign for BaseCounts {
    fn add_assign(&mut self, other: Self) {
        for (a, b) in self.counts.iter_mut().zip(other.counts.iter()) {
            *a += b
        }
    }
}

impl Serialize for BaseCounts {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let t = self.counts.iter().sum::<u64>() as f64;

        let mut map = serializer.serialize_map(Some(4))?;
        for (k, count) in BASE.iter().zip(self.counts.iter().copied()) {
            let proportion = if count > 0 { (count as f64) / t } else { 0.0 };
            let cp = CountProp { count, proportion };
            map.serialize_entry(k, &cp)?;
        }
        map.end()
    }
}

impl BaseCounts {
    pub fn incr_base(&mut self, b: u8) {
        self.counts[b as usize] += 1;
    }
}

#[derive(Debug, Serialize, Hash, Copy, Clone, Ord, PartialOrd, Eq, PartialEq)]
pub enum ReadType {
    Forwards,
    Reverse,
    ForwardsRead1,
    ForwardsRead2,
    ReverseRead1,
    ReverseRead2,
}

#[derive(Default, Debug, Serialize)]
pub struct Stats {
    #[serde(serialize_with = "serialize_counts")]
    counts: [usize; N_COUNTS],

    composition: BTreeMap<ReadType, BaseCounts>,

    #[serde(serialize_with = "serialize_matches")]
    mismatches: Matches,

    #[serde(serialize_with = "serialize_vec")]
    mapped_pctg: Vec<usize>,

    #[serde(serialize_with = "serialize_vec")]
    base_qual_hist: Vec<BaseCounts>,

    #[serde(serialize_with = "serialize_mm_vec")]
    #[serde(skip_serializing_if = "Vec::is_empty")]
    mismatch_pctg: Vec<usize>,

    #[serde(serialize_with = "serialize_vec")]
    indel_pctg: Vec<usize>,

    #[serde(serialize_with = "serialize_vec")]
    primary_mapq: Vec<usize>,

    read_len: BTreeMap<usize, usize>,

    n_splits: BTreeMap<usize, usize>,

    #[serde(serialize_with = "serialize_cov")]
    coverage: BTreeMap<usize, usize>,

    #[serde(skip_serializing_if = "BTreeMap::is_empty")]
    template_len: BTreeMap<usize, usize>,
}

fn add_vec<T: AddAssign + Copy>(v1: &mut [T], v2: &[T]) {
    for (ct, ct1) in v1.iter_mut().zip(v2.iter()) {
        *ct += *ct1
    }
}

fn add_btreemap<T: AddAssign + Copy + Default, K: Ord + Copy>(
    bt1: &mut BTreeMap<K, T>,
    vt2: &BTreeMap<K, T>,
) {
    for (x, y) in vt2.iter() {
        *bt1.entry(*x).or_insert(T::default()) += *y
    }
}

fn add_matches(x: &mut Matches, y: &Matches) {
    for (p, q) in x.iter_mut().zip(y.iter()) {
        add_vec(p, q)
    }
}

impl AddAssign for Stats {
    fn add_assign(&mut self, other: Self) {
        add_vec(&mut self.counts, &other.counts);
        add_vec(&mut self.mapped_pctg, &other.mapped_pctg);
        add_vec(&mut self.base_qual_hist, &other.base_qual_hist);
        add_vec(&mut self.mismatch_pctg, &other.mismatch_pctg);
        add_vec(&mut self.indel_pctg, &other.indel_pctg);
        add_vec(&mut self.primary_mapq, &other.primary_mapq);
        add_matches(&mut self.mismatches, &other.mismatches);
        add_btreemap(&mut self.read_len, &other.read_len);
        add_btreemap(&mut self.n_splits, &other.n_splits);
        add_btreemap(&mut self.coverage, &other.coverage);
        add_btreemap(&mut self.template_len, &other.template_len);
        add_btreemap(&mut self.composition, &other.composition);
    }
}

impl Stats {
    pub fn new() -> Self {
        let mut st = Self::default();
        st.primary_mapq.resize(256, 0);
        st.mapped_pctg.resize(101, 0);
        st.base_qual_hist.resize(101, BaseCounts::default());
        st.indel_pctg.resize(101, 0);
        st.mismatch_pctg.resize(1001, 0);
        st
    }

    pub fn incr(&mut self, ty: StatType) {
        self.counts[ty as usize] += 1
    }

    pub fn composition_get_mut(&mut self, bc: ReadType) -> &mut BaseCounts {
        self.composition
            .entry(bc)
            .or_insert_with(BaseCounts::default)
    }

    pub fn incr_n(&mut self, ty: StatType, n: usize) {
        self.counts[ty as usize] += n
    }

    pub fn incr_mapq(&mut self, q: u8) {
        self.primary_mapq[q as usize] += 1
    }

    pub fn incr_baseq_counts(&mut self, q: u8, b: usize) {
        self.base_qual_hist[q as usize].counts[b] += 1
    }

    pub fn update_readlen_stats(&mut self, rl: usize, used: usize) {
        assert!(used <= rl);
        let p = (100.0 * (used as f64) / (rl as f64)).round() as usize;
        *self.read_len.entry(rl).or_insert(0) += 1;
        self.counts[StatType::TotalBases as usize] += rl;
        self.counts[StatType::MappedBases as usize] += used;
        self.mapped_pctg[p] += 1
    }

    pub fn incr_mismatch_pctg(&mut self, p: f64) {
        assert!((0.0..=1.0).contains(&p), "Illegal mismatch percentage");
        self.mismatch_pctg[(p * 1000.0).round() as usize] += 1;
    }

    pub fn incr_matches(&mut self, m: &Matches) {
        add_matches(&mut self.mismatches, m)
    }

    pub fn incr_template_len(&mut self, x: usize) {
        *self.template_len.entry(x).or_insert(0) += 1;
    }

    pub fn incr_indel_pctg(&mut self, p: f64) {
        assert!((0.0..=1.0).contains(&p), "Illegal mismatch percentage");
        self.indel_pctg[(p * 100.0).round() as usize] += 1;
    }

    pub fn incr_coverage(&mut self, cov: usize) {
        *self.coverage.entry(cov).or_insert(0) += 1
    }

    pub fn incr_n_splits(&mut self, n_splits: usize) {
        *self.n_splits.entry(n_splits).or_insert(0) += 1;
    }

    pub fn incr_bisulfite_strand(&mut self, bs: BSStrand) {
        match bs {
            BSStrand::StrandC2T => self.incr(StatType::BisulfiteC2T),
            BSStrand::StrandG2A => self.incr(StatType::BisulfiteG2A),
        }
    }

    pub fn counts(&self, ty: StatType) -> usize {
        self.counts[ty as usize]
    }
}
