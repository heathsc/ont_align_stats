use std::{collections::BTreeMap, ops::AddAssign};

use serde::{
    self,
    ser::{SerializeMap, Serializer},
    Serialize,
};

pub enum StatType {
    Mappings = 0,
    Reads,
    Unmapped,
    Reversed,
    Secondary,
    QcFail,
    Duplicate,
    Paired,
    TotalBases,
    MappedBases,
}

const STAT_NAMES: [&str; 10] = [
    "Mappings",
    "Reads",
    "Unmapped",
    "Reversed",
    "Secondary",
    "QcFail",
    "Duplicate",
    "Paired",
    "TotalBases",
    "MappedBases",
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

fn serialize_vec<S>(ct: &[usize], serializer: S) -> Result<S::Ok, S::Error>
where
    S: Serializer,
{
    let l: usize = ct.iter().fold(0, |s, v| if *v > 0 { s + 1 } else { s });
    let mut map = serializer.serialize_map(Some(l))?;
    for (ix, val) in ct.iter().enumerate() {
        if *val > 0 {
            map.serialize_entry(&ix, &val)?;
        }
    }
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

#[derive(Default, Debug, Serialize)]
pub struct Stats {
    #[serde(serialize_with = "serialize_counts")]
    counts: [usize; 10],
    #[serde(serialize_with = "serialize_vec")]
    mapped_pctg: Vec<usize>,
    #[serde(serialize_with = "serialize_vec")]
    base_qual_pctg: Vec<usize>,
    #[serde(serialize_with = "serialize_mm_vec")]
    mismatch_pctg: Vec<usize>,
    #[serde(serialize_with = "serialize_vec")]
    indel_pctg: Vec<usize>,
    #[serde(serialize_with = "serialize_vec")]
    primary_mapq: Vec<usize>,
    read_len: BTreeMap<usize, usize>,
    n_splits: BTreeMap<usize, usize>,
    #[serde(serialize_with = "serialize_cov")]
    coverage: BTreeMap<usize, usize>,
}

fn add_vec<T: AddAssign + Copy>(v1: &mut [T], v2: &[T]) {
    for (ct, ct1) in v1.iter_mut().zip(v2.iter()) {
        *ct += *ct1
    }
}

fn add_btreemap<T: AddAssign + Copy + Default>(
    bt1: &mut BTreeMap<usize, T>,
    vt2: &BTreeMap<usize, T>,
) {
    for (x, y) in vt2.iter() {
        *bt1.entry(*x).or_insert(T::default()) += *y
    }
}

impl AddAssign for Stats {
    fn add_assign(&mut self, other: Self) {
        add_vec(&mut self.counts, &other.counts);
        add_vec(&mut self.mapped_pctg, &other.mapped_pctg);
        add_vec(&mut self.base_qual_pctg, &other.base_qual_pctg);
        add_vec(&mut self.mismatch_pctg, &other.mismatch_pctg);
        add_vec(&mut self.indel_pctg, &other.indel_pctg);
        add_vec(&mut self.primary_mapq, &other.primary_mapq);
        add_btreemap(&mut self.read_len, &other.read_len);
        add_btreemap(&mut self.n_splits, &other.n_splits);
        add_btreemap(&mut self.coverage, &other.coverage);
    }
}

impl Stats {
    pub fn new() -> Self {
        let mut st = Self::default();
        st.primary_mapq.resize(256, 0);
        st.mapped_pctg.resize(101, 0);
        st.base_qual_pctg.resize(101, 0);
        st.indel_pctg.resize(101, 0);
        st.mismatch_pctg.resize(1001, 0);
        st
    }

    pub fn incr(&mut self, ty: StatType) {
        self.counts[ty as usize] += 1
    }

    pub fn incr_n(&mut self, ty: StatType, n: usize) {
        self.counts[ty as usize] += n
    }

    pub fn incr_mapq(&mut self, q: u8) {
        self.primary_mapq[q as usize] += 1
    }

    pub fn incr_baseq_pctg(&mut self, q: u8) {
        self.base_qual_pctg[q as usize] += 1
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

    pub fn counts(&self, ty: StatType) -> usize {
        self.counts[ty as usize]
    }
}
