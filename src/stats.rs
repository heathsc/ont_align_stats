use std::{collections::BTreeMap, ops::AddAssign};

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

#[derive(Default, Debug)]
pub struct Stats {
    counts: [usize; 10],
    mapped_pctg: Vec<usize>,
    primary_mapq: Vec<usize>,
    read_len: BTreeMap<usize, usize>,
    n_splits: BTreeMap<usize, usize>,
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

    pub fn update_readlen_stats(&mut self, rl: usize, used: usize) {
        assert!(used <= rl);
        let p = (100.0 * (used as f64) / (rl as f64)).round() as usize;
        *self.read_len.entry(rl).or_insert(0) += 1;
        self.counts[StatType::TotalBases as usize] += rl;
        self.counts[StatType::MappedBases as usize] += used;
        self.mapped_pctg[p] += 1
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

    pub fn mapq(&self) -> &[usize] {
        &self.primary_mapq
    }

    pub fn read_len(&self) -> &BTreeMap<usize, usize> {
        &self.read_len
    }

    pub fn n_splits(&self) -> &BTreeMap<usize, usize> {
        &self.n_splits
    }

    pub fn mapped_pctg(&self) -> &[usize] {
        &self.mapped_pctg
    }
}
