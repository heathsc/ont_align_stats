pub enum StatType {
    Mappings = 0,
    Reads,
    Mapped,
    Unmapped,
    Reversed,
    Secondary,
    QcFail,
    Duplicate,
    Paired,
}

#[derive(Default)]
pub struct Stats {
    counts: [usize; 9],
    primary_mapq: Vec<usize>,
}

impl Stats {
    pub fn new() -> Self {
        let mut st = Self::default();
        st.primary_mapq.resize(256, 0);
        st
    }

    pub fn incr(&mut self, ty: StatType) {
        self.counts[ty as usize] += 1
    }

    pub fn incr_mapq(&mut self, q: u8) {
        self.primary_mapq[q as usize] += 1
    }

    pub fn counts(&self, ty: StatType) -> usize {
        self.counts[ty as usize]
    }

    pub fn mapq(&self) -> &[usize] {
        &self.primary_mapq
    }
}
