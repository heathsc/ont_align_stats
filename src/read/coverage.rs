use std::ops::AddAssign;

use m_htslib::{
    base::BaseQual,
    region::RegionCoords,
    sam::{
        CigarOp, SeqQualIter,
        record::{
            BamRec,
            bam1::{BAM_FPAIRED, BAM_FPROPER_PAIR},
        },
    },
};

use crate::{
    Config,
    stats::{BaseCounts, ReadType, StatType, Stats},
};

pub type Matches = [[u64; 5]; 4];

const MIN_PCTG_N: usize = 20;

const BASE_TAB: [u8; 256] = [
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
];

const MATCH_TAB: [[usize; 5]; 4] = [
    [0, 1, 1, 1, 2], // Observed A
    [1, 0, 1, 1, 2], // Observed C
    [1, 1, 0, 1, 2], // Observed G
    [1, 1, 1, 0, 2], // Observed T
];

#[derive(Default)]
pub(super) struct Coverage {
    start: usize,       // Start of region
    cov: Vec<u32>,      // Coverage 
    map: Option<Vec<u8>>, // Mappability bit map
    reference: Vec<u8>, // Reference sequence
    match_counts: [u32; 3],
}

impl Coverage {
    fn _calc_lens(start: usize, end: usize) -> (usize, usize) {
        assert!(end >= start);
        let l = end + 1 - start;
        (l, (l + 7) >> 3)
    }
    
    pub(super) fn new() -> Self {
        Self::default()
    }

    pub(super) fn start(&self) -> usize {
        self.start
    }

    pub(super) fn end(&self) -> usize {
        self.start + self.cov.len()
    }

    pub(super) fn has_reference(&self) -> bool {
        !self.reference.is_empty()
    }

    pub(super) fn clear_match_counts(&mut self) {
        self.match_counts = [0; 3];
    }

    pub(super) fn reset(&mut self, start: usize, end: usize, mappability: Option<&[RegionCoords]>, rf: Option<&[u8]>) {
        let (len, mlen) = Self::_calc_lens(start, end);
        self.start = start;
        self.cov.clear();
        self.cov.resize(len, 0);
        self.reference.clear();
        self.clear_match_counts();

        if let Some(p) = rf {
            self.reference.reserve(p.len());
            for c in p {
                self.reference.push(BASE_TAB[*c as usize])
            }
        }
        
        if let Some(map) = mappability {
            let m = if let Some(m) = self.map.as_mut() {
                m.clear();
                m.resize(mlen, 0);
                m
            } else {
                self.map = Some(vec![0; mlen]);
                self.map.as_mut().unwrap()
            };
            for mv in map.iter() {
                let s = (mv.start() as usize).max(start);
                let e = (mv.end().expect("Missing end value for mappability") as usize).min(end);
                for x in (s..e).map(|i| i - start) {
                    let ix = x >> 3;
                    let iy = x & 7;
                    m[ix] |= 1 << iy;
                }
            }
        } else {
            self.map = None
        }
    }

    pub(super) fn inc(
        &mut self,
        x: usize,
        z: BaseQual,
        thresh: u8,
        bc: &mut BaseCounts,
        st: &mut Stats,
    ) {
        if x >= self.start {
            let i = x - self.start;
            if let Some(c) = self.cov.get_mut(i) {
                let (b, q) = z.base_qual();
                if let Some(bs) = b.single_base() {
                    if q >= thresh {
                        if let Some(m) = self.map.as_mut() {
                            let ix = i >> 3;
                            let iy = i & 7;
                            if (m[ix] & (1 << iy)) != 0 {
                                *c += 1
                            }
                        } else {
                            *c += 1
                        }
                        bc.incr_base(bs);
                    }
                    st.incr_baseq_counts(q, bs as usize);
                }
            }
        }
    }

    pub(super) fn inc_with_mm(
        &mut self,
        x: usize,
        z: BaseQual,
        thresh: u8,
        bc: &mut BaseCounts,
        st: &mut Stats,
        matches: &mut Matches,
    ) {
        if x >= self.start {
            let i = x - self.start;
            if let Some(c) = self.cov.get_mut(i) {
                let (b, q) = z.base_qual();
                if let Some(bs) = b.single_base() {
                    if q >= thresh {
                        if let Some(m) = self.map.as_mut() {
                            let ix = i >> 3;
                            let iy = i & 7;
                            if (m[ix] & (1 << iy)) != 0 {
                                *c += 1
                            }
                        } else {
                            *c += 1
                        }
                        let rf = self.reference[i] as usize;
                        self.match_counts[MATCH_TAB[bs as usize][rf]] += 1;
                        matches[bs as usize][rf] += 1;
                        bc.incr_base(bs);
                    }
                    st.incr_baseq_counts(q, bs as usize);
                }
            }
        }
    }

    pub(super) fn mismatch_pctg(&self) -> Option<f64> {
        let t = self.match_counts[0] + self.match_counts[1];
        if t as usize >= MIN_PCTG_N {
            Some((self.match_counts[1] as f64) / (t as f64))
        } else {
            None
        }
    }

    pub(super) fn update_stats(&self, st: &mut Stats) {
        if let Some(map) = self.map.as_ref() {
            let mut it = self.cov.iter();
            for mut mask in map.iter().copied() {
                let mut n = 8;
                while mask != 0 {
                    if let Some(c) = it.next().map(|x| *x as usize) {
                        if (mask & 1) == 1 {
                            st.incr_coverage(c)
                        }
                        mask >>= 1;
                        n -= 1;
                    } else {
                        break;
                    }
                }
                if n > 0 {
                    it.nth(n - 1);
                }
            }
        } else {
            for c in self.cov.iter().map(|x| *x as usize) {
                st.incr_coverage(c)
            }
        }
    }
}

/// Collect coverage and other stats (mismatches, indel counts etc.) that require looking at each base.
/// Counts ae only performed for the bases that fall within the current region (as defined in cov)
/// to avoid double counting if a read spans multiple regions.
///
/// For paired reads where the proper pair flag is set, we identify the left most member of the pair (i.e., where
/// the template length is positive) and for this read we only collect coverage stats up to the start position
/// of the mate (on the reference).  Any bases that fall after that are counted as overlapping bases and do not
/// contribute to the coverage
pub(super) fn process_coverage(
    rec: &BamRec,
    mut seq_qual: SeqQualIter,
    cov: &mut Coverage,
    rd_type: ReadType,
    st: &mut Stats,
    cfg: &Config,
) {
    if let Some(cigar) = rec.cigar() {
        let slen = seq_qual.len();
        if let Some(l) = cfg.min_read_len()
            && slen < l
        {
            return;
        }
        if let Some(l) = cfg.max_read_len()
            && slen > l
        {
            return;
        }

        let min_qual = cfg.min_qual();
        let mut x = rec.pos().expect("No start position for mapped read");

        let mut indel_counts = [0; 2];
        let mut matches: Matches = [[0; 5]; 4];
        let (start, end) = (cov.start() as i64, cov.end() as i64);

        // Is this the left most member of a properly paired read pair?
        let end1 = if (rec.flag() & (BAM_FPAIRED | BAM_FPROPER_PAIR))
            == (BAM_FPAIRED | BAM_FPROPER_PAIR)
            && rec.template_len() > 0
        {
            rec.mpos()
                .expect("Missing mate position for correctly paired read")
                .min(end)
        } else {
            end
        };

        let mut bc = BaseCounts::default();
        for elem in cigar.iter() {
            // If we are already past the end of the region we can stop processing the cigar
            if x >= end {
                break;
            }

            // Get length of current cigar op
            let l = elem.op_len() as i64;
            assert!(l > 0, "Zero length cigar element");

            // Calculate how much of current op overlaps the region
            let l1 = if x >= end1 || x + l < start {
                0
            } else {
                let x1 = x.max(start);
                let y1 = (x + l).min(end1);
                y1 - x1
            };

            // Process cigar op
            match elem.op() {
                CigarOp::Match | CigarOp::Equal | CigarOp::Diff => {
                    if !cov.has_reference() {
                        for _ in 0..l {
                            let z = seq_qual
                                .next()
                                .expect("Mismatch between Cigar and sequence length");
                            if x < end1 {
                                cov.inc(x as usize, z, min_qual, &mut bc, st)
                            } else if x < end {
                                st.incr(StatType::OverlapBases)
                            }
                            x += 1;
                        }
                    } else {
                        for _ in 0..l {
                            let z = seq_qual
                                .next()
                                .expect("Mismatch between Cigar and sequence length");
                            if x < end1 {
                                cov.inc_with_mm(
                                    x as usize,
                                    z,
                                    min_qual,
                                    &mut bc,
                                    st,
                                    &mut matches,
                                )
                            } else if x < end {
                                st.incr(StatType::OverlapBases)
                            }
                            x += 1;
                        }
                    }
                    indel_counts[0] += l1;
                }
                CigarOp::SoftClip => {
                    seq_qual.nth(l as usize - 1);
                }
                CigarOp::Ins => {
                    seq_qual.nth(l as usize - 1);
                    indel_counts[1] += l1;
                }
                CigarOp::Del => {
                    x += l;
                    // For deletions we do not advance along the reference so if l1 > 0 then
                    // we score the entire deletion
                    if l1 > 0 {
                        indel_counts[1] += l
                    }
                }
                _ => (),
            }
        }
        st.composition_get_mut(rd_type).add_assign(bc);

        // Update stats
        if let Some(p) = cov.mismatch_pctg() {
            st.incr_mismatch_pctg(p);
        }
        cov.clear_match_counts();
        st.incr_matches(&matches);
        let t = indel_counts[0] + indel_counts[1];
        if t >= MIN_PCTG_N as i64 {
            let p = (indel_counts[1] as f64) / (t as f64);
            st.incr_indel_pctg(p);
        }
    }
}
