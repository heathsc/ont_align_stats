use r_htslib::{BamRec, CigarOp, SeqQual, BAM_FPAIRED, BAM_FPROPER_PAIR, BAM_FREVERSE};
use std::ops::AddAssign;

use crate::stats::{BaseComposition, ReadType, StatType, Stats};

use super::bisulfite::BSStrand;

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
    [0, 1, 1, 1, 1], // Observed A
    [1, 0, 1, 1, 1], // Observed C
    [1, 1, 0, 1, 1], // Observed G
    [1, 1, 1, 0, 1], // Observed T
];

const MATCH_TAB_C2T: [[usize; 5]; 4] = [
    [0, 1, 1, 1, 1], // Observed A
    [1, 0, 1, 1, 1], // Observed C
    [1, 1, 0, 1, 1], // Observed G
    [1, 0, 1, 0, 1], // Observed T matches C or T in reference
];

const MATCH_TAB_G2A: [[usize; 5]; 4] = [
    [0, 1, 0, 1, 1], // Observed A matches A or G in reference
    [1, 0, 1, 1, 1], // Observed C
    [1, 1, 0, 1, 1], // Observed G
    [1, 1, 1, 0, 1], // Observed T
];

#[derive(Default)]
pub(super) struct Coverage {
    start: usize,         // Start of region
    cov: Vec<u32>,        // Coverage
    map: Option<Vec<u8>>, // Mappability bit map
    reference: Vec<u8>,   // Reference sequence
    match_counts: [u32; 2],
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
        self.match_counts = [0; 2]
    }

    pub(super) fn reset(
        &mut self,
        start: usize,
        end: usize,
        mappability: Option<&Vec<[usize; 2]>>,
        rf: Option<&[u8]>,
    ) {
        let (len, mlen) = Self::_calc_lens(start, end);
        self.start = start;
        self.cov.clear();
        self.cov.resize(len, 0);
        self.reference.clear();
        self.match_counts = [0; 2];
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
                assert!(mv[0] >= start);
                for x in (mv[0]..=mv[1]).map(|i| i - start) {
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
        z: u8,
        thresh: u8,
        bc: &mut BaseComposition,
        st: &mut Stats,
    ) {
        if x >= self.start {
            let i = x - self.start;
            if let Some(c) = self.cov.get_mut(i) {
                let (q, b) = (z >> 2, z & 3);
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
                    bc.incr_base(b);
                }
                st.incr_baseq_pctg(q);
            }
        }
    }

    pub(super) fn inc_with_mm(
        &mut self,
        x: usize,
        z: u8,
        thresh: u8,
        bc: &mut BaseComposition,
        mm_tab: &[[usize; 5]; 4],
        st: &mut Stats,
    ) {
        if x >= self.start {
            let i = x - self.start;
            if let Some(c) = self.cov.get_mut(i) {
                let (q, b) = (z >> 2, z & 3);
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
                    self.match_counts[mm_tab[b as usize][rf]] += 1;
                    bc.incr_base(b);
                }
                st.incr_baseq_pctg(q);
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
    rec: &mut BamRec,
    seq_qual: &SeqQual,
    cov: &mut Coverage,
    min_qual: u8,
    rd_type: ReadType,
    bs_strand: Option<BSStrand>,
    st: &mut Stats,
) {
    if let Some(cigar) = rec.cigar() {
        let mut x = rec.pos().expect("No start position for mapped read") + 1;
        let mut sq = seq_qual.iter().copied();

        let mut indel_counts = [0; 2];
        let (start, end) = (cov.start(), cov.end());

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

        let mm_tab = match bs_strand {
            Some(BSStrand::StrandC2T) => MATCH_TAB_C2T,
            Some(BSStrand::StrandG2A) => MATCH_TAB_G2A,
            None => MATCH_TAB,
        };
        let mut bc = BaseComposition::default();
        for elem in cigar.iter() {
            // If we are already past the end of the region we can stop processing the cigar
            if x >= end {
                break;
            }

            // Get length of current cigar op
            let l = elem.op_len() as usize;
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
                            let z = sq
                                .next()
                                .expect("Mismatch between Cigar and sequence length");
                            if x < end1 {
                                cov.inc(x, z, min_qual, &mut bc, st)
                            } else if x < end {
                                st.incr(StatType::OverlapBases)
                            }
                            x += 1;
                        }
                    } else {
                        for _ in 0..l {
                            let z = sq
                                .next()
                                .expect("Mismatch between Cigar and sequence length");
                            if x < end1 {
                                cov.inc_with_mm(x, z, min_qual, &mut bc, &mm_tab, st)
                            } else if x < end {
                                st.incr(StatType::OverlapBases)
                            }
                            x += 1;
                        }
                    }
                    indel_counts[0] += l1;
                }
                CigarOp::SoftClip => {
                    sq.nth(l - 1);
                }
                CigarOp::Ins => {
                    sq.nth(l - 1);
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
        if (rec.flag() & BAM_FREVERSE) != 0 {
            bc.complement()
        }
        st.composition_get_mut(rd_type).add_assign(bc);

        // Update stats
        if let Some(p) = cov.mismatch_pctg() {
            st.incr_mismatch_pctg(p);
            cov.clear_match_counts();
        }
        let t = indel_counts[0] + indel_counts[1];
        if t >= MIN_PCTG_N {
            let p = (indel_counts[1] as f64) / (t as f64);
            st.incr_indel_pctg(p);
        }
    }
}
