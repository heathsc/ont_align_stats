use std::path::Path;

use crossbeam_channel::{Receiver, Sender};
use r_htslib::*;

use crate::{
    config::Config,
    input, regions,
    stats::{StatType, Stats},
};

fn get_start_end(cigar: &Cigar, rev: bool, read_len: usize) -> Option<(usize, usize)> {
    let mut first = 0;
    let mut last = 0;
    let mut started = false;
    let mut x = 0;
    for elem in cigar.iter() {
        match elem.op() {
            CigarOp::HardClip | CigarOp::SoftClip | CigarOp::Ins => x += elem.op_len() as usize,
            CigarOp::Match | CigarOp::Equal | CigarOp::Diff => {
                if !started {
                    started = true;
                    first = x;
                }
                x += elem.op_len() as usize;
                last = x;
            }
            _ => {}
        }
    }
    adjust_start_end(first, last, x, rev, read_len)
}

fn adjust_start_end(
    first: usize,
    last: usize,
    x: usize,
    rev: bool,
    read_len: usize,
) -> Option<(usize, usize)> {
    if read_len != x || last <= first {
        None
    } else if rev {
        Some((read_len - last, read_len - 1 - first))
    } else {
        Some((first, last - 1))
    }
}

fn sa_get_start_end(it: &mut SaTagIter, rev: bool, read_len: usize) -> Option<(usize, usize)> {
    let mut first = 0;
    let mut last = 0;
    let mut started = false;
    let mut x = 0;
    for (c, l) in it {
        match c {
            'H' | 'S' | 'I' => x += l,
            'M' | '=' | 'X' => {
                if !started {
                    started = true;
                    first = x;
                }
                x += l;
                last = x;
            }
            _ => {}
        }
    }
    adjust_start_end(first, last, x, rev, read_len)
}

struct SaTagIter<'a> {
    tag: &'a str,
    done: bool,
}

impl<'a> SaTagIter<'a> {
    fn new(tag: &'a str) -> Self {
        Self { tag, done: false }
    }
}

impl<'a> Iterator for SaTagIter<'a> {
    type Item = (char, usize);
    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }
        let mut it = self.tag.char_indices();
        let mut ix = None;
        let mut ch = None;
        for (iy, c) in it {
            if !c.is_numeric() {
                ch = Some((iy, c));
                break;
            }
            ix = Some(iy)
        }
        if let Some(ix) = ix {
            if let Ok(l) = self.tag[0..=ix].parse::<usize>() {
                if let Some((iy, c)) = ch {
                    if iy + 1 >= self.tag.len() {
                        self.done = true
                    }
                    self.tag = &self.tag[iy + 1..];
                    Some((c, l))
                } else {
                    None
                }
            } else {
                None
            }
        } else {
            None
        }
    }
}

fn process_primary_read(rec: &mut BamRec, st: &mut Stats, read_len: usize) {
    st.incr_mapq(rec.qual());
    if let Some(cigar) = rec.cigar() {
        // Find start and end points of aligned bases on read
        let reverse = (rec.flag() & BAM_FREVERSE) != 0;
        if let Some((start, stop)) = get_start_end(&cigar, reverse, read_len) {
            let mut v = vec![(start, stop)];
            // Look for supplementary mappings for this read
            if let Some(Ok(sa)) = rec.get_tag("SA", 'Z').map(|x| std::str::from_utf8(x)) {
                trace!("Found SA tag for read {}", rec.qname().unwrap());
                // Collect the start and end points of all mappings for this read
                for s in sa.split(';') {
                    if s.len() < 2 {
                        continue;
                    }
                    let fd: Vec<_> = s.split(',').collect();
                    if fd.len() == 6 {
                        if let Some(rev) = match fd[2] {
                            "+" => Some(false),
                            "-" => Some(true),
                            _ => {
                                warn!("Illegal SA tag {} - strand not + or -", s);
                                None
                            }
                        } {
                            let mut it = SaTagIter::new(fd[3]);
                            if let Some((a, b)) = sa_get_start_end(&mut it, rev, read_len) {
                                v.push((a, b))
                            } else {
                                warn!(
                                    "Illegal SA Tag {} for read {} (wrong size)",
                                    fd[3],
                                    rec.qname().unwrap()
                                )
                            }
                        }
                    } else {
                        warn!(
                            "Illegal SA tag {} (wrong number of fields {} instead of 6) >{}<",
                            s,
                            fd.len(),
                            sa
                        )
                    }
                }
                // Sort by starting point
                v.sort_unstable_by_key(|(x, _)| *x);
            }
            let n_splits = v.len();

            // Get total used bases excluding any overlaps
            let mut lim = v[0].1 + 1;
            let mut used = lim - v[0].0;
            for v1 in &v[1..] {
                if v1.1 >= lim {
                    used += v1.1 + 1 - lim;
                }
                lim = lim.max(v1.1 + 1);
            }
            trace!(
                "Read {}: len {}, n_splits {}, used {} ({}%)",
                rec.qname().unwrap(),
                read_len,
                n_splits,
                used,
                (used as f64) * 100.0 / (read_len as f64)
            );
            if used > read_len {
                for (a, b) in v.iter() {
                    println!(" -> {} {}", a, b);
                }
                panic!("{} {} {} {}", read_len, used, reverse, n_splits)
            }
            st.update_readlen_stats(read_len, used);
            st.incr_n_splits(n_splits);
        } else {
            warn!("Illegal CIGAR for read {}", rec.qname().unwrap())
        }
    }
}

fn process_coverage(
    rec: &mut BamRec,
    seq_qual: &SeqQual,
    cov: &mut [Option<usize>],
    begin: usize,
    min_qual: u8,
) {
    if let Some(cigar) = rec.cigar() {
        let mut x = rec.pos().expect("No start position for mapped read");
        let cov_len = cov.len();
        let mut sq = seq_qual.iter().map(|x| (*x & 3, *x >> 2));
        for elem in cigar.iter() {
            let l = elem.op_len() as usize;
            assert!(l > 0, "Zero length cigar element");
            match elem.op() {
                CigarOp::Match | CigarOp::Equal | CigarOp::Diff => {
                    for k in 0..l {
                        let (c, q) = sq
                            .next()
                            .expect("Mismatch between Cigar and sequence length");
                        if q >= min_qual && x >= begin {
                            if x - begin >= cov_len {
                                break;
                            }
                            if let Some(tp) = cov[x - begin].as_mut() {
                                *tp += 1
                            }
                        }
                        x += 1;
                    }
                }
                CigarOp::SoftClip | CigarOp::Ins => {
                    sq.nth(l - 1);
                }
                CigarOp::Del => x += l,
                _ => (),
            }
        }
    }
}

pub fn reader(
    cfg: &Config,
    ix: usize,
    in_file: &Path,
    tx: Sender<Stats>,
    rx: Receiver<(String, usize, &Vec<regions::Region>)>,
) {
    debug!("Starting reader thread {}", ix);

    let mut hts = input::open_input(
        in_file,
        true,
        cfg.reference(),
        cfg.threads_per_task().min(4),
    )
    .expect("Error opening input file in thread");
    let min_qual = cfg.min_qual().min(255) as u8;
    let min_mapq = cfg.min_mapq().min(255) as u8;

    let mut cov = Vec::new();
    let mut rec = BamRec::new().expect("Could not allocate new Bam Record");
    let mut st = Stats::new();
    let mut pair_warning = false;
    while let Ok((reg, ix, rvec)) = rx.recv() {
        let mappability = rvec[ix].mappability();
        let rlist = hts.make_region_list(&[&reg]);
        assert_eq!(rlist.len(), 1);
        let begin = rlist[0].begin() as usize;
        let end = rlist[0].end() as usize;
        let reg_len = end + 1 - begin;
        assert!(reg_len > 0);
        debug!(
            "Reader {} received {} {} - {} len {}",
            ix,
            reg,
            rlist[0].begin(),
            rlist[0].end(),
            reg_len
        );
        cov.clear();
        cov.reserve(reg_len as usize);

        if let Some(m) = mappability {
            if !m.is_empty() {
                let mut x = begin;
                for mv in m.iter() {
                    assert!(mv[0] >= x);
                    for _ in x..mv[0] {
                        cov.push(None)
                    }
                    for _ in mv[0]..=mv[1] {
                        cov.push(Some(0))
                    }
                    x = mv[1] + 1;
                }
                for _ in x..=rlist[0].end() as usize {
                    cov.push(None)
                }
            }
        } else {
            for _ in 0..reg_len {
                cov.push(Some(0))
            }
        };

        let mut rdr: HtsItrReader<BamRec> = hts.itr_reader(&rlist);
        while rdr.read(&mut rec).expect("Error reading from input file") {
            st.incr(StatType::Mappings);
            let flag = rec.flag();
            let chk_flg = |fg| (flag & fg) != 0;
            if chk_flg(BAM_FMUNMAP) {
                st.incr(StatType::Unmapped)
            } else {
                // Check if mapping could appear in another region
                let x = rec.pos().expect("Missing position for mapped read");
                let y = rec.endpos();
                // If mapping lies entirely within region then it can not appear in another region (as regions do not overlap)
                // If not then we check to see if this is the first region the mapping would appear in otherwise we skip
                if x < begin || y > end {
                    if let Err(i) = rvec.binary_search_by_key(&x, |r| r.start()) {
                        let j = if i > 0 {
                            // Check if the read overlaps previous region
                            if rvec[i - 1].end().map(|end| end >= x).unwrap_or(false) {
                                i - 1
                            } else {
                                i
                            }
                        } else {
                            i
                        };
                        assert!(j <= ix);
                        if j != ix {
                            continue;
                        }
                    }
                }
                if chk_flg(BAM_FPAIRED) {
                    warn!("Paired reads founds");
                    pair_warning = true;
                    st.incr(StatType::Paired)
                }

                st.incr(StatType::Mapped);
                let seq_qual = rec
                    .get_seq_qual()
                    .expect("Error getting sequence and qualities");
                let read_len = seq_qual.len();
                if chk_flg(BAM_FREVERSE) {
                    st.incr(StatType::Reversed)
                }
                if chk_flg(BAM_FSECONDARY) {
                    st.incr(StatType::Secondary)
                } else {
                    if !chk_flg(BAM_FSUPPLEMENTARY) {
                        //                       println!("{}\t{}\t{}", rec.qname().unwrap(), ix, flag);
                        st.incr(StatType::Reads)
                    }
                    if chk_flg(BAM_FDUP) {
                        st.incr(StatType::Duplicate)
                    } else if chk_flg(BAM_FQCFAIL) {
                        st.incr(StatType::QcFail)
                    } else {
                        if !chk_flg(BAM_FSUPPLEMENTARY) {
                            process_primary_read(&mut rec, &mut st, read_len);
                        }
                        if !cov.is_empty() && rec.qual() >= min_mapq {
                            process_coverage(&mut rec, &seq_qual, &mut cov, begin, min_qual)
                        }
                    }
                }
            }
        }
        for x in cov.iter().flatten() {
            st.incr_coverage(*x)
        }
        info!(
            "Reader {} read {} records {} reads from {}",
            ix,
            st.counts(StatType::Mappings),
            st.counts(StatType::Reads),
            reg
        );
    }
    tx.send(st).expect("Error sending Stats to collector");
    debug!("Terminating reader thread {}", ix);
}
