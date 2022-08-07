use std::path::Path;

use crossbeam_channel::{Receiver, Sender};
use r_htslib::*;

use crate::{
    config::Config,
    input,
    stats::{StatType, Stats},
};

fn get_start_end(cigar: &Cigar) -> (usize, usize) {
    let mut first = 0;
    let mut last = 0;
    let mut started = false;
    let mut x = 0;
    for elem in cigar.iter() {
        match elem.op() {
            CigarOp::HardClip | CigarOp::SoftClip | CigarOp::Ins => x += elem.op_len(),
            CigarOp::Match | CigarOp::Equal | CigarOp::Diff => {
                if !started {
                    started = true;
                    first = x;
                }
                x += elem.op_len();
                last = x;
            }
            _ => {}
        }
    }
    (first as usize, last as usize)
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
fn process_primary_read(rec: &mut BamRec, st: &mut Stats) {
    st.incr_mapq(rec.qual());
    if let Some(cigar) = rec.cigar() {
        // Find start and end points of aligned bases on read
        let reverse = (rec.flag() & BAM_FREVERSE) != 0;
        let (start, stop) = {
            let (x, y) = get_start_end(&cigar);
            if reverse {
                (y, x)
            } else {
                (x, y)
            }
        };

        // Look for supplementary mappings for this read
        if let Some(Ok(sa)) = rec.get_tag("SA", 'Z').map(|x| std::str::from_utf8(x)) {
            debug!("Found SA tag for read {}", rec.qname().unwrap());
            // Collect the start and end points of all mappings for this read
            let mut v = vec![(start, stop)];
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
                        println!("SA Tag: {}", fd[3]);
                        for (c, l) in it {
                            println!("{} {}", c, l);
                        }
                        println!();
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
        }
    }
}

pub fn reader(cfg: &Config, ix: usize, in_file: &Path, tx: Sender<usize>, rx: Receiver<String>) {
    debug!("Starting reader thread {}", ix);

    let mut hts = input::open_input(
        in_file,
        true,
        cfg.reference(),
        cfg.threads_per_task().min(4),
    )
    .expect("Error opening input file in thread");

    let mut rec = BamRec::new().expect("Could not allocate new Bam Record");
    let mut s = kstring_t::new();
    let mut st = Stats::new();
    let mut pair_warning = false;
    while let Ok(reg) = rx.recv() {
        debug!("Reader {} received {}", ix, reg);
        let rlist = hts.make_region_list(&[&reg]);
        let mut rdr: HtsItrReader<BamRec> = hts.itr_reader(&rlist);
        while rdr.read(&mut rec).expect("Error reading from input file") {
            //            rec.format(rdr.header().unwrap(), &mut s)?;
            // println!("Reader {} Read {}", ix, rec.qname().expect("No query"));
            st.incr(StatType::Mappings);
            let flag = rec.flag();
            let chk_flg = |fg| (flag & fg) != 0;
            if chk_flg(BAM_FPAIRED) {
                if !pair_warning {
                    warn!("Paired reads founds");
                    pair_warning = true;
                }
                st.incr(StatType::Paired)
            }
            if !chk_flg(BAM_FSUPPLEMENTARY) {
                st.incr(StatType::Reads)
            }
            if chk_flg(BAM_FMUNMAP) {
                st.incr(StatType::Unmapped)
            } else {
                if chk_flg(BAM_FREVERSE) {
                    st.incr(StatType::Reversed)
                }
                if chk_flg(BAM_FSECONDARY) {
                    st.incr(StatType::Secondary)
                } else if chk_flg(BAM_FDUP) {
                    st.incr(StatType::Duplicate)
                } else if chk_flg(BAM_FQCFAIL) {
                    st.incr(StatType::QcFail)
                } else if !chk_flg(BAM_FSUPPLEMENTARY) {
                    process_primary_read(&mut rec, &mut st);
                }
            }
        }
        info!(
            "Reader {} read {} records {} reads from {}",
            ix,
            st.counts(StatType::Mappings),
            st.counts(StatType::Reads),
            reg
        );
    }
    debug!("Terminating reader thread {}", ix);
}
