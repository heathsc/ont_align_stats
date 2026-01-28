pub mod coverage;
mod handle_sa_tag;
mod utils;

use std::{path::Path, sync::Arc};

use anyhow::Context;
use crossbeam_channel::{Receiver, Sender};
use m_htslib::{
    faidx::Sequence,
    hts::{HtsThreadPool, ReadRec},
    region::{Reg, RegCoords},
    sam::{
        BamRec, SamHdr, SamReader,
        record::bam1::{
            BAM_FDUP, BAM_FMUNMAP, BAM_FPAIRED, BAM_FPROPER_PAIR, BAM_FQCFAIL, BAM_FREAD1,
            BAM_FREAD2, BAM_FREVERSE, BAM_FSECONDARY, BAM_FSUPPLEMENTARY, BAM_FUNMAP,
        },
    },
};

use crate::{
    Config, input,
    stats::{ReadType, StatType, Stats},
};

use coverage::{Coverage, process_coverage};
use handle_sa_tag::*;
use utils::*;

fn get_read_type(rec: &BamRec) -> Option<ReadType> {
    let fg = rec.flag();
    if (fg & BAM_FUNMAP) == 0 {
        let reverse = (fg & BAM_FREVERSE) != 0;
        if (fg & BAM_FPAIRED) == 0 {
            // Non-paired
            Some(if reverse {
                ReadType::Reverse
            } else {
                ReadType::Forwards
            })
        } else {
            // paired record
            match fg & (BAM_FREAD1 | BAM_FREAD2) {
                BAM_FREAD1 => Some(if reverse {
                    ReadType::ReverseRead1
                } else {
                    ReadType::ForwardsRead1
                }),
                BAM_FREAD2 => Some(if reverse {
                    ReadType::ReverseRead2
                } else {
                    ReadType::ForwardsRead2
                }),
                _ => None,
            }
        }
    } else {
        None
    }
}

/// Collect read level stats that should only be performed on the primary mapping for a read
/// (i.e., not a secondary or supplementary mapping)
///
/// Because a read can overlap multiple regions, this function is only run on the first (lowest coordinate)
/// region overlapping this read.
///
/// We want to calculate the number of bases that are 'used' (i.e. that are not skips or insertions)
/// from the read. Since a read can have supplementary alignments, we lok for the SA tag to get the
/// same information on any other alignments for the read.  This requires making a list of all alignments
/// for the read, sorting on start order, and then removing any overlaps to ensure bases are only counted once
fn process_primary_read(
    cfg: &Config,
    rec: &mut BamRec,
    st: &mut Stats,
    read_len: usize,
) -> anyhow::Result<()> {
    let map_q = rec.mapq();
    let min_bq = cfg.min_qual();
    st.incr_mapq(map_q);
    let mapq_ok = map_q >= cfg.min_mapq();
    if let Some(cigar) = rec.cigar() {
        // Find start and end points of aligned bases on read
        let reverse = (rec.flag() & BAM_FREVERSE) != 0;
        if let Some((start, stop)) = get_start_end(cigar, reverse, read_len) {
            // Make list of all alignments from this read, starting with the primary alignment
            let mut v = vec![(start, stop)];

            // Look for supplementary mappings (int the SA tag) for this read
            // Add start,end points (on read) of supplementary mappings to v
            // On return v is sorted on starting position
            process_sa_tag(rec, read_len, &mut v)
                .with_context(|| "Error processing SA tag for read")?;

            // This is a least 1 as we include the primary alignment
            let n_splits = v.len();

            let get_used = |(p0, p1)| {
                if mapq_ok {
                    rec.qual()
                        .skip(p0)
                        .take(p1 + 1 - p0)
                        .filter(|q| *q >= min_bq)
                        .count()
                } else {
                    0
                }
            };

            // Get total used bases excluding any overlaps
            let mut lim = v[0].1 + 1;
            let mut mapped = lim - v[0].0;
            let mut used = get_used(v[0]);

            for v1 in &v[1..] {
                if v1.1 >= lim {
                    mapped += v1.1 + 1 - lim;
                    used += get_used((lim, v1.1))
                }
                lim = lim.max(v1.1 + 1);
            }
            let unique = if mapq_ok { mapped } else { 0 };

            trace!(
                "Read {:?}: len {}, n_splits {}, mapped {mapped}, unique {unique} used {used} ({}%)",
                rec.qname().unwrap(),
                read_len,
                n_splits,
                (used as f64) * 100.0 / (read_len as f64)
            );
            st.update_readlen_stats(read_len, mapped, unique, used);
            st.incr_n_splits(n_splits);
        } else {
            warn!("Illegal CIGAR for read {:?}", rec.qname().unwrap())
        }
    }
    Ok(())
}

fn update_flag_stats(flag: u16, pair_same_ctg: bool, tmpl_neg: bool, st: &mut Stats) {
    let chk_flg = |fg| (flag & fg) != 0;
    st.incr(StatType::Mappings);
    if chk_flg(BAM_FREVERSE) {
        st.incr(StatType::Reversed)
    }
    if chk_flg(BAM_FSECONDARY) {
        st.incr(StatType::Secondary)
    } else if !chk_flg(BAM_FSUPPLEMENTARY) {
        st.incr(StatType::Reads);
        if chk_flg(BAM_FDUP) {
            st.incr(StatType::Duplicate)
        }
        if chk_flg(BAM_FQCFAIL) {
            st.incr(StatType::QcFail)
        }
        if chk_flg(BAM_FPAIRED) {
            st.incr(StatType::TotalPairs);
            if chk_flg(BAM_FPROPER_PAIR) {
                st.incr(StatType::CorrectPairs);
            }
            if chk_flg(BAM_FMUNMAP) {
                st.incr(StatType::MateUnmapped)
            } else if pair_same_ctg {
                let orientation = flag & 0xf0;
                st.incr(match orientation {
                    0x50 | 0xa0 => StatType::OrientationRF,
                    0x60 | 0x90 => StatType::OrientationFR,
                    0x40 | 0x80 => StatType::OrientationFF,
                    0x70 | 0xb0 => StatType::OrientationRR,
                    _ => StatType::IllegalOrientation,
                });
                match orientation {
                    0x50 | 0x90 => {
                        if tmpl_neg {
                            if !chk_flg(BAM_FPROPER_PAIR) {
                                st.incr(StatType::BadTemplateLength)
                            }
                            st.incr(StatType::ConvergentPair)
                        } else {
                            st.incr(StatType::DivergentPair)
                        }
                    }
                    0x60 | 0xa0 => {
                        if tmpl_neg {
                            st.incr(StatType::DivergentPair)
                        } else {
                            if !chk_flg(BAM_FPROPER_PAIR) {
                                st.incr(StatType::BadTemplateLength)
                            }
                            st.incr(StatType::ConvergentPair)
                        }
                    }
                    _ => (),
                }
            } else {
                st.incr(StatType::DifferentContigs)
            }
        }
    }
}

// Read thread for indexed files
pub fn reader(
    cfg: &Config,
    ix: usize,
    in_file: &Path,
    tx: Sender<Stats>,
    rx: Receiver<(Reg, Option<Reg>, Option<Arc<Sequence>>)>,
    tpool: Option<&HtsThreadPool>,
) -> anyhow::Result<()> {
    debug!("Starting reader thread {}", ix);

    let mut hts = input::open_input(in_file, cfg.reference(), cfg.hts_threads(), tpool)
        .expect("Error opening input file in thread");

    if let Some(tp) = tpool.as_ref() {
        hts.set_thread_pool(tp)
            .with_context(|| "Error setting thread pool for input file")?;
    }

    let hdr = SamHdr::read(&mut hts).with_context(|| {
        format!(
            "Could not read SAM/BAM/CRAM header for input file {}",
            in_file.display()
        )
    })?;

    let min_mapq = cfg.min_mapq();

    let mut cov = Coverage::new();

    let mut rec = BamRec::new();
    let mut st = Stats::new();
    while let Ok((reg, prev_reg, seq)) = rx.recv() {
        assert!(!reg.is_all(), "Should not be getting the All Reg type here");

        debug!("Processing region {reg}. Sequence present: {}", seq.is_some());
        
        let rlist = if reg.is_unmapped() {
            None
        } else {
            let rlist = reg
                .make_htslib_region(&hdr)
                .with_context(|| format!("Bad region {reg}"))?;

            let begin = rlist.start() as usize;
            let end = rlist.end() as usize;
            let rf = seq.as_ref().map(|s| {
                s.get_seq(begin + 1, end + 1)
                    .expect("Error getting reference sequence for region")
            });
            cov.reset(begin, end, rf);
            Some(rlist)
        };

        let sam_reader = SamReader::new(&mut hts, &hdr);
        let mut rdr = sam_reader
            .region_iter(&reg)
            .with_context(|| format!("Can't make region iterator for {reg}"))?;

        while rdr
            .read_rec(&mut rec)
            .with_context(|| format!("Error reading from input file for region {reg}"))?
            .is_some()
        {
            let flag = rec.flag();
            let chk_flg = |fg| (flag & fg) != 0;

            let paired = chk_flg(BAM_FPAIRED);
            if chk_flg(BAM_FUNMAP) {
                st.incr(StatType::Unmapped);
                st.incr_n(StatType::TotalBases, rec.seq_len());
                st.incr(StatType::Reads);
                if paired {
                    st.incr(StatType::TotalPairs)
                }
            } else {
                let rl = rlist
                    .as_ref()
                    .expect("Should not be getting upmapped reads here");
                let x = rec.pos().expect("Missing position for mapped read");
                let y = rec.endpos();

                // If mapping lies entirely within region then it can not appear in another region (as regions do not overlap)
                // If not then we check to see if this is the first region the mapping would appear in otherwise we skip.
                // We just need to check if the read overlaps the previous region in the contig (if it exists)
                let primary_region = if (x < rl.start() || y > rl.end())
                    && let Some(preg) = prev_reg.as_ref()
                {
                    let prev_end = preg.coords().1.expect("Cannot have an open interval here!");
                    x > prev_end
                } else {
                    true
                };
                if primary_region {
                    let tid = rec.tid().expect("Invalid contig ID for mapped read");
                    let same_contig = rec.mtid().map(|x| x == tid).unwrap_or(false);
                    let tmpl_neg = if paired {
                        let tlen = rec.template_len();
                        if tlen >= 0 && chk_flg(BAM_FPROPER_PAIR) {
                            st.incr_template_len(tlen as usize)
                        }
                        tlen < 0
                    } else {
                        false
                    };
                    update_flag_stats(flag, same_contig, tmpl_neg, &mut st);
                    if !chk_flg(
                        BAM_FSECONDARY | BAM_FDUP | BAM_FQCFAIL | BAM_FUNMAP | BAM_FSUPPLEMENTARY,
                    ) {
                        let rl = rec.seq_len();
                        process_primary_read(cfg, &mut rec, &mut st, rl)
                            .with_context(|| format!("Error processing read {:?}", rec.qname()))?;
                    }
                }
                if rec.mapq() >= min_mapq
                    && !chk_flg(BAM_FSECONDARY | BAM_FDUP | BAM_FQCFAIL | BAM_FUNMAP)
                {
                    let sq = rec.seq_qual();
                    if let Some(rd_type) = get_read_type(&rec) {
                        process_coverage(&rec, sq, &mut cov, rd_type, &mut st, cfg);
                    } else {
                        warn!("Illegal read type")
                    }
                }
            }
        }
        cov.update_stats(&mut st);
        trace!(
            "Reader {} read {} records {} reads from {}",
            ix,
            st.counts(StatType::Mappings),
            st.counts(StatType::Reads),
            reg
        );
    }
    tx.send(st).expect("Error sending Stats to collector");
    debug!("Terminating reader thread {}", ix);
    Ok(())
}
