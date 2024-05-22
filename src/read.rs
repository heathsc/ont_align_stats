pub mod bisulfite;
pub mod coverage;
mod handle_sa_tag;
mod utils;

use std::{collections::HashMap, path::Path};

use crossbeam_channel::{Receiver, Sender, TryRecvError};
use r_htslib::{
    BamRec, HtsIterator, HtsItrReader, HtsRead, HtsThreadPool, Sequence, BAM_FDUP, BAM_FMUNMAP,
    BAM_FPAIRED, BAM_FPROPER_PAIR, BAM_FQCFAIL, BAM_FREAD1, BAM_FREAD2, BAM_FREVERSE,
    BAM_FSECONDARY, BAM_FSUPPLEMENTARY, BAM_FUNMAP,
};

use crate::{
    config::Config,
    input,
    regions::{self, Region, Regions, TaskRegion},
    stats::{ReadType, StatType, Stats},
};

use coverage::{process_coverage, Coverage};
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
fn process_primary_read(rec: &mut BamRec, st: &mut Stats, read_len: usize) {
    st.incr_mapq(rec.qual());
    if let Some(cigar) = rec.cigar() {
        // Find start and end points of aligned bases on read
        let reverse = (rec.flag() & BAM_FREVERSE) != 0;
        if let Some((start, stop)) = get_start_end(&cigar, reverse, read_len) {
            // Make list of all alignments from this read, starting with the primary alignment
            let mut v = vec![(start, stop)];

            // Look for supplementary mappings (int the SA tag) for this read
            // Add start,end points (on read) of supplementary mappings to v
            // On return v is sorted on starting position
            process_sa_tag(rec, read_len, &mut v);

            // This is a least 1 as we include the primary alignment
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
            st.update_readlen_stats(read_len, used);
            st.incr_n_splits(n_splits);
        } else {
            warn!("Illegal CIGAR for read {}", rec.qname().unwrap())
        }
    }
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
    rx: Receiver<(String, usize, Option<&Vec<Region>>, Option<Sequence>)>,
    tpool: Option<&HtsThreadPool>,
) {
    debug!("Starting reader thread {}", ix);

    let mut hts = input::open_input(in_file, false, cfg.reference(), cfg.hts_threads(), tpool)
        .expect("Error opening input file in thread");
    let min_mapq = cfg.min_mapq().min(255) as u8;

    let mut cov = Coverage::new();

    let mut rec = BamRec::new().expect("Could not allocate new Bam Record");
    let mut st = Stats::new();
    while let Ok((reg, reg_ix, rvec, seq)) = rx.recv() {
        let mappability = rvec.and_then(|v| v[reg_ix].mappability());
        let rlist = hts.make_region_list(&[&reg]);
        assert_eq!(rlist.len(), 1, "Empty region list for {}", reg);
        let begin = (rlist[0].begin() + 1) as usize;
        let end = rlist[0].end() as usize;
        if end >= begin {
            let rf = seq.as_ref().map(|s| {
                s.get_seq(begin, end)
                    .expect("Error getting reference sequence for region")
            });
            cov.reset(begin, end, mappability, rf);
        }
        let mut rdr: HtsItrReader<BamRec> = hts.itr_reader(&rlist);
        while rdr.read(&mut rec).expect("Error reading from input file") {
            let flag = rec.flag();
            let chk_flg = |fg| (flag & fg) != 0;

            let paired = chk_flg(BAM_FPAIRED);
            if chk_flg(BAM_FUNMAP) {
                st.incr(StatType::Unmapped);
                st.incr_n(StatType::TotalBases, rec.l_qseq() as usize);
                st.incr(StatType::Reads);
                if paired {
                    st.incr(StatType::TotalPairs)
                }
            } else {
                let rvec = rvec.expect("Empty region vec for mapped region");
                // Check if mapping could appear in another region
                let x = rec.pos().expect("Missing position for mapped read") + 1;
                let y = rec.endpos() + 1;
                // If mapping lies entirely within region then it can not appear in another region (as regions do not overlap)
                // If not then we check to see if this is the first region the mapping would appear in otherwise we skip
                let primary_region = if x < begin || y > end {
                    let v = regions::find_overlapping_regions(rvec, x, y);
                    assert!(!v.is_empty());
                    v[0] == reg_ix
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
                        let rl = rec.l_qseq() as usize;
                        process_primary_read(&mut rec, &mut st, rl);
                    }
                }
                if rec.qual() >= min_mapq
                    && !chk_flg(BAM_FSECONDARY | BAM_FDUP | BAM_FQCFAIL | BAM_FUNMAP)
                {
                    let sq = rec
                        .get_seq_qual()
                        .expect("Error getting sequence and qualities");
                    if let Some(rd_type) = get_read_type(&rec) {
                        let bs_strand = if cfg.bisulfite() {
                            bisulfite::get_bs_strand(&rec).map(|s| {
                                st.incr_bisulfite_strand(s);
                                s
                            })
                        } else {
                            None
                        };
                        process_coverage(&mut rec, &sq, &mut cov, rd_type, bs_strand, &mut st, cfg);
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
}

pub struct ReadHandlerData {
    rec: BamRec,
    reg_ix: usize, // Region index within task
    primary_region: bool,
    // Subsequent regions for this record
    task_info: Option<Vec<(usize, usize, Sender<Vec<ReadHandlerData>>)>>, // task index, region index within task, sender to task
}

// Read handler for non-indexed files
pub fn read_handler(
    cfg: &Config,
    ix: usize,
    region_list: &[TaskRegion],
    tx: Sender<Stats>,
    rx: Receiver<Vec<ReadHandlerData>>,
    bam_tx: Sender<Vec<BamRec>>,
) {
    debug!("Starting read handler thread {}", ix);

    let buf_size = cfg.bam_rec_thread_buffer();
    let mut brec_store = Vec::with_capacity(buf_size);
    let min_mapq = cfg.min_mapq().min(255) as u8;
    let mut st = Stats::new();
    let mut cov_vec = Vec::new();
    for tr in region_list.iter() {
        for (reg, _, seq) in tr.regions() {
            let rf = seq.as_ref().map(|s| {
                s.get_seq(reg.start().max(1), reg.end().unwrap())
                    .expect("Error getting reference sequence for region")
            });
            let mut cov = Coverage::new();
            cov.reset(reg.start(), reg.end().unwrap(), reg.mappability(), rf);
            cov_vec.push(cov);
        }
    }

    loop {
        let mut rd_blk = {
            match rx.try_recv() {
                Ok(rd) => rd,
                Err(TryRecvError::Empty) => {
                    if !brec_store.is_empty() {
                        let _ = bam_tx.send(brec_store);
                        brec_store = Vec::with_capacity(buf_size)
                    }
                    if let Ok(rd) = rx.recv() {
                        rd
                    } else {
                        break;
                    }
                }
                Err(_) => break,
            }
        };
        for rd in rd_blk.drain(..) {
            let ReadHandlerData {
                mut rec,
                mut reg_ix,
                primary_region,
                mut task_info,
            } = rd;
            let flag = rec.flag();
            let chk_flg = |fg| (flag & fg) != 0;

            if primary_region && !chk_flg(BAM_FSUPPLEMENTARY) {
                // Check primary read
                let rl = rec.l_qseq() as usize;
                process_primary_read(&mut rec, &mut st, rl);
            }
            let push_brec = |r: BamRec, mut br: Vec<BamRec>| {
                br.push(r);
                if br.len() >= buf_size {
                    let _ = bam_tx.send(br);
                    Vec::with_capacity(buf_size)
                } else {
                    br
                }
            };
            if rec.qual() >= min_mapq {
                let seq_qual = rec.get_seq_qual().expect("No sequence/quality data");
                if let Some(rd_type) = get_read_type(&rec) {
                    let bs_strand = if cfg.bisulfite() {
                        bisulfite::get_bs_strand(&rec).map(|s| {
                            st.incr_bisulfite_strand(s);
                            s
                        })
                    } else {
                        None
                    };
                    loop {
                        // Check coverage
                        process_coverage(
                            &mut rec,
                            &seq_qual,
                            &mut cov_vec[reg_ix],
                            rd_type,
                            bs_strand,
                            &mut st,
                            cfg,
                        );

                        // Remove next region if exists
                        if let Some((task_ix, r_ix, tx, t_info)) = task_info.map(|mut v| {
                            let (task_ix, r_ix, tx) = v.pop().unwrap();
                            if v.is_empty() {
                                (task_ix, r_ix, tx, None)
                            } else {
                                (task_ix, r_ix, tx, Some(v))
                            }
                        }) {
                            reg_ix = r_ix;
                            task_info = t_info;
                            if task_ix != ix {
                                let rd = ReadHandlerData {
                                    rec,
                                    reg_ix,
                                    primary_region: false,
                                    task_info,
                                };
                                tx.send(vec![rd])
                                    .expect("Error sending data to read handler");
                                break;
                            }
                        } else {
                            brec_store = push_brec(rec, brec_store);
                            break;
                        }
                    }
                } else {
                    brec_store = push_brec(rec, brec_store);
                }
            } else {
                brec_store = push_brec(rec, brec_store);
            }
        }
    }
    debug!("Read handle thread {} - updating stats", ix);
    if !brec_store.is_empty() {
        let _ = bam_tx.send(brec_store);
    }
    for cov in cov_vec {
        cov.update_stats(&mut st)
    }
    tx.send(st).expect("Error sending Stats to collector");
    debug!("Terminating read handler thread {}", ix);
}

// Read in BAM records from non-indexed file with multithreading
pub fn read_input_mt(
    cfg: &Config,
    in_file: &Path,
    regions: &Regions,
    task_regions: &[Vec<TaskRegion>],
    tx_vec: Vec<Sender<Vec<ReadHandlerData>>>,
    bam_rx: Receiver<Vec<BamRec>>,
    stats_tx: Sender<Stats>,
) -> anyhow::Result<()> {
    let nreg = regions.len();
    let n_tasks = cfg.n_tasks().min(nreg);
    let mut st = Stats::new();
    let n_brecs = n_tasks * 2 * cfg.bam_rec_thread_buffer();
    let mut brec_buf = Vec::with_capacity(n_brecs);
    for _ in 0..n_brecs {
        brec_buf.push(BamRec::new()?)
    }

    let block_size = cfg.non_index_buffer_size();
    let mut task_blocks: Vec<Vec<ReadHandlerData>> = (0..n_tasks)
        .map(|_| Vec::with_capacity(block_size))
        .collect();
    // Make hashtable so we can go from (ctg, idx) to (task_idx, task_region_idx, task_channel)
    let mut reg_task_hash = HashMap::new();
    for (task_ix, tv) in task_regions.iter().enumerate() {
        let mut tix: usize = 0;
        for tr in tv.iter() {
            for (_, reg_ix, _) in tr.regions() {
                let tx = tx_vec[task_ix].clone();
                reg_task_hash.insert((tr.ctg(), *reg_ix), (task_ix, tix, tx));
                tix += 1;
            }
        }
    }

    let mut pair_warning = true;
    let min_mapq = cfg.min_mapq().min(255) as u8;

    let mut hts = input::open_input(in_file, true, cfg.reference(), cfg.hts_threads(), None)?;

    loop {
        // If we have no empty bam rec then we wait until one appears
        if brec_buf.is_empty() {
            let mut v = bam_rx.recv()?;
            brec_buf.append(&mut v)
        }

        // Collect any other vectors of bam recs that are waiting
        for mut v in bam_rx.try_iter() {
            brec_buf.append(&mut v);
        }
        let mut rec = brec_buf.pop().unwrap();

        // Get next read if it exists
        if !rec.read(&mut hts)? {
            break;
        }

        let flag = rec.flag();
        let chk_flg = |fg| (flag & fg) != 0;

        let paired = chk_flg(BAM_FPAIRED);
        if chk_flg(BAM_FUNMAP) {
            st.incr(StatType::Unmapped);
            st.incr_n(StatType::TotalBases, rec.l_qseq() as usize);
            st.incr(StatType::Reads);
            if paired {
                st.incr(StatType::TotalPairs)
            }
            brec_buf.push(rec);
        } else if let Some(ctg) =
            regions.tid2ctg(rec.tid().expect("Missing contig for mapped read"))
        {
            let x = rec.pos().expect("Missing position for mapped read") + 1;
            let y = rec.endpos() + 1;
            let rvec = regions.ctg_regions(ctg).expect("Unknown contig");
            let regs = regions::find_overlapping_regions(rvec, x, y);
            if regs.is_empty() {
                // Read does not overlap requested areas
                brec_buf.push(rec);
            } else {
                if chk_flg(BAM_FPAIRED) && !pair_warning {
                    warn!("Paired reads founds");
                    pair_warning = true;
                }
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
                if !chk_flg(BAM_FSECONDARY | BAM_FDUP | BAM_FQCFAIL | BAM_FUNMAP) {
                    let mut v: Vec<_> = regs
                        .iter()
                        .rev()
                        .copied()
                        .map(|i| {
                            let (task_ix, reg_ix, tx) =
                                reg_task_hash.get(&(ctg, i)).expect("No task for region");
                            (*task_ix, *reg_ix, tx.clone())
                        })
                        .collect();

                    let (task_ix, task_reg, tx) = v.pop().unwrap();

                    let (task_info, sflag) = if rec.qual() >= min_mapq {
                        // Need to process all regions
                        if v.is_empty() {
                            (None, true)
                        } else {
                            (Some(v), true)
                        }
                    } else if !chk_flg(BAM_FSUPPLEMENTARY) {
                        // Only process the first region
                        (None, true)
                    } else {
                        // Do not process
                        (None, false)
                    };
                    if sflag {
                        let rd = ReadHandlerData {
                            rec,
                            reg_ix: task_reg,
                            primary_region: true,
                            task_info,
                        };
                        task_blocks[task_ix].push(rd);
                        if task_blocks[task_ix].len() >= block_size {
                            let tb = std::mem::replace(
                                &mut task_blocks[task_ix],
                                Vec::with_capacity(block_size),
                            );
                            tx.send(tb).expect("Error sending read to read handler");
                        }
                    } else {
                        brec_buf.push(rec)
                    }
                } else {
                    brec_buf.push(rec)
                }
            }
        } else {
            brec_buf.push(rec);
        }
    }
    for (ix, tb) in task_blocks.drain(..).enumerate() {
        if !tb.is_empty() {
            tx_vec[ix]
                .send(tb)
                .expect("Error sending read to read handler");
        }
    }
    stats_tx.send(st).expect("Error sending Stats to collector");
    Ok(())
}
