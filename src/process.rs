use std::{fmt::Display, fs, sync::Arc, thread};

use anyhow::Context;
use compress_io::compress::CompressIo;
use compress_io::compress_type::CompressType;
use crossbeam_channel::bounded;
use indexmap::IndexMap;
use m_htslib::{
    faidx::Faidx,
    hts::HtsThreadPool,
    region::{Reg, RegCtgName},
};

use serde::{
    self, Serialize,
    ser::{SerializeMap, Serializer},
};

use crate::{CompressOpt, Config, collect, metadata, read, stats::Stats};

fn serialize_im<S, K, V>(im: &IndexMap<K, V>, serializer: S) -> Result<S::Ok, S::Error>
where
    S: Serializer,
    K: Display + Serialize,
    V: Display + Serialize,
{
    let l: usize = im.len();
    let mut map = serializer.serialize_map(Some(l))?;
    for (key, val) in im.iter() {
        map.serialize_entry(key, val)?;
    }
    map.end()
}

#[derive(Serialize)]
struct Output {
    #[serde(serialize_with = "serialize_im")]
    metadata: IndexMap<&'static str, String>,
    stats: Stats,
}

pub fn process(cfg: Config) -> anyhow::Result<()> {
    let mut metadata = IndexMap::new();
    metadata::collect_starting_metadata(&mut metadata);

    let nreg = cfg.region_list().regions().count();

    let n_tasks = cfg.n_tasks().min(nreg);
    let nthr = cfg.hts_threads();
    debug!(
        "No. regions: {}, No. read tasks: {}, threads per task: {}",
        nreg, n_tasks, nthr
    );
    let cfg_ref = &cfg;
    let ifile = cfg.input();

    let tpool = HtsThreadPool::init(cfg.hts_threads());
    let tpool_ref = tpool.as_ref();

    let mut faidx = if let Some(p) = cfg.reference() {
        trace!("Try to open reference {}", p.display());
        let f = Faidx::load_or_create(p);
        if f.is_ok() {
            trace!("Reference index loaded successfully");
        } else {
            warn!("Could not load reference");
        }
        let mut fx = f.ok();
        if let Some(tp) = tpool_ref
            && let Some(fai) = fx.as_mut()
        {
            fai.set_thread_pool(tp);
        }
        fx
    } else {
        None
    };

    metadata.insert("n_tasks", format!("{}", n_tasks));
    metadata.insert("threads_per_reader", format!("{}", nthr));
    let min_qual = cfg.min_qual();
    let min_mapq = cfg.min_mapq();
    metadata.insert("min_qual", format!("{}", min_qual));
    metadata.insert("min_mapq", format!("{}", min_mapq));
    if let Some(l) = cfg.min_read_len() {
        metadata.insert("min_read_len", format!("{l}"));
    }
    if let Some(l) = cfg.max_read_len() {
        metadata.insert("max_read_len", format!("{l}"));
    }
    debug!("Processing input");

    let mut error = None;
    // Create thread scope so that we can share references across threads
    let st = thread::scope(|scope| {
        let (collector_tx, collector_rx) = bounded(n_tasks * 4);
        let collector_task = scope.spawn(|| collect::collector(collector_rx));

        // Spawn threads
        let (region_tx, region_rx) = bounded(n_tasks * 4);
        let mut tasks: Vec<_> = (0..n_tasks)
            .map(|ix| {
                let rx = region_rx.clone();
                let tx = collector_tx.clone();
                scope.spawn(move || read::reader(cfg_ref, ix, ifile, tx, rx, tpool_ref))
            })
            .collect();
        drop(collector_tx);

        let mut prev_reg: Option<Reg> = None;
        let mut seq = None;
        // Send mapped regions to readers
        for reg in cfg.region_list().regions() {
            trace!("Sending region {reg} for processing");

            let preg = if let Some(pctg) = prev_reg.and_then(|r| r.reg_contig())
                && let Some(ctg) = reg.reg_contig()
                && pctg == ctg
            {
                prev_reg
            } else {
                None
            };

            // Load the sequence for the contig if not already done so
            if prev_reg.is_none() {
                if reg.reg_contig().is_some() {
                    seq = None
                } else {
                    seq = if let Some(f) = faidx.as_mut() {
                        match f.fetch_seq(
                            reg.to_cstr().expect("Error getting contig name"),
                            0,
                            None,
                        ) {
                            Ok(sq) => {
                                trace!("Sequence for {} loaded successfully", reg.contig_name());
                                Some(Arc::new(sq))
                            }
                            Err(e) => {
                                warn!(
                                    "Error loading sequence for contig {}: {e}",
                                    reg.contig_name()
                                );
                                None
                            }
                        }
                    } else {
                        None
                    };
                }
            }

            prev_reg = Some(reg);
            let s = seq.as_ref().map(|sq| sq.clone());

            region_tx
                .send((reg, preg, s))
                .expect("Error sending region message");
        }
        drop(region_tx);

        // Wait for readers to finish
        for jh in tasks.drain(..) {
            if let Err(e) = jh.join().expect("Error joining read threads")
                && error.is_none()
            {
                error = Some(Box::new(e))
            }
        }
        // Wait for collector to finish
        collector_task.join()
    })
    .expect("Processing error");

    if let Some(e) = error {
        return Err(anyhow!("Processing error: {e:?}"))
    }
    // Create output dir if necessary
    if let Some(dir) = cfg.dir() {
        fs::create_dir_all(dir).with_context(|| "Could not create output directory")?;
        metadata.insert("output_dir", dir.to_string_lossy().to_string());
    }

    // Create results file name
    let fname = {
        let mut d = cfg.dir().map(|p| p.to_owned()).unwrap_or_default();
        d.push(cfg.prefix());
        d.set_extension("json");
        d
    };

    metadata.insert("json_output_name", fname.to_string_lossy().to_string());

    // Open output file

    // Set compress type
    let ct = match cfg.compress() {
        CompressOpt::None => CompressType::NoFilter,
        CompressOpt::Gzip => CompressType::Gzip,
        CompressOpt::Bzip2 => CompressType::Bzip2,
        CompressOpt::Zstd => CompressType::Zstd,
        CompressOpt::Xz => CompressType::Xz,
    };

    let mut wrt = CompressIo::new()
        .path(fname)
        .ctype(ct)
        .bufwriter()
        .with_context(|| "Could not open output file")?;

    //    let mut wrt = fs::File::create(&fname)
    //        .with_context(|| format!("Could not open output file {}", fname.display()))?;

    // Get elapsed time
    let mut t = cfg.elapsed().as_secs_f64();
    let mut s = String::new();
    if t > 24.0 * 3600.0 {
        let days = (t / (24.0 * 3600.0)).floor() as u32;
        s.push_str(format!("{}d", days).as_str());
        t -= (days as f64) * 24.0 * 3600.0;
    }
    if t > 3600.0 {
        let hours = (t / 3600.0).floor() as u32;
        s.push_str(format!("{}h", hours).as_str());
        t -= (hours as f64) * 3600.0;
    }
    if t > 60.0 {
        let mins = (t / 60.0).floor() as u32;
        s.push_str(format!("{}m", mins).as_str());
        t -= (mins as f64) * 60.0;
    }
    s.push_str(format!("{:.3}s", t).as_str());
    metadata.insert("elapsed_time", s);

    let out = Output {
        metadata,
        stats: st,
    };
    // Generate JSON
    serde_json::to_writer_pretty(&mut wrt, &out)
        .with_context(|| "Could not write out JSON stats file")
}
