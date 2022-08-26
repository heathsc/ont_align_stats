use std::{
    collections::HashMap,
    fmt::Display,
    fs,
    path::{Path, PathBuf},
    rc::Rc,
    thread,
};

use anyhow::Context;
use crossbeam_channel::{bounded, unbounded};
use indexmap::IndexMap;

use serde::{
    self,
    ser::{SerializeMap, Serializer},
    Serialize,
};

use crate::{
    collect,
    config::Config,
    read,
    regions::{Regions, TaskRegion},
    stats::Stats,
};

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

pub fn process(
    cfg: Config,
    input: PathBuf,
    regions: Regions,
    mut metadata: IndexMap<&'static str, String>,
) -> anyhow::Result<()> {
    let nreg = regions.len();
    let n_tasks = cfg.n_tasks().min(nreg);
    let nthr = cfg.threads_per_reader();
    debug!(
        "No. regions: {}, No. read tasks: {}, threads per task: {}",
        nreg, n_tasks, nthr
    );
    let cfg_ref = &cfg;
    let ifile: &Path = input.as_ref();

    metadata.insert("indexed", format!("{}", cfg.indexed()));
    metadata.insert("n_tasks", format!("{}", n_tasks));
    metadata.insert("threads_per_reader", format!("{}", nthr));
    let min_qual = cfg.min_qual().min(255) as u8;
    let min_mapq = cfg.min_mapq().min(255) as u8;
    metadata.insert("min_qual", format!("{}", min_qual));
    metadata.insert("min_mapq", format!("{}", min_mapq));

    let st = if cfg.indexed() {
        debug!("Processing input with index");
        // Create thread scope so that we can share references across threads
        thread::scope(|scope| {
            let (collector_tx, collector_rx) = bounded(n_tasks * 4);
            let collector_task = scope.spawn(|| collect::collector(collector_rx));

            // Spawn threads
            let (region_tx, region_rx) = bounded(n_tasks * 4);
            let mut tasks: Vec<_> = (0..n_tasks)
                .map(|ix| {
                    let rx = region_rx.clone();
                    let tx = collector_tx.clone();
                    scope.spawn(move || read::reader(cfg_ref, ix, ifile, tx, rx))
                })
                .collect();
            drop(collector_tx);

            // Send unmapped reads to readers
            region_tx
                .send((String::from('*'), 0, None))
                .expect("Error sending region message");

            // Send mapped regions to readers
            for (ctg, regv) in regions.iter() {
                let ctg_str = if ctg.contains(':') {
                    format!("{{{}}}", ctg)
                } else {
                    format!("{}", ctg)
                };
                for (ix, reg) in regv.iter().enumerate() {
                    let s = format!("{}{}", ctg_str, reg);
                    region_tx
                        .send((s, ix, Some(regv)))
                        .expect("Error sending region message");
                }
            }
            drop(region_tx);

            // Wait for readers to finish
            for jh in tasks.drain(..) {
                let _ = jh.join();
            }
            // Wait for collector to finish
            collector_task.join()
        })
    } else {
        debug!("Processing input without index using {} tasks", n_tasks);

        // First we assign regions to each task

        // Get total length of all regions
        let total_size: usize = regions
            .values()
            .map(|regv| {
                regv.iter()
                    .map(|r| r.length().expect("Open interval found"))
                    .sum::<usize>()
            })
            .sum();
        debug!("Total length of all regions: {}", total_size);

        // Allocate regions to tasks
        // Tasks are given contiguous sets of regions to decrease
        // the chance that a read spans regions handled by different tasks
        let mut region_hash = HashMap::new();
        let mut remainder = total_size;
        let mut task_ix = 0;
        let mut current_total = 0;
        let mut task_regions = Vec::with_capacity(n_tasks);
        let mut task_region_list = Vec::new();
        for (ctg, regv) in regions.iter() {
            let mut v = Vec::with_capacity(regv.len());
            let mut tv = Vec::new();
            for (reg_ix, reg) in regv.iter().enumerate() {
                let l = reg.length().unwrap();
                let avail = n_tasks - task_ix;
                if avail == 1 || current_total + l < remainder / avail {
                    current_total += l;
                    tv.push((reg, reg_ix))
                } else {
                    debug!("Task {}, total region length {}", task_ix, current_total);
                    if !tv.is_empty() {
                        let tr = TaskRegion::new(ctg.as_ref(), tv);
                        task_region_list.push(tr);
                        tv = Vec::new();
                    }
                    task_regions.push(task_region_list);
                    task_region_list = Vec::new();
                    task_ix += 1;
                    remainder -= current_total;
                    current_total = l;
                    tv.push((reg, reg_ix));
                }
                v.push(task_ix)
            }
            region_hash.insert(Rc::clone(ctg), v);
            if !tv.is_empty() {
                let tr = TaskRegion::new(ctg.as_ref(), tv);
                task_region_list.push(tr)
            }
        }
        task_regions.push(task_region_list);
        debug!("Task {}, total region length {}", task_ix, current_total);

        thread::scope(|scope| {
            // Spawn handling tasks

            let (collector_tx, collector_rx) = bounded(n_tasks * 4);
            let collector_task = scope.spawn(|| collect::collector(collector_rx));

            let (bam_tx, bam_rx) = unbounded();
            let mut task_tx = Vec::with_capacity(n_tasks);
            let mut tasks: Vec<_> = task_regions
                .iter()
                .enumerate()
                .map(|(ix, th)| {
                    let (tx, rx) = unbounded();
                    task_tx.push(tx);
                    let collect_tx = collector_tx.clone();
                    let b_tx = bam_tx.clone();
                    scope.spawn(move || read::read_handler(cfg_ref, ix, th, collect_tx, rx, b_tx))
                })
                .collect();
            drop(bam_tx);
            // Read from input
            read::read_input_mt(
                cfg_ref,
                ifile,
                &regions,
                &task_regions,
                task_tx,
                bam_rx,
                collector_tx,
            )
            .expect("Error opening input file");

            // Wait for read handlers to finish
            for jh in tasks.drain(..) {
                let _ = jh.join();
            }
            // Wait for collector to finish
            collector_task.join()
        })
    }
    .expect("Processing error");

    // Create output dir if necessary
    if let Some(dir) = cfg.dir() {
        fs::create_dir_all(dir).with_context(|| "Could not create output directory")?;
        metadata.insert("output_dir", dir.to_string_lossy().to_string());
    }

    // Create results file name
    let fname = {
        let mut d = cfg.dir().map(|p| p.to_owned()).unwrap_or_else(PathBuf::new);
        d.push(cfg.prefix().expect("Missing output prefix"));
        d.set_extension("json");
        d
    };

    metadata.insert("json_output_name", fname.to_string_lossy().to_string());

    // Open output file
    let mut wrt = fs::File::create(&fname)
        .with_context(|| format!("Could not open output file {}", fname.display()))?;

    // Get elapsed time
    if let Some(mut t) = cfg.elapsed().map(|d| d.as_secs_f64()) {
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
    }

    let out = Output {
        metadata,
        stats: st,
    };
    // Generate JSON
    serde_json::to_writer_pretty(&mut wrt, &out)
        .with_context(|| "Could not write out JSON stats file")
}
