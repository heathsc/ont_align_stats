use std::{
    fs,
    path::{Path, PathBuf},
    thread,
};

use anyhow::{Context, Error};
use crossbeam_channel::bounded;

use crate::{collect, config::Config, input, read, regions::Regions};

pub fn process(cfg: Config, input: PathBuf, regions: Regions) -> anyhow::Result<()> {
    debug!("Processing input with index");

    // Create output dir if necessary
    if let Some(dir) = cfg.dir() {
        fs::create_dir_all(dir).with_context(|| "Could not create output directory")?;
    }

    let nreg = regions.len();
    let n_tasks = cfg.n_tasks().min(nreg);
    let nthr = cfg.threads_per_task();
    debug!(
        "No. regions: {}, No. read tasks: {}, threads per task: {}",
        nreg, n_tasks, nthr
    );
    let cfg_ref = &cfg;
    let ifile: &Path = input.as_ref();
    // Create thread scope so that we can share references across threads
    thread::scope(|scope| {
        // Spawn threads
        let (collector_tx, collector_rx) = bounded(n_tasks * 4);
        let collector_task = scope.spawn(|| collect::collector(cfg_ref, collector_rx));
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

        // Send regions to readers
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
        let _ = collector_task.join();
    });
    Ok(())
}
