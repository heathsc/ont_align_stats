use crate::{config::Config, stats::Stats};
use crossbeam_channel::Receiver;

pub fn collector(cfg: &Config, rx: Receiver<Stats>) {
    debug!("Starting collector thread");
    let mut st = Stats::new();
    while let Ok(x) = rx.recv() {
        debug!("Collector received Stats");
        st += x;
    }
    println!("{:?}", st);
    debug!("Terminating collector thread");
}
