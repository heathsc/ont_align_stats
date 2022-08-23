use crate::stats::Stats;
use crossbeam_channel::Receiver;

pub fn collector(rx: Receiver<Stats>) -> Stats {
    debug!("Starting collector thread");
    let mut st = Stats::new();
    while let Ok(x) = rx.recv() {
        debug!("Collector received Stats");
        st += x;
    }
    debug!("Terminating collector thread");
    st
}
