use crate::config::Config;
use crossbeam_channel::Receiver;

pub fn collector(cfg: &Config, rx: Receiver<usize>) {
    debug!("Starting collector thread");
    while let Ok(x) = rx.recv() {
        debug!("Collector received {}", x)
    }
    debug!("Terminating collector thread");
}
