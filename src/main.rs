#[macro_use]
extern crate log;
#[macro_use]
extern crate anyhow;

mod cli;
mod collect;
mod config;
mod input;
mod mappability;
mod process_indexed;
mod process_unindexed;
mod read;
mod regions;
mod stats;

fn main() -> anyhow::Result<()> {
    let (cfg, input, regions) = cli::handle_cli()?;
    if cfg.indexed() {
        process_indexed::process(cfg, input, regions)
    } else {
        process_unindexed::process(cfg)
    }
}
