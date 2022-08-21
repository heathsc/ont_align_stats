#[macro_use]
extern crate log;
#[macro_use]
extern crate anyhow;

mod cli;
mod collect;
mod config;
mod input;
mod mappability;
mod process;
mod read;
mod regions;
mod stats;

fn main() -> anyhow::Result<()> {
    let (cfg, input, regions) = cli::handle_cli()?;
    process::process(cfg, input, regions)
}
