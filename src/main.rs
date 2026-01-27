#[macro_use]
extern crate log;
#[macro_use]
extern crate anyhow;

mod cli;
mod collect;
mod input;
mod metadata;
mod process;
mod read;
mod regions;
mod stats;

pub use cli::config::{Config, CompressOpt};

fn main() -> anyhow::Result<()> {
    let cfg = cli::handle_cli()?;
    process::process(cfg)
}
