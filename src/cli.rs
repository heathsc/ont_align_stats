mod cli_model;
pub mod config;
mod init_log;
mod log_level;

pub use log_level::LogLevel;

use config::Config;

pub fn handle_cli() -> anyhow::Result<Config> {
    let c = cli_model::cli_model();
    let m = c.get_matches();
    init_log::init_log(&m);
    Config::from_matches(&m)
}
