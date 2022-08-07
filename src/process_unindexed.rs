use crate::config::Config;

pub fn process(cfg: Config) -> anyhow::Result<()> {
    error!("Support for non-indexes files is not currently supported");

    Ok(())
}
