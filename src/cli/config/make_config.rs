use std::{
    path::{Path, PathBuf},
    time::Instant,
};

use anyhow::Context;
use clap::ArgMatches;
use m_htslib::{region::RegionList, sam::SamHdr};
use num_cpus::get_physical;

use super::{CompressOpt, Config};

use crate::{
    input::open_input,
    regions::{
        region_list_from_bed, region_list_from_sam_hdr, region_list_from_str_slice, split_regions,
    },
};

impl Config {
    pub fn from_matches(m: &ArgMatches) -> anyhow::Result<Self> {
        let start_time = Instant::now();
        let input = m
            .get_one::<PathBuf>("input")
            .expect("Missing input")
            .to_owned();

        let regions = make_region_list(m, &input)?;

        let reference = m.get_one::<PathBuf>("reference").map(|m| m.to_owned());

        let n_tasks = *m.get_one::<u64>("n_tasks").unwrap() as usize;
        let hts_threads = m.get_one::<u64>("hts_threads").map(|x| *x as usize).unwrap_or_else(get_physical);

        let min_mapq = *m.get_one::<u8>("min_mapq").unwrap();
        let min_qual = *m.get_one::<u8>("min_qual").unwrap();

        let min_read_len = m.get_one::<u64>("min_read_len").map(|x| *x as usize);
        let max_read_len = m.get_one::<u64>("max_read_len").map(|x| *x as usize);

        let prefix = m.get_one::<String>("prefix").unwrap().to_owned();
        let dir = m.get_one::<PathBuf>("dir").map(|p| p.to_owned());

        let compress = if m.get_flag("compress_gzip") {
            CompressOpt::Gzip
        } else if m.get_flag("compress_bzip2") {
            CompressOpt::Bzip2
        } else if m.get_flag("compress_zstd") {
            CompressOpt::Zstd
        } else if m.get_flag("compress_xz") {
            CompressOpt::Xz
        } else {
            CompressOpt::None
        };

        Ok(Self {
            input,
            region_list: regions,
            reference,
            n_tasks,
            hts_threads,
            min_mapq,
            min_qual,
            min_read_len,
            max_read_len,
            prefix,
            dir,
            compress,
            start_time,
        })
    }
}

fn make_region_list(m: &ArgMatches, input: &Path) -> anyhow::Result<RegionList> {
    let mut hts = open_input(input, None, 1, None)
        .with_context(|| format!("Error opening input file {}", input.display()))?;

    let hdr = SamHdr::read(&mut hts).with_context(|| {
        format!(
            "Could not read SAM/BAM/CRAM header for input file {}",
            input.display()
        )
    })?;

    // If region and/or regions options are present, process then both and merge.
    // Otherwise take all the contigs from the input file + unmapped reads
    let (mut rlist, chk_flag) = {
        let mut rl = None;
    
        if let Some(rl_file) = m.get_one::<PathBuf>("regions").map(|p| p.as_path()) {
            rl = Some(region_list_from_bed(rl_file, rl)?);
        }
        
        if let Some(rs) = m
            .get_many::<String>("region")
            .map(|v| v.map(|s| s.as_str()).collect::<Vec<&str>>())
        {
            rl = Some(region_list_from_str_slice(rs.as_ref(), rl)?);
        }
        
        // If All regions has been specified, set to None so that we can pull the ctgs from the sam header
        if rl.as_ref().map(|r| r.is_all_regions()).unwrap_or_default() {
            rl = None
        }
        
        if let Some(rlist) = rl.take() {
            (rlist, true)
        } else {
            (region_list_from_sam_hdr(&hdr)?, false)
        }
    };

    debug!("Initial region count: {}", rlist.regions().count());
    
    if chk_flag {
        for c in rlist.contigs() {
            if hdr.name2tid(c).is_err() {
                return Err(anyhow!("Contig {c:?} not found in SAM header"));
            }
        }
    }

    if let Some(map) = m.get_one::<PathBuf>("mappability").map(|p| p.as_path()) {
        debug!("Adding mappabililty information from {}", map.display());
        let map_list = region_list_from_bed(map, None)?;
        debug!("Read in {} regions from {}", map_list.regions().count(), map.display());
        rlist
            .intersect(&map_list)
            .with_context(|| "Error computing intersect between regions and mappability file")?;
        debug!("Regions after adding mappability: {}", rlist.regions().count());
    }

    let max_block_size = *m.get_one::<u64>("max_block_size").unwrap();
    split_regions(&mut rlist, max_block_size, &hdr);
    debug!("Final region count: {}", rlist.regions().count());
    Ok(rlist)
}
