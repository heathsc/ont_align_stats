use std::{io::BufRead, path::Path};

use anyhow::Context;
use compress_io::compress::CompressIo;

use m_htslib::{
    hts::HtsPos,
    region::{Reg, RegionCoords, RegionList},
    sam::SamHdr,
};

pub fn region_list_from_bed(
    p: &Path,
    mut prev_list: Option<RegionList>,
) -> anyhow::Result<RegionList> {
    let mut rdr = CompressIo::new()
        .path(p)
        .bufreader()
        .with_context(|| format!("Error opening bed file {}", p.display()))?;

    let mut s = String::new();
    let mut rlist = prev_list.take().unwrap_or_default();
    let mut line = 0;
    while rdr
        .read_line(&mut s)
        .with_context(|| format!("Error reading from bed file {}", p.display()))?
        > 0
    {
        line += 1;
        let s1 = s.trim();
        // Skip blanks lines and comments
        if s1.is_empty() || s1.starts_with('#') {
            continue;
        }

        // Try to parse line as bed
        let reg = match Reg::parse_bed_from_str(s1) {
            Err(e) => {
                // Check for special browser or track lines. We check now in case thre is a reference using one of these names
                if s1.starts_with("browser") || s1.starts_with("track") {
                    None
                } else {
                    return Err(anyhow!(
                        "Error parsing bed file {} at line {line}: {e}",
                        p.display()
                    ));
                }
            }
            Ok(r) => Some(r),
        };
        if let Some(r) = reg {
            rlist.add_reg(&r);
        }
        s.clear();
    }
    rlist.normalize();
    Ok(rlist)
}

pub fn region_list_from_str_slice(
    v: &[&str],
    mut prev_list: Option<RegionList>,
) -> anyhow::Result<RegionList> {
    let mut rlist = prev_list.take().unwrap_or_default();
    for s in v {
        let reg = Reg::try_from(s.trim()).with_context(|| format!("Error parsing region {s}"))?;
        rlist.add_reg(&reg);
    }
    rlist.normalize();
    Ok(rlist)
}

pub fn region_list_from_sam_hdr(hdr: &SamHdr) -> anyhow::Result<RegionList> {
    let mut rlist = RegionList::new();
    let nref = hdr.nref();
    for i in 0..nref {
        if let (Some(ctg), Some(l)) = (hdr.tid2name(i), hdr.tid2len(i))
            && l > 0
        {
            let reg =
                Reg::new_region(ctg, Some(&RegionCoords::new(0, Some(l as HtsPos)).unwrap()))?;
            rlist.add_reg(&reg);
        }
    }
    rlist.add_reg(&Reg::Unmapped);
    // DOn't need to normalize as we are only adding one region per contig
    Ok(rlist)
}

pub fn split_regions(rlist: &mut RegionList, max_block_size: u64, hdr: &SamHdr) {
    let max_block_size = max_block_size as i64;
    debug!("Splitting regions into block size of maximum {max_block_size}");

    for (ctg, crl) in rlist.contig_reg_lists_mut() {
        let tid = hdr.name2tid(ctg).expect("Unknown contig");
        let seq_len = hdr.tid2len(tid).expect("Bad ctg tid");

        let mut new_regions = None;
        if let Some(rs) = crl.regions() {
            let mut new_reg = Vec::new();
            for r in rs {
                let (start, end) = r.get_range(seq_len).expect("Invalid range");
                let l = end - start;
                if l <= max_block_size {
                    trace!("Keeping Region {:?}:{}-{}", ctg, start, end);
                    new_reg.push(*r)
                } else {
                    trace!("Splitting Region {:?}:{}-{}", ctg, start, end);
                    let mut ns = 1 + l / max_block_size;
                    let mut x = start;
                    while ns > 0 {
                        assert!(x < end);
                        let remainder = end - x;
                        let s = remainder / ns;
                        let reg1 = RegionCoords::new(x, Some(x + s)).unwrap();
                        trace!("Region {:?} (split)", reg1);
                        new_reg.push(reg1);
                        x += s;
                        ns -= 1;
                    }
                    assert_eq!(x, end);
                }
            }
            new_regions = Some(new_reg);
        }
        crl.set_regions(new_regions);
    }
}
