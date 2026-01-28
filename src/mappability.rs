use std::cmp::Ordering;

use m_htslib::region::{Reg, RegCoords, RegionCoords, RegionList};

pub struct Mappability {
    inner: RegionList,
}

impl Mappability {
    pub fn new(mut inner: RegionList) -> Self {
        inner.normalize();
        Self { inner }
    }

    pub fn region_intersect<'a>(&'a self, reg: &Reg) -> &'a [RegionCoords] {
        reg.reg_contig()
            .and_then(|ctg| self.inner.contig_reg_list(ctg))
            .and_then(|crl| crl.regions())
            .and_then(|rc| ctg_reg_intersect(rc, reg)).unwrap_or(&[])
    }
}

fn ctg_reg_intersect<'a>(rc: &'a [RegionCoords], reg: &Reg) -> Option<&'a [RegionCoords]> {
    let (start, end) = reg.coords();

    let a = start
        .map(|s| {
            let s1 = s + 1;
            match rc.binary_search_by(|r| r.end().map(|x| x.cmp(&s1)).unwrap_or(Ordering::Greater)) {
                Ok(i) => i,
                Err(i) => i,
            }
        })
        .unwrap_or_default();

    let l = rc.len();
    if a == l {
        None
    } else {
        let b = end
            .map(|e| match rc.binary_search_by_key(&e, |rc| rc.start() + 1) {
                Ok(i) => i + 1,
                Err(i) => i,
            })
            .unwrap_or(l);
        Some(&rc[a..b])
    }
}
