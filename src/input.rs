use std::{ffi::CString, path::Path};

use anyhow::Context;
use m_htslib::{
    hts::{HtsFile, HtsFormat, HtsThreadPool},
};

pub fn open_input<'a, P: AsRef<Path>>(
    name: P,
    reference: Option<&Path>,
    nthreads: usize,
    tpool: Option<&'a HtsThreadPool>,
) -> anyhow::Result<HtsFile<'a>> {
    let name = name.as_ref();
    debug!(
        "Try to open input file {} with reference {:?}",
        name.display(),
        reference,
    );
    // Set up htsFormat if we are specifying the reference
    let mut fmt = HtsFormat::default();
    let mut opt = None;
    if let Some(s) = reference {
        opt = Some(format!("reference={}", s.display()));
    }
    if nthreads > 0 && tpool.is_none() {
        debug!("Setting nthreads option to {} for input file", nthreads);
        let s = format!("nthreads={}", nthreads);
        if let Some(s1) = opt.as_mut() {
            s1.push(',');
            s1.push_str(s.as_str());
        } else {
            opt = Some(s);
        }
    }
    if let Some(s) = opt.as_ref() {
        let cs = CString::new(s.as_str()).expect("Error creating option");
        fmt.parse_opt_list(&cs)?;
    }
    // Try opening input file
    let mut hts = HtsFile::open_format(name, "r", &fmt)
        .with_context(|| format!("Failed to open input file {}", name.display()))?;

    // Try and attach thread pool if presemt
    if let Some(tp) = tpool {
        debug!("Try to attach thread pool to file");
        hts.set_thread_pool(tp)?
    }
    Ok(hts)
}
