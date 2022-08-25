use std::path::Path;

use anyhow::Context;
use r_htslib::*;

pub fn open_input<P: AsRef<Path>>(
    name: P,
    no_index: bool,
    reference: Option<&Path>,
    nthreads: usize,
) -> anyhow::Result<Hts> {
    let name = name.as_ref();
    debug!(
        "Try to open input file {} with threads {}, reference {:?} and no_index option {}",
        name.display(),
        nthreads,
        reference,
        no_index
    );
    // Set up htsFormat if we are specifying the reference
    let mut fmt = HtsFormat::default();
    if let Some(s) = reference {
        fmt.opt_add(format!("reference={}", s.display()))?
    }
    if nthreads > 1 {
        fmt.opt_add(format!("nthreads={}", nthreads))?
    }
    // Try opening input file
    let mut hts = Hts::open_format(name, "r", &fmt)
        .with_context(|| format!("Failed to open input file {}", name.display()))?;

    // Check that this is a SAM type file (SAM/BAM/CRAM)
    if !matches!(hts.rec_type(), Some(HtsRecType::Sam)) {
        Err(anyhow!(
            "Incorrect file format for input file {}",
            name.display()
        ))
    } else {
        // Try and load index
        if !no_index && hts.index_load().is_err() {
            warn!(
                "Warning - index not found for input file {}",
                name.display()
            );
        }
        Ok(hts)
    }
}
