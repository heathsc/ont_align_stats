use r_htslib::{BamAux, BamRec};

#[derive(Debug, Copy, Clone, PartialEq)]
pub enum BSStrand {
    StrandC2T,
    StrandG2A,
}

#[derive(PartialEq)]
enum Aligner {
    Unknown,
    GEM,
    Bowtie,
    Novoalign,
    BSMap,
    BWAMeth,
}

pub(super) fn get_bs_strand(b: &BamRec) -> Option<BSStrand> {
    let mut strand = None;
    if let Some(itr) = b.get_aux_iter() {
        for b in itr {
            let tag = b.data();
            let aligner = {
                if tag[0] == b'Z' {
                    if tag[1] == b'B' {
                        Aligner::Novoalign
                    } else if tag[1] == b'S' {
                        Aligner::BSMap
                    } else {
                        Aligner::Unknown
                    }
                } else if tag[0] == b'X' {
                    if tag[1] == b'G' {
                        Aligner::Bowtie
                    } else if tag[1] == b'B' {
                        Aligner::GEM
                    } else {
                        Aligner::Unknown
                    }
                } else if tag[0] == b'Y' && tag[1] == b'D' {
                    Aligner::BWAMeth
                } else {
                    Aligner::Unknown
                }
            };
            if aligner != Aligner::Unknown {
                match tag[2] {
                    b'A' if aligner == Aligner::GEM => {
                        if tag[3] == b'C' {
                            strand = Some(BSStrand::StrandC2T)
                        } else if tag[3] == b'G' {
                            strand = Some(BSStrand::StrandG2A)
                        }
                    }
                    b'Z' => match aligner {
                        Aligner::Bowtie | Aligner::Novoalign => {
                            if tag[3] == b'C' {
                                strand = Some(BSStrand::StrandC2T)
                            } else if tag[3] == b'G' {
                                strand = Some(BSStrand::StrandG2A)
                            }
                        }
                        Aligner::BSMap => {
                            if tag[3] == b'+' {
                                strand = Some(BSStrand::StrandC2T)
                            } else if tag[3] == b'-' {
                                strand = Some(BSStrand::StrandG2A)
                            }
                        }
                        Aligner::BWAMeth => {
                            if tag[3] == b'f' {
                                strand = Some(BSStrand::StrandC2T)
                            } else if tag[3] == b'r' {
                                strand = Some(BSStrand::StrandG2A)
                            }
                        }
                        _ => (),
                    },
                    _ => (),
                }
            }
            if strand.is_some() {
                break;
            }
        }
    }
    strand
}
