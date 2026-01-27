use m_htslib::sam::{Cigar, CigarOp};

pub(super) fn get_start_end(cigar: &Cigar, rev: bool, read_len: usize) -> Option<(usize, usize)> {
    let mut first = 0;
    let mut last = 0;
    let mut started = false;
    let mut x = 0;
    for elem in cigar.iter() {
        match elem.op() {
            CigarOp::HardClip | CigarOp::SoftClip | CigarOp::Ins => x += elem.op_len() as usize,
            CigarOp::Match | CigarOp::Equal | CigarOp::Diff => {
                if !started {
                    started = true;
                    first = x;
                }
                x += elem.op_len() as usize;
                last = x;
            }
            _ => {}
        }
    }
    adjust_start_end(first, last, x, rev, read_len)
}

pub(super) fn adjust_start_end(
    first: usize,
    last: usize,
    x: usize,
    rev: bool,
    read_len: usize,
) -> Option<(usize, usize)> {
    if read_len != x || last <= first {
        None
    } else if rev {
        Some((read_len - last, read_len - 1 - first))
    } else {
        Some((first, last - 1))
    }
}
