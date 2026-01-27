use m_htslib::sam::{BamAuxVal, BamRec};

use super::utils::adjust_start_end;

pub(super) fn sa_get_start_end(
    it: &mut SaTagIter,
    rev: bool,
    read_len: usize,
) -> Option<(usize, usize)> {
    let mut first = 0;
    let mut last = 0;
    let mut started = false;
    let mut x = 0;
    for (c, l) in it {
        match c {
            'H' | 'S' | 'I' => x += l,
            'M' | '=' | 'X' => {
                if !started {
                    started = true;
                    first = x;
                }
                x += l;
                last = x;
            }
            _ => {}
        }
    }
    adjust_start_end(first, last, x, rev, read_len)
}

pub(super) struct SaTagIter<'a> {
    tag: &'a str,
    done: bool,
}

impl<'a> SaTagIter<'a> {
    pub(super) fn new(tag: &'a str) -> Self {
        Self { tag, done: false }
    }
}

impl<'a> Iterator for SaTagIter<'a> {
    type Item = (char, usize);
    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }
        let it = self.tag.char_indices();
        let mut ix = None;
        let mut ch = None;
        for (iy, c) in it {
            if !c.is_numeric() {
                ch = Some((iy, c));
                break;
            }
            ix = Some(iy)
        }
        if let Some(ix) = ix {
            if let Ok(l) = self.tag[0..=ix].parse::<usize>() {
                if let Some((iy, c)) = ch {
                    if iy + 1 >= self.tag.len() {
                        self.done = true
                    }
                    self.tag = &self.tag[iy + 1..];
                    Some((c, l))
                } else {
                    None
                }
            } else {
                None
            }
        } else {
            None
        }
    }
}

pub(super) fn process_sa_tag(
    rec: &mut BamRec,
    read_len: usize,
    v: &mut Vec<(usize, usize)>,
) -> anyhow::Result<()> {
    if let Some(b) = rec.get_tag("SA")?
        && let Ok(bval) = b.get_val()
        && let BamAuxVal::String(cs) = bval
        && let Ok(sa) = cs.to_str()
    {
        trace!("Found SA tag {sa:?} for read {:?}", rec.qname().unwrap());
        // Collect the start and end points of all supplementary mappings for this read
        for s in sa.split(';') {
            if s.len() < 11 {
                continue;
            }
            let fd: Vec<_> = s.split(',').collect();
            if fd.len() == 6 {
                if let Some(rev) = match fd[2] {
                    "+" => Some(false),
                    "-" => Some(true),
                    _ => {
                        warn!("Illegal SA tag {} - strand not + or -", s);
                        None
                    }
                } {
                    let mut it = SaTagIter::new(fd[3]);
                    if let Some((a, b)) = sa_get_start_end(&mut it, rev, read_len) {
                        v.push((a, b))
                    } else {
                        warn!(
                            "Illegal SA Tag {} for read {:?} (wrong size)",
                            fd[3],
                            rec.qname().unwrap()
                        )
                    }
                }
            } else {
                warn!(
                    "Illegal SA tag {} (wrong number of fields {} instead of 6) >{}<",
                    s,
                    fd.len(),
                    sa
                )
            }
        }

        // Sort mappings  by starting point (on read)
        v.sort_unstable_by_key(|(x, _)| *x);
    }
    Ok(())
}
