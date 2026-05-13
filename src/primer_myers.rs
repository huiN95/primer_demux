use crate::bam_record_extention::ReadRecord;
use crate::core_context::PrimerCandidate;
use bio::alignment::Alignment;
use bio::pattern_matching::myers::long::Myers as MyersLong;
use bio::pattern_matching::myers::Myers as Myers64;
use bio::pattern_matching::myers::MyersBuilder;
use std::collections::HashMap;
use std::sync::Arc;

pub enum MayersPattern {
    Myers64 {
        myers: Myers64<u64>,
        pattern: String,
    },
    MyersLong {
        myers: MyersLong<u8>,
        pattern: String,
    },
}

impl MayersPattern {
    #[inline]
    pub fn pattern(&self) -> &str {
        match self {
            MayersPattern::Myers64 { pattern, .. } => pattern,
            MayersPattern::MyersLong { pattern, .. } => pattern,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]

pub enum Direction {
    Forward,
    Reverse,
}

const AMBIGS: &[(u8, &[u8])] = &[
    (b'M', b"AC"),
    (b'R', b"AG"),
    (b'W', b"AT"),
    (b'S', b"CG"),
    (b'Y', b"CT"),
    (b'K', b"GT"),
    (b'V', b"ACG"),
    (b'H', b"ACT"),
    (b'D', b"AGT"),
    (b'B', b"CGT"),
    (b'N', b"ACGT"),
];

/// Convert a batch of `ReadRecord` into `HashMap<id, MyersType>`
// pub fn build_map(records: &[ReadRecord], builder: &MyersBuilder) -> HashMap<Arc<str>, MayersType> {
//     records
//         .iter()
//         .map(|rec| {
//             let myers = if rec.sequence.len() <= 64 {
//                 MayersType::Myers64(builder.build_64(&rec.sequence))
//             } else {
//                 MayersType::MyersLong(builder.build_long(&rec.sequence))
//             };
//             (rec.id.clone(), myers)
//         })
//         .collect()
// }

// pub fn get_myers_from_read_record_with_additional_equalities(
//     patterns: &[ReadRecord],
// ) -> HashMap<Arc<str>, MayersType> {
//     let mut builder = MyersBuilder::new();
//     for &(base, equivalents) in AMBIGS {
//         builder.ambig(base, equivalents);
//     }
//     let patterns_myers = build_map(patterns, &builder);
//     patterns_myers
// }

pub fn get_alignments_from_myers(
    name: Arc<str>,
    myers: &mut MayersPattern,
    target: &[u8],
    max_distance: u8,
    candidates: &mut Vec<PrimerCandidate>,
) {
    let mut aln = Alignment::default();
    match myers {
        MayersPattern::Myers64 {
            myers,
            pattern: _pat,
        } => {
            let mut matches = myers.find_all(target, max_distance.into());
            // println!("current myers name={}", &*name);
            // println!("current pattern {}", _pat);

            while matches.next_alignment(&mut aln) {
                // println!("current myers name={}", &*name);
                // println!("current pattern {}", _pat);
                candidates.push(PrimerCandidate {
                    start: aln.ystart,
                    end: aln.yend,
                    distance: aln.score,
                    name: name.clone(),
                });
            }
        }
        MayersPattern::MyersLong { myers, pattern: _ } => {
            let mut matches = myers.find_all(target, max_distance.into());
            while matches.next_alignment(&mut aln) {
                candidates.push(PrimerCandidate {
                    start: aln.ystart,
                    end: aln.yend,
                    distance: aln.score,
                    name: name.clone(),
                });
            }
        }
    }
}

// pub fn get_seq_reverse_complement(seq: &[u8]) -> Vec<u8> {
//     seq.iter()
//         .rev()
//         .map(|&base| match base {
//             b'A' => b'T',
//             b'T' => b'A',
//             b'C' => b'G',
//             b'G' => b'C',
//             _ => base, // Keep other characters unchanged
//         })
//         .collect()
// }

/// Wrap builder results into MayersPattern
#[inline]
fn build_myers_enum(builder: &MyersBuilder, pat: &[u8]) -> MayersPattern {
    // Rule: use Myers64 for length up to 64, MyersLong for longer
    let pattern = String::from_utf8_lossy(pat).into_owned();

    if pat.len() <= 64 {
        let m64 = builder.build_64(pat); // ← modify according to actual API name
        MayersPattern::Myers64 {
            myers: m64,
            pattern,
        }
    } else {
        let ml = builder.build_long(pat); // ← modify according to actual API name
        MayersPattern::MyersLong { myers: ml, pattern }
    }
}

fn compl_iupac(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' | b'U' => b'A',
        b'R' => b'Y',
        b'Y' => b'R',
        b'S' => b'S',
        b'W' => b'W',
        b'K' => b'M',
        b'M' => b'K',
        b'B' => b'V',
        b'V' => b'B',
        b'D' => b'H',
        b'H' => b'D',
        b'N' => b'N',
        _ => b,
    }
}

fn revcomp_iupac(seq: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(seq.len());
    for &b in seq.iter().rev() {
        out.push(compl_iupac(b.to_ascii_uppercase()));
    }
    out
}

pub fn get_myers_from_primers(
    patterns: &[ReadRecord],
    direction: Direction,
    cut_off_bases: u8,
) -> HashMap<Arc<str>, MayersPattern> {
    let mut builder = MyersBuilder::new();
    for &(base, equivalents) in AMBIGS {
        builder.ambig(base, equivalents);
    }

    patterns
        .iter()
        .filter_map(|r| {
            let id = r.id.as_ref();
            let mut pattern = r.sequence.as_slice();
            if cut_off_bases > 0 {
                if id.ends_with("_L") {
                    pattern = &r.sequence[cut_off_bases as usize..r.sequence.len()];
                } else {
                    pattern = &r.sequence[0..r.sequence.len() - cut_off_bases as usize];
                }
            }
            // println!("{}", String::from_utf8_lossy(pattern));

            // Only used to determine _L/_T, no longer used to generate key
            let is_l = if id.ends_with("_L") {
                true
            } else if id.ends_with("_T") {
                false
            } else {
                // Should panic here, this pattern does not match the preset
                return None;
            };

            // Decide whether to use revcomp based on direction + suffix
            let need_revcomp = match (direction, is_l) {
                (Direction::Forward, true) => false,
                (Direction::Forward, false) => false, // _L: Forward original sequence
                (Direction::Reverse, true) => true,   // _L: Reverse revcomp
                (Direction::Reverse, false) => true,  // _T: Forward revcomp
            };

            let myers = if need_revcomp {
                let rc = revcomp_iupac(pattern);
                build_myers_enum(&builder, &rc)
            } else {
                build_myers_enum(&builder, pattern)
            };

            // ✅ use original id (including _L/_T) as key to avoid name collisions
            Some((Arc::<str>::from(id), myers))
        })
        .collect()
}
