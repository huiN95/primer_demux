use crate::{bam_record_extention::ReadRecord, core_context::PrimerCandidate};
use metrics::counter;
use std::collections::BTreeSet;
// use std::sync::Arc;

pub fn merge_non_overlapping_no_copy(
    candidates: &mut Vec<PrimerCandidate>,
) -> Vec<PrimerCandidate> {
    if candidates.is_empty() {
        return Vec::new();
    }

    // 1) Sort by start first (in-place sorting, no cloning)
    candidates.sort_by_key(|c| c.start);

    // 2) Move elements from candidates, merging into out while iterating
    let mut out: Vec<PrimerCandidate> = Vec::with_capacity(candidates.len());

    for cand in candidates.drain(..) {
        if out.is_empty() {
            out.push(cand);
            continue;
        }

        let last = out.last_mut().unwrap();

        if cand.start >= last.end {
            // No overlap, add directly
            out.push(cand);
        } else {
            // Overlap: choose which one to keep based on rules
            counter!("primer_overlap_cnts").increment(1);

            let better = (cand.distance < last.distance)
                || (cand.distance == last.distance
                    && (cand.end - cand.start) > (last.end - last.start));

            if better {
                *last = cand; // This is a move: overwrite last with cand (no cloning)
            }
            // Otherwise discard cand (automatically dropped)
        }
    }

    // At this point, candidates has been cleared by drain
    out
}

pub fn get_pattern_keys(pattern_records: &Vec<ReadRecord>) -> Vec<String> {
    let mut set: BTreeSet<String> = pattern_records
        .iter()
        .map(|r| {
            // Assume field name is r.id: Arc<str> (modify according to actual field name)
            let s: &str = r.id.as_ref();

            // Remove last 2 characters (if less than 2, becomes empty string)
            let trimmed = if s.len() >= 2 { &s[..s.len() - 2] } else { "" };

            trimmed.to_string()
        })
        .collect();
    set.insert("uncertain".into());
    set.into_iter().collect()
}
// pub fn merge_non_overlapping_no_copy(
//     candidates: &mut Vec<PrimerCandidate>,
// ) -> Vec<PrimerCandidate> {
//     if candidates.is_empty() {
//         return Vec::new();
//     }

//     // 1. Sort by start first
//     candidates.sort_by_key(|c| c.start);

//     // 2. Use index i to manage the end of the "merged interval"
//     let mut i = 0;
//     for j in 1..candidates.len() {
//         if candidates[j].start > candidates[i].end {
//             // No overlap, increment index i
//             i += 1;
//             candidates[i] = candidates[j].clone();
//         } else {
//             counter!("primer overlap cnts").increment(1);
//             // Overlap, replace if distance is smaller
//             if candidates[j].distance < candidates[i].distance {
//                 candidates[i] = candidates[j].clone();
//             } else if candidates[j].distance == candidates[i].distance {
//                 if candidates[j].end - candidates[j].start > candidates[i].end - candidates[i].start
//                 {
//                     // If distance is the same and new interval is longer, replace
//                     candidates[i] = candidates[j];
//                 }
//             }
//         }
//     }

//     // 3. truncate to keep the first i+1 elements, which are the merge result
//     candidates.truncate(i + 1);

//     // 4. Use std::mem::take to move these i+1 elements, avoiding cloning
//     std::mem::take(candidates)
// }

#[cfg(test)]
mod test {
    use bio::alignment::Alignment;
    use bio::pattern_matching::myers::MyersBuilder;

    #[test]
    fn test_ambig_pattern() {
        let text: &'static [u8; 15] = b"GGATGAGCGCCATAG";
        let pattern = b"GAGGC";

        let mut myers = MyersBuilder::new().ambig(b'N', b"ACGT").build_64(pattern);
        let mut reuslt = myers.find_all(text, 2);
        let mut aln = Alignment::default();
        while reuslt.next_alignment(&mut aln) {
            println!(
                "start: {}, end: {}, distance: {}",
                aln.ystart, aln.yend, aln.score
            );
        }
    }
}
