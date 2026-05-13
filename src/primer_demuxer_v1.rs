use crate::bam_record_extention::ReadRecord;
use crate::core_context::PrimerCandidate;
use crate::find_pattern::merge_non_overlapping_no_copy;
use std::sync::Arc;

use crate::primer_myers::{get_alignments_from_myers, MayersPattern};

use metrics::{self, counter, histogram};

use std::collections::HashMap;

use crate::core_context::PrimerPair;
// use tracing::{debug, info, level_filters};

pub fn primer_demuxer_v1(
    min_read_length: usize,
    max_primer_distance: u8,

    primer_l_pattern: &mut HashMap<Arc<str>, MayersPattern>,
    primer_t_pattern: &mut HashMap<Arc<str>, MayersPattern>,
    target: &ReadRecord,
) -> Result<Vec<PrimerPair>, Box<dyn std::error::Error>> {
    let mut result = Vec::<PrimerPair>::new();
    // If sequence length is <= 2 * min_read_length, return empty result
    if target.sequence.len() <= 2 * min_read_length {
        counter!("filtered_reads_too_short").increment(1 as u64);
        return Ok(result);
    }

    // let mut all_primer_pos = IndexMap::<Arc<str>, Vec<PrimerCandidate>>::new();
    // let mut leading_primer_pos = Vec::<PrimerCandidate>::new();

    // loop all the patterns, to see which one is the best match
    let mut leading_keys: Vec<_> = primer_l_pattern.keys().cloned().collect();
    leading_keys.sort();

    // println!("primer demux: {:?} patterns found, ", keys,);
    // for (name, myers) in inside_patterns_myers.iter_mut() {
    let mut candidates = Vec::new();
    for name in leading_keys.iter() {
        // println!("primer name {name}");
        let myers = primer_l_pattern.get_mut(name).unwrap();

        get_alignments_from_myers(
            name.clone(),
            myers,
            &target.sequence,
            max_primer_distance,
            &mut candidates,
        );
    }
    let leading_primer_pos = merge_non_overlapping_no_copy(&mut candidates);

    // for primer in leading_primer_pos.iter() {
    //     println!(
    //         "leading primer {} has distance {} start with {}",
    //         primer.name, primer.distance, primer.start
    //     )
    // }
    let mut trailing_keys: Vec<_> = primer_t_pattern.keys().cloned().collect();
    trailing_keys.sort();
    // let mut trailing_primer_pos = Vec::<PrimerCandidate>::new();
    let mut candidates = Vec::new();

    for name in trailing_keys.iter() {
        let myers: &mut MayersPattern = primer_t_pattern.get_mut(name).unwrap();

        get_alignments_from_myers(
            name.clone(),
            myers,
            &target.sequence,
            max_primer_distance,
            &mut candidates,
        );
        // Theoretically, there should be no overlaps
    }
    let trailing_primer_pos = merge_non_overlapping_no_copy(&mut candidates);
    // for primer in trailing_primer_pos.iter() {
    //     println!(
    //         "trailing primer {} has distance {} start with {}",
    //         primer.name, primer.distance, primer.start
    //     )
    // }
    // Leading and trailing primer sequences are now sorted.
    // Both should not be empty
    if trailing_primer_pos.is_empty() || leading_primer_pos.is_empty() {
        return Ok(result);
    } else {
        // for
        // let start_primer_max_idx = leading_primer_pos.len() - 1;
        // let end_primer_max_idx = trailing_primer_pos.len() - 1;
        let mut leading_idx = 0;
        let mut trailing_idx = 0;
        let mut pre_pair_end_postion = 0;
        while leading_idx < leading_primer_pos.len() && trailing_idx < trailing_primer_pos.len() {
            // let start_primer = &leading_primer_pos[leading_idx];
            let mut left_move_step = 0;
            let mut right_move_step = 0;
            if leading_primer_pos[leading_idx].start < pre_pair_end_postion {
                left_move_step = 1;
            } else {
                (left_move_step, right_move_step) = primer_pair_status(
                    &leading_primer_pos[leading_idx],
                    &trailing_primer_pos[trailing_idx],
                );
                if left_move_step == 1 && right_move_step == 1 {
                    let pair_id = &leading_primer_pos[leading_idx].name
                        [0..&leading_primer_pos[leading_idx].name.len() - 2];
                    let left_distance = leading_primer_pos[leading_idx].distance;
                    let right_distance = trailing_primer_pos[trailing_idx].distance;

                    let distance = (left_distance, right_distance);
                    // Construct as left-closed and right-open.
                    let outter_position = (
                        leading_primer_pos[leading_idx].start,
                        trailing_primer_pos[trailing_idx].end,
                    );
                    let inner_position = (
                        leading_primer_pos[leading_idx].end,
                        trailing_primer_pos[trailing_idx].start,
                    );
                    let primer_pair = PrimerPair {
                        primer_id: pair_id.into(),
                        distance: distance,
                        outter_position: outter_position,
                        inner_position: inner_position,
                        single_end: false,
                    };
                    result.push(primer_pair);
                    pre_pair_end_postion = outter_position.1;
                }
            }

            leading_idx += left_move_step;
            trailing_idx += right_move_step;
        }
    }
    // if primer
    histogram!("adapter_match_cnts_per_read").record(result.len() as f64);

    // println!("primer pair length  {}", result.len());
    Ok(result)
}

fn primer_pair_status(
    left_primer: &PrimerCandidate,
    right_primer: &PrimerCandidate,
) -> (usize, usize) {
    // let name_len = left_primer.name.len() - 2;
    if (left_primer.name[0..left_primer.name.len() - 2]
        == right_primer.name[0..right_primer.name.len() - 2])
        && (left_primer.name[left_primer.name.len() - 2..left_primer.name.len()]
            != right_primer.name[right_primer.name.len() - 2..right_primer.name.len()])
        && (left_primer.end < right_primer.start)
    {
        return (1, 1);
    } else {
        if left_primer.end > right_primer.start {
            return (0, 1);
        } else {
            return (1, 0);
        }
    }
}
