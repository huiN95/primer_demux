use crate::primer_demuxer_v1::primer_demuxer_v1;
use crate::primer_myers::{get_myers_from_primers, Direction, MayersPattern};

use crate::bam_record_extention::ReadRecord;
// use crate::core_context::{Annotated, PrimerCandidate, PrimerMeta};
// use crate::find_pattern::merge_non_overlapping_no_copy;
use crate::get_demuxed_reads::{prepare_record_to_writer, RecordType};
use crossbeam::channel::{Receiver, Sender};

use metrics::{self, counter, histogram};

use crate::core_context::PrimerPair;

use tracing::{debug, info, level_filters};

pub fn demux_reads_by_primer(
    patterns: &[ReadRecord],
    max_distance: u8,
    receiver: Receiver<ReadRecord>,
    sender: Sender<(Vec<RecordType>, Vec<RecordType>, Vec<PrimerPair>)>,
    ouput_folder: &str,
    min_subread_len: usize,
    keep_primer_flag: bool,
    tail_cutoff: u8,
    output_format: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut primer_f_myers = get_myers_from_primers(patterns, Direction::Forward, tail_cutoff);
    let mut primer_r_myers = get_myers_from_primers(patterns, Direction::Reverse, tail_cutoff);
    // println!("{:?}", barcode_f_myers.keys());
    for current_record in receiver {
        // let mut current_primer_distance: u8 = max_primer_distance;
        // // TODO Limit read length here
        //     counter!("filtered_reads_too_short").increment(1 as u64);

        // print!("max primer tolerance: {}", max_primer_distance);
        let primer_info = primer_demuxer_v1(
            min_subread_len,
            max_distance,
            &mut primer_f_myers,
            &mut primer_r_myers,
            &current_record,
        )
        .unwrap();
        let demuxed_reads = prepare_record_to_writer(
            &current_record,
            primer_info,
            min_subread_len,
            ouput_folder,
            keep_primer_flag,
            output_format,
        )?;
        // println!(demuxed_reads);
        if let Err(e) = sender.send(demuxed_reads) {
            eprintln!("Failed to send demux: {}", e);
            // Optionally return early or break here
        }
        // }
    }
    Ok(())
}
