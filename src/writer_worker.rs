use crate::core_context::PrimerPair;
use crate::get_demuxed_reads::RecordType;
use crate::io_utils::{flush_sink, make_writers, open_line_sink, LineSink};

use crate::pbar::{get_spin_pb, DEFAULT_INTERVAL};
// use clap::error::ContextKind;
use crossbeam::channel::Receiver;
use tracing::error;
// use flate2::write;
use metrics::counter;
// use std::error::Error;
use std::path::{Path, PathBuf};

pub fn write_primer_demux_result_with_seperated_reads(
    input_file: &str,
    output_format: &str,
    output_folder: &str,
    primer_demux_info: Receiver<(Vec<RecordType>, Vec<RecordType>, Vec<PrimerPair>)>,
    // q_threshold: Option<u8>,
    save_primer: bool,
    primer_names: Vec<String>,
) -> anyhow::Result<()> {
    print!("enter write process");
    let path: &Path = output_folder.as_ref();

    let _stem = path.with_extension(""); // Remove .fasta / .fastq / .bam
                                        // let ext = input_path
                                        //     .extension()
                                        //     .and_then(|s| s.to_str())
                                        //     .unwrap_or("None")
                                        //     .to_ascii_lowercase();

    // Open JSONL file once only if primers need to be saved:
    // You can change the suffix to ".jsonl.gz" to enable compression
    let primer_path: Option<PathBuf> = if save_primer {
        Some(_stem.with_extension("primer.tsv"))
        // Or: Some(_stem.with_extension("barcode.jsonl.gz"))
    } else {
        None
    };
    let mut _line_sink: Option<LineSink> = if let Some(ref p) = primer_path {
        Some(open_line_sink(p)?)
    } else {
        None
    };
    let _distance_scale: u16 = 10;
    static _EPS: f32 = 0.001;

    let pb = get_spin_pb("Writing demuxed reads".to_string(), DEFAULT_INTERVAL);
    let mut writers = make_writers(output_folder, primer_names, input_file, output_format)?;
    match output_format {
        "fa" | "fasta" | "fastq" | "fq" => {
            while let Ok(demuxed_record) = primer_demux_info.recv() {
                counter!("total_channels").increment(1 as u64); // Skip if both are empty

                let (demuxed_reads, unused_reads, primer_pair) = demuxed_record;
                if demuxed_reads.is_empty() && unused_reads.is_empty() {
                    counter!("no_reads_channels").increment(1 as u64); // Skip if both are empty

                    continue;
                } else {
                    if !unused_reads.is_empty() {
                        // writers.1.write_prepared_record(subreads);
                        if let Err(e) = writers
                            .get_mut("uncertain")
                            .unwrap()
                            .write_prepared_record(unused_reads.as_slice())
                        {
                            error!(?e, "Failed to write invalid subreads file");
                        }
                    }
                    if !demuxed_reads.is_empty() {
                        for (idx, read) in demuxed_reads.iter().enumerate() {
                            if let Err(e) = writers
                                .get_mut(&primer_pair[idx].primer_id)
                                .unwrap()
                                .write_primer_record(read)
                            {
                                error!(?e, "Failed to write valid subreads file");
                            }
                            counter!("valid_reads").increment(1 as u64);
                        }
                    }
                }
                pb.inc(1);

                // pb.finish_with_message("Finished writing demuxed reads");
            }
        }
        "bam" => {
            while let Ok(demuxed_record) = primer_demux_info.recv() {
                let (demuxed_reads, unused_reads, primer_pair) = demuxed_record;
                if demuxed_reads.is_empty() && unused_reads.is_empty() {
                    continue; // Skip if both are empty
                } else {
                    if !unused_reads.is_empty() {
                        // writers.1.write_prepared_record(subreads);
                        if let Err(e) = writers
                            .get_mut("uncertain")
                            .unwrap()
                            .write_prepared_record(unused_reads.as_slice())
                        {
                            error!(?e, "Failed to write valid subreads BAM");
                        }
                    }
                    if !demuxed_reads.is_empty() {
                        for (idx, read) in demuxed_reads.iter().enumerate() {
                            let _ = writers
                                .get_mut(&primer_pair[idx].primer_id)
                                .unwrap()
                                .write_primer_record(read);
                        }
                    }
                }
                pb.inc(1);
            }
        }
        _ => anyhow::bail!("Unsupported write format: {output_format}"),
    }

    Ok(())
}
