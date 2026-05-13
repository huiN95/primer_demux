use crate::cli::Cli;

use crate::core_context::PrimerPair;
use crate::demux_primer::demux_reads_by_primer;
use crate::find_pattern::get_pattern_keys;
use crate::get_demuxed_reads::RecordType;
use crate::writer_worker::write_primer_demux_result_with_seperated_reads;
use crossbeam::channel::{Receiver, Sender};
use std::error::Error;

use std::thread;
use tracing::info;

// use tracing_appender;
// use tracing_subscriber::{fmt, EnvFilter};
use crate::io_utils::{ensure_output_dir, read_sequences};
use crate::reader_worker::read_sequences_to_queue;
use crossbeam::channel::bounded;

pub fn demux_pipeline_v1(cli: &Cli) -> Result<(), Box<dyn Error>> {
    ensure_output_dir(&cli.output_folder)?;
    let patterns = read_sequences(&cli.primer)?;

    let output_names = get_pattern_keys(&patterns);

    type WriteMsg = (Vec<RecordType>, Vec<RecordType>, Vec<PrimerPair>);
    let queue_len = 100000;
    std::thread::scope(|scope| {
        let (seq_sender, seq_receiver) = bounded(queue_len);
        // let (primer_sender, primer_receiver) = bounded(100000);
        let (primer_sender, primer_receiver): (Sender<WriteMsg>, Receiver<WriteMsg>) =
            bounded(queue_len);
        // let
        scope.spawn(move || {
            if let Err(e) = read_sequences_to_queue(&cli.input_file, seq_sender) {
                eprintln!("[Producer] Error: {}", e);
            }
            println!("[Producer] Completed sending.");
        });

        // b) Multiple demux threads
        // ---------------------------
        // let worker_count = 20;

        let worker_count: usize = match cli.threads {
            0 => std::thread::available_parallelism()
                .map(|n| n.get())
                .unwrap_or(1)
                .saturating_sub(2)
                .max(1),
            n => n as usize,
        };

        println!("[Producer] Starting {} demux threads...\n", worker_count);
        let patterns = &patterns;

        for _ in 0..worker_count {
            // If patterns are large, Arc::clone is needed; patterns.clone() is also fine

            let seq_receiver_clone = seq_receiver.clone();
            let primer_sender_clone = primer_sender.clone();
            let max_distance = cli.max_distance;

            scope.spawn(
                move || {
                    println!("Demuxing with primers...");
                    if let Err(e) = demux_reads_by_primer(
                        patterns,
                        max_distance,
                        seq_receiver_clone,
                        primer_sender_clone,
                        &cli.output_folder,
                        cli.min_subread_len,
                        cli.keep_primer,
                        cli.tail_cutoff,
                        &cli.output_format,
                    ) {
                        eprintln!("[Demux] error: {e}");
                    }

                    // println!("[Demux] finished");
                }, // Similarly, if the writer needs to know no more data is available, 
                   // drop(primer_sender_clone) in the last thread before all demux ends.
            );
        }
        drop(seq_receiver);
        drop(primer_sender);
        // print!(
        //     "save barcode flag {:?} after multiprcoessing",
        //     cli.keep_primer
        // );

        write_primer_demux_result_with_seperated_reads(
            &cli.input_file,
            &cli.output_format,
            &cli.output_folder,
            primer_receiver,
            // Some(cli.q_threshold),
            cli.keep_primer,
            output_names,
        )
        .unwrap();
    });

    Ok(())
}
