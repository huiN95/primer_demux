// mod alignment_statistics;
mod bam_record_extention;
mod cli;
mod find_pattern;
mod get_demuxed_reads;
mod io_utils;
mod pbar;
mod reader_worker;
mod run_logger;
// mod statitcs;
mod core_context;
mod demux_pipeline_v1;
mod demux_primer;
mod primer_demuxer_v1;
mod primer_myers;
mod writer_worker;
use clap::Parser;
use cli::Cli;

use run_logger::{init_metrics, init_tracing_log};

use crate::demux_pipeline_v1::demux_pipeline_v1;
// use crate::io_utils::ensure_output_dir;
use tracing::info;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // 1. Parse command line arguments

    let cli = Cli::parse();
    let log_path = &cli.log_folder;

    let _log_writer = init_tracing_log(&cli);
    // let metric_writer = start_metrics_file_writer(&log_path);
    let _metric_writer = init_metrics(&log_path);
    info!("Start processing");
    info!(?cli, "parsed CLI");
    match cli.pipeline_version.as_str() {
        "1" => {
            info!("Using demux pipeline version 1");
            demux_pipeline_v1(&cli).unwrap();
        }
        other => panic!("Unsupported pipeline version: {}", other),
    }
    info!("End processing");
    Ok(())
}

#[cfg(test)]
mod test {
    use clap::Parser;

    use super::*;

    #[test]
    fn test_fasta_format() {
        let args = vec![
            "my_program",
            "--input_file",
            "/mnt/data_7t/adam/primer_demux/simulated_seq/ground_truth_mg1655.fasta",
            "--primer",
            "/mnt/data_7t/adam/adapter_bcs/primers_40.fasta",
            "--output_folder",
            "/mnt/data_7t/adam/primer_demux/simulated_seq/tmp",
            "--log_folder",
            "/mnt/data_7t/adam/primer_demux/simulated_seq/tmp",
            "--tail_cutoff",
            "0",
            "--max_distance",
            "3",
        ];
        let matches = Cli::parse_from(args);
        let result = demux_pipeline_v1(&matches);
        assert!(result.is_ok(), "Pipeline failed: {:?}", result.err());
    }

    #[test]
    fn test_fastq_format() {
        let args = vec![
            "my_program",
            "--input_file",
            "/mnt/data_7t/adam/primer_demux/temp1.fastq",
            "--primer",
            "/mnt/data_7t/adam/adapter_bcs/primer40_v2.fasta",
            "--output_folder",
            "/mnt/data_7t/adam/primer_demux/simulated_seq/tmp",
            "--log_folder",
            "/mnt/data_7t/adam/primer_demux/simulated_seq/tmp",
            "--tail_cutoff",
            "3",
            "--max_distance",
            "2",
            "--output_format",
            "fastq",
            "--keep_primer",
        ];
        let matches = Cli::parse_from(args);
        let result = demux_pipeline_v1(&matches);
        assert!(result.is_ok(), "Pipeline failed: {:?}", result.err());
    }
}
