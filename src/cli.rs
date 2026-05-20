
use clap::Parser;

/// Define CLI structure
#[derive(Parser, Debug)]
#[command(
    name = "primer demux",
    version = "v0.0.4",
    about = "Approximate matching of given patterns in sequence files"
)]
pub struct Cli {

    #[arg(long = "pipeline-version", alias = "pipeline_version", default_value ="1")]
    pub pipeline_version: String,

    /// Specify primer sequence file (FASTA format)
    #[arg(short = 'p', long = "primer")]
    pub primer: String,


    // Maximum error rate allowed when matching patterns
    #[arg(short = 'd', long = "max-distance", alias = "max_distance",
          value_parser = clap::value_parser!(u8),
          default_value_t = 1)]
    pub max_distance: u8,


    /// Specify input sequence file (FASTA/FASTQ/BAM)
    #[arg(short = 'i', long = "input-file", alias = "input_file")]
    pub input_file: String,

    /// Set log level (optional), e.g., "info", "debug", "trace"
    #[arg(long, env = "RUST_LOG")]
    pub log: Option<String>,

    #[arg(short = 'o', long = "output-folder", alias = "output_folder")]
    pub output_folder: String,
    
    #[arg(long = "log-folder", alias = "log_folder")]
    pub log_folder: String,

    // Whether to keep primer sequences in reads, default is false
    #[arg(long = "keep-primer", alias = "keep_primer", default_value_t = false)]
    pub keep_primer: bool,

    #[arg(long = "tail-cutoff", alias = "tail_cutoff",
          value_parser = clap::value_parser!(u8).range(0..=6),
          default_value_t = 2)]
    pub tail_cutoff: u8,

    #[arg(long = "threads", alias = "threads",
          value_parser = clap::value_parser!(u8).range(0..=6),
          default_value_t = 2)]
    pub threads:u8,
    
    #[arg(long = "output-format", alias = "output_format", default_value ="fasta" ,       
        value_parser = ["fasta", "fastq", "bam","fq", "fa"])]
    pub output_format: String,
   
    /// Minimum subread length, default 50. If primers are kept, they are included in this length, so increase this parameter accordingly.
    #[arg(short = 'l', 
    long = "min-subread-len", 
    alias = "min_subread_len",
    default_value_t = 50,
    value_parser = clap::value_parser!(usize),)]
    pub min_subread_len: usize,

    /// Only process reads exceeding this threshold for demux, otherwise pass. Not yet implemented.
    #[arg(long = "min-q", alias = "min_q",
          value_parser = clap::value_parser!(u8).range(0..=60),
          default_value_t = 20)]
    pub min_q: u8,

}
