# Primer Demux (v0.0.4)

`primer_demux` is a high-performance tool written in Rust designed for demultiplexing sequencing reads based on primer sequences (including connected reads). It supports multiple input and output formats (FASTA, FASTQ, BAM) and leverages the Myers bit-parallel string matching algorithm for fast approximate matching.

## Key Features

- **Approximate Matching**: Uses the Myers bit-parallel algorithm to handle mismatches and indels in primer sequences, and support connected reads.
- **Multi-format Support**: Process sequences in FASTA, FASTQ, and BAM formats.

- **Flexible Configuration**: Control edit distance, primer retention, tail cutoffs, and minimum subread lengths.
- **Logging & Metrics**: Integrated tracing logs and metrics for monitoring and debugging.

## Installation

### Prerequisites

- Rust (MSRV 1.70+)
- `libhts` (for BAM support, usually provided by `rust-htslib`)

### Build from Source

```bash
git clone <repository_url>
cd primer_demux
cargo build --release
```

The binary will be available at `target/release/primer_demux`.

## Usage

### Basic Command

```bash
primer_demux -i input.bam -p primers.fasta -o output_dir --log_folder logs --keep_primer
```

### CLI Arguments

| Argument | Short | Default | Description |
| :--- | :--- | :--- | :--- |
| `--input_file` | `-i` | (Required) | Path to input sequencing file (FASTA/FASTQ/BAM). |
| `--primer` | `-p` | (Required) | Path to primer sequence file (FASTA). |
| `--output_folder` | `-o` | (Required) | Directory to save demultiplexed files. |
| `--log_folder` | | (Required) | Directory to store logs and metrics. |
| `--max_distance` | `-d` | `1` | Maximum edit distance allowed for primer matching. |
| `--output_format` | | `fasta` | Output format: `fasta`, `fastq`, `bam`, `fq`, `fa`. |
| `--keep_primer` | | `false` | Whether to keep the primer sequence in the output reads. |
| `--min_subread_len`| `-l` | `50` | Minimum length of the subread after primer removal. |
| `--tail_cutoff` | | `2` | Cutoff for searching primers near the ends of sequences, used for connected reads when the primers might lose some bp when they are connected. |
| `--min_q` | | `20` | Minimum average quality threshold (implementation pending). |
| `--threads`| | `0` | CPU threads used for demux. |
| `--pipeline_version` | | `1` | Pipeline logic version. |

## How it Works

1. **Primer Indexing**: Reads the primer FASTA file and builds Myers matching patterns.
2. **Approximate Search**: For each read, the tool searches for the specified primers within the search bounds, allowing for mismatches up to the specified `max_distance`.
3. **Demultiplexing**: Reads are classified based on the matched primers, if the reads contain multiple connected reads, it will also demux them.
4. **Subread Extraction**: Extracts the sequence between identified primers (or from the primer to the end), optionally keeping the primer sequence itself.

## Example
1. **Primer fasta**: "_L" and "_T" are the key words to point out the leading and trailing primer sequences.See [primer.fasta](examples/primers.fasta).


## License

This project is licensed under the [MIT License](LICENSE).

## Contact

[Your Name/Team] - [Email/GitLab Profile]
