use crate::bam_record_extention::{AuxExt, ReadRecord};
// use crate::demux_worker::PrimerPosition;
use crate::get_demuxed_reads::RecordType;
use crate::reader_worker::collect_bam_files;
use anyhow::{bail, Context, Result};
use bio::io::fasta::Reader as FastaReader;
use bio::io::fastq::Reader as FastqReader;
use bio::io::{fasta, fastq};
use rust_htslib::bam::Reader as BAMReader;
use rust_htslib::bam::{self, Read};
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufReader};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::Arc;

use serde::Serialize;

// pub type SharedWriter = Arc<Mutex<Box<dyn ReadWriter + Send>>>;

// main writer is exclusive; uncertain is shared
// pub type DoubleWriter = (Box<dyn ReadWriter + Send>, SharedWriter);

pub fn read_sequences(input_path: &str) -> Result<Vec<ReadRecord>, anyhow::Error> {
    let path: &Path = input_path.as_ref();
    // let stem = path.with_extension(""); // Remove .fasta / .fastq / .bam
    let ext = path
        .extension()
        .and_then(|s| s.to_str())
        .expect("not supported input format")
        .to_ascii_lowercase();

    match ext.as_str() {
        "fa" | "fasta" => read_fasta_sequences(input_path),
        "fq" | "fastq" => read_fastq_sequences(input_path),
        "bam" => read_bam_sequences(input_path),
        _ => bail!("Unsupported input format: {ext}"),
    }
}

pub fn read_fasta_sequences(input_path: &str) -> Result<Vec<ReadRecord>> {
    let file = File::open(input_path)?;
    let reader = FastaReader::new(BufReader::new(file));
    let mut records = Vec::new();

    for result in reader.records() {
        let record = result.map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        let mut cur_record = ReadRecord::default();
        cur_record.id = record.id().into();
        cur_record.sequence = record.seq().to_vec();
        records.push(cur_record);
    }

    Ok(records)
}

pub fn read_fastq_sequences(input_path: &str) -> Result<Vec<ReadRecord>, anyhow::Error> {
    // let file = File::open(input_path).unwrap();
    let file = File::open(input_path).with_context(|| format!("open fastq: {input_path}"))?;
    let reader = FastqReader::new(BufReader::new(file));
    let mut records = Vec::new();

    for result in reader.records() {
        let record = result.map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        let mut cur_record = ReadRecord::default();
        cur_record.id = record.id().into();
        cur_record.sequence = record.seq().to_vec();
        cur_record.quality = Some(record.qual().to_vec());
        records.push(cur_record);
    }

    Ok(records)
}

pub fn read_bam_sequences(input_path: &str) -> Result<Vec<ReadRecord>, anyhow::Error> {
    let mut bam_reader = BAMReader::from_path(input_path).unwrap();
    let mut records = Vec::new();

    // Iterate over each record in the BAM file
    for record in bam_reader.records() {
        let record = record.map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        // println!("record {:?}", record);
        let mut cur_record = ReadRecord::default();
        cur_record.id = Arc::<str>::from(String::from_utf8(record.qname().to_vec()).unwrap());
        cur_record.sequence = record.seq().as_bytes().to_vec();
        let qual_ascii: Vec<u8> = record // &bam::Record
            .qual() // &[u8] (naked PHRED)
            .to_vec();
        cur_record.quality = Some(qual_ascii);
        if let Some(dw) = record.array_u8(b"dw") {
            cur_record.dw = Some(dw);
        }
        if let Some(cr) = record.array_u8(b"cr") {
            cur_record.cr = Some(cr);
        }
        if let Some(rq) = record.float(b"rq") {
            cur_record.rq = Some(rq);
        }
        if let Some(np) = record.i32(b"np") {
            cur_record.np = Some(np);
        }
        if let Some(cx) = record.i32(b"cx") {
            cur_record.cx = Some(cx);
        }
        if let Some(ch) = record.i32(b"ch") {
            cur_record.ch = Some(ch);
        }
        if let Some(sn) = record.array_f32(b"sn") {
            cur_record.sn = Some(sn);
        }
        if let Some(rg) = record.string(b"RG") {
            cur_record.RG = Some(rg.to_string());
        }
        if let Some(be) = record.array_u32(b"be") {
            cur_record.be = Some(be);
        }
        if let Some(cq) = record.float(b"cq") {
            cur_record.cq = Some(cq);
        }
        if let Some(ar) = record.array_u32(b"ar") {
            cur_record.ar = Some(ar);
        }
        // println!("record dw {:?}", cur_record);
        // break;
        records.push(cur_record);
    }
    Ok(records)
}

pub trait ReadWriter {
    fn write_prepared_record(&mut self, record: &[RecordType]) -> io::Result<()>;
    fn write_primer_record(&mut self, record: &RecordType) -> io::Result<()>;

    /// Some formats need flush/close (bam closes on Drop),
    /// optional operations can be reserved in the trait.
    fn finish(&mut self) -> io::Result<()> {
        Ok(())
    }
}

// type DoubleWriter = (
//     Box<dyn ReadWriter>, // normal reads
//     Box<dyn ReadWriter>, // uncertain
// );
pub struct FastaWriter {
    inner: fasta::Writer<std::fs::File>,
}

impl FastaWriter {
    pub fn new<P: AsRef<Path>>(p: P) -> anyhow::Result<Self> {
        Ok(Self {
            inner: fasta::Writer::to_file(p)?,
        })
    }
}

impl Drop for FastaWriter {
    fn drop(&mut self) {
        let _ = self.inner.flush();
    }
}
pub struct FastqWriter {
    inner: fastq::Writer<std::fs::File>,
}
impl FastqWriter {
    pub fn new<P: AsRef<Path>>(p: P) -> anyhow::Result<Self> {
        Ok(Self {
            inner: fastq::Writer::to_file(p)?,
        })
    }
}

impl ReadWriter for FastqWriter {
    fn finish(&mut self) -> io::Result<()> {
        self.inner.flush()?;
        Ok(())
    }

    fn write_prepared_record(&mut self, records: &[RecordType]) -> io::Result<()> {
        for r in records {
            match r {
                RecordType::Fastq(rec) => {
                    // Method ①: most concise - lend the whole record to the writer
                    self.inner.write_record(&rec)?;
                }
                _other => {
                    // Enable during debugging
                    eprintln!("FastqWriter got non-fastq record:");
                }
            }
        }

        Ok(())
    }
    fn write_primer_record(&mut self, record: &RecordType) -> io::Result<()> {
        if let RecordType::Fastq(rec) = record {
            // Method ①: most concise - lend the whole record to the writer
            self.inner.write_record(&rec)?;
        }
        Ok(())
    }
}

impl ReadWriter for FastaWriter {
    fn write_prepared_record(&mut self, records: &[RecordType]) -> io::Result<()> {
        for r in records {
            if let RecordType::Fasta(rec) = r {
                // Method ①: most concise - lend the whole record to the writer
                self.inner.write_record(&rec)?;
            }
        }
        Ok(())
    }
    fn write_primer_record(&mut self, record: &RecordType) -> io::Result<()> {
        if let RecordType::Fasta(rec) = record {
            // Method ①: most concise - lend the whole record to the writer
            self.inner.write_record(&rec)?;
        }
        Ok(())
    }
}

impl ReadWriter for bam::Writer {
    fn write_prepared_record(&mut self, record: &[RecordType]) -> io::Result<()> {
        for rec in record {
            if let RecordType::BAM(bam_rec) = rec {
                // Method ①: most concise - lend the whole record to the writer
                self.write(&bam_rec)
                    .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;
            }
        }

        Ok(())
    }
    fn write_primer_record(&mut self, record: &RecordType) -> io::Result<()> {
        if let RecordType::BAM(rec) = record {
            // Method ①: most concise - lend the whole record to the writer
            self.write(rec)
                .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;
        }
        Ok(())
    }
}

fn sanitize_filename(s: &str) -> String {
    // Make input_name safe for use as a filename (enhance as needed)
    s.replace('/', "_").replace('\\', "_").replace(' ', "_")
}

fn _make_main_path_for_name(stem: &Path, name: &str, out_ext: &str) -> PathBuf {
    let stem_s = stem.to_string_lossy();
    let name = sanitize_filename(name);
    let base = format!("{stem_s}_{name}");
    PathBuf::from(format!("{base}.{out_ext}"))
}

pub fn make_writers(
    output_folder: &str,
    input_names: Vec<String>,
    input_file: &str,
    output_format: &str,
) -> anyhow::Result<HashMap<String, Box<dyn ReadWriter + Send>>> {
    let out_ext = match output_format {
        "fa" | "fasta" => "fa",
        "fq" | "fastq" => "fq",
        "bam" => "bam",
        other => return Err(anyhow::anyhow!("Unsupported output format: {other}")),
    };

    // Output directory
    let out_dir = Path::new(output_folder);
    std::fs::create_dir_all(out_dir)?; // Ensure existence

    let mut map: HashMap<String, Box<dyn ReadWriter + Send>> =
        HashMap::with_capacity(input_names.len());

    for name in input_names {
        // Suggest sanitizing to prevent primer names with / or spaces from causing path issues
        let safe_name = sanitize_filename(&name); // Use if already in your project; otherwise implement yourself
        let main_path = out_dir.join(format!("{safe_name}.{out_ext}"));

        let main_writer: Box<dyn ReadWriter + Send> = match out_ext {
            "fa" => Box::new(FastaWriter::new(&main_path)?),
            "fq" => Box::new(FastqWriter::new(&main_path)?),
            "bam" => {
                let header = make_header(input_file);
                let mut w = rust_htslib::bam::Writer::from_path(
                    &main_path,
                    &header,
                    rust_htslib::bam::Format::Bam,
                )?;
                w.set_threads(4)?;
                Box::new(w)
            }
            _ => unreachable!(),
        };

        // key still uses original name (doesn't affect map lookup), filename uses safe_name
        map.insert(name, main_writer);
    }

    Ok(map)
}

fn make_header<P: AsRef<Path>>(input: P) -> bam::Header {
    let path = input.as_ref();
    let cmdline: String = std::env::args().collect::<Vec<_>>().join(" ");

    let mut header: bam::Header = if path.is_dir() {
        // Recursively / non-recursively collect .bam
        let bam_files: Vec<PathBuf> = collect_bam_files(path);
        if bam_files.is_empty() {
            println!("bam path: {:?}", path);
            println!("bam files: {:?}", bam_files);
            panic!("No BAM files found");
        }
        let bam_reader = BAMReader::from_path(&bam_files[0]).expect("open first bam in folder");
        bam::Header::from_template(bam_reader.header())
    } else {
        let bam_reader = BAMReader::from_path(path).expect("open bam file");
        bam::Header::from_template(bam_reader.header())
    };

    let mut hd = bam::header::HeaderRecord::new(b"PG");
    hd.push_tag(b"PN", &"primer_demux");
    hd.push_tag(b"ID", &"v0.0.4");
    hd.push_tag(b"VN", &"v0.0.4");
    hd.push_tag(b"CL", &cmdline);
    header.push_record(&hd);

    header
}

pub enum OutputFormat {
    Jsonl,
    Tsv,
}
pub fn detect_format_from_path(path: &Path) -> io::Result<OutputFormat> {
    let fname = path
        .file_name()
        .and_then(|s| s.to_str())
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "invalid filename"))?
        .to_ascii_lowercase();

    if fname.ends_with(".jsonl") {
        Ok(OutputFormat::Jsonl)
    } else if fname.ends_with(".tsv") {
        Ok(OutputFormat::Tsv)
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "filename must end with .jsonl or .tsv",
        ))
    }
}

/// Line writer (plain text, no compression)
pub enum LineSink {
    Jsonl {
        inner: BufWriter<File>,
    },
    Tsv {
        inner: BufWriter<File>,
        wrote_header: bool,
    },
}

pub fn open_line_sink(path: &Path) -> io::Result<LineSink> {
    let fmt = detect_format_from_path(path)?;
    let f = File::create(path)?;
    let inner = BufWriter::new(f);
    Ok(match fmt {
        OutputFormat::Jsonl => LineSink::Jsonl { inner },
        OutputFormat::Tsv => LineSink::Tsv {
            inner,
            wrote_header: false,
        },
    })
}
/// TSV first line header (write once)
fn write_header_if_needed(sink: &mut LineSink) -> io::Result<()> {
    if let LineSink::Tsv {
        inner,
        wrote_header,
    } = sink
    {
        if !*wrote_header {
            inner.write_all(b"primer_id\tchannel_idx\tprimer_distance\tsingle_end\n")?;
            *wrote_header = true;
        }
    }
    Ok(())
}

#[derive(Serialize)]
pub struct CompactRec {
    pub b: String, // barcode_index
    pub c: i32,    // channel_idx
    pub d: u16,    // distance_scaled (e.g., 1.5 -> 15)
    pub s: u8,     // single_end: 1/0
}
#[inline]
pub fn to_compact_rec(
    primer_idx: String,
    channel_idx: i32,
    distance: f32,
    single_end: bool,
    distance_scale: u16, // Usually 10
) -> CompactRec {
    let d_scaled = (distance * distance_scale as f32).round() as u16;
    CompactRec {
        b: primer_idx,
        c: channel_idx,
        d: d_scaled,
        s: if single_end { 1 } else { 0 },
    }
}

/// Write a record
pub fn write_record(sink: &mut LineSink, rec: &CompactRec) -> io::Result<()> {
    match sink {
        LineSink::Jsonl { inner } => {
            serde_json::to_writer(&mut *inner, rec)?;
            inner.write_all(b"\n")?;
        }
        LineSink::Tsv { .. } => {
            write_header_if_needed(sink)?;
            if let LineSink::Tsv { inner, .. } = sink {
                // b  c  d  s
                write!(inner, "{}\t{}\t{}\t{}\n", rec.b, rec.c, rec.d, rec.s)?;
            }
        }
    }
    Ok(())
}

/// flush (optional explicit call)
pub fn flush_sink(sink: &mut LineSink) -> io::Result<()> {
    match sink {
        LineSink::Jsonl { inner } => inner.flush(),
        LineSink::Tsv { inner, .. } => inner.flush(),
    }
}

pub fn ensure_output_dir(output_folder: &str) -> anyhow::Result<()> {
    let p = Path::new(output_folder);

    if p.exists() {
        if !p.is_dir() {
            bail!(
                "output_folder exists but is not a directory: {}",
                p.display()
            );
        }
        return Ok(());
    }

    std::fs::create_dir_all(p)
        .with_context(|| format!("failed to create output_folder: {}", p.display()))?;
    Ok(())
}
