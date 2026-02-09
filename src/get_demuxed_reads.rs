use crate::bam_record_extention::ReadRecord;
use crate::core_context::{PrimerCandidate, PrimerPair};

use bio::io::{fasta, fastq};

use rust_htslib::bam::{
    self,
    record::{Aux, Record},
};
// use std::cmp::min;
use std::{ops::Deref, path::Path};
// use tracing_subscriber::layer::SubscriberExt;
pub enum RecordType {
    Fastq(fastq::Record),
    Fasta(fasta::Record),
    BAM(Record),
}

pub fn prepare_record_to_writer(
    rec: &ReadRecord,
    primer_paris: Vec<PrimerPair>,
    min_subread_len: usize,
    output_folder: &str,
    keep_primer_flag: bool,
    output_format: &str,
) -> anyhow::Result<(Vec<RecordType>, Vec<RecordType>, Vec<PrimerPair>)> {
    // let primer_records = get_primer_record(rec, &primer_paris, output_file)?;
    let subread_records = get_subread_record_with_keep_primer_option(
        rec,
        &primer_paris,
        min_subread_len,
        // None, // 质量阈值
        output_folder,
        keep_primer_flag,
        output_format,
    )?;

    Ok((subread_records.0, subread_records.1, primer_paris))
    // Ok((primer_records, subread_records))
}

fn make_name(rec_id: &str, primer_id: &str, s: usize, e: usize) -> String {
    // 预估容量：id+primer+数字
    let mut name = String::with_capacity(rec_id.len() + primer_id.len() + 32);
    use std::fmt::Write;
    write!(&mut name, "{}:{}/{}-{}", rec_id, primer_id, s, e).unwrap();
    name
}

fn get_subread_record_with_keep_primer_option(
    rec: &ReadRecord,
    primer_position: &Vec<PrimerPair>,
    min_subread_len: usize,
    // q_threshold: Option<u8>,
    output_file: &str,
    keep_primer: bool,
    output_format: &str,
) -> anyhow::Result<(Vec<RecordType>, Vec<RecordType>)> {
    // let path: &Path = output_file.as_ref();
    // let _stem = path.with_extension(""); // 去掉 .fasta / .fastq / .bam
    // let ext = path
    //     .extension()
    //     .and_then(|s| s.to_str())
    //     .unwrap_or("fasta")
    //     .to_ascii_lowercase();
    let n = primer_position.len();

    let mut demuxed_reads = Vec::with_capacity(n + 1);
    let mut unused_reads = Vec::with_capacity(n + 1);

    match output_format {
        // ---------- FASTA ----------
        "fa" | "fasta" => {
            let mut pre_end = 0;
            for (idx, pos) in primer_position.iter().enumerate() {
                let start = pos.outter_position.0;
                let end = pos.outter_position.1;
                if pre_end + min_subread_len < start {
                    // let name = format!("{}:{}/{}-{}", rec.id, pos.primer_id, pre_end, start);
                    let name = make_name(&rec.id, &pos.primer_id, pre_end, start);
                    let fasta_rec = fasta::Record::with_attrs(
                        &name,                         // &str
                        None,                          // description/comment
                        &rec.sequence[pre_end..start], // &[u8]
                    );
                    unused_reads.push(RecordType::Fasta(fasta_rec));
                }
                if keep_primer {
                    let name = make_name(&rec.id, &pos.primer_id, start, end);

                    let fasta_rec = fasta::Record::with_attrs(
                        &name,                     // &str
                        None,                      // description/comment
                        &rec.sequence[start..end], // &[u8]
                    );
                    demuxed_reads.push(RecordType::Fasta(fasta_rec));
                } else {
                    // 读长不够
                    let name = format!(
                        "{}:{}/{}-{}",
                        rec.id, pos.primer_id, pos.inner_position.0, pos.inner_position.1
                    );
                    let fasta_rec = fasta::Record::with_attrs(
                        &name,                                                     // &str
                        None, // description/comment
                        &rec.sequence[pos.inner_position.0..pos.inner_position.1], // &[u8]
                    );
                    demuxed_reads.push(RecordType::Fasta(fasta_rec));
                }
                if idx == primer_position.len() - 1 {
                    if rec.sequence.len() > end + min_subread_len {
                        let name = format!(
                            "{}:{}/{}-{}",
                            rec.id,
                            pos.primer_id,
                            pos.outter_position.1,
                            &rec.sequence.len()
                        );
                        let fasta_rec = fasta::Record::with_attrs(
                            &name,                                                    // &str
                            None, // description/comment
                            &rec.sequence[pos.outter_position.1..rec.sequence.len()], // &[u8]
                        );
                        unused_reads.push(RecordType::Fasta(fasta_rec));
                    }
                }
                pre_end = end;
            }
        }
        // ---------- FASTQ ----------
        "fq" | "fastq" => {
            let mut pre_end = 0;
            for (idx, pos) in primer_position.iter().enumerate() {
                let start = pos.outter_position.0;
                let end = pos.outter_position.1;
                // 头部的数据够长才行
                if pre_end + min_subread_len < start {
                    let name = format!("{}:{}/{}-{}", rec.id, pos.primer_id, pre_end, start);

                    let fastq_rec = fastq::Record::with_attrs(
                        &name,                         // &str
                        None,                          // description/comment
                        &rec.sequence[pre_end..start], // &[u8]
                        rec.quality
                            .as_deref() // Option<&[u8]>
                            .expect("missing quality")
                            .get(pre_end..start) // Option<&[u8]> 取子片段
                            .ok_or_else(|| anyhow::anyhow!("bad range: {}..{}", pre_end, start))?, // quality: Option<&[u8]>
                                                                                                   // .as_ref(),
                    );

                    unused_reads.push(RecordType::Fastq(fastq_rec));
                }
                if keep_primer {
                    let name = format!("{}:{}/{}-{}", rec.id, pos.primer_id, start, end);

                    let fastq_rec = fastq::Record::with_attrs(
                        &name,                     // &str
                        None,                      // description/comment
                        &rec.sequence[start..end], // &[u8]
                        rec.quality
                            .as_deref() // Option<&[u8]>
                            .expect("missing quality")
                            .get(start..end) // Option<&[u8]> 取子片段
                            .expect("bad range"), // quality: Option<&[u8]>
                    );

                    demuxed_reads.push(RecordType::Fastq(fastq_rec));
                } else {
                    let name = format!(
                        "{}:{}/{}-{}",
                        rec.id, pos.primer_id, pos.inner_position.0, pos.inner_position.1
                    );
                    let fastq_rec = fastq::Record::with_attrs(
                        &name, // &str
                        None,  // description/comment
                        &rec.sequence[pos.inner_position.0..pos.inner_position.1],
                        rec.quality
                            .as_deref() // Option<&[u8]>
                            .expect("missing quality")
                            .get(pos.inner_position.0..pos.inner_position.1) // Option<&[u8]> 取子片段
                            .ok_or_else(|| {
                                anyhow::anyhow!(
                                    "bad range: {}..{}",
                                    pos.inner_position.0,
                                    pos.inner_position.1
                                )
                            })?, // quality: Option<&[u8]>
                    ); // &[u8]

                    demuxed_reads.push(RecordType::Fastq(fastq_rec));
                }
                if idx == primer_position.len() - 1 {
                    // 最后一段够长才行，否则就不要了
                    if rec.sequence.len() > pos.outter_position.1 + min_subread_len {
                        let name = format!(
                            "{}:{}/{}-{}",
                            rec.id,
                            pos.primer_id,
                            pos.outter_position.1,
                            rec.sequence.len()
                        );
                        let fastq_rec = fastq::Record::with_attrs(
                            &name, // &str
                            None,  // description/comment
                            &rec.sequence[pos.outter_position.1..rec.sequence.len()],
                            rec.quality
                                .as_deref() // Option<&[u8]>
                                .expect("missing quality")
                                .get(pos.outter_position.1..rec.sequence.len()) // Option<&[u8]> 取子片段
                                .ok_or_else(|| {
                                    anyhow::anyhow!(
                                        "bad range: {}..{}",
                                        pos.outter_position.1,
                                        rec.sequence.len()
                                    )
                                })?, // quality: Option<&[u8]> // &[u8]
                        );
                        unused_reads.push(RecordType::Fastq(fastq_rec));
                    }
                }
                pre_end = end;
            }
        }
        // ---------- BAM ----------
        "bam" => {
            // let channel_num = rec.id.split('_').nth(1).unwrap_or("0").parse().unwrap_or(0);
            let channel_num = rec.ch.unwrap_or(0);

            let mut pre_end = 0;
            // let mut cnt: i32 = 0;

            // 头部的数据够长才行
            for (idx, pos) in primer_position.iter().enumerate() {
                let start = pos.outter_position.0;
                let end = pos.outter_position.1;
                if pre_end + min_subread_len < start {
                    let name = format!("{}/{}/{}/{}", rec.id, channel_num, "subread", idx);
                    let mut bam_rec = bam::Record::new();

                    bam_rec.set(
                        name.as_bytes(),
                        None,
                        &rec.sequence[pre_end..start],
                        rec.quality // &Option<Vec<u8>>
                            .as_ref() // Option<&Vec<u8>>
                            .expect("missing quality")[pre_end..start]
                            .as_ref(),
                    ); // 第 2 个参数是 CIGAR 向量
                    bam_rec.set_flags(0x4);
                    // bam_rec.push_aux(b"np", Aux::I8(1)).unwrap();
                    // bam_rec.push_aux(b"cx", Aux::I8(3)).unwrap();
                    bam_rec.push_aux(b"ch", Aux::I32(channel_num)).unwrap();
                    bam_rec.push_aux(b"RG", Aux::String("0425")).unwrap();
                    bam_rec.push_aux(b"rq", Aux::Float(0.8f32)).unwrap();
                    bam_rec
                        .push_aux(b"cq", Aux::Float(rec.cq.unwrap_or(0.8f32)))
                        .unwrap();
                    // bam_rec
                    //     .push_aux(
                    //         b"sn",
                    //         Aux::ArrayFloat((&rec.sn.as_ref().unwrap()[..]).into()),
                    //     )
                    //     .unwrap();
                    // bam_rec
                    //     .push_aux(b"cr", Aux::ArrayU8((&vec![0u8; start - end][..]).into()))
                    //     .unwrap();
                    // bam_rec
                    //     .push_aux(b"dw", Aux::ArrayU8((&dw[start..end]).into()))
                    //     .unwrap();

                    let start_end_position = vec![start as u32, end as u32];

                    bam_rec
                        .push_aux(b"be", Aux::ArrayU32((&start_end_position[..]).into()))
                        .unwrap();
                    unused_reads.push(RecordType::BAM(bam_rec));
                }
                if keep_primer {
                    let name = format!("{}:{}/{}-{}", rec.id, pos.primer_id, start, end);

                    let mut bam_rec = bam::Record::new();

                    bam_rec.set(
                        name.as_bytes(),
                        None,
                        &rec.sequence[start..end],
                        rec.quality // &Option<Vec<u8>>
                            .as_ref() // Option<&Vec<u8>>
                            .expect("missing quality")[start..end]
                            .as_ref(),
                    ); // 第 2 个参数是 CIGAR 向量
                    bam_rec.set_flags(0x4);
                    // bam_rec.push_aux(b"np", Aux::I8(1)).unwrap();
                    // bam_rec.push_aux(b"cx", Aux::I8(3)).unwrap();
                    bam_rec.push_aux(b"ch", Aux::I32(channel_num)).unwrap();
                    bam_rec.push_aux(b"RG", Aux::String("0425")).unwrap();
                    bam_rec.push_aux(b"rq", Aux::Float(0.8f32)).unwrap();
                    bam_rec
                        .push_aux(b"cq", Aux::Float(rec.cq.unwrap_or(0.8f32)))
                        .unwrap();
                    bam_rec
                        .push_aux(
                            b"sn",
                            Aux::ArrayFloat((&rec.sn.as_ref().unwrap()[..]).into()),
                        )
                        .unwrap();
                    // bam_rec
                    //     .push_aux(b"cr", Aux::ArrayU8((&vec![0u8; start - end][..]).into()))
                    //     .unwrap();
                    // bam_rec
                    //     .push_aux(b"dw", Aux::ArrayU8((&dw[start..end]).into()))
                    //     .unwrap();

                    // let start_end_position = vec![start as u32, end as u32];
                    let be = [start as u32, end as u32];

                    bam_rec
                        .push_aux(b"be", Aux::ArrayU32((&be[..]).into()))
                        .unwrap();
                    demuxed_reads.push(RecordType::BAM(bam_rec));
                } else {
                    let name = format!(
                        "{}:{}/{}-{}",
                        rec.id, pos.primer_id, pos.inner_position.0, pos.inner_position.1
                    );
                    let mut bam_rec = bam::Record::new();
                    bam_rec.set(
                        name.as_bytes(),
                        None,
                        &rec.sequence[pos.inner_position.0..pos.inner_position.1],
                        rec.quality // &Option<Vec<u8>>
                            .as_ref() // Option<&Vec<u8>>
                            .expect("missing quality")[pos.inner_position.0..pos.inner_position.1]
                            .as_ref(),
                    ); // 第 2 个参数是 CIGAR 向量

                    bam_rec.set_flags(0x4);
                    // bam_rec.push_aux(b"np", Aux::I8(1)).unwrap();

                    // bam_rec.push_aux(b"cx", Aux::I8(3)).unwrap();
                    bam_rec.push_aux(b"ch", Aux::I32(channel_num)).unwrap();
                    bam_rec.push_aux(b"RG", Aux::String("0425")).unwrap();
                    // bam_rec
                    //     .push_aux(b"rq", Aux::Float(rec.rq.unwrap_or(0.8f32)))
                    //     .unwrap();
                    bam_rec.push_aux(b"rq", Aux::Float(0.8f32)).unwrap();
                    bam_rec
                        .push_aux(b"cq", Aux::Float(rec.cq.unwrap_or(0.8f32)))
                        .unwrap();

                    // let start_end_position = vec![start as u32, end as u32];
                    let be = vec![start as u32, end as u32];

                    bam_rec
                        .push_aux(b"be", Aux::ArrayU32((&be[..]).into()))
                        .unwrap();

                    demuxed_reads.push(RecordType::BAM(bam_rec));
                }
                if idx == primer_position.len() - 1 {
                    let tail_start = pos.outter_position.1;
                    if rec.sequence.len() < tail_start + min_subread_len {
                        // 最后一段够长才行，否则就不要了
                        continue;
                    }
                    let name = format!(
                        "{}:{}/{}-{}",
                        rec.id,
                        pos.primer_id,
                        pos.outter_position.1,
                        rec.sequence.len()
                    );
                    let mut bam_rec = bam::Record::new();
                    bam_rec.set(
                        name.as_bytes(),
                        None,
                        &rec.sequence[pos.outter_position.1..rec.sequence.len()],
                        rec.quality // &Option<Vec<u8>>
                            .as_ref() // Option<&Vec<u8>>
                            .expect("missing quality")[pos.outter_position.1..rec.sequence.len()]
                            .as_ref(),
                    ); // 第 2 个参数是 CIGAR 向量

                    bam_rec.set_flags(0x4);
                    // bam_rec.push_aux(b"np", Aux::I8(1)).unwrap();

                    // bam_rec.push_aux(b"cx", Aux::I8(3)).unwrap();
                    bam_rec.push_aux(b"ch", Aux::I32(channel_num)).unwrap();
                    bam_rec.push_aux(b"RG", Aux::String("0425")).unwrap();
                    // bam_rec
                    //     .push_aux(b"rq", Aux::Float(rec.rq.unwrap_or(0.8f32)))
                    //     .unwrap();
                    bam_rec.push_aux(b"rq", Aux::Float(0.8f32)).unwrap();
                    bam_rec
                        .push_aux(b"cq", Aux::Float(rec.cq.unwrap_or(0.8f32)))
                        .unwrap();
                    bam_rec
                        .push_aux(
                            b"sn",
                            Aux::ArrayFloat((&rec.sn.as_ref().unwrap()[..]).into()),
                        )
                        .unwrap();

                    // let start_end_position =
                    //     vec![pos.outter_position.1 as u32, rec.sequence.len() as u32];
                    let be: [u32; 2] = [pos.outter_position.1 as u32, rec.sequence.len() as u32];

                    bam_rec
                        .push_aux(b"be", Aux::ArrayU32((&be[..]).into()))
                        .unwrap();
                    unused_reads.push(RecordType::BAM(bam_rec));
                }
                pre_end = end;
            }
        }

        // ---------- 其他格式 ----------
        _ => anyhow::bail!("不支持的输出格式: {output_format}"),
    }
    Ok((demuxed_reads, unused_reads))
}

fn slice_seq_qual<'a>(
    rec: &'a ReadRecord,
    s: usize,
    e: usize,
) -> anyhow::Result<(&'a [u8], &'a [u8])> {
    let seq = rec
        .sequence
        .get(s..e)
        .ok_or_else(|| anyhow::anyhow!("bad seq range {s}..{e}"))?;
    let qual_all = rec
        .quality
        .as_deref()
        .ok_or_else(|| anyhow::anyhow!("missing quality"))?;
    let qual = qual_all
        .get(s..e)
        .ok_or_else(|| anyhow::anyhow!("bad qual range {s}..{e}"))?;
    anyhow::ensure!(
        seq.len() == qual.len(),
        "seq/qual length mismatch at {s}..{e}"
    );
    Ok((seq, qual))
}
