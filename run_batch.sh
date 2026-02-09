#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

PRIMER="/mnt/data_7t/adam/adapter_bcs/primer40_v2.fasta"
IN_DIR="/mnt/data_7t/adam/primer_demux/barcode_filter_reads"
OUT_BASE="/mnt/data_7t/adam/primer_demux/primer_demuxed"   # 输出放在同一个父目录下

MAX_DISTANCE=2
TAIL_CUTOFF=0
KEEP_PRIMER="--keep_primer"
OUTPUT_FORMAT="fq"

# 可选：如果你想跳过已经完成的样本，可以用一个 done 标记文件
# DONE_FILE=".done"

# 支持 fastq / fq / fastq.gz / fq.gz
files=( "$IN_DIR"/*.fastq "$IN_DIR"/*.fq "$IN_DIR"/*.fastq.gz "$IN_DIR"/*.fq.gz )

for f in "${files[@]}"; do
  # 文件名（去掉路径）
  base="$(basename "$f")"
  # 去掉压缩后缀
  base="${base%.gz}"
  # 去掉 fastq/fq 后缀
  name="${base%.fastq}"
  name="${name%.fq}"

  out_dir="$OUT_BASE/$name"
  log_dir="$out_dir"   # 你原来 log_folder 和 output_folder 一样，这里也保持一致

  mkdir -p "$out_dir"

  echo "==> Processing: $f"
  cargo run --release -- \
    --primer "$PRIMER" \
    --max_distance "$MAX_DISTANCE" \
    --input_file "$f" \
    --output_folder "$out_dir" \
    --log_folder "$log_dir" \
    --tail_cutoff "$TAIL_CUTOFF" \
    $KEEP_PRIMER \
    --output_format "$OUTPUT_FORMAT"
done
