#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

IN_BASE="/mnt/data_7t/adam/primer_demux/primer_demuxed"
OUT_DIR="/mnt/data_7t/adam/primer_demux/merged_fq"   # Change to your desired output directory

mkdir -p "$OUT_DIR"

# 1) Collect all fq "filenames" (excluding path), and remove duplicates
mapfile -t names < <(find "$IN_BASE" -type f -name '*.fq' -printf '%f\n' | sort -u)

# 2) For each filename, concatenate all files with the same name sorted by path
for name in "${names[@]}"; do
  out="$OUT_DIR/$name"
  : > "$out"  # Clear/create output file

  # Find all files with the same name (sorted by path for reproducibility)
  mapfile -t files < <(find "$IN_BASE" -type f -name "$name" | sort)

  echo "Merging ${#files[@]} files -> $out"
  for f in "${files[@]}"; do
    cat "$f" >> "$out"
  done
done
