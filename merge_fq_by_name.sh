#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

IN_BASE="/mnt/data_7t/adam/primer_demux/primer_demuxed"
OUT_DIR="/mnt/data_7t/adam/primer_demux/merged_fq"   # 改成你想输出的目录

mkdir -p "$OUT_DIR"

# 1) 收集所有 fq 的“文件名”（不含路径），去重
mapfile -t names < <(find "$IN_BASE" -type f -name '*.fq' -printf '%f\n' | sort -u)

# 2) 对每个文件名，把所有同名文件按路径排序后依次拼接
for name in "${names[@]}"; do
  out="$OUT_DIR/$name"
  : > "$out"  # 清空/创建输出文件

  # 找到所有同名文件（按路径排序保证可重复）
  mapfile -t files < <(find "$IN_BASE" -type f -name "$name" | sort)

  echo "Merging ${#files[@]} files -> $out"
  for f in "${files[@]}"; do
    cat "$f" >> "$out"
  done
done
