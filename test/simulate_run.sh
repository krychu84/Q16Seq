#!/usr/bin/env bash
# Simulate a sequencing run: copy fastq.gz files into a watch folder in batches,
# the way MinKNOW does (write to a hidden temp name, then rename).
#
# usage: simulate_run.sh <source_dir_with_fastq.gz> <watch_dir> [files_per_batch=10] [pause_seconds=120]
set -euo pipefail
src=${1:?source folder}; dst=${2:?watch folder}; per=${3:-10}; pause=${4:-120}
mkdir -p "$dst"
i=0
for f in $(ls "$src"/*.fastq.gz | sort -V); do
    b=$(basename "$f")
    cp -L "$f" "$dst/.$b.tmp" && mv "$dst/.$b.tmp" "$dst/$b"
    i=$((i + 1))
    if (( i % per == 0 )); then
        echo "$(date +%T) dropped $i files"
        sleep "$pause"
    fi
done
echo "$(date +%T) simulation done ($i files)"
