#!/bin/bash

count_alignments_bam() {
    local threads="$1"; local bam="$2"

    [[ ! -f "$bam" ]] && {
      echo "[ERROR] BAM not found: $bam" >&2
      echo 0; return 1
    }

    # Exclude unmapped, secondary, supplementary; then pick only the same FLAGs 
    local count=$( \
      samtools view -@ "$threads" -F 0x904 "$bam" \
      | awk '$2==99   || $2==1123 \
            || $2==163  || $2==1187 \
            || $2==0    || $2==16   \
            || $2==1024 || $2==1040' \
      | wc -l \
    )

    echo "$count"
}

