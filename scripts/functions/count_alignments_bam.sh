#!/bin/bash

# Function: count_alignments_bam
# Usage: count_alignments_bam <threads> <bam_file>
# Returns: number of **mapped reads**, assuming paired-end BAMs

count_alignments_bam() {
    local threads="$1"
    local bam="$2"

    if [[ ! -f "$bam" ]]; then
        echo "[ERROR] BAM file not found: $bam" >&2
        echo "0"
        return 1
    fi

    # Count number of **mapped reads**, assuming paired-end
    # Exclude unmapped reads (flag 0x4), count only properly aligned reads
    local count
    count=$(samtools view -@ "${threads}" -f 0x2 -F 0x4 -c "$bam")

    echo "$count"
}