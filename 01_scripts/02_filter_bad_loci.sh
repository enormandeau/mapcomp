#!/bin/bash
# Identify the loci that map uniquely to the reference with a MAPQ score >= 10

# Usage:
# ./01_scripts/get_wanted_loci_ids.sh bamfile

# Global variables
readonly BAMFILE="$1"
readonly DIRECTORY="$(dirname $1)"

# Filter bam file
echo "  Filtering non-unique or bad quality alignments..."
samtools view ${BAMFILE} | \
    awk '$5 >= 10 {print $1}' | \
    sort | \
    uniq -c | \
    sort -nr | \
    awk '$1 == 1 {print $2}' > ${DIRECTORY}/wanted_loci.ids
