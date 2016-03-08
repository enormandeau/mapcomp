#!/bin/bash

## AIMS  : Identify the loci that map uniquely to the reference with a MAPQ score >= 10
## USAGE : 02_filter_bad_loci.sh bamfile outfile

# Global variables
readonly BAMFILE="$1"
readonly DIRECTORY="$(dirname $1)"

# default value set for retro-compatibility reason
readonly OUTFILE="${2-${DIRECTORY}/wanted_loci.ids}"

# Filter bam file
echo "  Filtering bad alignments..."
samtools view ${BAMFILE} | \
    awk '$5 >= 10 {print $1}' | \
    sort | \
    uniq -c | \
    sort -nr | \
    awk '$1 == 1 {print $2}' > ${DIRECTORY}/wanted_loci.ids
