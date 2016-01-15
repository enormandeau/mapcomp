#!/bin/bash

# Global variables
readonly USER_INPUT=${1}
readonly INPUT_CSV_FILE="input_markers.csv"
readonly OUTPUT_CSV_FILE="input_markers_with_total_potision.csv"
readonly INPUT_FASTA="02_raw_data/markers.fasta"

# Function to print script usage
usage () {
cat << EOF
Usage:
    ./01_scripts/00_prepare_input_fasta_file_from_csv.sh CSV_FILE

Where:
    CSV_FILE is the input data for MapComp in the exact format described in the
    README.md file
EOF
exit
}

test -z ${USER_INPUT} && usage

# Copy input csv file to expected file
cp ${USER_INPUT} ${INPUT_CSV_FILE}

# Source R script to find total positions
echo "Formating data for MapComp..."
R --slave -q -e 'source("01_scripts/utility_scripts/total_linkage_group_position.R")'

# Format intermediate output into fasta file
awk -F, 'BEGIN{OFS="";} {print ">"$1"_"$2"_"$3"_"$4"_"$5"\n"$6}' \
    ${OUTPUT_CSV_FILE} > ${INPUT_FASTA}

# Cleanup
rm ${INPUT_CSV_FILE}
rm ${OUTPUT_CSV_FILE}

echo "You can now run './mapcomp'"
