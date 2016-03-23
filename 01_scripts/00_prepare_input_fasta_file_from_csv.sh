#!/bin/bash

## USAGE : 00_prepare_input_fasta_file_from_csv.sh CSV_INFILE CSV_OUTFILE FASTA_OUTFILE

set -o errexit
set -o nounset

# Global variables
readonly USER_INPUT="${1}"
readonly OUTPUT_CSV_FILE="${2-02_data/markers_with_total_potision.csv}"
readonly OUTPUT_FASTA="${3-02_data/markers.fasta}"

readonly TMP_CSV_FILE=$( mktemp markers_XXXX ) #"02_data/.temp_input_markers.csv"

readonly PROGDIR=$( readlink -e $( dirname $0 ) )

# Function to print script usage
usage () {
cat << EOF
Usage:
    00_prepare_input_fasta_file_from_csv.sh CSV_INFILE CSV_OUTFILE FASTA_OUTFILE

Where:
    CSV_FILE is the input data for MapComp in the exact format described in the
    README.md file
EOF
exit
}

test -z ${USER_INPUT} && usage

# Copy input csv file to expected file
# with no hard coded values this may be useless.
#   R script should be able to work on USER_INPUT
cp ${USER_INPUT} ${TMP_CSV_FILE}

# Source R script to find total positions
echo "Formating data for MapComp..."
Rscript --vanilla ${PROGDIR}/utility_scripts/total_linkage_group_position.R \
    "${TMP_CSV_FILE}" "${OUTPUT_CSV_FILE}"

# Format intermediate output into fasta file
awk -F, 'BEGIN{OFS="";} {print ">"$1"_"$2"_"$3"_"$4"_"$5"\n"$6}' \
    ${OUTPUT_CSV_FILE} > ${OUTPUT_FASTA}

echo "You can now run './mapcomp'"

# Cleanup
rm ${TMP_CSV_FILE}
