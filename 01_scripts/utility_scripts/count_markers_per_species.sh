#!/bin/bash

# Count the proportion of markers that mapped on the genome

# Global variables
SPECIES_NAMES=species_names.txt
FASTA_FILE=02_raw_data/markers.fasta
SAM_FILE=03_mapped/wanted_loci.sam
INFO_FILE=03_mapped/wanted_loci.info

# Find species names
awk '{ print $1 }' ${SAM_FILE} | awk -F'_' '{print $1}' | sort | uniq > ${SPECIES_NAMES}

# Count number of available and mapping markers per species and
echo "Proportion of markers mapping on the genome:"
cat ./${SPECIES_NAMES} |
    sort |
    while read i
    do
        mapped=$(grep -E ^"$i"_ ${SAM_FILE} | awk '{ print $1 }' - | uniq | wc -l)
        total=$(grep -E "$i"_ ${FASTA_FILE} | wc -l)
        percent=$(echo "100.0 * ${mapped} / ${total}" | bc)
        echo -e "  ${i} ${mapped} / ${total} (${percent}%)" 
    done

# Cleanup
rm ${SPECIES_NAMES}
