#!/bin/bash

# This script is to automate the counting of the markers to be mapped, and that do map during MapComp, lauch from the main directory

# Global variables
SPECIESDESIGNATION=species_designations.txt

# Find out designations for species in the comparison
awk '{ print $1 }' 03_mapped/wanted_loci.sam | awk -F'_' '{print $1}' | uniq > $SPECIESDESIGNATION

# Count how many markers are in each of those designations
# and how many successful mappings against the reference genome there are
cat ./$SPECIESDESIGNATION |
    sort |
    while read i
    do
        echo "Number of markers available: " $i
        # TODO Remove explicit file names and replace by variables
        grep -E "$i"_ 02_raw_data/all_species-30-11-15.fasta | wc -l
        
        echo "Number of mappings to genome: " $i
        grep -E ^"$i"_ 03_mapped/wanted_loci.sam | awk '{ print $1 }' - | uniq | wc -l
        echo "Number of pairs between $i and Sfon"
        grep -E "^$i\t" 03_mapped/wanted_loci.info | grep -E 'Sfon\t' | wc -l
    done


