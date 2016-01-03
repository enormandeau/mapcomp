#!/bin/bash

# This script is to automate the counting of the markers to be mapped, and that do map during MapComp, lauch from the main directory

# first find out what designations you have for species in the comparison
awk '{ print $1 }' 03_mapped/wanted_loci.sam | awk -F'_' '{print $1}' | uniq > species_designations.txt

# then find out how many markers are in each of those designations
# then find out how many successful mappings against the reference genome
cat ./species_designations.txt |
    sort |
    while read i
    do
        echo "Number of markers avail: " $i
        grep -E "$i"_ 02_raw_data/all_species-30-11-15.fasta | wc -l
        
        echo "Number of mappings to genome: " $i
        grep -E ^"$i"_ 03_mapped/wanted_loci.sam | awk '{ print $1 }' - | uniq | wc -l
    done


