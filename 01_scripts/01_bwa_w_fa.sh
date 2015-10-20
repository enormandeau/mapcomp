#!/bin/bash
# use bwa mem to map map loci in fasta format against reference genome

# global variables (note: point to REFERENCE)
INPUT_FOLDER="02_raw_data"
MAPPED_FOLDER="03_mapped"
REFERENCE="/project/lbernatchez/drobo/users/bensuth/00_resources/Ssa_ASM_3.6.fasta"

# Notes:
# Assumes reference is already indexed
# Assumes adapters have been removed from reads
# We keep only mapped reads using the -F 4 flag in samtools view

# Map reads using bwa mem 
ls -1 $INPUT_FOLDER/*.fasta |
    sort -u |
    while read i
    do
        echo $i
        bwa mem -t 2 $REFERENCE $i > $i.sam
        samtools view -Sb -F 4 $i.sam > $i.mapped_only.bam
        samtools sort $i.mapped_only.bam $i.sorted
        samtools index $i.sorted.bam
    done

# clean up space
rm ./$INPUT_FOLDER/*.mapped_only.bam
mv ./$INPUT_FOLDER/*.sam ./$INPUT_FOLDER/*.sorted.bam ./$INPUT_FOLDER/*.sorted.bam.bai ./$MAPPED_FOLDER/
