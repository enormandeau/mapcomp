#!/bin/bash
# use bwa mem to map map loci in fasta format against reference genome

# global variables (note: point to REFERENCEGENOME)
INPUT_FOLDER="02_raw_data"
INPUT_FASTA="markers.fasta"
MAPPED_FOLDER="03_mapped"
REFERENCEGENOME="02_raw_data/genome/genome.fasta"
NUMPROCESSORS=16

# Notes:
# Assumes reference is already indexed
# Otherwise, index it with:
#   bwa index /path/to/genome
# Assumes adapters have been removed from reads

# Extract informations for figures later
echo -e "Species\tLG\tTotalPosition" > $INPUT_FOLDER/$INPUT_FASTA.temp
grep ">" $INPUT_FOLDER/$INPUT_FASTA |
    perl -pe 's/>//; s/_/\t/g' |
    awk 'BEGIN {OFS="\t"} {print $1,$2,$4}' >> $INPUT_FOLDER/$INPUT_FASTA.temp

./01_scripts/utility_scripts/find_min_max_total_positions.py $INPUT_FOLDER/$INPUT_FASTA.temp $INPUT_FOLDER/$INPUT_FASTA.info

rm $INPUT_FOLDER/$INPUT_FASTA.temp

# Map reads using bwa mem 
ls -1 $INPUT_FOLDER/markers.fasta |
    sort -u |
    while read i
    do
        echo $i

        # Align reads
        bwa mem -t $NUMPROCESSORS $REFERENCEGENOME $i > $i.sam

        # keep only mapped reads (-F 4)
        samtools view -Sb -F 4 $i.sam > $i.mapped_only.bam

        # Sort reads
        samtools sort $i.mapped_only.bam $i.sorted

        # Index bam file
        samtools index $i.sorted.bam
    done

# clean temporary files
rm ./$INPUT_FOLDER/*.mapped_only.bam
mv ./$INPUT_FOLDER/*.sam ./$INPUT_FOLDER/*.sorted.bam ./$INPUT_FOLDER/*.sorted.bam.bai ./$MAPPED_FOLDER/
