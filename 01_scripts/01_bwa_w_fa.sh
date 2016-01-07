#!/bin/bash
# use bwa mem to map map loci in fasta format against reference genome

# global variables (note: point to REFERENCEGENOME)
INPUT_FASTA="markers.fasta"
REFERENCEGENOME="02_raw_data/genome/genome.fasta"
NUMPROCESSORS=16
INPUT_FOLDER="02_raw_data"
MAPPED_FOLDER="03_mapped"

# Notes:
# Assumes reference is already indexed
# Otherwise, index it with:
#   bwa index /path/to/genome
# Assumes adapters have been removed from reads

echo "Running MapComp"

# Extract informations for figures later
echo -e "Species\tLG\tTotalPosition" > $INPUT_FOLDER/$INPUT_FASTA.temp
grep ">" $INPUT_FOLDER/$INPUT_FASTA |
    perl -pe 's/>//; s/_/\t/g' |
    awk 'BEGIN {OFS="\t"} {print $1,$2,$4}' >> $INPUT_FOLDER/$INPUT_FASTA.temp

./01_scripts/utility_scripts/find_min_max_total_positions.py $INPUT_FOLDER/$INPUT_FASTA.temp $INPUT_FOLDER/$INPUT_FASTA.info

rm $INPUT_FOLDER/$INPUT_FASTA.temp

# Map reads using bwa mem 
input_file=$INPUT_FOLDER/$INPUT_FASTA

# Align reads
echo "  Aligning reads to genome..."
bwa mem -t $NUMPROCESSORS $REFERENCEGENOME $input_file > $input_file.sam 2> /dev/null

# keep only mapped reads (-F 4)
echo "  Filering non-mapped reads..."
samtools view -Sb -F 4 $input_file.sam > $input_file.mapped_only.bam

# Sort reads
echo "  Sorting bam file..."
samtools sort $input_file.mapped_only.bam $input_file.sorted

# Index bam file
echo "  Indexing bam file..."
#samtools index $input_file.sorted.bam

# clean temporary files
rm $input_file.mapped_only.bam
mv ./$INPUT_FOLDER/*.sam ./$INPUT_FOLDER/*.sorted.bam ./$MAPPED_FOLDER/
