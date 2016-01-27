# MapComp

A Genetic Map Comparison Pipeline

# Introduction

MapComp facilitates the visual comparison among maps of similar taxons in order
to assess their quality and makes exploring the genetic evolution of these
taxons simpler. The novelty of the approach lies in the use of a reference
genome to maximize the number of marker pairs that can be compared among maps,
even when completely different library preparation protocols have been used to
create the map markers. As such, it requires the existence of an assembled
genome for a species that is phylogenetically close to the species with the
maps that are being compared.

# Using MapComp

The main steps in using MapComp are:

- Get and prepare markers data of different taxons
- Map marker sequences on genome scaffolds
- Filter out non-unique and bad quality alignments
- Keep only the best marker pairs
- Create figures

# Tutorial
```
# Prepare fasta file
./01_scripts/00_prepare_input_fasta_file_from_csv.sh 02_raw_data/tutorial_markers.csv

# Rename and index genome
mv 02_raw_data/genome/tutorial_genome.fasta 02_raw_data/genome/genome.fasta
bwa index 02_raw_data/genome/genome.fasta

# Run mapcomp
./mapcomp

# Look at figures!
```

# TODO

## Data preparation

Note: Importance of respecting the EXACT format

- Collect data
- Download reference genome
- Create marker fasta file
- Genome -> `./02_raw_data/genome/genome.fasta`
- Markers -> `./02_raw_data/markers.fasta`

## Collecting the data
# TODO

## Creating the marker fasta file
# TODO

## Running MapComp in one step
# TODO

## Running MapComp step by step

### Mapping the markers on the genome

- Index the genome if needed with

```
bwa index 02_raw_data/genome/genome.fasta
```

- Launch `01_scripts/01_bwa_w_fa.sh`

```
./01_scripts/01_bwa_w_fa.sh
```

### Extracting meaningful hits

```
./01_scripts/02_get_wanted_loci_ids.sh \
    03_mapped/markers.fasta.sorted.bam

./01_scripts/03_get_wanted_loci_sam.py \
    03_mapped/markers.fasta.sam \
    03_mapped/wanted_loci.ids \
    03_mapped/wanted_loci.sam
```

### Keeping only the best marker pairs

```
./01_scripts/04_pair_markers_by_target.py \
    03_mapped/wanted_loci.sam \
    03_mapped/wanted_loci.ids \
    03_mapped/wanted_loci.info \
    10000000
```

### Creating MapComp figures

```
R -q -e 'source("01_scripts/05_create_figures.R")'
```

# Citing
# TODO
If you use MapComp in your research, please cite:

Sutherland et al. 2016...

# License

MapComp is licensed under the GNU General Public Licence version 3 (GPL3). See
the LICENCE file for more details.
