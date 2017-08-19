# MapComp

Genetic Map Comparison

## Introduction

MapComp facilitates visual comparisons among linkage maps of closely-related
species in order to assess their quality and to simplify the exploration of
their chromosomal differences. The novelty of the approach lies in the use of a
reference genome in order to maximize the number of comparable marker pairs
among pairs of maps, even when completely different library preparation
protocols have been used to generate the markers. As such, MapComp requires a
reference genome, at least a contig-level genome assembly, for a species that
is phylogenetically close to the target species.

## Using MapComp

The main steps in using MapComp are:

- Get a reference genome and put here: `02_data/genome/genome.fasta`
- Index the reference genome (`bwa index 02_data/genome/genome.fasta`)
- Get marker data from two or more taxa
- Prepare .csv marker file (see `02_data/tutorial_markers.csv` for exact format)
- Prepare markers fasta file automatically from .csv file
- Run mapcomp, which will:
  - Map marker sequences on reference genome scaffolds
  - Filter out non-unique and bad quality alignments
  - Keep only the best marker pairs
  - Create figures

## Dependencies

In order to use MapComp, you will need the following:

- Linux or MacOS
- Python 2.7
- numpy (Python library)
- bwa
- samtools (1.x release) 
- The R statistical language

If you are using a Debian derived Linux distribution, for example Ubuntu or
Linux Mint, you can install all the required tools with the following command:

```
sudo apt-get install bwa samtools r-base-core
```

## Tutorial

A tutorial data set of markers for two species and a reference genome are
included in MapComp. Both the genome and marker data used for the tutorial were
created *in silico*. As a result, the figures will look really perfect.
However, the goal of the tutorial to run a full MapComp analysis once to learn
how to use it with your real data. Additionally, the tutorial .csv data file
serves as an example of the exact format required for the marker .csv file,
which contains the marker information for the analyzed species.

Once you have produced the figures from the tutorial data, then using MapComp
on your data will be as easy as preparing the .csv file, automatically creating
the markers fasta file, getting and indexing the reference genome and running
`./mapcomp`.

### Tutorial run

```
# Rename and index genome
cp 02_data/genome/tutorial_genome.fasta 02_data/genome/genome.fasta
bwa index 02_data/genome/genome.fasta

# Prepare fasta file
./01_scripts/00_prepare_input_fasta_file_from_csv.sh 02_data/tutorial_markers.csv

# Run mapcomp
./mapcomp
```

You can now look at the figures in the `04_figures` folder and at the linkage
group correspondance among the species in the `05_results` folder.

## Data preparation

In order to compare linkage maps, you will need to collect the following
information about each marker:

- Species name (eg: hsapiens)
- Linkage Group number (eg: 1, 2, 3...)
- Position in centi Morgans, or cM (eg: 0, 5.32, 22.8)
- Marker Identifier (eg: marker0001)
- Marker Nucleotide Sequence (60 base pairs of more)

Once you have all this information about the markers, you will need to create a
.csv file containing these informations. The .csv file will feature one extra
column containing zeroes and be in the following format:

```
SpeciesName,LG,Position,Zeroes,markerName,markerSequence
```

Here is what the .csv file may look like:

```
hsapiens,1,0.58,0,marker0001,CGGCACCTCCACTGCGGCACGAAGAGTTAGGCCCCGTGCTTTGCGG
hsapiens,1,5.74,0,marker0002,CGGCACCTCCACTGCGGCACGAAGAGTTAGGCCCCGTGCTTTGCGG
...
hsapiens,1,122.39,0,marker0227,CGGCACCTCCACTGCGGCACGAAGAGTTAGGCCCCGTGCTTTGCGG
```

Use the `02_data/tutorial_markers.csv` file as a template for your own .csv
file.

Note that:

- There is no header line in the .csv file
- There are 6 columns of information
- The different columns are separated by a comma (`,`)
- The fourth column is filled with zeroes (`0`)
- You need more than one map in the .csv file
- You should avoid special characters, including underscores (`_`) in the marker names
- You must use the period (`.`) as the decimal separator (no comma (`,`))

## Automatically creating the markers fasta file

The .csv file will be used to create a fasta file using the following script:

```
./01_scripts/00_prepare_input_fasta_file_from_csv.sh <your_file.csv>
```

This will produce a file named `02_data/marker.fasta`.

## Preparing the reference genome

Once you have a reference genome in fasta format, copy it here:
`02_data/genome/genome.fasta` and index it with bwa:

```
bwa index 02_data/genome/genome.fasta
```

## Running MapComp

Once your data has been prepared and your reference genome is indexed, running
mapcomp is as easy launching the following command:

```
./mapcomp
```

## Exploring Results

After MapComp finishes, visual plots comparing the different linkage maps will be 
found in `04_figures` and a summary of the results in `05_results`. For more detailed
results, one can inspect the `03_mapped/wanted_loci.info` file. This file
contains the details of the marker pairs found for each species pair, and can be
useful to obtain exact mapping locations of markers on the reference genome.

Example output image from the tutorial markers and genome:
![Alt text](00_archive/tutorial_figure.png?raw=true  "Output Example")

## Citing

If you use MapComp in your research, please cite:

Sutherland BJG, Gosselin T, Normandeau E, Lamothe M, Isabel N, Bernatchez L.
Salmonid Chromosome Evolution as Revealed by a Novel Method for Comparing
RADseq Linkage Maps. Genome Biol Evol (2016) 8 (12): 3600-3617.
DOI: https://doi.org/10.1093/gbe/evw262

(preprint version: bioRxiv. 2016: 1–44. doi:10.1101/039164)

## Troubleshooting

A Google Group for MapComp is available at: https://groups.google.com/forum/#!forum/mapcomp

## License

MapComp is licensed under the GNU General Public Licence version 3 (GPL3). See
the LICENCE file for more details.
