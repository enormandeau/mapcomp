# MapComp

A Genetic Map Comparison Pipeline

# Introduction

MapComp facilitates the visual comparison among linkage maps of similar taxons
in order to assess their quality and makes exploring the genetic evolution of
these taxons simpler. The novelty of the approach lies in the use of a
reference genome to maximize the number of marker pairs that can be compared
among maps, even when completely different library preparation protocols have
been used to create the map markers. As such, it requires the existence of an
assembled genome for a species that is phylogenetically close to the species
with the maps that are being compared.

# Using MapComp

The main steps in using MapComp are:

- Get reference genome and put here: `02_data/genome/genome.fasta`
- Index reference genome
- Get marker data from different taxa
- Prepare .csv file (see `02_data/tutorial_markers.csv` for exact format)
- Prepare fasta file of marker automatically from .csv file
- Run mapcomp, which will:
  - Map marker sequences on genome scaffolds
  - Filter out non-unique and bad quality alignments
  - Keep only the best marker pairs
  - Create figures

# Dependencies

In order to use MapComp, you will need the following:

- Linux or MacOS
- Python 2.7
- bwa
- samtools
- The R statistical language

If you are using a Debian derived Linux distribution for example Ubuntu or
Linux Mint, you can install all the required tools with the following command:

```
sudo apt-get install bwa samtools r-base-core
```

# Tutorial

A tutorial data set of markers for two species and a reference genome are
included in MapComp. Both the genome and marker data used for the tutorial were
created *in silico*. As a result, the figures will look too perfect. However,
the goal of the tutorial is for you to run a full MapComp analysis once so that
you then know how to use it on your real data. Additionally, the tutorial data
serves as an example of the exact format required for the .csv file containing
the information about the markers of the different linkage map.

Once you have produced the figures from the tutorial data, then using MapComp
on your data will be as easy as preparing the .csv file, automatically creating
the markers fasta file, getting and indexing the reference genome and running
`./mapcomp`.

## Tutorial run

```
# Rename and index genome
cp 02_data/genome/tutorial_genome.fasta 02_data/genome/genome.fasta
bwa index 02_data/genome/genome.fasta

# Prepare fasta file
./01_scripts/00_prepare_input_fasta_file_from_csv.sh 02_data/tutorial_markers.csv

# Run mapcomp
./mapcomp
```

You can now look at the figures in the `04_figures` folder.

# Data preparation

In order to compare linkage maps, you will need to collect the following
information about each marker:

- Species name (eg: hsapiens)
- Linkage Group number
- Position (eg: in centi Morgans, or cM)
- Marker Identifier (eg: marker_0001)
- Marker Nucleotide Sequence (60 base pairs of more)

Once you have all this information about the markers, you will need to create a
.csv file containing these informations, plus one extra column containing
zeroes, in the following format:

```
SpeciesName,LG,Position,Zeroes,markerName,markerSequence
```

Here is what the .csv file may look like:

```
hsapiens,1,0.58,0,marker_0001,CGGCACCTCCACTGCGGCACGAAGAGTTAGGCCCCGTGCTTTGCGG
hsapiens,1,5.74,0,marker_0002,CGGCACCTCCACTGCGGCACGAAGAGTTAGGCCCCGTGCTTTGCGG
...
hsapiens,1,122.39,0,marker_0227,CGGCACCTCCACTGCGGCACGAAGAGTTAGGCCCCGTGCTTTGCGG
```

Use the `02_data/tutorial_markers.csv` file as a template for your own .csv
file.

Note that:

- There is no header line in the .csv file
- There are 6 columns of information
- The different columns are separated by a comma (`,`)
- The fourth column is filled with zeroes (`0`)
- You need more than one map in the .csv file

# Automatically creating the markers fasta file

The .csv file will be used to create a fasta file using the following script:

```
./01_scripts/00_prepare_input_fasta_file_from_csv.sh <your_file.csv>
```

This will produce a file named `02_data/marker.fasta`.

# Preparing the reference genome

Once you have a reference genome in fasta format, copy it here:
`02_data/genome/genome.fasta` and index it with bwa:

```
bwa index 02_data/genome/genome.fasta
```

# Running MapComp

Once your data has been prepared and your reference genome is indexed, running
mapcomp is as easy launching the following command:

```
./mapcomp
```

# Citing
If you use MapComp in your research, please cite:

Sutherland et al. 2016...

# License

MapComp is licensed under the GNU General Public Licence version 3 (GPL3). See
the LICENCE file for more details.
