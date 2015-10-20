# Map Comparison
Genetic Map Comparison Pipeline

### Preparation
- Collect and format data
- Put fasta markers in the `02_markers` folder
- Copy path to bwa indexed reference genome in `01_scripts/01_bwa_w_fa.sh`

### Mapping the markers
- Launch `01_scripts/01_bwa_w_fa.sh`

```
./01_scripts/01_bwa_w_fa.sh
```

### Extracting meaningful hits
A) Extract unique hits
```
./01_scripts/02_get_wanted_loci_ids.sh <bamfile>
```

- Produces `wanted_loci.ids` in the folder containing the bamfile

B) Obtain the records for the wanted loci from the samfile

```
./01_scripts/03_get_wanted_loci_sam.py <input_sam_file> <wanted_loci.ids> <output (e.g. wanted_loci.sam)>
```

C) Pair by species and contig

```
./01_scripts/04_pair_markers_by_target.py <wanted_loci.sam> <wanted_loci.ids> <output file (e.g. wanted_loci.info)> <max_dist_bw_pairs>
```

 ### Creating figures 
- Prepare data with Python script (not existing yet)
- Create figures (need new R script)

