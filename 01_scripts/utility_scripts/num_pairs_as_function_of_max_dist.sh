#!/bin/bash
for i in 1 2 5 10 20 50 100 200 500 1000 2000 5000 10000 20000 50000 100000 200000 500000 1000000 2000000 5000000 10000000
do
    ./01_scripts/04_pair_markers_by_target.py 03_mapped/wanted_loci.sam 03_mapped/wanted_loci.ids 03_mapped/wanted_loci.info $i; echo $i $(wc -l 03_mapped/wanted_loci.info)
done | tee num_pairs_as_function_of_max_dist.csv
