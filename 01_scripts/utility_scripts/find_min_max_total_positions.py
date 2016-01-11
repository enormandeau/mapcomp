#!/usr/bin/env python
"""Format information about min and max total position for each linkage group

Usage:
    ./01_scripts/utility_scripts/find_min_max_total_positions.py input_file output_file
"""

# Modules
from collections import defaultdict
import sys

# Parsing user input
try:
    input_file = sys.argv[1]
    output_file = sys.argv[2]
except:
    print __doc__
    sys.exit(1)

# Create dictionary
minimums = defaultdict(lambda: defaultdict(dict))
maximums = defaultdict(lambda: defaultdict(dict))

# Read and treat input file
with open(input_file) as ifile:
    for line in ifile:
        if line.startswith("Species"):
            continue

        species, lg, total_position = line.strip().split()
        lg = float(lg)
        total_position = float(total_position)

        if lg not in minimums[species] or minimums[species][lg] > total_position:
            minimums[species][lg] = total_position

        if lg not in maximums[species] or maximums[species][lg] < total_position:
            maximums[species][lg] = total_position

# Write output
with open(output_file, "w") as ofile:
    for sp in sorted(minimums):
        for lg in minimums[sp]:
            ofile.write("\t".join([sp, str(lg), str(minimums[sp][lg]), str(maximums[sp][lg])]) + "\n")
