#!/usr/bin/env python
"""Find pairs of inter-species markers that map to the same target

Usage:
    python extract_wanted_loci.py sam_file wanted_loci output_file max_dist
"""

# Modules
from collections import defaultdict
import numpy as np
import warnings
import random
import math
import sys

# Classes
class Locus(object):
    """One locus of one individual and its hit information
    """

    def __init__(self, info):

        info = info.split("\t")
        species_info = info[0].split("_")

        self.species = species_info[0]
        self.linkage_group = species_info[1]
        self.position = species_info[2]
        self.total_position = species_info[3]
        self.name = species_info[4]

        self.target = info[2]
        self.hit_position = info[3]

    def __repr__(self):
        return "\t".join([self.species,
                          self.linkage_group,
                          self.position,
                          self.total_position,
                          self.name,
                          self.target,
                          self.hit_position])

# Functions
def get_best_pairs(distance_matrix, max_dist):
    """Find best pairs of markers

    Each marker is used a maximum of one time
    """

    # Compute row and column minimums
    row_minimums = np.nanmin(distance_matrix, axis=1)
    col_minimums = np.nanmin(distance_matrix, axis=0)

    # Find best pairs
    best_pairs = []
    nrow, ncol = distance_matrix.shape
    already_found = set()

    for row in xrange(nrow):
        for col in xrange(ncol):
            value = distance_matrix[row, col]
            if value <= max_dist and value == row_minimums[row] and value == col_minimums[col] and value not in already_found:
                already_found.add(value)
                best_pairs.append((row, col))

    return best_pairs

# Main
if __name__ == '__main__':
    # Suppress numpy nan means warnings
    warnings.filterwarnings("ignore")

    # Parsing user input
    try:
        sam_file = sys.argv[1]
        wanted_loci = sys.argv[2]
        output_file = sys.argv[3]
        max_dist = float(sys.argv[4])
    except:
        print __doc__
        sys.exit(1)

    # Building set of wanted loci
    wanted = set()

    with open(wanted_loci) as wfile:
        for line in wfile:
            wanted.add(line.strip())

    # Build dictionary of targets (keys) and a all the loci who mapped to them (values)
    targets = defaultdict(list)

    with open(sam_file) as sfile:
        for line in sfile:
            locus = Locus(line)
            targets[locus.target].append(locus)

    # Go through targets and find pairs of cross species markers 
    print("  Find best pairs of markers...")

    with open(output_file, "w") as outfile:
        for target in targets:
            loci = targets[target]

            # Use only targets where there is more than one species
            species = set([x.species for x in loci])

            if len(species) > 1:
                for sp1 in species:
                    sp1_loci = [x for x in loci if x.species == sp1]

                    for sp2 in [x for x in species if x != sp1]:
                        sp2_loci = [x for x in loci if x.species == sp2]

                        # Create distance matrix
                        nrow = len(sp1_loci)
                        ncol = len(sp2_loci)
                        distance_matrix = np.zeros((nrow, ncol))

                        for l1 in range(len(sp1_loci)):
                            for l2 in range(len(sp2_loci)):
                                abs_dist = abs(int(sp1_loci[l1].hit_position) -
                                          int(sp2_loci[l2].hit_position))

                                distance_matrix[l1, l2] = abs_dist

                        # Find best pairs
                        best_pairs = get_best_pairs(distance_matrix, max_dist)
                        
                        for pair in best_pairs:
                            i1, i2 = pair
                            l1 = sp1_loci[i1]
                            l2 = sp2_loci[i2]
                            outfile.write("\t".join([str(l1), str(l2)]) + "\n")
