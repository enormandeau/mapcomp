#!/usr/bin/env python
"""Find pairs of inter-species markers that map to the same target

Usage:
    python extract_wanted_loci.py sam_file wanted_loci output_file max_dist
"""

# Modules
from collections import defaultdict
import numpy as np
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
def find_closest_locus(self_locus, other_loci):
    """For one locus in 'self' species, find closest one in 'other' species
    In case of ties, send only the first one
    """

    self_position = float(self_locus.position)
    other_positions = [abs( self_position - float(x.position)) for x in other_loci]
    result = other_loci[other_positions.index(min(other_positions))]
    return result

def find_min_coord(distance_matrix, max_dist):
    """Find coordinates of smallest distance(s) in a distance matrix
    """

    coord = []

    flat = [item for sublist in distance_matrix for item in sublist if item >= 0]
    minimum = min([x for x in flat])

    for i in range(len(distance_matrix)):
        for j in range(len(distance_matrix[0])):
            if distance_matrix[i][j] == minimum:
                if minimum >= max_dist:
                    cleanup_matrix(distance_matrix, (i, j))
                    return coord

                row_min = min([abs(x) for x in distance_matrix[i]])
                col = [x[j] for x in distance_matrix]
                col_min = min([abs(x) for x in col])

                if distance_matrix[i][j] > min([row_min, col_min]):
                    cleanup_matrix(distance_matrix, (i, j))
                    return coord

                else:
                    coord.append((i, j))

    return coord

def cleanup_matrix(distance_matrix, coord):
    """Remove other possible pairs for items that were paired
    """

    x, y = coord
    for i in range(len(distance_matrix)):
        for j in range(len(distance_matrix[0])):
            if i == x or j == y:
                # Remove other possibilities
                if distance_matrix[i][j] > 0:
                    distance_matrix[i][j] *= -1

    return distance_matrix

def is_finished(distance_matrix):
    """Determine if all best locus pairs have been found
    """

    flat = [item for sublist in distance_matrix for item in sublist if item >= 0]
    return set(flat) == set([])

def find_best_pairs(distance_matrix, max_dist):
    """Find best pairs of markers

    Each marker is used a maximum of one time
    """

    finished = False
    good_pairs = []

    while not finished:
        coord = find_min_coord(distance_matrix, max_dist)

        if len(coord) > 0:
            distance_matrix = cleanup_matrix(distance_matrix, coord[0])
            good_pairs.append(coord[0])

        finished = is_finished(distance_matrix)

    return good_pairs

# Main
if __name__ == '__main__':
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
    with open(output_file, "w") as outfile:
        #count = 1
        #max_count = 1000000000
        for target in targets:
            loci = targets[target]

            # Use only targets where there is more than one species
            species = set([x.species for x in loci])

            if len(species) > 1:
                for sp1 in species:
                    sp1_loci = [x for x in loci if x.species == sp1]
                    for sp2 in [x for x in species if x != sp1]:

                        #if sp1 != sp2: # TODO redundant?
                        sp2_loci = [x for x in loci if x.species == sp2]

                        #if count > max_count:
                        #    sys.exit(0)
                        #else:
                        #    count += 1

                        # Create distance matrix
                        distance_matrix = []
                        for l1 in range(len(sp1_loci)):
                            distance_matrix.append([])
                            for l2 in range(len(sp2_loci)):
                                abs_dist = 0.1 + abs(int(sp1_loci[l1].hit_position) -
                                          int(sp2_loci[l2].hit_position))

                                distance_matrix[l1].append(abs_dist)

                        # Find best pairs
                        best_pairs = find_best_pairs(distance_matrix, max_dist)

                        for pair in best_pairs:
                            i1, i2 = pair
                            l1 = sp1_loci[i1]
                            l2 = sp2_loci[i2]

                            outfile.write("\t".join([str(l1), str(l2)]) + "\n")
