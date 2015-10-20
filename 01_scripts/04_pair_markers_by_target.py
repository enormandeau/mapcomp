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

def find_min_coord(m, max_dist):
    coord = []

    flat = [item for sublist in m for item in sublist if item >= 0]
    minimum = min([x for x in flat])
    #print("Minimum: {}".format(minimum))

    for i in range(len(m)):
        for j in range(len(m[0])):
            if m[i][j] == minimum:
                if minimum >= max_dist:
                    #print("Minimum is too high")
                    cleanup_matrix(m, (i, j))
                    return coord

                row_min = min([abs(x) for x in m[i]])
                col = [x[j] for x in m]
                col_min = min([abs(x) for x in col])
                if m[i][j] > min([row_min, col_min]):
                    #print("Suboptimal pair removed")
                    cleanup_matrix(m, (i, j))
                    return coord

                else:
                    coord.append((i, j))

    return coord

def cleanup_matrix(m, coord):
    x, y = coord
    for i in range(len(m)):
        for j in range(len(m[0])):

            # Remove other possibilities
            if i == x or j == y:
                if m[i][j] > 0:
                    m[i][j] *= -1

    return m

def is_finished(m):
    flat = [item for sublist in m for item in sublist if item >= 0]
    return set(flat) == set([])

def find_best_pairs(m, max_dist):
    finished = False
    good_pairs = []

    while not finished:
        coord = find_min_coord(m, max_dist)

        if len(coord) > 0:
            m = cleanup_matrix(m, coord[0])
            good_pairs.append(coord[0])

        #for i in range(len(m)):
        #    print(m[i])

        #print("Good pairs: {}".format(good_pairs))
        #print("")

        finished = is_finished(m)

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
        count = 1
        max_count = 1000000000
        for target in targets:
            loci = targets[target]

            # Use only targets where there is more than one species
            # find species (continue if only one species present)
            species = set([x.species for x in loci])

            if len(species) > 1:
                for sp1 in species:
                    sp1_loci = [x for x in loci if x.species == sp1]
                    for sp2 in [x for x in species if x != sp1]:
                        if sp1 != sp2:
                            sp2_loci = [x for x in loci if x.species == sp2]

                            if count > max_count:
                                sys.exit(0)
                            else:
                                count += 1

                            #print len(sp1_loci), len(sp2_loci)
                            #print sp1_loci[0].hit_position, sp1_loci[0]

                            # Create distance matrix
                            m = []
                            for l1 in range(len(sp1_loci)):
                                m.append([])
                                for l2 in range(len(sp2_loci)):
                                    abs_dist = 0.1 + abs(int(sp1_loci[l1].hit_position) -
                                              int(sp2_loci[l2].hit_position))

                                    m[l1].append(abs_dist)

                            # Find best pairs
                            #print len(sp1_loci), len(sp2_loci)
                            best_pairs = find_best_pairs(m, max_dist)
                            #print "Num best pairs:", len(best_pairs)

                            for pair in best_pairs:
                                i1, i2 = pair
                                #print i1, len(sp1_loci) - 1
                                #print i2, len(sp2_loci) - 1

                                l1 = sp1_loci[i1]
                                l2 = sp2_loci[i2]

                                outfile.write("\t".join([str(l1), str(l2)]) + "\n")

                            #for l in sp1_loci:
                            #    #for ll in sp2_loci:
                            #    #    outfile.write(str(l) + "\t" + str(ll) + "\n")

                            #    ## Write to file
                            #    closest_locus = find_closest_locus(l, sp2_loci)
                            #    outfile.write("\t".join([str(l), str(closest_locus)]) + "\n")
                            #    
                            #for l in sp2_loci:
                            #    #for ll in sp1_loci:
                            #    #    outfile.write(str(l) + "\t" + str(ll) + "\n")

                            #    ## Write to file
                            #    closest_locus = find_closest_locus(l, sp1_loci)
                            #    #outfile.write("\t".join([str(l), str(closest_locus)]) + "\n")

