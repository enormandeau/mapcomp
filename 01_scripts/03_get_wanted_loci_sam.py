#!/usr/bin/env python
"""Subsest sam file to keep wanted loci

Usage:
    ./01_scripts/03_get_wanted_loci_sam.py input_sam_file wanted_loci output_sam_file
"""

# Modules
import sys

# Classes

# Functions

# Main
if __name__ == '__main__':
    try:
        input_file = sys.argv[1]
        wanted_file = sys.argv[2]
        output_file = sys.argv[3]
    except:
        print __doc__
        sys.exit(1)

    infile = open(input_file)
    wfile = open(wanted_file)
    outfile = open(output_file, "w")

    wanted_ids = set()
    for line in wfile:
        wanted_ids.add(line.strip())

    for line in infile:
        if line.startswith("@"):
            continue

        l = line.strip().split()
        if l[0] in wanted_ids:
            outfile.write(line)

    infile.close()
    wfile.close()
    outfile.close()

