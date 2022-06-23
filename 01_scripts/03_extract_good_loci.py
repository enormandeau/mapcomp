#!/usr/bin/env python3
"""Subsest sam file to keep wanted loci

Usage:
    ./01_scripts/03_get_wanted_loci_sam.py input_sam_file wanted_loci output_sam_file
"""

# Modules
import sys

# Main
if __name__ == '__main__':
    # Parse user input
    try:
        input_file = sys.argv[1]
        wanted_file = sys.argv[2]
        output_file = sys.argv[3]
    except:
        print(__doc__)
        sys.exit(1)

    # Open file handles
    infile = open(input_file)
    wfile = open(wanted_file)
    outfile = open(output_file, "w")

    # Get set of wanted ids
    wanted_ids = set()
    for line in wfile:
        wanted_ids.add(line.strip())

    # Iterate over input_sam_file
    print("  Extracting good reads...")
    for line in infile:
        if line.startswith("@"):
            continue

        l = line.strip().split()

        # Write line if mapping quality is sufficient
        if l[0] in wanted_ids and int(l[4]) >= 10:
            outfile.write(line)

    # Close file handles
    infile.close()
    wfile.close()
    outfile.close()
