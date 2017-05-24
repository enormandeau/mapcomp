# Pair using MapComp, then remove the anonymous paired markers from the input fasta file and redo the pairings again
# NOTE: this is currently destructive to the markers.fasta file in 02_data/ so always have a second version saved elsewere that can be easily re-copied back into this repo. 

# Dependencies:
# Requires biopython
# Requires 01_scripts/utility_scripts/fasta_remove.py
#   The fasta_remove.py script was developed by Eric Normandeau
#   https://github.com/enormandeau/Scripts

# Setup:
# Set the species name of the anon markers
ANON="Salp.anon"
# this will correspond to the first field in the fasta accession name (e.g. >Salp.anon_1_0_0_29)

# Run all from the main directory:
# ./01_scripts/utility_scripts/remove_paired_anon_and_pair_again.sh

# TODO Fix using the .pre to replace fasta at the end as this can result in problems if one stops the process mid-process



# Standard files
WANTED_LOCI="03_mapped/wanted_loci.info"
PAIRINGS_OUT="03_mapped/pairings_out.txt"
PAIRINGS_NAMES="03_mapped/paired_anon_names.txt"

# Prevent carry-over from previous rounds
rm 02_data/num_markers_assigned.txt 2>/dev/null
rm $PAIRINGS_OUT 2>/dev/null

# Backup the markers file as it will be affected during the run
cp 02_data/markers.fasta 02_data/markers.fasta.pre

for i in $(seq 1 10 )
do
    echo "## MapComp to assign position to $ANON markers"
    ./mapcomp
    
    # Retain paired output
    echo "## Number of paired $ANON markers identified this round = "
    grep -E "^$ANON" 03_mapped/wanted_loci.info | wc -l
    grep -E "^$ANON" 03_mapped/wanted_loci.info >> 03_mapped/pairings_out.txt 
    echo "## Number of paired $ANON markers identified in total = "
    wc -l $PAIRINGS_OUT

    # Obtain names of anonymous markers that were paired this round
    awk -v var=$ANON '{ print var "_1_0_0_" $5 }' $PAIRINGS_OUT > $PAIRINGS_NAMES
    PAIRINGS_NUM=`wc -l $PAIRINGS_NAMES`
    echo $PAIRINGS_NUM >> 02_data/num_markers_assigned.txt
    echo "## Number of $ANON markers assigned positions this round"
    echo $PAIRINGS_NUM
    
    # Remove paired markers from fasta
    ./01_scripts/utility_scripts/fasta_remove.py ./02_data/markers.fasta $PAIRINGS_NAMES 02_data/markers_remaining.fasta
     
    # Replace original markers file with the removed markers file
    mv 02_data/markers_remaining.fasta 02_data/markers.fasta
    
    echo "## Number of $ANON markers remaining = "
    grep -c ">$ANON" 02_data/markers.fasta 
done

# Replace the marker removed markers.fasta with the original one saved at the start 
echo "## Replacing the reduced fasta file with the original fasta file"
mv 02_data/markers.fasta.pre 02_data/markers.fasta
wc -l 02_data/markers.fasta

