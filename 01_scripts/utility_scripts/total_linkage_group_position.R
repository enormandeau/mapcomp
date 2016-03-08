## AIMS  : Add marker total positions to a CSV input file prepared as described in README.md
## USAGE : Rscript --vanilla total_linkage_group_position.R input.csv output.csv

# Example format of CSV input file below
#SpeciesName,LG,Pos,TotPos,MarkerName,Sequence

#species_name,1,0,0,63776,GTATGAGGTTTGTCTTTAACAAAGGTCTCCAGTCAGAAACAGAGATGATGTGTCTTTAACCCTCCAGT
#species_name,1,6.135,0,64642,CATCAAGTTATAAAAGTAAATCAAGTTGACATGTTAATGTACACCTCAAAACAGCTCTTTTGATTCAG
#species_name,1,7.349,0,2835,TGACTGATACTATAGGAAACTAAGGTAGTACCGATTTAGGTCAGTTGATTTGTGTCCACCATTTCCCC
#species_name,1,8.381,0,97480,TCGTATCCCTCTAGGTCCATGTTCAAACACAGCTCCATCCTCCCTACGAGTTGAGCTCAGCCAGCCTC
#species_name,1,13.62,0,61360,TCAGCCGTGTGCATTACTCGTAAAATCCCATTTCTTCGTGAACAGGCCCCACGGTTTCCAGCCCGGGA
#species_name,1,14.347,0,100001,CTGTGTAAGGGCTATTTAACCAAGAAGGAGAGTGATGGATTGCTGCATCAGATGACCTGGCCTCCACA
#species_name,1,15.925,0,19208,TTTTCTTCAGAGGGGATCTGCTGAGAGAGAAGCCCCTCCTACCAGGGGAGGAGATGACCATGCCACGC
#species_name,1,16.856,0,44753,TGAAAAAGACAGATGTGAAGGTCCTGGACTGGCATGGTTACACATGGTCTGCGGTTGTGAGGCCGGTT
#species_name,1,18.172,0,29838,ACACGAGAGGGACGGTGTTGGTGACATGATGTTAGCTGACAGGCAGGAAAACTGATGACTTTTTCATG
#species_name,1,18.582,0,37275,CCTTAGAGCTAGGCTACAGTACCTCATACAATAATTTATTTGCTTTGTATGATTGAGTCAGAGGTTGT

# Clear workspace
rm(list=ls())

args = commandArgs( TRUE )

# Global variables / default values
input.csv  = ifelse( exists( args[1] ), args[1], "02_data/.temp_input_markers.csv" )
output.csv = ifelse( exists( args[2] ), args[2], "02_data/markers_with_total_potision.csv" )

# Load data
data = read.csv(input.csv, header = F, col.names = c("sp","lg","pos","totpos","mname","seq"))

# Initialize variables
data.new = NULL
data$totpos = 0

# TODO make variable names more explicit
# Iterate over the species
for (species in levels(data$sp)) {
    cat("  Treating:", species, "\n")
    spec.sp.data = data[data$sp == species, ]
    lg.unique = unique(spec.sp.data$lg)

    if (length(lg.unique) > 1) {
        # Iterate over the LGs
        for (i in lg.unique) {
            maximum = max(spec.sp.data$pos[as.integer(spec.sp.data$lg) == i])

            for (j in lg.unique) {
                if (j > i) {
                    spec.sp.data$totpos[as.integer(spec.sp.data$lg) == j] = 
                        spec.sp.data$totpos[as.integer(spec.sp.data$lg) == j] + maximum
                }
            }
        }
    }

    data.new = rbind(data.new, spec.sp.data)    
}

data.new$totpos = data.new$totpos + data.new$pos

write.table(data.new, output.csv, row.names =F, quote = F, sep = ",", col.names = F)
