# it takes its arguements as a dataframe with two columns,
#   the LG, the position on the LG.
# it generates a third column, the totpos

rm(list=ls())

#setwd("~/Documents/bernatchez/01_Sfon_projects/01_SfQTL/05_mapping_against_genomes/01_raw_materials/merged_cleaned_species_data")
data <- read.csv(file = 
                   "all_species.csv",
                 header = F, col.names = c("sp","lg","pos","totpos","mname","seq"))
head(data$pos)
levels(data$sp)
str(data) # make sure lg is integer
str(data$lg) # needed to make these all parallel so I can make a for loop to generate absolute position for each species
#data$lg <- as.numeric(data$lg)


## 

lg.unique = NULL
data.new = NULL
spec.sp.data = NULL
data$totpos = 0

for (species in levels(data$sp)) {
    print(species)
    spec.sp.data <- data[data$sp == species, ]
    #print(paste("test",species))
    lg.unique = unique(spec.sp.data$lg)
    count = 0
    
    for (i in 1:(length(lg.unique) - 1)) {
      print(i)
      maximum = max(spec.sp.data$pos[as.integer(spec.sp.data$lg) == i])
      print(maximum)
      for (j in (i+1):length(lg.unique)) {
         spec.sp.data$totpos[as.integer(spec.sp.data$lg) == j] = 
            spec.sp.data$totpos[as.integer(spec.sp.data$lg) == j] + maximum
      }
    }
    data.new <- rbind(data.new, spec.sp.data)    
}

data.new$totpos = data.new$totpos + data.new$pos

head(data.new, n = 1000)[,1:4] # just to check

write.table(x = data.new, file = "all_species_w_totpos.csv", row.names =F, quote = F, sep = ",", col.names = F)

# then do this in bash:
# awk -F, 'BEGIN{OFS="";} {print $1"_"$2"_"$3"_"$4"_"$5"\n"$6}' *totpos.csv > all_sp_RAD_maps-08-17-15-mapping_ready.fasta
# fasta record name will be in this order: sp_lg_pos_totpos_mname

# this file is the input for the map_comparison repo, so copy it to the raw data folder

