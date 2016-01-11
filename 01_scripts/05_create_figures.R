# Create map to map comparison figures

# Empty workspace
rm(list=ls())

# Load data.pairs
data.pairs = read.table("03_mapped/wanted_loci.info")
data.loci = read.table("02_raw_data/markers.fasta.info")

# Global variables
minimum_number_of_points = 5 # Minimum number hits between 2 LGs to display in dark
figure_folder = "04_figures"
sp1_col=1
sp2_col=8

# Loop over species pairs
cat("  Creating pairwise map comparison figures...\n")

for (sp1 in levels(data.pairs[,sp1_col])) {
    cat(paste("    ", sp1, "\tVS\t", sep=""))

    for (sp2 in levels(data.pairs[,sp2_col])) {
        if (sp1 != sp2) {
            cat(paste(sp2, " ", sep=""))

            # Subset data.pairs to sp1 and sp2
            # Treat only half the sp1 / sp2 pairwise comparisons when sp1 != sp2
            d = data.pairs[data.pairs[,sp1_col] == sp1 & data.pairs[,sp2_col] == sp2,]

            # Create new figure
            #figure_name = paste(sp2, "_", sp1, ".pdf", sep="")
            #pdf(paste(figure_folder, figure_name, sep="/"), width=12, height=12)
            figure_name = paste(sp2, "_", sp1, ".png", sep="")
            png(paste(figure_folder, figure_name, sep="/"), width=1100, height=1100)

                # Create empty plot of the appropriate dimensions
                plot(d[,(sp1_col+3)], d[,(sp2_col+3)],
                     main=paste("Map comparison for species", sp2, "and", sp1),
                     xlab=paste(sp1, "map"),
                     ylab=paste(sp2, "map"),
                     pch=19, col="#00000088", cex=0.5, type="n")

                max.x = max(d[,(sp1_col+3)])
                max.y = max(d[,(sp2_col+3)])

                # Adding linkage group names for species 1 on x axis
                sp1.lgs = data.loci[data.loci[,1] == sp1, 2]
                sp2.lgs = data.loci[data.loci[,1] == sp2, 2]
                sp1.data = data.loci[data.loci[,1] == sp1, ]
                sp2.data = data.loci[data.loci[,1] == sp2, ]

                #for (lg in sort(unique(d[,(sp1_col+1)]))){
                for (lg in sp1.lgs){
                    minimum = sp1.data[sp1.data[,2] == lg, 3]
                    maximum = sp1.data[sp1.data[,2] == lg, 4]
                    text((maximum + minimum) / 2, -30, lg, cex=0.8)
                }

                # Adding linkage group names for species 2 on y axis
                #for (lg in sort(unique(d[,(sp2_col+1)]))){
                for (lg in sp2.lgs){
                    minimum = sp2.data[sp2.data[,2] == lg, 3]
                    maximum = sp2.data[sp2.data[,2] == lg, 4]
                    text(-30, (maximum + minimum) / 2, lg, cex=0.8)
                }

                # Adding linkage group rectangles for species 1 on x axis
                for (lg1 in sp1.lgs){
                    for (lg2 in sp2.lgs){

                        # New (using all markers to have proper LG boundaries)
                        sp1.lg1.data = data.loci[data.loci[,1] == sp1 & data.loci[,2] == lg1, ]
                        sp2.lg2.data = data.loci[data.loci[,1] == sp2 & data.loci[,2] == lg2, ]

                        min.x = sp1.lg1.data[,3]
                        max.x = sp1.lg1.data[,4]
                        min.y = sp2.lg2.data[,3]
                        max.y = sp2.lg2.data[,4]
                        rect(min.x, min.y, max.x, max.y, col="#00000022", border=F)
                    }
                }

                # Iterate through linkage group pairs
                for (lg1 in sort(unique(d[,(sp1_col+1)]))) {
                    for (lg2 in sort(unique(d[,(sp2_col+1)]))) {

                        # Subset data for that linkage group pair
                        dd = d[d[,(sp1_col+1)] == lg1 & d[,(sp2_col+1)] == lg2, ]

                        # Skip if there are no markers for this LG pair
                        if (nrow(dd) == 0) {
                            next
                        }

                        # Color points black or red as a function of whether
                        # they are in a quadrant with enough data points or not
                        if (nrow(dd) >= minimum_number_of_points) {
                            color = "black"
                        } else {
                            color = "red"
                        }

                        # Adding the points
                        points(dd[,(sp1_col+3)], dd[,(sp2_col+3)],
                               xlab=paste(sp1, "map"),
                               ylab=paste(sp2, "map"),
                               pch=19, col=color, cex=0.1)
                    }
                }
            dev.off()
        }
    }
    cat("\n")
}
