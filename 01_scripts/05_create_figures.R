## AIMS  : Create map to map comparison figures
## USAGE : Rscript --vanilla 05_create_figures.R pairs_file loci_file output_lg_correspondance figure_folder

# Empty workspace
rm(list=ls())

args = commandArgs( TRUE )

# Global variables
# default value set for retro-compatibility reason
x = args[1]
pairs_file = ifelse( exists( x ), x, "03_mapped/wanted_loci.info" )
x = args[2]
loci_file  = ifelse( exists( x ), x, "02_data/markers.fasta.info" )
x <- args[3]
output_lg_correspondance = ifelse( exists( x ), x, "05_results/linkage_group_correspondance.csv" )
x <- args[4]
figure_folder = ifelse( exists( x ), x, "04_figures" )

minimum_number_of_points = 5 # Minimum number hits between 2 LGs to display in dark
sp1_col=1
sp2_col=8

# Load data.pairs
# Note : var name or file name are not well chosen, imho
data.pairs = read.table( pairs_file )
data.loci = read.table( loci_file )

# Loop over species pairs
cat("  Creating map comparison figures...\n")

lg.correspondance = NULL

for (sp1 in levels(data.pairs[,sp1_col])) {
    cat(paste("    ", sp1, " VS ", sep=""))

    for (sp2 in levels(data.pairs[,sp2_col])) {
        if (sp1 != sp2) {
            cat(paste(sp2, " ", sep=""))

            # Subset data.pairs to sp1 and sp2
            data.sp = data.pairs[data.pairs[,sp1_col] == sp1 & data.pairs[,sp2_col] == sp2,]

            # PDF output
            #figure_name = paste(sp2, "_", sp1, ".pdf", sep="")
            #pdf(paste(figure_folder, figure_name, sep="/"), width=12, height=12)

            # png output
            figure_name = paste(sp2, "_", sp1, ".png", sep="")
            png(paste(figure_folder, figure_name, sep="/"), width=1000, height=1000)

            # Create empty figure of the appropriate dimensions
            plot(data.sp[,(sp1_col+3)], data.sp[,(sp2_col+3)],
                 main=paste("Map comparison for species", sp2, "and", sp1),
                 xlab=paste(sp1, "map (distances in CM)"),
                 ylab=paste(sp2, "map (distances in CM)"),
                 pch=19,
                 col="#00000088",
                 cex=0.5,
                 type="n",
                 bty="n",
                 xaxt="n",
                 yaxt="n")

            # Prepare data minimum and maximum positions
            sp1.lgs = data.loci[data.loci[,1] == sp1, 2]
            sp2.lgs = data.loci[data.loci[,1] == sp2, 2]
            sp1.data = data.loci[data.loci[,1] == sp1, ]
            sp2.data = data.loci[data.loci[,1] == sp2, ]

            # Creating the axes
            axis.increment = 500 # in CM

            # X axis
            maximum.axis.position.x = max(sp1.data[, 4])
            axis.positions.x = seq(0, maximum.axis.position.x, by=axis.increment)
            axis.text.x = as.integer(axis.positions.x)
            axis(1, at=axis.positions.x, labels=axis.text.x, las=1, cex.axis=0.8)

            # Y axis
            maximum.axis.position.y = max(sp2.data[, 4])
            axis.positions.y = seq(0, maximum.axis.position.y, by=axis.increment)
            axis.text.y = as.integer(axis.positions.y)
            axis(2, at=axis.positions.y, labels=axis.text.y, las=1, cex.axis=0.8)

            # Adding linkage group names for species 1 on x axis
            for (lg in sp1.lgs){
                minimum = sp1.data[sp1.data[,2] == lg, 3]
                maximum = sp1.data[sp1.data[,2] == lg, 4]
                text((maximum + minimum) / 2, -30, lg, cex=0.8)
            }

            # Adding linkage group names for species 2 on y axis
            for (lg in sp2.lgs){
                minimum = sp2.data[sp2.data[,2] == lg, 3]
                maximum = sp2.data[sp2.data[,2] == lg, 4]
                text(-30, (maximum + minimum) / 2, lg, cex=0.8)
            }

            # Iterate through linkage group pairs
            for (lg1 in sp1.lgs){
                for (lg2 in sp2.lgs){

                    # New (using all markers to have proper LG boundaries)
                    sp1.lg1.data = data.loci[data.loci[,1] == sp1 & data.loci[,2] == lg1, ]
                    sp2.lg2.data = data.loci[data.loci[,1] == sp2 & data.loci[,2] == lg2, ]

                    min.x = sp1.lg1.data[,3]
                    max.x = sp1.lg1.data[,4]
                    min.y = sp2.lg2.data[,3]
                    max.y = sp2.lg2.data[,4]

                    # Subset data for that linkage group pair
                    dd = data.sp[data.sp[,(sp1_col+1)] == lg1 & data.sp[,(sp2_col+1)] == lg2, ]

                    # Color quadrants lighter or darker grey depending on weather
                    # they are in a quadrant with enough data points or not
                    if (nrow(dd) >= minimum_number_of_points) {
                        rect(min.x, min.y, max.x, max.y, col="#00000022", border=F)
                        lg.correspondance = rbind(lg.correspondance, c(sp1, lg1, sp2, lg2))
                    } else {
                        rect(min.x, min.y, max.x, max.y, col="#00000012", border=F)
                    }

                    # Adding the points
                    points(dd[,(sp1_col+3)], dd[,(sp2_col+3)],
                           xlab=paste(sp1, "map"),
                           ylab=paste(sp2, "map"),
                           pch=19, col="black", cex=0.1)
                }
            }
            dev.off()
        }
    }
    cat("\n")
}

write.table(lg.correspondance, output_lg_correspondance, sep="\t", row.names=F, col.names=F, quote=F)
