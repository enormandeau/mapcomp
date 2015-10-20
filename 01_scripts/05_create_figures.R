# Create map to map comparison figures

# Empty workspace
rm(list=ls())

# Load data
data = read.table("03_mapped/wanted_loci.info")

# Global variables
sp1_col=1
sp2_col=8
minimum_number_of_points = 5
figure_folder = "04_figures"

# Loop over species pairs
for (sp1 in levels(data[,sp1_col])) {
    print(sp1)
    #started = F
    for (sp2 in levels(data[,sp2_col])) {
        print(paste("  ", sp2, sep=""))

        # Number of points in quadrants with more than minimum_number_of_points
        good = 0
        bad = 0

        # Treat only half the sp1 / sp2 pairwise comparisons when sp1 != sp2
        #if (started) {
        if (sp1 != sp2) {

            # Subset data to sp1 and sp2
            d = data[data[,sp1_col] == sp1 & data[,sp2_col] == sp2,]

            # Move to next dataset if there is no data
            if (nrow(d) == 0) {next}

            # Maximum dimensions for rectangles
            bottom = -1000
            top = 10000
            left = -1000
            right = 10000

            # Create new pdf figure
            figure_name = paste(sp1, "_", sp2, ".pdf", sep="")
            #png(paste(figure_folder, figure_name, sep="/"), width=1000, height=1000)
            pdf(paste(figure_folder, figure_name, sep="/"), width=12, height=12)

                # Create empty plot of the appropriate dimensions
                plot(d[,(sp1_col+3)], d[,(sp2_col+3)], main=paste(sp1, "vs.", sp2),
                     xlab=paste(sp1, "map"),
                     ylab=paste(sp2, "map"),
                     pch=19, col="#00000088", cex=0.5, type="n")

                max.x = max(d[,(sp1_col+3)])
                max.y = max(d[,(sp2_col+3)])

                # Adding linkage group rectangles for species 1 on x axis
                for (lg in sort(unique(d[,(sp1_col+1)]))){
                    minimum = min(d[d[,(sp1_col+1)] == lg, (sp1_col+3)])
                    maximum = max(d[d[,(sp1_col+1)] == lg, (sp1_col+3)])
                    #rect(minimum, bottom, maximum, top, col="#00000022", border=F)

                    # Add LG number
                    text((maximum + minimum) / 2, 0, lg, cex=0.8)
                }

                # Adding linkage group rectangles for species 2 on y axis
                for (lg in sort(unique(d[,(sp2_col+1)]))){
                    minimum = min(d[d[,(sp2_col+1)] == lg, (sp2_col+3)])
                    maximum = max(d[d[,(sp2_col+1)] == lg, (sp2_col+3)])
                    #rect(left, minimum, right, maximum, col="#00000022", border=F)

                    # Add LG number
                    text(0, (maximum + minimum) / 2, lg, cex=0.8)
                }

                # Adding linkage group rectangles for species 1 on x axis
                for (lg1 in sort(unique(d[,(sp1_col+1)]))){
                    for (lg2 in sort(unique(d[,(sp2_col+1)]))){
                        min.x = min(d[d[,(sp1_col+1)] == lg1, (sp1_col+3)])
                        max.x = max(d[d[,(sp1_col+1)] == lg1, (sp1_col+3)])
                        min.y = min(d[d[,(sp2_col+1)] == lg2, (sp2_col+3)])
                        max.y = max(d[d[,(sp2_col+1)] == lg2, (sp2_col+3)])
                        rect(min.x, min.y, max.x, max.y, col="#00000022", border=F)
                    }
                }

                # Iterate through linkage group pairs
                for (lg1 in sort(unique(d[,(sp1_col+1)]))) {
                    for (lg2 in sort(unique(d[,(sp2_col+1)]))) {

                        # Subset data for that linkage group pair
                        dd = d[d[,(sp1_col+1)] == lg1 & d[,(sp2_col+1)] == lg2, ]

                        # Color points black or red as a function of whether
                        # they are in a quadrant with enough data points
                        if (nrow(dd) >= minimum_number_of_points) {
                            color = "#00000066"
                            good = good + nrow(dd)
                        } else {
                            color = "#ff000033"
                            bad = bad + nrow(dd)
                        }

                        # Adding the points
                        points(dd[,(sp1_col+3)], dd[,(sp2_col+3)], main=paste(sp1, "vs.", sp2),
                               xlab=paste(sp1, "map"),
                               ylab=paste(sp2, "map"),
                               pch=19, col=color, cex=0.5)
                    }
                }

                # Add proportion of good data to the plot
                #text(0, 0, paste("Proportion good data:", round(good / (good + bad),2)), adj=c(0,1))

            dev.off()

        }

        # Criteria to start creating the graph:
        # Treat only half the sp1 / sp2 pairwise comparisons when sp1 != sp2
        #if (sp1 == sp2) {
        #    started = T
        #}
    }
}
