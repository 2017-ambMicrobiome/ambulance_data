# Import metadata
div <- read.csv("~/Dropbox/Jakes_Docs/masonLab/projects/Ambulance_Study/cleaningStatus_stackedBarGraph/diversity/div_pairs_relAbun_0.1")
# Preview 1:10
div2 <- div[1:2, ]
head(div)

# Calculate number of spp present based on total counds and total reads as compared to surface
y_values <- t(as.matrix((div$num_spp/div$tot_counts)/(div$tot_reads/1E6)))
x_values <- t(as.matrix(div$cleanSt_surface))


#Create Plot
par(las=2) # make label text perp
par(mar=c(20,15,2,2)) # increase y-axis margin
barplot(y_values, 
        horiz=TRUE,
        names.arg = x_values,
        cex.names=1)

title_main <- "Number of Species Present"
xlabel <- "(# spp/total counts used)/(total reads/1E6 per sample)"

title(main = title_main)
title(main = "Number of Species Present", col.main="black",
      xlab = xlabel,
      col.lab = "black",
      cex.lab=1.2,
      line = 5)
