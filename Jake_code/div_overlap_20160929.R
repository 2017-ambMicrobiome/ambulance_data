library("vegan")
library(stringr)
library(plyr)
library(xtable)

# import data, this tsv file has header and firs # removed IOT for importation and is also removed taxonomy data (now stored as OTU)
amb_data <- read.table("/home/jake/ambulance/overlap_diversity/spp_only_amb.biom.tsv", sep = "\t", header = TRUE, numerals = c("allow.loss"), colClasses = "numeric")
x1 <- data.frame(t(amb_data))

x1_dim <- matrix(dim(x1))

# Reformats each value in the matrix to an integer
for (i in 1:x1_dim[1,1])
{ 
  for (j in 1:x1_dim[2,1])
  {
    x1[i,j] <- as.integer(x1[i,j])
  }
}

# Diversity shannon-weaver index
div_amb_simp <- matrix(diversity(x1, index = "simpson"))
div_amb_shan <- matrix(diversity(x1, index = "shannon"))
div_amb_stat <- matrix(c(mean(div_amb_simp),sd(div_amb_simp), mean(div_amb_shan),sd(div_amb_shan)), nrow = 2, ncol = 2)

# Creating table
div_amb_stat <- round(div_amb_stat, 3)
rownames(div_amb_stat) <- c("Mean", "SD")
colnames(div_amb_stat) <- c("Simpsons", "Shannon")
div_amb_stat

sink("/home/jake/amb_stat_overlap.tex", append = TRUE, split = TRUE)
xtable(div_amb_stat)

# Create csv table with shannon and simpsons diversity indices
div_amb_simp_r <- round(div_amb_simp, 3)
div_amb_shan_r <- round(div_amb_shan, 3)
x2 <- as.matrix(row.names(x1))
x3 <- cbind(x2[2:398,], div_amb_simp_r[2:398,], div_amb_shan_r[2:398,])
colnames(x3) <- c("Sample.ID", "Simpsons.index", "Shannon.index")

write.csv(x3, file = "Diversity.indices.overlap.csv", quote = FALSE, row.names = FALSE)
