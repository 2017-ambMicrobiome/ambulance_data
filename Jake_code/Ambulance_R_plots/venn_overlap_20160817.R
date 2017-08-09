source("/home/cem2009/R/cem_plot2.R") 
library("plotrix")
setwd("/zenodotus/masonlab/rotcollab_scratch/har2011/ambulance/docs/clark_overlap/")

# Create venn diagrams
# Import data
counts_ovr_spp <- data.frame(read.csv(file="overlap_counts_spp.csv", head=TRUE, sep=","))
names_spp <- counts_ovr_spp$Sample
df_spp <- t(counts_ovr_spp[,-1])
colnames(df_spp) <- names_spp

counts_ovr_gen <- data.frame(read.csv(file="overlap_counts_gen.csv", head=TRUE, sep=","))
names_gen <- counts_ovr_gen$Sample
df_gen <- t(counts_ovr_gen[,-1])
colnames(df_gen) <- names_gen


# Creating matrix for Cem venn diagram format
m_rows <- as.numeric(df_spp["Clark_spp_count","Total/Avg", drop=TRUE] +
                       df_spp["Metaphlan_spp_count", "Total/Avg", drop=TRUE] - 
                       df_spp["Overlap_spp_count", "Total/Avg", drop=TRUE])
                     
m <- matrix(0, ncol = 2,nrow = m_rows)
for (i in 1:df_spp["Clark_spp_count","Total/Avg", drop=TRUE]){
  m[i,1] <- 1
}

for (i in (df_spp["Clark_spp_count","Total/Avg", drop=TRUE] - df_spp["Overlap_spp_count", "Total/Avg", drop=TRUE] + 1):m_rows){
  m[i,2] <- 1
}

df <- data.frame(m)

m_rows2 <- as.numeric(df_gen["Clark_gen_count","Total/Avg", drop=TRUE] +
                       df_gen["Metaphlan_gen_count", "Total/Avg", drop=TRUE] - 
                       df_gen["Overlap_gen_count", "Total/Avg", drop=TRUE])

m2 <- matrix(0, ncol = 2,nrow = m_rows2)
for (i in 1:df_gen["Clark_gen_count","Total/Avg", drop=TRUE]){
  m2[i,1] <- 1
}

for (i in (df_gen["Clark_gen_count","Total/Avg", drop=TRUE] - df_gen["Overlap_gen_count", "Total/Avg", drop=TRUE] + 1):m_rows2){
  m2[i,2] <- 1
}

df2 <- data.frame(m2)


# Output Venns
sample = data.frame(CLARK= df[,1], Metaphlan=df[,2])
venn.js(sample, "Clark and Metaphlan Species Overlap", "Clark_Meta_overlap_spp.html") 

sample2 = data.frame(CLARK=df2[,1], Metaphlan=df2[,2])
venn.js(sample2, "Clark and Metaphlan Genus Overlap", "Clark_Meta_overlap_gen.html") 

# General stats for species and genus by method
clark_spp_mean <- mean(counts_ovr_spp$Clark_spp_count[1:399], na.rm = TRUE)
clark_spp_stdE <- std.error(counts_ovr_spp$Clark_spp_count[1:399], na.rm)

meta_spp_mean <- mean(counts_ovr_spp$Metaphlan_spp_count[1:399], na.rm =TRUE)
meta_spp_stdE <- std.error(counts_ovr_spp$Metaphlan_spp_count[1:399], na.rm = TRUE)

clark_gen_mean <- mean(counts_ovr_gen$Clark_gen_count[1:399], na.rm = TRUE)
clark_gen_stdE <- std.error(counts_ovr_gen$Clark_gen_count[1:399], na.rm)

meta_gen_mean <- mean(counts_ovr_gen$Metaphlan_gen_count[1:399], na.rm =TRUE)
meta_gen_stdE <- std.error(counts_ovr_gen$Metaphlan_gen_count[1:399], na.rm = TRUE)

ovr_spp_mean <- mean(counts_ovr_spp$Overlap_spp_count[1:399], na.rm =TRUE)
ovr_spp_stdE <- std.error(counts_ovr_spp$Overlap_spp_count[1:399], na.rm = TRUE)

ovr_gen_mean <- mean(counts_ovr_gen$Overlap_gen_count[1:399], na.rm =TRUE)
ovr_gen_stdE <- std.error(counts_ovr_gen$Overlap_gen_count[1:399], na.rm = TRUE)


means <- c(clark_spp_mean, clark_gen_mean, meta_spp_mean, meta_gen_mean, ovr_spp_mean, ovr_gen_mean)
stdE <- c(clark_spp_stdE, clark_gen_stdE, meta_spp_stdE, meta_gen_stdE, ovr_spp_stdE, ovr_gen_stdE)

basicStat_df <-data.frame(means, stdE) 
rownames(basicStat_df) <- c("Clark Species", "Clark Genera", "Metaphlan Species", "Metaphlan Genera", "Overlap Species", "Overlap Genera")
colnames(basicStat_df) <- c("Mean", "Standard Error")

sink("Mean_stdError_Overlap.csv")
basicStat_df
sink()
