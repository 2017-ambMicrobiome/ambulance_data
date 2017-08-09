library(ggplot2)
library(devtools)
source_gist("524eade46135f6348140")
library(wesanderson)

# Working dir
setwd("/home/jake/ambulance/results/CARD_db_abundance/31kmer_150score")

# Import daat
df_overlap <- data.frame(read.csv("overlap.csv"))
df_eval_tot_hits <- data.frame(read.csv("tot_hits_eval.csv"))
colnames(df_eval_tot_hits) <- c("Sample", "hits")
df_vir_tot_hits <- data.frame(read.csv("tot_hits_31kmer_card.csv"))
colnames(df_vir_tot_hits) <- c("Sample", "hits")
df_eval_tot_hits[is.na(df_eval_tot_hits)] <- 0
df_vir_tot_hits[is.na(df_vir_tot_hits)] <- 0

# Determine means
eval_mean <- mean(df_eval_tot_hits$hits, na.rm = TRUE)
vir_mean <- mean(df_vir_tot_hits$hits, na.rm = TRUE)
df_overlap <- cbind(df_overlap, df_overlap$Overlap./df_overlap$Vir_marker.)
df_overlap$align <- 'align'
colnames(df_overlap) <- c("Sample", "Virulence Spp. #", "Number of Spp. Overlap", "Percent Overlap", "align")
df_overlap[is.na(df_overlap)] <- 0
df_overlap_mean_per <- mean(df_overlap$`Percent Overlap`)
hist_bub <- hist(df_overlap$`Percent Overlap`)
hist_bub$align <- 'align'
hist_bub_df <- cbind.data.frame(hist_bub$breaks[2:11], hist_bub$counts, hist_bub$density, hist_bub$align)
colnames(hist_bub_df) <- c("Breaks", "Counts", "Density", "Align")

# Scaling for histogram
df_eval_tot_hits$scale <- scale(df_eval_tot_hits$hits)
df_vir_tot_hits$scale <- scale(df_vir_tot_hits$hits)

# Processing for line Graph
df_combined_ln <- cbind.data.frame(df_vir_tot_hits$Sample, df_vir_tot_hits$hits, df_eval_tot_hits$hits)
colnames(df_combined_ln) <- c("Sample", "Vir.hits", "Bac.hits")

# Labeling types of databases in database column and combining
df_eval_tot_hits$database <- 'bacterial'
df_vir_tot_hits$database <- 'virulence'
df_com_hist_hits <- rbind.data.frame(df_eval_tot_hits$scale, df_vir_tot_hits$scale)
df_com_hits_scale_ln <- cbind.data.frame(df_eval_tot_hits$scale, df_vir_tot_hits$scale)
colnames(df_com_hits_scale_ln) <- c("Bacs.scale", "Vir.scale")
df_com_hist_db <- rbind.data.frame(df_eval_tot_hits$database, df_vir_tot_hits$database)
df_com_hist <- cbind(df_com_hist_hits, df_com_hist_db)
colnames(df_com_hist) <- c("hits.scale", "database")

# Plotting
ggplot(df_combined_ln, aes(x = Bac.hits, y = Vir.hits)) + 
  geom_point(shape=1) +
  geom_smooth(method = "lm") +
  labs(x = "Bacteria DB hits", y = "Virulence DB hits", title = "Scatter Plot and Regression Comparing Databases") +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))

ggplot(df_com_hist, aes(x=hits.scale, fill=database, color=database)) + 
  geom_histogram(position = "identity", alpha=0.2, binwidth = 0.3) +
  labs(x = "Z-Score Scaled Counts", y = "Frequency", title = "Comparing Distribution of Counts between DBs") +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))

ggplot(hist_bub_df, aes(x = Align, y = Breaks, size = Counts, colour = Density)) + 
  geom_point(alpha=0.4) +
  scale_size(range = c(1,30)) +
  labs(x = "Bubbles", y = "Percent Overlap", title = "Percent of Overlapping Classifications") +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))

