library("pheatmap")
library(vegan)
library(Rtsne)
library(ggplot2)
source("/home/jake/R_scripts/functions/multiplot.R")

# Import data
setwd("/home/jake/Dropbox/ambulance/clark_abundance/CARD_db_class/31kmer/score_50_all")
filenames <- list.files(path=getwd())
df_list = lapply(filenames, read.table, sep = ",", fill = TRUE)
filenames <- noquote(as.character(filenames))
df_meta <- data.frame(read.csv("/home/jake/Dropbox/Ambulance_Study/data_sample_info/metadata_compiled/ambulance_metaphlan_hais_diversity_june_2016.csv"))

#Only Pertinent Meta
df_meta <- df_meta[df_meta$sequenced == 1,1:20]

# Setup matrix
n.rows <- length(df_list[[1]][,1])
n.cols <- 0
data_mat <- matrix(nrow = n.rows, ncol = n.cols)

# Column bind all documents
for (column in 1:length(filenames)){
  data_mat_combined <- cbind(data_mat,df_list[[column]][,2])
  data_mat <- data_mat_combined
}

# Create shortened row names
# noquote(strsplit(colnames(data_mat), "_")[[3]][5])
lrowNM <- length(as.character(df_list[[1]][,1]))
rowNM <- noquote(strsplit(as.character(df_list[[1]][1,1]), "_")[[1]][5])
for (rowName in 2:lrowNM){
  rowNM_com <- rbind(rowNM, noquote(strsplit(as.character(df_list[[1]][rowName,1]), "_")[[1]][5]))
  rowNM <- rowNM_com
}
rowNM <- noquote(rowNM)

len_st_NM <- length(df_meta[,9])
st_NM <- sapply(strsplit(as.character(df_meta[1,9]), "_"), tail, 1)
for (stateName in 2:len_st_NM){
  st_NM_com <- rbind(st_NM, sapply(strsplit(as.character(df_meta[stateName,9]), "_"), tail, 1))
  st_NM <- st_NM_com
}

rownames(data_mat) <- rowNM
colnames(data_mat) <- df_meta[,1]
#metaphlan_m_colLT20k <- metaphlan_m[,colSums(metaphlan_m) < 12000]

t_data_mat <- t(data_mat)
t_data_mat_colSumGT0 <- t_data_mat[,colSums(t_data_mat) > 0]
t_data_mat_nz <- t_data_mat_colSumGT0
t_data_mat_colSumGT0LT30k <- t_data_mat[,colSums(t_data_mat_colSumGT0) < 30000]

# Log transform data
log_mat <- log10(t_data_mat+1)
log2_mat <- log2(t_data_mat+1)
log_nz_mat <- log10(t_data_mat_colSumGT0+1)

# Run rtsne for each
rtsne_data_mat <- Rtsne(log_mat, max_iter = 100000, check_duplicates = FALSE , perplexity = 100, theta = 0, initial_dims = 100, 
                        verbose = TRUE)

rtsne_data_nz_mat <- Rtsne(log_nz_mat, max_iter = 100000, check_duplicates = FALSE , perplexity = 100, theta = 0, initial_dims = 100, 
                                             verbose = TRUE)

# Produce plots
setwd("/home/jake/Dropbox/ambulance/results/CARD_db_abundance/31kmer_50score")
jpeg(filename = "ClarkCard_Ambulance_allSamples_byCity.jpeg",
     width = 15000, height = 10000, units = "px", pointsize = 12)
pheatmap(log_mat)
dev.off()

df_plot <- data.frame( st_NM, df_meta[,9], df_meta[,12], df_meta[,3], df_meta[,7], df_meta[,10], df_meta[,11], df_meta[,19], 
                       rtsne_data_mat$Y[,1], rtsne_data_mat$Y[,2])
colnames(df_plot) <- c("state", "city", "region", "surface", "date", "lat", "long", "temp", "x", "y")

city_all <- ggplot(df_plot, aes(x=x, y=y, colour=city)) +
  geom_point(size=2) +
  labs(title="t-SNE CLARK CARD by City") + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave(city_all, filename = "tsne_city_all_clark_card.png", path = "/home/jake/Dropbox/ambulance/results/CARD_db_abundance/31kmer_50score", 
       width = 8, height = 8, units = "in", dpi = 300, 
       limitsize = FALSE, device = "png")

region_all <-ggplot(df_plot, aes(x=x, y=y, colour=region)) +
  geom_point(size=2) +
  labs(title="t-SNE CLARK CARD by Region") + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave(region_all, filename = "tsne_region_all_clark_card.png", path = "/home/jake/Dropbox/ambulance/results/CARD_db_abundance/31kmer_50score", 
       width = 8, height = 8, units = "in", dpi = 300, 
       limitsize = FALSE, device = "png")

surface_all <- ggplot(df_plot, aes(x=x, y=y, colour=surface)) +
  geom_point(size=2) + 
  labs(title="t-SNE CLARK CARD by Surface") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(surface_all, filename = "tsne_surface_all_clark_card.png", path = "/home/jake/Dropbox/ambulance/results/CARD_db_abundance/31kmer_50score", 
       width = 8, height = 8, units = "in", dpi = 300, 
       limitsize = FALSE, device = "png")

state_all <- ggplot(df_plot, aes(x=x, y=y, colour=state)) +
  geom_point(size=2) + 
  labs(title="t-SNE CLARK CARD by State") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(state_all, filename = "tsne_state_all_clark_card.png", path = "/home/jake/Dropbox/ambulance/results/CARD_db_abundance/31kmer_50score", 
       width = 8, height = 8, units = "in", dpi = 300, 
       limitsize = FALSE, device = "png")

# Produce non-zero plots

setwd("/home/jake/Dropbox/ambulance/results/CARD_db_abundance/31kmer_50score")
jpeg(filename = "ClarkCard_Ambulance_allSamples_byCity_nonzero.jpeg",
     width = 3000, height = 5000, units = "px", pointsize = 12)
pheatmap(log_nz_mat)
dev.off()

df_plot_nz <- data.frame( st_NM, df_meta[,9], df_meta[,12], df_meta[,3], df_meta[,7], df_meta[,10], df_meta[,11], df_meta[,19], 
                       rtsne_data_nz_mat$Y[,1], rtsne_data_nz_mat$Y[,2])
colnames(df_plot_nz) <- c("state", "city", "region", "surface", "date", "lat", "long", "temp", "x", "y")

city_all <- ggplot(df_plot_nz, aes(x=x, y=y, colour=city)) +
  geom_point(size=2) +
  labs(title="t-SNE CLARK CARD by City") + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave(city_all, filename = "nonzero_tsne_city_all_clark_card.png", path = "/home/jake/Dropbox/ambulance/results/CARD_db_abundance/31kmer_50score", 
       width = 8, height = 8, units = "in", dpi = 300, 
       limitsize = FALSE, device = "png")

region_all <-ggplot(df_plot_nz, aes(x=x, y=y, colour=region)) +
  geom_point(size=2) +
  labs(title="t-SNE CLARK CARD by Region") + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave(region_all, filename = "nonzero_tsne_region_all_clark_card.png", path = "/home/jake/Dropbox/ambulance/results/CARD_db_abundance/31kmer_50score", 
       width = 8, height = 8, units = "in", dpi = 300, 
       limitsize = FALSE, device = "png")

surface_all <- ggplot(df_plot_nz, aes(x=x, y=y, colour=surface)) +
  geom_point(size=2) + 
  labs(title="t-SNE CLARK CARD by Surface") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(surface_all, filename = "nonzero_tsne_surface_all_clark_card.png", path = "/home/jake/Dropbox/ambulance/results/CARD_db_abundance/31kmer_50score", 
       width = 8, height = 8, units = "in", dpi = 300, 
       limitsize = FALSE, device = "png")

state_all <- ggplot(df_plot_nz, aes(x=x, y=y, colour=state)) +
  geom_point(size=2) + 
  labs(title="t-SNE CLARK CARD by State") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(state_all, filename = "nonzero_tsne_state_all_clark_card.png", path = "/home/jake/Dropbox/ambulance/results/CARD_db_abundance/31kmer_50score", 
       width = 8, height = 8, units = "in", dpi = 300, 
       limitsize = FALSE, device = "png")

# Working on combined
rownames(t_data_mat_nz) <- df_meta[,12]
# Coasts loop to sum columns
coasts.list <- c("s_w_w_coast", "e_coast", "^w_coast", "^w$", "s_e")
coasts <- list()
colSum_coasts <- list()
for (x in 1:length(coasts.list)){
  coasts[[x]] <- t_data_mat_nz[grep(coasts.list[x], rownames(t_data_mat_nz)),]
  colSum_coasts[[x]] <- colSums(coasts[[x]])
}
regionCombine <- rbind(colSum_coasts[[1]], colSum_coasts[[2]], colSum_coasts[[3]], colSum_coasts[[4]], colSum_coasts[[5]])
rownames(regionCombine) <- c("South West West Coast", "East Coast", "West Coast", "West", "South East")

# Log transform and plot tsne and heatmap
regionCombine_log <- log10(regionCombine+1)

setwd("/home/jake/Dropbox/ambulance/results/CARD_db_abundance/31kmer_50score")
jpeg(filename = "ClarkCard_Ambulance_allSamples_heatmap_byRegion.jpeg",
     width = 1600, height = 500, units = "px", pointsize = 12)
pheatmap(regionCombine_log)
dev.off()


regionCombine_rtsne <- Rtsne(regionCombine_log, max_iter = 100000, check_duplicates = FALSE , perplexity = 1, theta = 0, initial_dims = 100, 
                             verbose = TRUE)
df_combinePlot <- data.frame(rownames(regionCombine), regionCombine_rtsne$Y[,1], regionCombine_rtsne$Y[,2])
colnames(df_combinePlot) <- c("region", "x", "y")
plot(regionCombine_rtsne$Y)

combine_region <- ggplot(df_combinePlot, aes(x=x, y=y, colour=region)) +
  geom_point(size=3) + 
  labs(title="t-SNE CLARK CARD Combined Region") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(combine_region, filename = "nonzero_tsne_combinedRegion_all_clark_card.png", path = "/home/jake/Dropbox/ambulance/results/CARD_db_abundance/31kmer_50score", 
       width = 6, height = 6, units = "in", dpi = 300, 
       limitsize = FALSE, device = "png")



#simplify, sum all cities and look for region. sum all cities and look for state etc
# City list and combine
rownames(t_data_mat_nz) <- df_meta[,9]
city.list <- c("Brooklyn", "Broomfield", "Davie", "Lakewood", "Mesa", "Miami", "Plantation", "Reno", "San_Diego", "San_Leandro",
               "Scottsdale", "Sierra_Vista", "Sun_City_West", "Tempe", "Tuscon", "Vale", "Yuma")
cities <- list()
sumcolumns_cities <-list()
for (x in 1:length(city.list)){
  cities[[x]] <- t_data_mat_nz[grep(city.list[x], rownames(t_data_mat_nz)),]
  sumcolumns_cities[[x]] <- colSums(cities[[x]])
}
cityCombine <- rbind(sumcolumns_cities[[1]], sumcolumns_cities[[2]], sumcolumns_cities[[3]], sumcolumns_cities[[4]],
                     sumcolumns_cities[[5]], sumcolumns_cities[[6]], sumcolumns_cities[[7]], sumcolumns_cities[[8]],
                     sumcolumns_cities[[9]], sumcolumns_cities[[10]], sumcolumns_cities[[11]], sumcolumns_cities[[12]],
                     sumcolumns_cities[[13]], sumcolumns_cities[[14]], sumcolumns_cities[[15]], sumcolumns_cities[[16]],
                     sumcolumns_cities[[17]])
rownames(cityCombine) <- city.list

cityCombine_log <- log10(cityCombine+1)

pheatmap(cityCombine_log)

cityCombine_log_rtsne <- Rtsne(cityCombine_log, max_iter = 100000, check_duplicates = FALSE , perplexity = 5, theta = 0, initial_dims = 100, 
                               verbose = TRUE)

#Create list to add to data frame column
df_combinePlot <- data.frame(rownames(cityCombine), cityCombine_log_rtsne$Y[,1], cityCombine_log_rtsne$Y[,2])
df_combinePlot["region"] <- c("east coast", "west", "south east", "west", "south west west coast", "south east", 
                              "south east", "west coast", "south west west coast", "west coast", "south west west coast", 
                              "south west west coast", "south west west coast", "south west west coast", "south west west coast",
                              "south west west coast", "south west west coast")
df_combinePlot["state"] <- c("NY", "CO", "FL", "CO", "AZ", "FL", "FL", "NV", "CA", "CA", "AZ", "AZ", "AZ", "AZ", "AZ", "AZ", 
                             "AZ")


colnames(df_combinePlot) <- c("city", "x", "y", "region", "state")
plot(cityCombine_log_rtsne$Y)

combine_regionState <- ggplot(df_combinePlot, aes(x=x, y=y, colour=region, shape=state)) +
  geom_point(size=3) + 
  labs(title="t-SNE CLARK CARD Combined Region and State") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(combine_regionState, filename = "nonzero_tsne_combinedRegionandState_all_clark_card.png", path = "/home/jake/Dropbox/ambulance/results/CARD_db_abundance/31kmer_50score", 
       width = 6, height = 6, units = "in", dpi = 300, 
       limitsize = FALSE, device = "png")

ggplot(df_combinePlot, aes(x=x, y=y, colour=region)) + 
  geom_point(size=3)