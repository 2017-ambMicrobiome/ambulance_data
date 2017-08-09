# Declare libraries 
library(vegan)
library(Rtsne)
library(ggplot2)
source("/home/jake/R_scripts/functions/multiplot.R")

# Import data
# Meta data and overlap diversity
df_meta <- data.frame(read.csv("/home/jake/Dropbox/Ambulance_Study/metadata_compiled/ambulance_metaphlan_hais_diversity_june_2016_mod.csv"))
df_div <- data.frame(read.csv("/home/jake/Dropbox/Ambulance_Study/results/meta_clark_overlap_diversity/Diversity.indices.overlap.csv"))
colnames(df_div) <- (c("sample.ID", "Simpsons.ovr.index", "Shannon.ovr.index"))

# Metaphlan results
df_metaphlan_data <- data.frame(read.table("~/ambulance_cleanedup_summary.txt", sep ="\t", header=TRUE))
id <- df_metaphlan_data$ID
df_metaphlan_data <- as.data.frame(t(df_metaphlan_data[,-1]))
colnames(df_metaphlan_data) <- id
df_meta2 <- data.frame(read.csv("/home/jake/Dropbox/Ambulance_Study/metadata_compiled/ambulance_metaphlan_hais_diversity_june_2016.csv"))

# Overlap data
spp_amb_overlap <- read.table("/home/jake/ambulance/overlap_diversity/spp_only_amb.biom.tsv", sep = "\t", header = TRUE, numerals = c("allow.loss"), colClasses = "numeric")
id <- spp_amb_overlap$OTU.ID
spp_amb_overlap <- as.data.frame(t(spp_amb_overlap[,-1]))
colnames(spp_amb_overlap) <- id

# Overlap functional data
overlap_functional <- read.table("/home/jake/Dropbox/feature_selection_ambulance/all399_allnodes_seed_name_2.txt", sep = "\t", header = TRUE, numerals = c("allow.loss"), colClasses = c("character","numeric"))
# Removing SEED row
overlap_functional <- overlap_functional[2:2519,]
id <- overlap_functional$Datasets
overlap_functional <- as.data.frame(t(overlap_functional[,-1]))
colnames(overlap_functional) <- id
  #Remove null from dataset and AS1510
 # overlap_functional <- overlap_functional[rownames(overlap_functional) != "null",]
  #overlap_functional <- overlap_functional[rownames(overlap_functional) != "AS1510",]
  # Replace variable metabolic pathways with numbers
  for ( i in 1:2518){
    colnames(overlap_functional)[i] <- i
  }
    
# Keep only pertinent meta
df_meta <- df_meta[df_meta$sequenced == 1,1:20]
df_meta2 <- df_meta2[df_meta2$sequenced == 1,1:20]

# Merge dfs
df_merge <- merge(df_meta, df_div, by.x=c("sample.ID"))

# Convert to matrices
#m_shannon <- as.matrix(df_merge[,21])
#row.names(m_shannon) <- df_merge[,1]
#m_shannon <- unique(m_shannon)
#m_simpsons <- as.matrix(df_merge[,22])
#row.names(df_merge) <- df_merge[,1]
#m_simpsons <- unique(m_simpsons)

combined_m <- as.matrix(df_merge[,c(21,22)])
row.names(combined_m) <- df_merge[,1]
#combined_m <- unique(combined_m)


# Principal component analyses
#p_shannon <- prcomp(m_shannon)
#p_simpsons <- prcomp(m_simpsons)
#plot(p_shannon$x)
#ordiplot(p_shannon$x)

combined_m_log <- log10(combined_m+1)

# Run Rtsne
combined_rtsne <- Rtsne(combined_m_log, check_duplicates = FALSE)

# Combine rtsne output with meta data
df_merge_out <- df_merge
df_merge_out["tsne_x"] <- combined_rtsne$Y[,2]
df_merge_out["tsne_y"] <- combined_rtsne$Y[,1]

# Plots
div_sur <- ggplot(df_merge_out, aes(x=tsne_x, y=tsne_y, colour=sample.Surface )) +
  geom_point(size=3) +
  labs(title="t-SNE of Diversity Colored by Sample Surface")

div_serviceType <- ggplot(df_merge_out, aes(x=tsne_x, y=tsne_y, colour=sample.Service_type )) +
  geom_point(size=3) +
  labs(title="t-SNE of Diversity Colored by Service Type")

div_region <- ggplot(df_merge_out, aes(x=tsne_x, y=tsne_y, colour=region )) +
  geom_point(size=3) +
  labs(title="t-SNE of Diversity Colored by Region")

div_city <- ggplot(df_merge_out, aes(x=tsne_x, y=tsne_y, colour=city )) +
  geom_point(size=3) +
  labs(title="t-SNE of Diversity Colored by City")

div_temp <- ggplot(df_merge_out, aes(x=tsne_x, y=tsne_y, colour=mean_temp_F )) +
  geom_point(size=3) +
  labs(title="t-SNE of Diversity Colored by Temperature")

multiplot(div_sur, div_serviceType, div_region, div_city, div_temp, cols=2)

# Run Rtsne

df_metaphlan_data_log <- log10(df_metaphlan_data+1)

metaphlan_rtsne <- Rtsne(df_metaphlan_data_log, check_duplicates = FALSE)

df_metaphlan_out <- df_meta2
df_metaphlan_out["tsne_x"] <- metaphlan_rtsne$Y[,1]
df_metaphlan_out["tsne_y"] <- metaphlan_rtsne$Y[,2]

# Plots
orig_plot <- plot(metaphlan_rtsne$Y)

meta_sur <- ggplot(df_metaphlan_out, aes(x=tsne_x, y=tsne_y, colour=sample.Surface )) +
  geom_point(size=2) +
  labs(title="t-SNE of Metaphlan Results Colored by Surface")

meta_service <- ggplot(df_metaphlan_out, aes(x=tsne_x, y=tsne_y, colour=sample.Service_type )) +
  geom_point(size=2) +
  labs(title="t-SNE of Metaphlan Results Colored by Service Type")

meta_region <- ggplot(df_metaphlan_out, aes(x=tsne_x, y=tsne_y, colour=region )) +
  geom_point(size=2) +
  labs(title="t-SNE of Metaphlan Results Colored by Region")

meta_city <- ggplot(df_metaphlan_out, aes(x=tsne_x, y=tsne_y, colour=city )) +
  geom_point(size=2) +
  labs(title="t-SNE of Metaphlan Results Colored by City")

meta_temp <- ggplot(df_metaphlan_out, aes(x=tsne_x, y=tsne_y, colour=mean_temp_F )) +
  geom_point(size=2) +
  labs(title="t-SNE of Metaphlan Results Colored by Temperature")

multiplot(meta_sur, meta_service, meta_region, meta_city, meta_temp, cols = 2)

# Run Rtsne

spp_amb_overlap_log <- log10(spp_amb_overlap+1)

overlap_rtsne <- Rtsne(spp_amb_overlap_log, check_duplicates = FALSE)

# Combine rtsne output with meta data
df_overlap_out <- df_merge
df_overlap_out["tsne_x"] <- overlap_rtsne$Y[,1]
df_overlap_out["tsne_y"] <- overlap_rtsne$Y[,2]

#Plots
plot(overlap_rtsne$Y)

ovr_sur <- ggplot(df_overlap_out, aes(x=tsne_x, y=tsne_y, colour=sample.Surface )) +
  geom_point(size=2) +
  labs(title="t-SNE of Overlap Results Colored by Surface")

ovr_service <- ggplot(df_overlap_out, aes(x=tsne_x, y=tsne_y, colour=sample.Service_type )) +
  geom_point(size=2) +
  labs(title="t-SNE of Overlap Results Colored by Service Type")

ovr_region <- ggplot(df_overlap_out, aes(x=tsne_x, y=tsne_y, colour=region )) +
  geom_point(size=2) +
  labs(title="t-SNE of Overlap Results Colored by Region")

ovr_city <- ggplot(df_overlap_out, aes(x=tsne_x, y=tsne_y, colour=city )) +
  geom_point(size=2) +
  labs(title="t-SNE of Overlap Results Colored by City")

ovr_temp <- ggplot(df_overlap_out, aes(x=tsne_x, y=tsne_y, colour=mean_temp_F )) +
  geom_point(size=2) +
  labs(title="t-SNE of Overlap Results Colored by Temperature")

multiplot(ovr_sur, ovr_service, ovr_region, ovr_city, ovr_temp, cols = 2)

overlap_functional_rtsne <- Rtsne(overlap_functional)
