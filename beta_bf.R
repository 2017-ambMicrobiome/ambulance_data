source('http://bioconductor.org/biocLite.R')
biocLite('phyloseq')
require(phyloseq)
require(vegan)
require(dplyr)

#top 20
metaphlanoverlap_heatmap_top20 <- read.csv("~/Downloads/metaphlanoverlap_heatmap_top20.csv")
spp.top20 <- metaphlanoverlap_heatmap_top20[,2:ncol(metaphlanoverlap_heatmap_top20)]
colnames(metaphlanoverlap_heatmap_top20)[1] <- "Surface"
env.top20 <- metaphlanoverlap_heatmap_top20
env.top20$front[env.top20$Surface == "SteeringWheel_DriverControls" | env.top20$Surface == "Computer" | env.top20$Surface == "FrontHandles"] <- "Front"
env.top20$front <- ifelse(env.top20$front == "NA", "Back", "Front")

#group by surface then randomly sample without replacement to balance 
surf <- env.top20 %>% group_by(Surface)  
sum.surf <- summarise(surf, samp.n = n())

surfs <- filter(sum.surf, samp.n > 20)

surf.names <- as.character(surfs$Surface)
env.top20.2 <- env.top20[env.top20$Surface %in% surf.names,]
sampled.top20 <-env.top20.2 %>% group_by(Surface) %>% sample_n(size = 25, replace=FALSE) 
sampled.spp <- sampled.top20[,2:56]

#get bray curtis distances
top20.bc.dist.2 <- vegdist(sampled.spp) #default dist is bray curt
top20.bc.dist.2[is.na(top20.bc.dist.2)] <- 0 # change NA to 0

#get PCs 
mds.top20 <- dudi.pco(top20.bc.dist.2, scannf=F) 
VariationExplainedPC1 <- mds.top20$eig[1]/sum(mds.top20$eig) # compute explained var
VariationExplainedPC2 <- mds.top20$eig[2]/sum(mds.top20$eig)

#permanova for top 20 (using unbalanced)
unsamp <- vegdist(spp.top20)
unsamp[is.na(unsamp)] <- 0
adonis(unsamp~Surface, data = env.top20)

#set up plot
cols <- c("#336000", "#0000FF", "#00FF00", "#FFFF00", "Dark Red", "Orange", "#CC9966", "#00FFCC", "#0066CC", "#FF3300")

ppp <- ggplot() + coord_fixed() + 
  labs(x="Comp1, Axis1", y="Comp2, Axis2") +
  geom_hline(yintercept=0, col="darkgrey") + 
  geom_vline(xintercept=0, col="darkgrey")
# make the scree plot in a viewport
myscree <- function(eigs, x=0.8, y=0.1, just=c("right","bottom")){
  vp <- viewport(x=x, y=y, width=0.2, height=0.2, just=just)
  sp <- qplot(factor(1:length(eigs)), eigs, 
              geom="bar", stat="identity") +  
    labs(x = NULL, y = NULL)
  print(sp, vp=vp)
}

#colored by surfaces
a <- ppp + geom_point(data=data.frame(mds.top20$li, sample_data(data.frame(sampled.top20))), 
                 aes(x=A1, y=A2, col=sampled.top20$Surface), size = 2, alpha=.7) +
  labs(title="PCoA: Top 20 Species") + scale_colour_manual(values = cols) + xlab("PC1, 11.9% Variation Explained")  + ylab("PC2, 10.8% Variation Explained") + labs(colour = "Surfaces")
#colored by front v back
b <- ppp + geom_point(data=data.frame(mds.top20$li, sample_data(data.frame(sampled.top20))), 
                 aes(x=A1, y=A2, col=sampled.top20$front), size = 3, alpha=.7) +
  labs(title="PCoA: Top 20 Species") + scale_fill_hue(c=45, l=80) + xlab("PC1, 11.9% Variation Explained")  + ylab("PC2, 10.8% Variation Explained") + labs(colour = "Front Vs Back")
VariationExplainedPC1 #total variation explained for PC 1
VariationExplainedPC2 #total variation explained for PC 2

###top 50
metaphlanoverlap_heatmap_top50 <- read.csv("~/Downloads/metaphlanoverlap_heatmap_top50.csv")
rownames(metaphlanoverlap_heatmap_top50) <- metaphlanoverlap_heatmap_top50[,1]
metaphlanoverlap_heatmap_top50 <- metaphlanoverlap_heatmap_top50[,-1]
spp.top50 <- metaphlanoverlap_heatmap_top50
env.top50<- metaphlanoverlap_heatmap_top50
env.top50$Surface <- metaphlanoverlap_heatmap_top20$Surface
env.top50$front[env.top50$Surface == "SteeringWheel_DriverControls" | env.top50$Surface == "Computer" | env.top50$Surface == "FrontHandles"] <- "Front"
env.top50$front <- ifelse(env.top50$front == "NA", "Back", "Front")

#group by surface then randomly sample to balance 
env.top50.2 <- env.top50[env.top50$Surface %in% surf.names,]
sampled.top50.2 <-env.top50.2 %>% group_by(Surface) %>% sample_n(size = 25, replace=FALSE) 
spp.top50.2 <- sampled.top50.2[,1:90]
hits.top50 <- ifelse(spp.top50.2==0,0,1)
rownames(hits.top50) <- rownames(spp.top50.2)

#get bray curtis dists 
top50.bc.dist <- vegdist(spp.top50.2) # default dist is bray curt
top50.bc.dist[is.na(top50.bc.dist)] <- 0 # change NA to 0

#get pcs
mds.top50 <- dudi.pco(top50.bc.dist, scannf=F) #PCs
VariationExplainedPC1 <- mds.top50$eig[1]/sum(mds.top50$eig) #compute variation explained 
VariationExplainedPC2 <- mds.top50$eig[2]/sum(mds.top50$eig)


### permanova for full dataset
unsamp.50 <- vegdist(spp.top50)
unsamp.50[is.na(unsamp.50)] <- 0
adonis(unsamp.50~Surface, data = env.top50)

ppp + geom_point(data=data.frame(mds.top50$li, sample_data(data.frame(sampled.top50.2))), 
                 aes(x=A1, y=A2, col=sampled.top50.2$Surface), size = 2, alpha=.7) +
  labs(title="PCoA: Top 50 Species") +  scale_colour_manual(values = cols) + xlab("PC1, 13.5% Variation Explained")  + ylab("PC2, 10.7% Variation Explained") + labs(colour = "Surfaces")

d <- ppp + geom_point(data=data.frame(mds.top50$li, sample_data(data.frame(sampled.top50.2))), 
                 aes(x=A1, y=A2, col=sampled.top50.2$front), size = 3, alpha=.4) +
  labs(title="PCoA: Top 50 Species") + scale_fill_brewer(palette = "") + xlab("PC1, 13.5% Variation Explained")  + ylab("PC2, 10.7% Variation Explained") + labs(colour = "Front Vs Back")

VariationExplainedPC1 #total variation explained for PC 1
VariationExplainedPC2 #total variation explained for PC 2


 







 
