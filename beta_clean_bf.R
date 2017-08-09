require(tidyr)
require(vegan)
require(reshape)
require(dplyr)
require(plyr)
require(ape)
require(ade4)
require(ggplot2)
require(ape)
require(ecodist)
require(ade4)
require(phyloseq)
library("devtools")
require(phyloseq)
require(PerformanceAnalytics)
install.packages("corrplot")
require(corrplot)
require(Hmisc)
require(readr)
library(RColorBrewer)

## Data prep/clean
ambdata_jun16 <- read_csv("~/Dropbox/Ambulance_Study/data_sample_info/metadata_compiled/ambulance_metaphlan_hais_diversity_june_2016.csv")
amb.env <- ambdata_jun16[,1:20]
rownames(amb.env) <- amb.env$sample.ID

# move overlap files (in folder spp_outlap_files_meta_clark_v1) to local wd or set wd to dropbox: ~/Dropbox/Ambulance_Study/results/overlap_meta_clark/v1/spp_outlap_files_meta_clark_v1

files <- list.files(pattern="*_spp")
olf <- NULL
for (f in files) {
  dat <- read.csv(f, header=F, sep=",", na.strings="", colClasses="character")
  dat$file <- unlist(strsplit(f,split=".",fixed=T))[1]
  olf <- rbind(olf, dat)
}

colnames(olf) <- c("Species","Clark","Meta", "File")
meta.olf <- olf[,-2] #remove clark
# make wide instead of long
meta.olf.2 <- reshape(meta.olf, timevar = "Species", idvar = "File", direction = "wide")
rownames(meta.olf.2) <- meta.olf.2$File

meta.olf.3 <- meta.olf.2[,-2]

colnames(meta.olf.3)[1] <- "sample.ID"
meta.olf.3$sample.ID <- gsub("_clean", "", meta.olf.3$sample.ID)

rownames(amb.env) <- amb.env[,1]

## merge the env data from full dataset with the olf to get env data for overlap
olf.env <- merge(meta.olf.3, amb.env, by = "sample.ID", all.x = FALSE)
rownames(olf.env) <- olf.env[,1]

meta.olf.4 <- olf.env[,2:128]
meta.olf.4[is.na(meta.olf.4)] <- 0 #change NAs to O
#
meta.olf.5 <- apply(meta.olf.4, 2, as.numeric) #change to numeric data -- might not be necessary
rownames(meta.olf.5) <- rownames(olf.env)
rownames(meta.olf.5) <- gsub("_clean", "", rownames(meta.olf.5))
colnames(meta.olf.5) <- gsub("Meta.", "", colnames(meta.olf.5))
colnames(meta.olf.5) <- gsub("\\s+","_",colnames(meta.olf.5))
meta.full <- cbind(meta.olf.5, olf.env[,129:ncol(olf.env)]) #combine env DF with spp DF
#add front/back (location of sampled surface) column
meta.full$front <- NA # add column for front/back; front = 1, back = 0
meta.full$front[meta.full$sample.Surface == "SteeringWheel_DriverControls" | meta.full$sample.Surface == "Computer" | meta.full$sample.Surface == "FrontHandles"] <- 1
meta.full$front[is.na(meta.full$front)] <- 0
meta.full$fact.front <- factor(meta.full$front, labels = c("back", "front"))

# remove the "after" cleaning status
meta.full.1 <- meta.full[-grep("AFTER",meta.full$sample.Clean_status),]
meta.spp.1 <- meta.full.1[,1:127]

##########

overlap <- read_csv("~/Dropbox/Ambulance_Study/results/overlap_meta_clark/overlap_v1_metaphlan.csv")
spp <- overlap[,2:128]
sim <- with(overlap, simper(spp, sample.Surface))
# summary(sim)

par(mar=c(2,2,2,2))

spp <- data.frame(spp)
spp.l <- apply(spp,1,log)
spp.l <- t(spp.l)

spp.2 <- spp[,colSums(spp)>10]

# chart.Correlation(spp.2, histogram=TRUE, pch=19)

# pearson's  and spearman's

pcor <- cor(spp.2, method = "pearson")
scor <-  cor(spp.2, method = "spearman")


# get pvals, plot (green=NA)
# pearson, log-transform 
spp.3 <- log(spp.2+1)
pclog <- rcorr(as.matrix(spp.3), type = "pearson")
rlog <- pclog$r
plog <- pclog$P
corrplot(rlog, type = "lower", order = "hclust", p.mat = plog, 
         sig.level = 0.05, insig = "blank",  
         tl.col="black", tl.cex = .6, tl.srt =45, 
         col = brewer.pal(n = 9, name = "PuOr"), bg = "darkgreen")

#pearson, no transform
pcp <- rcorr(as.matrix(spp.2), type = "pearson")
M <- pcp$r
p_mat <- pcp$P

corrplot(M, type = "lower", order = "hclust", p.mat = p_mat, sig.level = 0.05, insig = "blank",  
         tl.col="black", tl.cex = .6, tl.srt =45, col = brewer.pal(n = 9, name = "PuOr"), bg = "darkgreen")

#spearman - best option for these data 

scp <- rcorr(as.matrix(spp.2), type = "spearman")
R <- scp$r
p_matp <- scp$P

corrplot(R, type = "lower", order = "hclust", p.mat = p_matp, sig.level = 0.05, insig = "blank",  
        tl.col="black", tl.cex = .6, tl.srt =45, col = brewer.pal(n = 9, name = "PuOr"), bg = "darkgreen")


corrplot(R, type = "upper", order = "hclust", p.mat = p_matp, sig.level = 0.05, insig = "blank",  
         tl.col="black", tl.cex = .6, tl.srt =45, col = brewer.pal(n = 9, name = "PuOr"), bg = "darkgreen")

## flatten output
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

cms <- flattenCorrMatrix(scp$r, scp$P)
cms2 <- filter(cms, abs(cor > .5), p <= .05)
lowcor <- filter(cms, abs(cor<.11),  p <= .05)
View(cms2)
View(lowcor)
summary(cms2)
summary(lowcor)

#check out some plots: s>p for most 

strep <- spp$Streptococcus_mitis
acid <- spp$Acidovorax_ebreus
plot(strep, acid)

cmp <- flattenCorrMatrix(pcp$r, pcp$P)
cmp2 <- filter(cmp, abs(cor > .5), p < .05)
View(cmp2)
prop <- spp$Propionibacterium_acnes
geo <- spp$Geodermatophilus_obscurus
roth <- spp$Rothia_mucilaginosa
strep <- spp$Streptococcus_mitis
plot(prop, roth)
plot(prop, geo)
plot(strep, roth)
cor(roth, geo)
cor(strep, roth)


##### PCoA -- balance data, remove ec, etc

##Region
#see which regions have small n 
reg <- overlap  %>% group_by(region) 
sum.reg <- summarise(reg, samp.n = n())

##remove eastcoast 
no.ec<- overlap[-grep("e_coast",overlap$region ),]


## sample from regions n=36
sampled.region <- no.ec %>% group_by(region) %>% sample_n(size = 36, replace=FALSE) 
sampled.spp <- sampled.region[,2:128]
sampled.spp <- data.frame(sampled.spp)
spp.hel <- decostand(sampled.spp, "range", MARGIN = 1) # standardize rows [0,1]

volf.sn <- vegdist(spp.hel)
volf.sn[is.na(volf.sn)] <- 0 

mds.volf.sn <- dudi.pco(volf.sn, scannf=F)
VariationExplainedPC1 <- mds.volf.sn$eig[1]/sum(mds.volf.sn$eig)
VariationExplainedPC2 <- mds.volf.sn$eig[2]/sum(mds.volf.sn$eig)

#set up plot

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

ppp + geom_point(data=data.frame(mds.volf.sn$li, sample_data(data.frame(sampled.region))), 
                      aes(x=A1, y=A2, col=sampled.region$region), size = 2, alpha=.6) +
  labs(title="PCoA: Regions") + scale_fill_hue(c=45, l=80) + xlab("PC1, 22.9% Variation Explained")  + ylab("PC2, 17.5% Variation Explained") + labs(colour = "Regions")

VariationExplainedPC1
VariationExplainedPC2

##Surfaces
#see what surfaces have low n
surf <- overlap %>% group_by(sample.Surface)  
sum.surf <- summarise(surf, samp.n = n())

surfs <- filter(sum.surf, samp.n > 20)

surf.names <- as.character(surfs$sample.Surface)
overlap.2 <- overlap[overlap$sample.Surface %in% surf.names,]

## sample from surfaces n=27 

sampled.ol <- overlap.2 %>% group_by(sample.Surface) %>% sample_n(size = 27, replace=FALSE) 
sampled.spp <- sampled.ol[,2:128]

## cor by surface 

sp <- sampled.spp[colSums(sampled.spp)>10]
sp$sample.Surface <- sampled.ol$sample.Surface
rho <- dlply(sp, .(sample.Surface), function(x) rcorr(as.matrix(x[,1:33]), type="spearman")$r)
P <- dlply(sp, .(sample.Surface), function(x) rcorr(as.matrix(x[,1:33]), type="spearman")$P)

nanfun <- function(x){
  xx <- unlist(x)
  xx[is.nan(xx)] <- 0
  return(xx)
}

R <- lapply(rho, nanfun)
Pvals <- lapply(P, nanfun)


Corfun <- function(cormat) {
  cormat <- as.data.frame(cormat)
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut]
  )
}

Pfun <- function(pv) {
  pv <- as.data.frame(pv)
  ut <- upper.tri(pv)
  data.frame(p = pv[ut])
}

cors <- lapply(R, Corfun)
pvs <- lapply(Pvals, Pfun)

comp <- data.frame(cors$Computer, pvs$Computer)
c <- with(comp, comp[order(cor, decreasing = TRUE),])[1:3,]
c$surface <- rep("computer", nrow(c))

fh <- data.frame(cors$FrontHandles, pvs$FrontHandles)
f <- with(fh, fh[order(cor, decreasing = TRUE),])[1:3,]
f$surface <- rep("front handles", nrow(f))

steth <- data.frame(cors$Stethoscope, pvs$Stethoscope)
ste <- with(steth, steth[order(cor, decreasing = TRUE),])[1:3,]
ste$surface <- rep("stethoscope", nrow(ste))

stretch <- data.frame(cors$Stretcher, pvs$Stretcher)
str <- with(stretch, stretch[order(cor, decreasing = TRUE),])[1:3,]
str$surface <- rep("stretcher", nrow(str))


rbs <- data.frame(cors$RearBench_seats, pvs$RearBench_seats)
rb <- with(rbs, rbs[order(cor, decreasing = TRUE),])[1:3,]
rb$surface <- rep("rear bench seats", nrow(rb))

rcc <- data.frame(cors$RearCabinets_counters, pvs$RearCabinets_counters)
rc <- with(rcc, rcc[order(cor, decreasing = TRUE),])[1:3,]
rc$surface <- rep("rear cabinet controls", nrow(rc))

rhr <- data.frame(cors$RearHandles_rails, pvs$RearHandles_rails)
rh <- with(rhr, rhr[order(cor, decreasing = TRUE),])[1:3,]
rh$surface <- rep("rear handles rails", nrow(rh))

rlcp <- data.frame(cors$RearLights_controlPanel, pvs$RearLights_controlPanel)
rlc <- with(rlcp, rlcp[order(cor, decreasing = TRUE),])[1:3,]
rlc$surface <- rep("rear lights control", nrow(rlc))

swdc <- data.frame(cors$SteeringWheel_DriverControls, pvs$SteeringWheel_DriverControls)
swd <- with(swdc, swdc[order(cor, decreasing = TRUE),])[1:3,]
swd$surface <- rep("steering wheel driver controls", nrow(swd))

suct <- data.frame(cors$Suction_O2, pvs$Suction_O2)
suc <- with(suct, suct[order(cor, decreasing = TRUE),])[1:3,]
suc$surface <- rep("steering wheel driver controls", nrow(suc))

top3 <- rbind(c, f, ste, str, rb, rc, rh, rlc, swd, suc)






corrplot(as.matrix(comp.R), type = "lower", order = "hclust", p.mat = as.matrix(comp.P), sig.level = 0.05, insig = "blank",  
         tl.col="black", tl.cex = .6, tl.srt =45, col = brewer.pal(n = 9, name = "PuOr"), bg = "darkgreen")

comp.cm <- flattenCorrMatrix(R$Computer, Pvals$Computer)
comp.highcor <- filter(comp.cm, abs(cor > .5), p <= .05)
comp.nocor <- filter(comp.cm, cor == 0,  p <= .05)




fh <- flattenCorrMatrix(R$FrontHandles, Pvals$FrontHandles)
fh$surface <- rep("fh", nrow(fh))
steth <- flattenCorrMatrix(R$Stethoscope, Pvals$Stethoscope)
steth$surface <- rep("steth", nrow(steth))
stretch <- flattenCorrMatrix(R$Stretcher, Pvals$Stretcher)
stretch$surface <- rep("stretch", nrow(stretch))
rbs <- flattenCorrMatrix(R$RearBench_seats, Pvals$RearBench_seats)
rbs$surface <- rep("rbs", nrow(rbs))
rcc <- flattenCorrMatrix(R$RearCabinets_counters, Pvals$RearCabinets_counters)
rcc$surface <- rep("rcc", nrow(rcc))







## standardize 

spp.3 <- decostand(sampled.spp, "range", MARGIN = 1) #standardized row values [0,1]

#####

volf.sn.2 <- vegdist(spp.3)
volf.sn.2[is.na(volf.sn.2)] <- 0 

mds.volf.sn.2 <- dudi.pco(volf.sn.2, scannf=F)
VariationExplainedPC1 <- mds.volf.sn.2$eig[1]/sum(mds.volf.sn.2$eig)
VariationExplainedPC2 <- mds.volf.sn.2$eig[2]/sum(mds.volf.sn.2$eig)
VariationExplainedPC1
VariationExplainedPC2

ppp + geom_point(data=data.frame(mds.volf.sn.2$li, sample_data(data.frame(sampled.ol))), 
                      aes(x=A1, y=A2, col=sampled.ol$sample.Surface), size = 2, alpha=.6) +
  labs(title="PCoA: Surfaces") + scale_fill_hue(c=45, l=80) + xlab("PC1, 23.7% Variation Explained") + ylab("PC2, 17.1% Variation Explained") + labs(colour = "Surfaces")

#front/back
ppp + geom_point(data=data.frame(mds.volf.sn.2$li, sample_data(data.frame(sampled.meta.2))), 
                 aes(x=A1, y=A2, col=sampled.meta.2$fact.front), size = 2, alpha=.6) +
  labs(title="PCoA: Front vs Back of Ambulances") + scale_fill_hue(c=45, l=80)


#### permanova 

a <- vegdist(spp.2)
a[is.na(a)] <- 0

perma.1 <- adonis(a~region, data = overlap, permutations = 2000)
perma.2 <- adonis(a~sample.Surface, data = overlap, permutations = 2000)
perma.3 <- adonis(a~region*sample.Surface, data = overlap, permutations = 2000)
round(perma.2$aov.tab[1:6], 4)
round(perma.1$aov.tab[1:6], 4)













