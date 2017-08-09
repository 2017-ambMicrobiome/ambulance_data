#!/usr/bin/env Rscript

# load packages
lapply(c('ggplot2','stringr','dplyr','purrr','lubridate','readr','tidyr','scales','gridExtra','ggrepel',
         'biom','phyloseq','vegan','RColorBrewer','lme4','DESeq2',
         'glmnet','caret','snow','doSNOW','randomForest',
         'pROC','e1071','fastICA','Rtsne','MLmetrics'),
       require,character.only=TRUE)

work_os <- ifelse(Sys.info()['sysname'] == 'Linux','~','D:')
work_path <- sprintf('%s/Dropbox/Hospital',work_os)
out_path <- file.path(work_path,'Out')
source(file.path(work_path,'Code/hospital_fxns.R'))


# prep data
# hospital meta data
# create groups for front/rear and interior/exterior
META <- read_delim(file.path(work_path,'Data/ambulance_metadata_june_2016_new.txt'),'\t') %>%
  dplyr::rename(SampleID=`#SampleID`) %>%
  as.data.frame() %>%
  dplyr::mutate(FR=ifelse(sample.Surface %in% 
                            c('SteeringWheel_DriverControls','FrontHandles','Computer'),
                          'Front','Rear'),
                IE=ifelse(sample.Surface %in% 
                            c('RearHandles_rails','FrontHandles'),'Exterior','Interior'))
rownames(META) <- META$SampleID

# load overlap data
OTU_OVERLAP <- read_csv(file.path(work_path,'Data/metaphlan_overlap.csv'))
TAXA_OVERLAP <- OTU_OVERLAP[,1]
OTU_OVERLAP <- t(OTU_OVERLAP[,-1])
colnames(OTU_OVERLAP) <- unlist(TAXA_OVERLAP)
OTU_OVERLAP <- OTU_OVERLAP[rowSums(OTU_OVERLAP) > 0, colSums(OTU_OVERLAP) > 0]
TAXA_OVERLAP <- colnames(OTU_OVERLAP)
META_OVERLAP <- META[rownames(OTU_OVERLAP),]

# load metaphlan data
OTU_METAPHLAN <- read_delim(file.path(work_path,'Data/ambulance_cleanedup_summary.txt'),delim='\t')[-1,]
TAXA_METAPHLAN <- unlist(OTU_METAPHLAN[,1])
OTU_METAPHLAN <- t(OTU_METAPHLAN[,-1])
class(OTU_METAPHLAN) <- 'numeric'
colnames(OTU_METAPHLAN) <- TAXA_METAPHLAN
rownames(OTU_METAPHLAN) <- gsub('_clean.*$','',rownames(OTU_METAPHLAN))
SPECIES_METAPHLAN <- grepl('s__',TAXA_METAPHLAN)
SAMPLES_METAPHLAN <- intersect(rownames(OTU_METAPHLAN),rownames(META))
META_METAPHLAN <- META[SAMPLES_METAPHLAN,]
OTU_METAPHLAN <- OTU_METAPHLAN[SAMPLES_METAPHLAN,]
OTU_METAPHLAN <- OTU_METAPHLAN[rowSums(OTU_METAPHLAN) > 0, colSums(OTU_METAPHLAN) > 0]
TAXA_METAPHLAN <- colnames(OTU_METAPHLAN)
META_METAPHLAN <- META_METAPHLAN[rownames(OTU_METAPHLAN),]

# load megan data
MEG <- read_delim(file.path(work_path,'Data/all399_allnodes_seed_name_2.txt'),'\t') #MEGAN
colnames(MEG) <- str_trim(colnames(MEG))
MEG <- MEG %>%   
  dplyr::rename(Datasets=`#Datasets`) %>%
  filter(!(Datasets %in% c('SEED','No hits','Not assigned'))) %>% 
  select(-null)
KEGG <- MEG[,1]
MEG <- MEG %>% select(-Datasets) %>% dplyr::mutate_all(funs(as.numeric)) %>% t()
colnames(MEG) <- unlist(KEGG)

# convert megan to relative abundances
RAMEG <- MEG/rowSums(MEG)
META_RAMEG <- META[rownames(RAMEG),]

# load humann data
HUM <- as.data.frame(read_delim(file.path(work_path,'Data/complete_ambulance_humann2_output.renormcpm.merged.nostrat.tsv'),'\t')) #HUMAN
rownames(HUM) <- HUM[,1]
HUM <- HUM[,-1]
colnames(HUM) <- gsub('_.*','',colnames(HUM))
HUM <- t(HUM)
META_HUM <- META[rownames(HUM),]

HUMPW <- readr::read_csv(file.path(work_path,'Data/humann_pw.csv'))
HUMPW <- HUMPW[,-c(2,4)]
HUMPW <- HUMPW[!is.na(HUMPW[,1]),]
HUMPW <- as.data.frame(apply(HUMPW,2,stringr::str_trim,side='both'),stringsAsFactors=FALSE)
colnames(HUMPW) <- c('class','superclass1','superclass2','superclass3')
HUMPW <- rbind(HUMPW,data.frame(class=colnames(HUM)[!(colnames(HUM) %in% HUMPW$class)],superclass1=NA,superclass2=NA,superclass3=NA))
rownames(HUMPW) <- HUMPW$class

# convert humann to relative abundances
RAHUM <- HUM/rowSums(HUM)
SAMPELS_HUM_OVERLAP <- intersect(rownames(OTU_OVERLAP),rownames(RAHUM))
OTU_HUM_OVERLAP <- cbind(OTU_OVERLAP[SAMPELS_HUM_OVERLAP,],RAHUM[SAMPELS_HUM_OVERLAP,])
OTU_HUM_OVERLAP <- OTU_HUM_OVERLAP[,colSums(OTU_HUM_OVERLAP)>0]
OTU_HUM_OVERLAP <- OTU_HUM_OVERLAP[rowSums(OTU_HUM_OVERLAP)>0,]
META_HUM_OVERLAP <- META[rownames(OTU_HUM_OVERLAP),]

# tsne ordination plots for exploration

# center and scale ovrerlap features
zOTU <- apply(OTU_METAPHLAN,2,function(x) znormalize(x))
tsne <- Rtsne(zOTU,3,theta=0)
df <- data.frame(tsne$Y,META_METAPHLAN[rownames(zOTU),]) %>% 
  filter(!is.na(sample.Yield)) %>%
  group_by(sample.Surface) %>%
  filter(n() >= 5)
p1 <-  ggplot() +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  geom_point(data=df,aes(x=X1,y=X2,colour=FR)) +
  scale_colour_brewer(type='qual',palette=3) +
  facet_wrap(~sample.Surface,nrow=1) +
  theme(aspect.ratio=1,legend.position='none') +
  labs(x='Axis 1',y='Axis 2') 
p2 <-  ggplot() +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  geom_point(data=df,aes(x=X1,y=X3,colour=FR)) +
  scale_colour_brewer(type='qual',palette=3) +
  facet_wrap(~sample.Surface,nrow=1) +
  theme(aspect.ratio=1,legend.position='none') +
  labs(x='Axis 1',y='Axis 3') 
p3 <-  ggplot() +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  geom_point(data=df,aes(x=X2,y=X3,colour=FR)) +
  scale_colour_brewer(type='qual',palette=3) +
  facet_wrap(~sample.Surface,nrow=1) +
  theme(aspect.ratio=1,legend.position='none') +
  labs(x='Axis 2',y='Axis 3') 

zOTU <- apply(OTU_OVERLAP,2,function(x) znormalize(x))
tsne <- Rtsne(zOTU,3,theta=0)
df <- data.frame(tsne$Y,META_OVERLAP[rownames(zOTU),]) %>% 
  filter(!is.na(sample.Yield)) %>%
  group_by(sample.Surface) %>%
  filter(n() >= 5)
p4 <-  ggplot() +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  geom_point(data=df,aes(x=X1,y=X2,colour=FR)) +
  scale_colour_brewer(type='qual',palette=3) +
  facet_wrap(~sample.Surface,nrow=1) +
  theme(aspect.ratio=1,legend.position='none') +
  labs(x='Axis 1',y='Axis 2') 
p5 <-  ggplot() +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  geom_point(data=df,aes(x=X1,y=X3,colour=FR)) +
  scale_colour_brewer(type='qual',palette=3) +
  facet_wrap(~sample.Surface,nrow=1) +
  theme(aspect.ratio=1,legend.position='none') +
  labs(x='Axis 1',y='Axis 3') 
p6 <-  ggplot() +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  geom_point(data=df,aes(x=X2,y=X3,colour=FR)) +
  scale_colour_brewer(type='qual',palette=3) +
  facet_wrap(~sample.Surface,nrow=1) +
  theme(aspect.ratio=1,legend.position='none') +
  labs(x='Axis 2',y='Axis 3') 

pdf(file.path(out_path,'metaphlan_vs_overlap_tsne.pdf'),height=20,width=20)
grid.arrange(p1,p2,p3,p4,p5,p6,ncol=1,top='Metaphlan vs. Overlap')
dev.off()


# convert overlap to presence/absence then tsne
zOTU <- t(apply(OTU_METAPHLAN,1,function(x) ifelse(x>0,1,0)))
zOTU <- unique(zOTU)
tsne <- Rtsne(zOTU,3,theta=0)
df <- data.frame(tsne$Y,META_METAPHLAN[rownames(zOTU),]) %>% 
  filter(!is.na(sample.Yield)) %>%
  group_by(sample.Surface) %>%
  filter(n() >= 5)
p1 <-  ggplot() +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  geom_point(data=df,aes(x=X1,y=X2,colour=FR)) +
  scale_colour_brewer(type='qual',palette=3) +
  facet_wrap(~sample.Surface,nrow=1) +
  theme(aspect.ratio=1,legend.position='none') +
  labs(x='Axis 1',y='Axis 2') 
p2 <-  ggplot() +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  geom_point(data=df,aes(x=X1,y=X3,colour=FR)) +
  scale_colour_brewer(type='qual',palette=3) +
  facet_wrap(~sample.Surface,nrow=1) +
  theme(aspect.ratio=1,legend.position='none') +
  labs(x='Axis 1',y='Axis 3') 
p3 <-  ggplot() +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  geom_point(data=df,aes(x=X2,y=X3,colour=FR)) +
  scale_colour_brewer(type='qual',palette=3) +
  facet_wrap(~sample.Surface,nrow=1) +
  theme(aspect.ratio=1,legend.position='none') +
  labs(x='Axis 2',y='Axis 3') 

zOTU <- t(apply(OTU_OVERLAP,1,function(x) ifelse(x>0,1,0)))
zOTU <- unique(zOTU)
tsne <- Rtsne(zOTU,3,theta=0)
df <- data.frame(tsne$Y,META_OVERLAP[rownames(zOTU),]) %>% 
  filter(!is.na(sample.Yield)) %>%
  group_by(sample.Surface) %>%
  filter(n() >= 5)
p4 <-  ggplot() +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  geom_point(data=df,aes(x=X1,y=X2,colour=FR)) +
  scale_colour_brewer(type='qual',palette=3) +
  facet_wrap(~sample.Surface,nrow=1) +
  theme(aspect.ratio=1,legend.position='none') +
  labs(x='Axis 1',y='Axis 2') 
p5 <-  ggplot() +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  geom_point(data=df,aes(x=X1,y=X3,colour=FR)) +
  scale_colour_brewer(type='qual',palette=3) +
  facet_wrap(~sample.Surface,nrow=1) +
  theme(aspect.ratio=1,legend.position='none') +
  labs(x='Axis 1',y='Axis 3') 
p6 <-  ggplot() +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  geom_point(data=df,aes(x=X2,y=X3,colour=FR)) +
  scale_colour_brewer(type='qual',palette=3) +
  facet_wrap(~sample.Surface,nrow=1) +
  theme(aspect.ratio=1,legend.position='none') +
  labs(x='Axis 2',y='Axis 3') 

pdf(file.path(out_path,'metaphlan_vs_overlap_pres_aps_tsne.pdf'),height=20,width=20)
grid.arrange(p1,p2,p3,p4,p5,p6,ncol=1,top='Metaphlan vs. Overlap')
dev.off()


# megan
# center and scale
zOTU <- apply(MEG,2,function(x) znormalize(x))
tsne <- Rtsne(zOTU,3,theta=0)
df <- data.frame(tsne$Y,META_RAMEG[rownames(zOTU),],
                 Coverage=log10(rowSums(MEG)+1),
                 Shannon=diversity(MEG,index='shannon',MARGIN=1),
                 Prevelance=rowSums(MEG>0)) %>% 
  filter(!is.na(sample.Yield)) %>%
  group_by(sample.Surface) 
pc1 <- ggplot() +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  geom_point(data=df,aes(x=X1,y=X2,colour=Coverage),size=4,alpha=.7) +
  scale_colour_distiller(palette = 'Spectral') +
  theme(aspect.ratio=1) +
  labs(x='Axis 1',y='Axis 2') +
  ggtitle('Raw Megan')
ps1 <- ggplot() +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  geom_point(data=df,aes(x=X1,y=X2,colour=Shannon),size=4,alpha=.7) +
  scale_colour_distiller(palette = 'Spectral') +
  theme(aspect.ratio=1) +
  labs(x='Axis 1',y='Axis 2') +
  ggtitle('Raw Megan')
pp1 <- ggplot() +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  geom_point(data=df,aes(x=X1,y=X2,colour=Prevelance),size=4,alpha=.7) +
  scale_colour_distiller(palette = 'Spectral') +
  theme(aspect.ratio=1) +
  labs(x='Axis 1',y='Axis 2') +
  ggtitle('Raw Megan')
pss1 <- ggplot() +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  geom_point(data=df,aes(x=X1,y=X2,colour=sample.Surface),size=4,alpha=.7) +
  theme(aspect.ratio=1) +
  labs(x='Axis 1',y='Axis 2') +
  ggtitle('Raw Megan')
psr1 <- ggplot() +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  geom_point(data=df,aes(x=X1,y=X2,colour=region),size=4,alpha=.7) +
  theme(aspect.ratio=1) +
  labs(x='Axis 1',y='Axis 2') +
  ggtitle('Raw Megan')
psy1 <- ggplot() +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  geom_point(data=df,aes(x=X1,y=X2,colour=log10(sample.Yield)),size=4,alpha=.7) +
  scale_colour_distiller(palette = 'Spectral') +
  theme(aspect.ratio=1) +
  labs(x='Axis 1',y='Axis 2',colour='Biomass') +
  ggtitle('Raw Megan')

# center and scale megan relative abundances
zOTU <- apply(RAMEG,2,function(x) znormalize(x))
tsne <- Rtsne(zOTU,3,theta=0)
df <- data.frame(tsne$Y,META_RAMEG[rownames(zOTU),],
                 Coverage=log10(rowSums(MEG)+1),
                 Shannon=diversity(MEG,index='shannon',MARGIN=1),
                 Prevelance=rowSums(MEG>0)) %>% 
  filter(!is.na(sample.Yield)) %>%
  group_by(sample.Surface) 
pc2 <- ggplot() +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  geom_point(data=df,aes(x=X1,y=X2,colour=Coverage),size=4,alpha=.7) +
  scale_colour_distiller(palette = 'Spectral') +
  theme(aspect.ratio=1) +
  labs(x='Axis 1',y='Axis 2') +
  ggtitle('RA Megan')
ps2 <- ggplot() +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  geom_point(data=df,aes(x=X1,y=X2,colour=Shannon),size=4,alpha=.7) +
  scale_colour_distiller(palette = 'Spectral') +
  theme(aspect.ratio=1) +
  labs(x='Axis 1',y='Axis 2') +
  ggtitle('RA Megan')
pp2 <- ggplot() +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  geom_point(data=df,aes(x=X1,y=X2,colour=Prevelance),size=4,alpha=.7) +
  scale_colour_distiller(palette = 'Spectral') +
  theme(aspect.ratio=1) +
  labs(x='Axis 1',y='Axis 2') +
  ggtitle('RA Megan')
pss2 <- ggplot() +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  geom_point(data=df,aes(x=X1,y=X2,colour=sample.Surface),size=4,alpha=.7) +
  theme(aspect.ratio=1) +
  labs(x='Axis 1',y='Axis 2') +
  ggtitle('RA Megan')
psr2 <- ggplot() +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  geom_point(data=df,aes(x=X1,y=X2,colour=region),size=4,alpha=.7) +
  theme(aspect.ratio=1) +
  labs(x='Axis 1',y='Axis 2') +
  ggtitle('RA Megan')
psy2 <- ggplot() +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  geom_point(data=df,aes(x=X1,y=X2,colour=log10(sample.Yield)),size=4,alpha=.7) +
  scale_colour_distiller(palette = 'Spectral') +
  theme(aspect.ratio=1) +
  labs(x='Axis 1',y='Axis 2',colour='Biomass') +
  ggtitle('RA Megan')

# center and scale humann
zOTU <- apply(HUM,2,function(x) znormalize(x))
tsne <- Rtsne(zOTU,3,theta=0)
df <- data.frame(tsne$Y,META_RAMEG[rownames(zOTU),],
                 Coverage=log10(rowSums(HUM)+1),
                 Shannon=diversity(HUM,index='shannon',MARGIN=1),
                 Prevelance=rowSums(HUM>0)) %>% 
  filter(!is.na(sample.Yield)) %>%
  group_by(sample.Surface) 
pc3 <- ggplot() +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  geom_point(data=df,aes(x=X1,y=X2,colour=Coverage),size=4,alpha=.7) +
  scale_colour_distiller(palette = 'Spectral') +
  theme(aspect.ratio=1) +
  labs(x='Axis 1',y='Axis 2') +
  ggtitle('Raw Humann')
ps3 <- ggplot() +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  geom_point(data=df,aes(x=X1,y=X2,colour=Shannon),size=4,alpha=.7) +
  scale_colour_distiller(palette = 'Spectral') +
  theme(aspect.ratio=1) +
  labs(x='Axis 1',y='Axis 2') +
  ggtitle('Raw Humann')
pp3 <- ggplot() +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  geom_point(data=df,aes(x=X1,y=X2,colour=Prevelance),size=4,alpha=.7) +
  scale_colour_distiller(palette = 'Spectral') +
  theme(aspect.ratio=1) +
  labs(x='Axis 1',y='Axis 2') +
  ggtitle('Raw Humann')
pss3 <- ggplot() +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  geom_point(data=df,aes(x=X1,y=X2,colour=sample.Surface),size=4,alpha=.7) +
  theme(aspect.ratio=1) +
  labs(x='Axis 1',y='Axis 2') +
  ggtitle('Raw Humann')
psr3 <- ggplot() +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  geom_point(data=df,aes(x=X1,y=X2,colour=region),size=4,alpha=.7) +
  theme(aspect.ratio=1) +
  labs(x='Axis 1',y='Axis 2') +
  ggtitle('Raw Humann')
psy3 <- ggplot() +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  geom_point(data=df,aes(x=X1,y=X2,colour=log10(sample.Yield)),size=4,alpha=.7) +
  scale_colour_distiller(palette = 'Spectral') +
  theme(aspect.ratio=1) +
  labs(x='Axis 1',y='Axis 2',colour='Biomass') +
  ggtitle('Raw Humann')
zOTU <- apply(RAHUM,2,function(x) znormalize(x))
tsne <- Rtsne(zOTU,3,theta=0)
df <- data.frame(tsne$Y,META_RAMEG[rownames(zOTU),],
                 Coverage=log10(rowSums(HUM)+1),
                 Shannon=diversity(HUM,index='shannon',MARGIN=1),
                 Prevelance=rowSums(HUM>0)) %>% 
  filter(!is.na(sample.Yield)) %>%
  group_by(sample.Surface) 
pc4 <- ggplot() +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  geom_point(data=df,aes(x=X1,y=X2,colour=Coverage),size=4,alpha=.7) +
  scale_colour_distiller(palette = 'Spectral') +
  theme(aspect.ratio=1) +
  labs(x='Axis 1',y='Axis 2') +
  ggtitle('RA Humann')
ps4 <- ggplot() +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  geom_point(data=df,aes(x=X1,y=X2,colour=Shannon),size=4,alpha=.7) +
  scale_colour_distiller(palette = 'Spectral') +
  theme(aspect.ratio=1) +
  labs(x='Axis 1',y='Axis 2') +
  ggtitle('RA Humann')
pp4 <- ggplot() +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  geom_point(data=df,aes(x=X1,y=X2,colour=Prevelance),size=4,alpha=.7) +
  scale_colour_distiller(palette = 'Spectral') +
  theme(aspect.ratio=1) +
  labs(x='Axis 1',y='Axis 2') +
  ggtitle('RA Humann')
pss4 <- ggplot() +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  geom_point(data=df,aes(x=X1,y=X2,colour=sample.Surface),size=4,alpha=.7) +
  theme(aspect.ratio=1) +
  labs(x='Axis 1',y='Axis 2') +
  ggtitle('RA Humann')
psr4 <- ggplot() +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  geom_point(data=df,aes(x=X1,y=X2,colour=region),size=4,alpha=.7) +
  theme(aspect.ratio=1) +
  labs(x='Axis 1',y='Axis 2') +
  ggtitle('RA Humann')
psy4 <- ggplot() +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  geom_point(data=df,aes(x=X1,y=X2,colour=log10(sample.Yield)),size=4,alpha=.7) +
  scale_colour_distiller(palette = 'Spectral') +
  theme(aspect.ratio=1) +
  labs(x='Axis 1',y='Axis 2',colour='Biomass') +
  ggtitle('RA Humann')

pdf(file.path(out_path,'megan_vs_human_coverage_tsne.pdf'),height=25,width=15)
grid.arrange(pc1,pc2,ps1,ps2,pp1,pp2,pss1,pss2,psr1,psr2,psy1,psy2,ncol=2)
grid.arrange(pc3,pc4,ps3,ps4,pp3,pp4,pss3,pss4,psr3,psr4,psy3,psy4,ncol=2)
dev.off()



# make abundance prevelance plots to assess what samples/features to filter

# megan
pm1 <- qplot(x=colSums(MEG),y=colSums(MEG>0)/nrow(MEG),geom='point',alpha=.7) +
  guides(alpha='none',size='none') +
  scale_x_log10() + labs(x='Abundance',y='Prevelance') + ggtitle('Megan Pathways')
pm2 <- qplot(x=rowSums(MEG),y=rowSums(MEG>0)/ncol(MEG),geom='point',colour=META_RAMEG[rownames(MEG),'sample.Surface'],size=4,alpha=.7) +
  guides(alpha='none',size='none') +
  scale_x_log10() + labs(x='Abundance',y='Prevalence',colour='Surface') + ggtitle('Megan Samples')
pm3 <- qplot(x=rowSums(MEG),y=rowSums(MEG>0)/ncol(MEG),geom='point',colour=META_RAMEG[rownames(MEG),'region'],size=4,alpha=.7) +
  guides(alpha='none',size='none') +
  scale_x_log10() + labs(x='Abundance',y='Prevalence',colour='Region') + ggtitle('Megan Samples')
pm4 <- qplot(x=rowSums(MEG),y=rowSums(MEG>0)/ncol(MEG),geom='point',colour=META_RAMEG[rownames(MEG),'FR'],size=4,alpha=.7) +
  guides(alpha='none',size='none') +
  scale_x_log10() + labs(x='Abundance',y='Prevalence',colour='Front/Rear') + ggtitle('Megan Samples')
pm5 <- qplot(x=rowSums(MEG),y=rowSums(MEG>0)/ncol(MEG),geom='point',colour=log10(META_RAMEG[rownames(MEG),'sample.Yield']+1),size=4,alpha=.7) +
  scale_colour_distiller(palette = 'Spectral') +
  guides(alpha='none',size='none') +
  scale_x_log10() + labs(x='Abundance',y='Prevalence',colour='Biomass') + ggtitle('Megan Samples') +
  theme_dark()
pm6 <- qplot(x=log10(META_RAMEG[rownames(MEG),'sample.Yield']+1),
             y=rowSums(MEG>0)/ncol(MEG),geom='point',
             colour=log10(META_HUM[rownames(MEG),'sample.Yield']+1),size=4,alpha=.9) +
  scale_colour_distiller(palette = 'Spectral') +
  guides(alpha='none',size='none') +
  scale_x_log10() + labs(x='Biomass',y='Prevalence',colour='Biomass') + ggtitle('Megan Samples') +
  theme_dark()

# humann
ph1 <- qplot(x=colSums(HUM),y=colSums(HUM>0)/nrow(HUM),geom='point',alpha=.7) + 
  guides(alpha='none',size='none') +
  scale_x_log10() + labs(x='Abundance',y='Prevelance') + ggtitle('Humann Pathways')
ph2 <- qplot(x=rowSums(HUM),y=rowSums(HUM>0)/ncol(HUM),geom='point',colour=META_HUM[rownames(HUM),'sample.Surface'],size=4,alpha=.7) +
  guides(alpha='none',size='none') +
  scale_x_log10() + labs(x='Abundance',y='Prevalence',colour='Surface') + ggtitle('Humann Samples')
ph3 <- qplot(x=rowSums(HUM),y=rowSums(HUM>0)/ncol(HUM),geom='point',colour=META_HUM[rownames(HUM),'region'],size=4,alpha=.7) +
  guides(alpha='none',size='none') +
  scale_x_log10() + labs(x='Abundance',y='Prevalence',colour='Region') + ggtitle('Humann Samples')
ph4 <- qplot(x=rowSums(HUM),y=rowSums(HUM>0)/ncol(HUM),geom='point',colour=META_HUM[rownames(HUM),'FR'],size=4,alpha=.7) +
  guides(alpha='none',size='none') +
  scale_x_log10() + labs(x='Abundance',y='Prevalence',colour='Front/Rear') + ggtitle('Humann Samples')
ph5 <- qplot(x=rowSums(HUM),y=rowSums(HUM>0)/ncol(HUM),geom='point',colour=log10(META_HUM[rownames(HUM),'sample.Yield']+1),size=4,alpha=.7) +
  scale_colour_distiller(palette = 'Spectral') +
  guides(alpha='none',size='none') +
  scale_x_log10() + labs(x='Abundance',y='Prevalence',colour='Biomass') + ggtitle('Humann Samples') +
  theme_dark()
ph6 <- qplot(x=log10(META_HUM[rownames(HUM),'sample.Yield']+1),
             y=rowSums(HUM>0)/ncol(HUM),geom='point',
             colour=log10(META_HUM[rownames(HUM),'sample.Yield']+1),size=4,alpha=.9) +
  scale_colour_distiller(palette = 'Spectral') +
  guides(alpha='none',size='none') +
  scale_x_log10() + labs(x='Biomass',y='Prevalence',colour='Biomass') + ggtitle('Humann Samples') +
  theme_dark()

pmp1 <- qplot(x=colSums(OTU_METAPHLAN),y=colSums(OTU_METAPHLAN>0)/nrow(OTU_METAPHLAN),geom='point',alpha=.7) + 
  guides(alpha='none',size='none') +
  scale_x_log10() + labs(x='Abundance',y='Prevelance') + ggtitle('Metaphlan Taxa')
pmp2 <- qplot(x=rowSums(OTU_METAPHLAN),y=rowSums(OTU_METAPHLAN>0)/ncol(OTU_METAPHLAN),geom='point',colour=META_METAPHLAN[rownames(OTU_METAPHLAN),'sample.Surface'],size=4,alpha=.7) +
  guides(alpha='none',size='none') +
  scale_x_log10() + labs(x='Abundance',y='Prevalence',colour='Surface') + ggtitle('Metaphlan Samples')
pmp3 <- qplot(x=rowSums(OTU_METAPHLAN),y=rowSums(OTU_METAPHLAN>0)/ncol(OTU_METAPHLAN),geom='point',colour=META_METAPHLAN[rownames(OTU_METAPHLAN),'region'],size=4,alpha=.7) +
  guides(alpha='none',size='none') +
  scale_x_log10() + labs(x='Abundance',y='Prevalence',colour='Region') + ggtitle('Metaphlan Samples')
pmp4 <- qplot(x=rowSums(OTU_METAPHLAN),y=rowSums(OTU_METAPHLAN>0)/ncol(OTU_METAPHLAN),geom='point',colour=META_METAPHLAN[rownames(OTU_METAPHLAN),'FR'],size=4,alpha=.7) +
  guides(alpha='none',size='none') +
  scale_x_log10() + labs(x='Abundance',y='Prevalence',colour='Front/Rear') + ggtitle('Metaphlan Samples')
pmp5 <- qplot(x=rowSums(OTU_METAPHLAN),y=rowSums(OTU_METAPHLAN>0)/ncol(OTU_METAPHLAN),geom='point',colour=log10(META_METAPHLAN[rownames(OTU_METAPHLAN),'sample.Yield']+1),size=4,alpha=.7) +
  scale_colour_distiller(palette = 'Spectral') +
  guides(alpha='none',size='none') +
  scale_x_log10() + labs(x='Abundance',y='Prevalence',colour='Biomass') + ggtitle('Metaphlan Samples') +
  theme_dark()
pmp6 <- qplot(x=log10(META_METAPHLAN[rownames(OTU_METAPHLAN),'sample.Yield']+1),
              y=rowSums(OTU_METAPHLAN>0)/ncol(MEG),geom='point',
              colour=log10(META_METAPHLAN[rownames(OTU_METAPHLAN),'sample.Yield']+1),size=4,alpha=.9) +
  scale_colour_distiller(palette = 'Spectral') +
  guides(alpha='none',size='none') +
  scale_x_log10() + labs(x='Biomass',y='Prevalence',colour='Biomass') + ggtitle('Metaphlan Samples') +
  theme_dark()

# ensure metaphlan and humann plots show same samples
SAMPLES_MHM <- intersect(intersect(rownames(OTU_METAPHLAN),rownames(HUM)),rownames(MEG))

X_MH <- (rowSums(OTU_METAPHLAN>0)/ncol(OTU_METAPHLAN))[SAMPLES_MHM]
Y_MH <- (rowSums(MEG>0)/ncol(MEG))[SAMPLES_MHM]
pmm1 <- qplot(x=X_MH,y=Y_MH,geom='point',
              colour=META_METAPHLAN[SAMPLES_MHM,'sample.Surface'] ,size=4,alpha=.7) +
  guides(alpha='none',size='none') +
  scale_x_log10() + labs(x='Metaphlan Prevalence',y='Megan Prevalence',colour='Surface') + ggtitle('Samples')
pmm2 <- qplot(x=X_MH,y=Y_MH,geom='point',
              colour=META_METAPHLAN[SAMPLES_MHM,'region'] ,size=4,alpha=.7) +
  guides(alpha='none',size='none') +
  scale_x_log10() + labs(x='Metaphlan Prevalence',y='Megan Prevalence',colour='Region') + ggtitle('Samples')
pmm3 <- qplot(x=X_MH,y=Y_MH,geom='point',
              colour=META_METAPHLAN[SAMPLES_MHM,'FR'] ,size=4,alpha=.7) +
  guides(alpha='none',size='none') +
  scale_x_log10() + labs(x='Metaphlan Prevalence',y='Megan Prevalence',colour='Front/Rear') + ggtitle('Samples')
pmm4 <- qplot(x=X_MH,y=Y_MH,geom='point',
              colour=log10(META_METAPHLAN[SAMPLES_MHM,'sample.Yield']+1),size=4,alpha=.7) +
  scale_colour_distiller(palette = 'Spectral') +
  guides(alpha='none',size='none') +
  labs(x='Metaphlan Prevalence',y='Megan Prevalence',colour='Biomass') + ggtitle('Samples') +
  theme_dark()

X_MM <- (rowSums(OTU_METAPHLAN>0)/ncol(OTU_METAPHLAN))[SAMPLES_MHM]
Y_MM <- (rowSums(HUM>0)/ncol(HUM))[SAMPLES_MHM]
pmh1 <- qplot(x=X_MM,y=Y_MM,geom='point',
              colour=META_METAPHLAN[SAMPLES_MHM,'sample.Surface'] ,size=4,alpha=.7) +
  guides(alpha='none',size='none') +
  scale_x_log10() + labs(x='Metaphlan Prevalence',y='Humann Prevalence',colour='Surface') + ggtitle('Samples')
pmh2 <- qplot(x=X_MM,y=Y_MM,geom='point',
              colour=META_METAPHLAN[SAMPLES_MHM,'region'] ,size=4,alpha=.7) +
  guides(alpha='none',size='none') +
  scale_x_log10() + labs(x='Metaphlan Prevalence',y='Humann Prevalence',colour='Region') + ggtitle('Samples')
pmh3 <- qplot(x=X_MM,y=Y_MM,geom='point',
              colour=META_METAPHLAN[SAMPLES_MHM,'FR'] ,size=4,alpha=.7) +
  guides(alpha='none',size='none') +
  scale_x_log10() + labs(x='Metaphlan Prevalence',y='Humann Prevalence',colour='Front/Rear') + ggtitle('Samples')
pmh4 <- qplot(x=X_MM,y=Y_MM,geom='point',
              colour=log10(META_METAPHLAN[SAMPLES_MHM,'sample.Yield']+1),size=4,alpha=.7) +
  scale_colour_distiller(palette = 'Spectral') +
  guides(alpha='none',size='none') +
  labs(x='Metaphlan Prevalence',y='Humann Prevalence',colour='Biomass') + ggtitle('Samples') +
  theme_dark()


pdf(file.path(out_path,'megan_vs_human_checks.pdf'),height=15,width=25)
grid.arrange(pm1,ph1,pmp1,ncol=3)
grid.arrange(pm2,ph2,pmp2,
             pm3,ph3,pmp3,
             pm4,ph4,pmp4,
             pm5,ph5,pmp5,
             pm6,ph6,pmp6,ncol=3)
grid.arrange(pmm1,pmh1,
             pmm2,pmh2,
             pmm3,pmh3,
             pmm4,pmh4,ncol=2)
dev.off()

# tsne ordination plots for centered/scaled megan relative abundances

zOTU <- apply(RAMEG,2,function(x) znormalize(x))
tsne <- Rtsne(zOTU,3,theta=0)
df <- data.frame(tsne$Y,META_METAPHLAN[rownames(zOTU),]) %>% 
  filter(!is.na(sample.Yield)) %>%
  group_by(sample.Surface) %>%
  filter(n() >= 5)
p1 <-  ggplot() +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  geom_point(data=df,aes(x=X1,y=X2,colour=FR)) +
  scale_colour_brewer(type='qual',palette=3) +
  facet_wrap(~sample.Surface,nrow=1) +
  theme(aspect.ratio=1,legend.position='none') +
  labs(x='Axis 1',y='Axis 2') 
p2 <-  ggplot() +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  geom_point(data=df,aes(x=X1,y=X3,colour=FR)) +
  scale_colour_brewer(type='qual',palette=3) +
  facet_wrap(~sample.Surface,nrow=1) +
  theme(aspect.ratio=1,legend.position='none') +
  labs(x='Axis 1',y='Axis 3') 
p3 <-  ggplot() +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  geom_point(data=df,aes(x=X2,y=X3,colour=FR)) +
  scale_colour_brewer(type='qual',palette=3) +
  facet_wrap(~sample.Surface,nrow=1) +
  theme(aspect.ratio=1,legend.position='none') +
  labs(x='Axis 2',y='Axis 3') 

zOTU <- apply(RAHUM,2,function(x) znormalize(x))
tsne <- Rtsne(zOTU,3,theta=0)
df <- data.frame(tsne$Y,META_OVERLAP[rownames(zOTU),]) %>% 
  filter(!is.na(sample.Yield)) %>%
  group_by(sample.Surface) %>%
  filter(n() >= 5)
p4 <-  ggplot() +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  geom_point(data=df,aes(x=X1,y=X2,colour=FR)) +
  scale_colour_brewer(type='qual',palette=3) +
  facet_wrap(~sample.Surface,nrow=1) +
  theme(aspect.ratio=1,legend.position='none') +
  labs(x='Axis 1',y='Axis 2') 
p5 <-  ggplot() +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  geom_point(data=df,aes(x=X1,y=X3,colour=FR)) +
  scale_colour_brewer(type='qual',palette=3) +
  facet_wrap(~sample.Surface,nrow=1) +
  theme(aspect.ratio=1,legend.position='none') +
  labs(x='Axis 1',y='Axis 3') 
p6 <-  ggplot() +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  geom_point(data=df,aes(x=X2,y=X3,colour=FR)) +
  scale_colour_brewer(type='qual',palette=3) +
  facet_wrap(~sample.Surface,nrow=1) +
  theme(aspect.ratio=1,legend.position='none') +
  labs(x='Axis 2',y='Axis 3') 

pdf(file.path(out_path,'megan_vs_humann_tsne.pdf'),height=20,width=20)
grid.arrange(p1,p2,p3,p4,p5,p6,ncol=1,top='Megan (RA) vs. Humann (RA)')
dev.off()

# classification workflow
# metaphlan
# multiclass: surface
# assess classifiers on training data repeated CV

# 60 cores to parallelize
ncores <- 60
cl <- makeCluster(ncores,type='SOCK')
registerDoSNOW(cl)

# metaphlan data
OTU <- OTU_METAPHLAN
META <- META_METAPHLAN

# filter surfaces with N<=20
sMETA <- META[rownames(OTU),] %>%
  dplyr::select(id=SampleID,surface=sample.Surface,region=region,FR=FR,IE=IE) %>%
  group_by(surface) %>%
  filter(n() > 20) 

sOTU <- OTU[sMETA$id,]
FILT <- (colSums(sOTU>0) < 3) # so no empties in test set
sOTU <- sOTU[,!FILT]
sOTU <- sOTU[rowSums(sOTU) > 0,]
sMETA <- sMETA %>% filter(id %in% rownames(sOTU))
Y <- as.factor(sMETA$surface)

# quantile normalized features
qOTU <- apply(sOTU,2,qnormalize)
# center/scale features
zOTU <- apply(sOTU,2,znormalize)

Xz <- zOTU[sMETA$id,]
Xq <- qOTU[sMETA$id,]

names(Y) <- sMETA$id

set.seed(78) 

# make 80/20 training/testing datasets
p <- .8
train_idx <- createDataPartition(Y,times=1,p=p)$Resample1
train_names <- names(Y)[train_idx]
train_Y <- Y[names(Y) %in% train_names]
test_Y <- Y[!(names(Y) %in% train_names)]

X <- Xz #Xz
train_X <- X[rownames(X) %in% train_names,]
test_X <- X[!(rownames(X) %in% train_names),]

# test classifiers on training set

# 10-fold cv, repeated 10 times with down sampling to overcome class imbalance
nmin <- table(train_Y) 
tr_ctrl <- trainControl(method='repeatedCV', number=10, repeats=10, 
                        classProbs=TRUE,
                        summaryFunction=multiClassSummary,
                        sampling='down',
                        allowParallel=TRUE)

# rf: sweep number of sampled features at each split
# 128 trees, maximize mean roc of multiclass classification
param_sweep <- expand.grid(mtry=floor(seq(ncol(train_X)^.25,ncol(train_X)^.75,length=5)))
out_rf <- train(x=train_X,
                y=train_Y,
                method='rf',
                ntree=128,
                tuneGrid=param_sweep,
                metric='Mean_ROC', 
                maximize=TRUE, 
                trControl=tr_ctrl)

# en: sweep en and sparsity parameters
param_sweep <- expand.grid(alpha=c(.1,.5,.9),lambda=c(.0005,.005,.05))
out_en <- train(x=train_X,
                y=train_Y,
                method='glmnet',
                tuneGrid=param_sweep,
                maxit=1000000,
                metric='Mean_ROC', 
                maximize=TRUE, #FALSE,
                trControl=tr_ctrl)

# regularized rf
param_sweep <- expand.grid(mtry=floor(seq(ncol(train_X)^.25,ncol(train_X)^.75,length=3)),coefReg=c(.1,.5,.9))
out_rrf <- train(x=train_X,
                 y=train_Y,
                 method='RRFglobal',
                 ntree=128,
                 tuneGrid=param_sweep,
                 metric='Mean_ROC', 
                 maximize=TRUE, #FALSE,
                 trControl=tr_ctrl)

# linear svm
out_svmlinear <- train(x=train_X,
                       y=train_Y,
                       method='svmLinear',
                       metric='Mean_ROC', 
                       maximize=TRUE, #FALSE,
                       trControl=tr_ctrl)

# rbf svm
out_svmrbf <- train(x=train_X,
                    y=train_Y,
                    method='svmRadial',
                    metric='Mean_ROC', 
                    maximize=TRUE, #FALSE,
                    trControl=tr_ctrl)

# gradient boosting
out_gbm <- train(x=train_X,
                 y=train_Y,
                 method='gbm',
                 metric='Mean_ROC', 
                 maximize=TRUE, #FALSE,
                 trControl=tr_ctrl)

# partial least squares
out_pls <- train(x=train_X,
                 y=train_Y,
                 method='pls',
                 metric='Mean_ROC', 
                 maximize=TRUE, #FALSE,
                 trControl=tr_ctrl)

# k nearest neighbors
out_knn <- train(x=train_X,
                 y=train_Y,
                 method='kknn',
                 metric='Mean_ROC', 
                 maximize=TRUE, #FALSE,
                 trControl=tr_ctrl)

# c5.0 decision tree
out_c50 <- train(x=train_X,
                 y=train_Y,
                 method='C5.0Tree',
                 metric='Mean_ROC', 
                 maximize=TRUE, #FALSE,
                 trControl=tr_ctrl)

stopCluster(cl)


# assess model performance on training/validation set CV
models <- list(rf=out_rf,rrf=out_rrf,en=out_en,svmrbf=out_svmrbf,svmlinear=out_svmlinear,
               gbm=out_gbm,pls=out_pls,knn=out_knn,c50=out_c50)
resample_models <- resamples(models)
summary(resample_models,metric=c('Kappa','Mean_Balanced_Accuracy'))

saveRDS(list(models=models,train_X=train_X,train_Y=train_Y,test_X=test_X,test_Y=test_Y),file.path(out_path,'ml_metaphlan.rds'))

pdf(file.path(out_path,'ml_metaphlan.pdf'),height=6,width=10)
bwplot(resample_models,metric=c('Kappa','Mean_Balanced_Accuracy'))
dev.off()

rm(list=ls()[grepl('out_',ls()) & ls() != 'out_path'])

# classification workflow
# overlap
# multiclass: surface
# assess classifiers on training data repeated CV

ncores <- 60
cl <- makeCluster(ncores,type='SOCK')
registerDoSNOW(cl)

OTU <- OTU_OVERLAP
META <- META_OVERLAP

sMETA <- META[rownames(OTU),] %>%
  dplyr::select(id=SampleID,surface=sample.Surface,region=region,FR=FR,IE=IE) %>%
  group_by(surface) %>%
  filter(n() > 20) 

sOTU <- OTU[sMETA$id,]
FILT <- (colSums(sOTU>0) < 3) # so no empties in test set
sOTU <- sOTU[,!FILT]
sOTU <- sOTU[rowSums(sOTU) > 0,]
sMETA <- sMETA %>% filter(id %in% rownames(sOTU))
Y <- as.factor(sMETA$surface)

qOTU <- apply(sOTU,2,qnormalize)
zOTU <- apply(sOTU,2,znormalize)

Xz <- zOTU[sMETA$id,]
Xq <- qOTU[sMETA$id,]

names(Y) <- sMETA$id

set.seed(78) 

p <- .8
train_idx <- createDataPartition(Y,times=1,p=p)$Resample1
train_names <- names(Y)[train_idx]
train_Y <- Y[names(Y) %in% train_names]
test_Y <- Y[!(names(Y) %in% train_names)]

X <- Xz #Xz
train_X <- X[rownames(X) %in% train_names,]
test_X <- X[!(rownames(X) %in% train_names),]

nmin <- table(train_Y) 
tr_ctrl <- trainControl(method='repeatedCV', number=10, repeats=10, 
                        classProbs=TRUE,
                        summaryFunction=multiClassSummary,
                        sampling='down',
                        allowParallel=TRUE)

param_sweep <- expand.grid(mtry=floor(seq(ncol(train_X)^.25,ncol(train_X)^.75,length=5)))
out_rf <- train(x=train_X,
                y=train_Y,
                method='rf',
                ntree=128,
                tuneGrid=param_sweep,
                metric='Mean_ROC', 
                maximize=TRUE, #FALSE,
                trControl=tr_ctrl)

# save rf seed to reproduce down sampling and CV for future runs
rf_seeds <- out_rf$control$seeds

param_sweep <- expand.grid(alpha=c(.1,.5,.9),lambda=c(.0005,.005,.05))
out_en <- train(x=train_X,
                y=train_Y,
                method='glmnet',
                tuneGrid=param_sweep,
                maxit=1000000,
                metric='Mean_ROC', 
                maximize=TRUE, #FALSE,
                trControl=tr_ctrl)

param_sweep <- expand.grid(mtry=floor(seq(ncol(train_X)^.25,ncol(train_X)^.75,length=3)),coefReg=c(.1,.5,.9))
out_rrf <- train(x=train_X,
                 y=train_Y,
                 method='RRFglobal',
                 ntree=128,
                 tuneGrid=param_sweep,
                 metric='Mean_ROC', 
                 maximize=TRUE, #FALSE,
                 trControl=tr_ctrl)

out_svmlinear <- train(x=train_X,
                       y=train_Y,
                       method='svmLinear',
                       metric='Mean_ROC', 
                       maximize=TRUE, #FALSE,
                       trControl=tr_ctrl)

out_svmrbf <- train(x=train_X,
                    y=train_Y,
                    method='svmRadial',
                    metric='Mean_ROC', 
                    maximize=TRUE, #FALSE,
                    trControl=tr_ctrl)

# polynomial SVM
out_svmpoly <- train(x=train_X,
                     y=train_Y,
                     method='svmPoly',
                     metric='Mean_ROC', 
                     maximize=TRUE, #FALSE,
                     trControl=tr_ctrl)

out_gbm <- train(x=train_X,
                 y=train_Y,
                 method='gbm',
                 metric='Mean_ROC', 
                 maximize=TRUE, #FALSE,
                 trControl=tr_ctrl)

out_pls <- train(x=train_X,
                 y=train_Y,
                 method='pls',
                 metric='Mean_ROC', 
                 maximize=TRUE, #FALSE,
                 trControl=tr_ctrl)

out_knn <- train(x=train_X,
                 y=train_Y,
                 method='kknn',
                 metric='Mean_ROC', 
                 maximize=TRUE, #FALSE,
                 trControl=tr_ctrl)

out_c50 <- train(x=train_X,
                 y=train_Y,
                 method='C5.0Tree',
                 metric='Mean_ROC', 
                 maximize=TRUE, #FALSE,
                 trControl=tr_ctrl)

stopCluster(cl)


models <- list(rf=out_rf,rrf=out_rrf,en=out_en,svmpoly=out_svmpoly,svmrbf=out_svmrbf,svmlinear=out_svmlinear,
               gbm=out_gbm,pls=out_pls,knn=out_knn,c50=out_c50)
resample_models <- resamples(models)
summary(resample_models,metric=c('Kappa','Mean_Balanced_Accuracy'))

saveRDS(list(models=models,train_X=train_X,train_Y=train_Y,test_X=test_X,test_Y=test_Y),file.path(out_path,'ml_overlap.rds'))

pdf(file.path(out_path,'ml_overlap.pdf'),height=6,width=10)
bwplot(resample_models,metric=c('Kappa','Mean_Balanced_Accuracy'))
dev.off()

rm(list=ls()[grepl('out_',ls()) & ls() != 'out_path'])

# classification workflow
# multiclass: surface
# assess datasets on training data repeated CV, using RF

ncores <- 60
cl <- makeCluster(ncores,type='SOCK')
registerDoSNOW(cl)

# overlap
OTU <- OTU_OVERLAP
META <- META_OVERLAP

sMETA <- META[rownames(OTU),] %>%
  dplyr::select(id=SampleID,surface=sample.Surface,region=region,FR=FR,IE=IE) %>%
  group_by(surface) %>%
  filter(n() > 20) 

sOTU <- OTU[sMETA$id,]
FILT <- (colSums(sOTU>0) < 3) # so no empties in test set
sOTU <- sOTU[,!FILT]
sOTU <- sOTU[rowSums(sOTU) > 0,]
sMETA <- sMETA %>% filter(id %in% rownames(sOTU))
Y <- as.factor(sMETA$surface)

qOTU <- apply(sOTU,2,qnormalize)
zOTU <- apply(sOTU,2,znormalize)

Xz <- zOTU[sMETA$id,]
Xq <- qOTU[sMETA$id,]

names(Y) <- sMETA$id

set.seed(78) 

train_Y <- Y[names(Y) %in% train_names]
test_Y <- Y[!(names(Y) %in% train_names)]

X <- Xz
train_X <- X[rownames(X) %in% train_names,]
test_X <- X[!(rownames(X) %in% train_names),]

tr_ctrl <- trainControl(method='repeatedCV', number=10, repeats=10, 
                        classProbs=TRUE,
                        summaryFunction=multiClassSummary,
                        sampling='down',
                        seeds=rf_seeds,
                        allowParallel=TRUE)

param_sweep <- expand.grid(mtry=floor(seq(ncol(train_X)^.25,ncol(train_X)^.75,length=5)))
out_rf_overlap <- train(x=train_X,
                        y=train_Y,
                        method='rf',
                        ntree=128,
                        tuneGrid=param_sweep,
                        metric='Mean_ROC', 
                        maximize=TRUE, #FALSE,
                        trControl=tr_ctrl)

# metaphlan
OTU <- OTU_METAPHLAN
META <- META_METAPHLAN

sMETA <- META[rownames(OTU),] %>%
  dplyr::select(id=SampleID,surface=sample.Surface,region=region,FR=FR,IE=IE) %>%
  group_by(surface) %>%
  filter(n() > 20) 

sOTU <- OTU[sMETA$id,]
FILT <- (colSums(sOTU>0) < 3) # so no empties in test set
sOTU <- sOTU[,!FILT]
sOTU <- sOTU[rowSums(sOTU) > 0,]
sMETA <- sMETA %>% filter(id %in% rownames(sOTU))
Y <- as.factor(sMETA$surface)


qOTU <- apply(sOTU,2,qnormalize)
zOTU <- apply(sOTU,2,znormalize)

Xz <- zOTU[sMETA$id,]
Xq <- qOTU[sMETA$id,]

names(Y) <- sMETA$id

set.seed(78) 

train_Y <- Y[names(Y) %in% train_names]
test_Y <- Y[!(names(Y) %in% train_names)]

X <- Xz
train_X <- X[rownames(X) %in% train_names,]
test_X <- X[!(rownames(X) %in% train_names),]

tr_ctrl <- trainControl(method='repeatedCV', number=10, repeats=10, 
                        classProbs=TRUE,
                        summaryFunction=multiClassSummary,
                        sampling='down',
                        seeds=rf_seeds,
                        allowParallel=TRUE)

param_sweep <- expand.grid(mtry=floor(seq(ncol(train_X)^.25,ncol(train_X)^.75,length=5)))
out_rf_metaphlan <- train(x=train_X,
                          y=train_Y,
                          method='rf',
                          ntree=128,
                          tuneGrid=param_sweep,
                          metric='Mean_ROC', 
                          maximize=TRUE, #FALSE,
                          trControl=tr_ctrl)

# humann
OTU <- HUM
META <- META_HUM

sMETA <- META[rownames(OTU),] %>%
  dplyr::select(id=SampleID,surface=sample.Surface,region=region,FR=FR,IE=IE) %>%
  group_by(surface) %>%
  filter(n() > 20) 

sOTU <- OTU[sMETA$id,]
FILT <- (colSums(sOTU>0) < 3) # so no empties in test set
sOTU <- sOTU[,!FILT]
sOTU <- sOTU[rowSums(sOTU) > 0,]
sMETA <- sMETA %>% filter(id %in% rownames(sOTU))
Y <- as.factor(sMETA$surface)


qOTU <- apply(sOTU,2,qnormalize)
zOTU <- apply(sOTU,2,znormalize)

Xz <- zOTU[sMETA$id,]
Xq <- qOTU[sMETA$id,]

names(Y) <- sMETA$id

set.seed(78) 

train_Y <- Y[names(Y) %in% train_names]
test_Y <- Y[!(names(Y) %in% train_names)]

X <- Xz
train_X <- X[rownames(X) %in% train_names,]
test_X <- X[!(rownames(X) %in% train_names),]

tr_ctrl <- trainControl(method='repeatedCV', number=10, repeats=10, 
                        classProbs=TRUE,
                        summaryFunction=multiClassSummary,
                        sampling='down',
                        seeds=rf_seeds,
                        allowParallel=TRUE)

param_sweep <- expand.grid(mtry=floor(seq(ncol(train_X)^.25,ncol(train_X)^.75,length=5)))
out_rf_humann <- train(x=train_X,
                       y=train_Y,
                       method='rf',
                       ntree=128,
                       tuneGrid=param_sweep,
                       metric='Mean_ROC', 
                       maximize=TRUE, #FALSE,
                       trControl=tr_ctrl)

# megan
OTU <- MEG
META <- META_RAMEG

sMETA <- META[rownames(OTU),] %>%
  dplyr::select(id=SampleID,surface=sample.Surface,region=region,FR=FR,IE=IE) %>%
  group_by(surface) %>%
  filter(n() > 20) 

sOTU <- OTU[sMETA$id,]
FILT <- (colSums(sOTU>0) < 3) # so no empties in test set
sOTU <- sOTU[,!FILT]
sOTU <- sOTU[rowSums(sOTU) > 0,]
sMETA <- sMETA %>% filter(id %in% rownames(sOTU))
Y <- as.factor(sMETA$surface)


qOTU <- apply(sOTU,2,qnormalize)
zOTU <- apply(sOTU,2,znormalize)

Xz <- zOTU[sMETA$id,]
Xq <- qOTU[sMETA$id,]

names(Y) <- sMETA$id

set.seed(78) 

#p <- .8
#train_idx <- createDataPartition(Y,times=1,p=p)$Resample1
#train_names <- names(Y)[train_idx]
train_Y <- Y[names(Y) %in% train_names]
test_Y <- Y[!(names(Y) %in% train_names)]

X <- Xz
train_X <- X[rownames(X) %in% train_names,]
test_X <- X[!(rownames(X) %in% train_names),]

tr_ctrl <- trainControl(method='repeatedCV', number=10, repeats=10, 
                        classProbs=TRUE,
                        summaryFunction=multiClassSummary,
                        sampling='down',
                        seeds=rf_seeds,
                        allowParallel=TRUE)

param_sweep <- expand.grid(mtry=floor(seq(ncol(train_X)^.25,ncol(train_X)^.75,length=5)))
out_rf_megan <- train(x=train_X,
                      y=train_Y,
                      method='rf',
                      ntree=128,
                      tuneGrid=param_sweep,
                      metric='Mean_ROC', 
                      maximize=TRUE, #FALSE,
                      trControl=tr_ctrl)


stopCluster(cl)


# assess performance of datasets using rf
models <- list(rf_overlap=out_rf_overlap,rf_metaphlan=out_rf_metaphlan,rf_megan=out_rf_megan,rf_humann=out_rf_humann)
resample_models <- resamples(models)
summary(resample_models,metric=c('Kappa','Mean_Balanced_Accuracy'))

saveRDS(list(models=models,train_X=train_X,train_Y=train_Y,test_X=test_X,test_Y=test_Y),file.path(out_path,'ml_datasets.rds'))

pdf(file.path(out_path,'ml_datasets.pdf'),height=6,width=10)
bwplot(resample_models,metric=c('Kappa','Mean_Balanced_Accuracy'))
dev.off()

rm(list=ls()[grepl('out_',ls()) & ls() != 'out_path'])


### assess generalization error using overlap dataset and rf
### fit on training/validation data, then assess error on test set

ncores <- 60
cl <- makeCluster(ncores,type='SOCK')
registerDoSNOW(cl)

OTU <- OTU_OVERLAP
META <- META_OVERLAP

sMETA <- META[rownames(OTU),] %>%
  dplyr::select(id=SampleID,surface=sample.Surface,region=region,FR=FR,IE=IE) %>%
  group_by(surface) %>%
  filter(n() > 20) 

sOTU <- OTU[sMETA$id,]
FILT <- (colSums(sOTU>0) < 3) # so no empties in test set
sOTU <- sOTU[,!FILT]
sOTU <- sOTU[rowSums(sOTU) > 0,]
sMETA <- sMETA %>% filter(id %in% rownames(sOTU))
Y <- as.factor(sMETA$surface)

qOTU <- apply(sOTU,2,qnormalize)
zOTU <- apply(sOTU,2,znormalize)

Xz <- zOTU[sMETA$id,]
Xq <- qOTU[sMETA$id,]

names(Y) <- sMETA$id

set.seed(78) 

train_Y <- Y[names(Y) %in% train_names]
test_Y <- Y[!(names(Y) %in% train_names)]

X <- Xz
train_X <- X[rownames(X) %in% train_names,]
test_X <- X[!(rownames(X) %in% train_names),]

tr_ctrl <- trainControl(method='repeatedCV', number=10, repeats=10, 
                        classProbs=TRUE,
                        sampling='down',
                        seeds=rf_seeds,
                        summaryFunction=multiClassSummary,
                        allowParallel=TRUE)

param_sweep <- expand.grid(mtry=floor(seq(ncol(train_X)^.25,ncol(train_X)^.75,length=5)))
out_rf <- train(x=train_X,
                y=train_Y,
                method='rf',
                ntree=128,
                tuneGrid=param_sweep,
                importance=TRUE,
                proximity=TRUE,
                metric='Mean_ROC', 
                maximize=TRUE, #FALSE,
                trControl=tr_ctrl)

stopCluster(cl)

out_rf

pred_p <- predict(out_rf,test_X,type='prob')
pred_Y <- factor(colnames(pred_p)[apply(pred_p,1,which.max)],levels=colnames(pred_p))
pred_out <- data.frame(obs=test_Y,
                       pred=pred_Y,
                       pred_p)
multiClassSummary(pred_out,lev=levels(pred_out$pred))
caret::confusionMatrix(pred_Y,test_Y)


fm <- out_rf$finalModel
pred_p <- fm$votes
pred_Y <- fm$predicted
pred_out <- data.frame(obs=fm$y,
                       pred=pred_Y,
                       pred_p)
multiClassSummary(pred_out,lev=levels(pred_out$pred))
caret::confusionMatrix(pred_Y,fm$y)


### plot most important features from rf on best performing classes
surfaces <- c('Stethoscope','RearLights_controlPanel','RearBench_seats')
top_taxa <- sapply(seq_along(surfaces), function(i) fm$xNames[order(varImp(fm)[,surfaces[i]],decreasing=TRUE)[1]])
names(top_taxa) <- surfaces
p1 <- data.frame(OTU_OVERLAP[,top_taxa],
                 Surface=META_OVERLAP[rownames(OTU_OVERLAP),'sample.Surface']) %>%
  group_by(Surface) %>%
  filter(Surface %in% surfaces) %>%
  ungroup() %>%
  gather(Top,Abundance,-Surface) %>%
  group_by(Surface) %>%
  dplyr::mutate(N=length(Surface)/length(surfaces)) %>%
  ungroup() %>%
  dplyr::mutate(Bin=cut(Abundance,breaks=c(-1,seq(min(OTU_OVERLAP[OTU_OVERLAP!=0])/10,max(OTU_OVERLAP),length=10)),ordered_result=TRUE)) %>%
  group_by(Surface,Top,Bin) %>%
  dplyr::summarise(Samples=length(Bin)/unique(N)) %>%
  ungroup() %>%
  dplyr::mutate(Pair=ifelse(Top == top_taxa[as.character(Surface)],'A','B'))
levels(p1$Bin)[1] <- '0'
levels(p1$Bin)[2] <- gsub('.*,(.*)\\]$','(0,\\1]',levels(p1$Bin)[2])
p1 <- p1 %>%  ggplot() + geom_bar(aes(x=Bin,y=Samples,fill=Pair),colour='black',stat='identity') + facet_grid(Surface ~ Top) +
  theme(legend.position='none', axis.text.x = element_text(angle=45,hjust=1)) + xlab('RPK Bin')

pdf(file.path(out_path,'rf_final_import_abund.pdf'),height=8,width=8)
p1
dev.off()

# make proximity plots
rf_prox <- cmdscale(1 - out_rf$finalModel$proximity)
p1 <- data.frame(rf_prox,Surface=out_rf$finalModel$y) %>%
  ggplot(aes(x = X1, y = X2, colour = Surface)) +
  geom_point(size=4,alpha=.7) +
  facet_wrap(~ Surface) +
  geom_hline(yintercept=0) + geom_vline(xintercept=0) +
  labs(x = "Axis 1", y = "Axis 2")

pdf(file.path(out_path,'rf_final_prox.pdf'),height=10,width=10)
p1
dev.off()

saveRDS(list(model=out_rf,test_X=test_X,test_Y=test_Y),file.path(out_path,'rf_final.rds'))

pdf(file.path(out_path,'rf_final_features.pdf'),height=15,width=10)
plot(varImp(out_rf),top=25)
dev.off()

rm(list=ls()[grepl('out_',ls()) & ls() != 'out_path'])

# classification workflow
# metaphlan
# multiclass: surface
# assess classifiers on training data repeated CV

ncores <- 60
cl <- makeCluster(ncores,type='SOCK')
registerDoSNOW(cl)

OTU <- OTU_METAPHLAN
META <- META_METAPHLAN

sMETA <- META[rownames(OTU),] %>%
  dplyr::select(id=SampleID,surface=sample.Surface,region=region,FR=FR,IE=IE) %>%
  group_by(region) %>%
  filter(n() > 10) 

sOTU <- OTU[sMETA$id,]
FILT <- (colSums(sOTU>0) < 3) # so no empties in test set
sOTU <- sOTU[,!FILT]
sOTU <- sOTU[rowSums(sOTU) > 0,]
sMETA <- sMETA %>% filter(id %in% rownames(sOTU))
Y <- as.factor(sMETA$region)

qOTU <- apply(sOTU,2,qnormalize)
zOTU <- apply(sOTU,2,znormalize)

Xz <- zOTU[sMETA$id,]
Xq <- qOTU[sMETA$id,]

names(Y) <- sMETA$id

set.seed(78) 

p <- .8
train_idx <- createDataPartition(Y,times=1,p=p)$Resample1
train_names <- names(Y)[train_idx]
train_Y <- Y[names(Y) %in% train_names]
test_Y <- Y[!(names(Y) %in% train_names)]

X <- Xz #Xz
train_X <- X[rownames(X) %in% train_names,]
test_X <- X[!(rownames(X) %in% train_names),]

nmin <- table(train_Y) 
tr_ctrl <- trainControl(method='repeatedCV', number=10, repeats=10, 
                        classProbs=TRUE,
                        summaryFunction=multiClassSummary,
                        sampling='up',
                        allowParallel=TRUE)

param_sweep <- expand.grid(mtry=floor(seq(ncol(train_X)^.25,ncol(train_X)^.75,length=5)))
out_rf <- train(x=train_X,
                y=train_Y,
                method='rf',
                ntree=128,
                tuneGrid=param_sweep,
                metric='Mean_ROC', 
                maximize=TRUE, #FALSE,
                trControl=tr_ctrl)

param_sweep <- expand.grid(alpha=c(.1,.5,.9),lambda=c(.0005,.005,.05))
out_en <- train(x=train_X,
                y=train_Y,
                method='glmnet',
                tuneGrid=param_sweep,
                maxit=1000000,
                metric='Mean_ROC', 
                maximize=TRUE, #FALSE,
                trControl=tr_ctrl)

param_sweep <- expand.grid(mtry=floor(seq(ncol(train_X)^.25,ncol(train_X)^.75,length=3)),coefReg=c(.1,.5,.9))
out_rrf <- train(x=train_X,
                 y=train_Y,
                 method='RRFglobal',
                 ntree=128,
                 tuneGrid=param_sweep,
                 metric='Mean_ROC', 
                 maximize=TRUE, #FALSE,
                 trControl=tr_ctrl)

out_svmlinear <- train(x=train_X,
                       y=train_Y,
                       method='svmLinear',
                       metric='Mean_ROC', 
                       maximize=TRUE, #FALSE,
                       trControl=tr_ctrl)

out_svmrbf <- train(x=train_X,
                    y=train_Y,
                    method='svmRadial',
                    metric='Mean_ROC', 
                    maximize=TRUE, #FALSE,
                    trControl=tr_ctrl)

out_gbm <- train(x=train_X,
                 y=train_Y,
                 method='gbm',
                 metric='Mean_ROC', 
                 maximize=TRUE, #FALSE,
                 trControl=tr_ctrl)

out_pls <- train(x=train_X,
                 y=train_Y,
                 method='pls',
                 metric='Mean_ROC', 
                 maximize=TRUE, #FALSE,
                 trControl=tr_ctrl)

out_knn <- train(x=train_X,
                 y=train_Y,
                 method='kknn',
                 metric='Mean_ROC', 
                 maximize=TRUE, #FALSE,
                 trControl=tr_ctrl)

out_c50 <- train(x=train_X,
                 y=train_Y,
                 method='C5.0Tree',
                 metric='Mean_ROC', 
                 maximize=TRUE, #FALSE,
                 trControl=tr_ctrl)

stopCluster(cl)


models <- list(rf=out_rf,rrf=out_rrf,en=out_en,svmrbf=out_svmrbf,svmlinear=out_svmlinear,
               gbm=out_gbm,pls=out_pls,knn=out_knn,c50=out_c50)
resample_models <- resamples(models)
summary(resample_models,metric=c('Kappa','Mean_Balanced_Accuracy'))

saveRDS(list(models=models,train_X=train_X,train_Y=train_Y,test_X=test_X,test_Y=test_Y),file.path(out_path,'ml_region_metaphlan.rds'))

pdf(file.path(out_path,'ml_region_metaphlan.pdf'),height=6,width=10)
bwplot(resample_models,metric=c('Kappa','Mean_Balanced_Accuracy'))
dev.off()

rm(list=ls()[grepl('out_',ls()) & ls() != 'out_path'])

# classification workflow
# overlap
# multiclass: region
# assess classifiers on training data repeated CV

ncores <- 60
cl <- makeCluster(ncores,type='SOCK')
registerDoSNOW(cl)

OTU <- OTU_OVERLAP
META <- META_OVERLAP

sMETA <- META[rownames(OTU),] %>%
  dplyr::select(id=SampleID,surface=sample.Surface,region=region,FR=FR,IE=IE) %>%
  group_by(region) %>%
  filter(n() > 10) 

sOTU <- OTU[sMETA$id,]
FILT <- (colSums(sOTU>0) < 3) # so no empties in test set
sOTU <- sOTU[,!FILT]
sOTU <- sOTU[rowSums(sOTU) > 0,]
sMETA <- sMETA %>% filter(id %in% rownames(sOTU))
Y <- as.factor(sMETA$region)

qOTU <- apply(sOTU,2,qnormalize)
zOTU <- apply(sOTU,2,znormalize)

Xz <- zOTU[sMETA$id,]
Xq <- qOTU[sMETA$id,]

names(Y) <- sMETA$id

set.seed(78) 

p <- .8
train_idx <- createDataPartition(Y,times=1,p=p)$Resample1
train_names <- names(Y)[train_idx]
train_Y <- Y[names(Y) %in% train_names]
test_Y <- Y[!(names(Y) %in% train_names)]

X <- Xz #Xz
train_X <- X[rownames(X) %in% train_names,]
test_X <- X[!(rownames(X) %in% train_names),]

### performing upsampling here due to one class being underrepresented
nmin <- table(train_Y) 
tr_ctrl <- trainControl(method='repeatedCV', number=10, repeats=10, 
                        classProbs=TRUE,
                        summaryFunction=multiClassSummary,
                        sampling='up',
                        allowParallel=TRUE)

param_sweep <- expand.grid(mtry=floor(seq(ncol(train_X)^.25,ncol(train_X)^.75,length=5)))
out_rf <- train(x=train_X,
                y=train_Y,
                method='rf',
                ntree=128,
                tuneGrid=param_sweep,
                metric='Mean_ROC', 
                maximize=TRUE, #FALSE,
                trControl=tr_ctrl)

# preserve seeds for rf for subsequent runs
rf_seeds <- out_rf$control$seeds

param_sweep <- expand.grid(alpha=c(.1,.5,.9),lambda=c(.0005,.005,.05))
out_en <- train(x=train_X,
                y=train_Y,
                method='glmnet',
                tuneGrid=param_sweep,
                maxit=1000000,
                metric='Mean_ROC', 
                maximize=TRUE, #FALSE,
                trControl=tr_ctrl)

param_sweep <- expand.grid(mtry=floor(seq(ncol(train_X)^.25,ncol(train_X)^.75,length=3)),coefReg=c(.1,.5,.9))
out_rrf <- train(x=train_X,
                 y=train_Y,
                 method='RRFglobal',
                 ntree=128,
                 tuneGrid=param_sweep,
                 metric='Mean_ROC', 
                 maximize=TRUE, #FALSE,
                 trControl=tr_ctrl)

out_svmlinear <- train(x=train_X,
                       y=train_Y,
                       method='svmLinear',
                       metric='Mean_ROC', 
                       maximize=TRUE, #FALSE,
                       trControl=tr_ctrl)

out_svmrbf <- train(x=train_X,
                    y=train_Y,
                    method='svmRadial',
                    metric='Mean_ROC', 
                    maximize=TRUE, #FALSE,
                    trControl=tr_ctrl)

out_svmpoly <- train(x=train_X,
                     y=train_Y,
                     method='svmPoly',
                     metric='Mean_ROC', 
                     maximize=TRUE, #FALSE,
                     trControl=tr_ctrl)

out_gbm <- train(x=train_X,
                 y=train_Y,
                 method='gbm',
                 metric='Mean_ROC', 
                 maximize=TRUE, #FALSE,
                 trControl=tr_ctrl)

out_pls <- train(x=train_X,
                 y=train_Y,
                 method='pls',
                 metric='Mean_ROC', 
                 maximize=TRUE, #FALSE,
                 trControl=tr_ctrl)

out_knn <- train(x=train_X,
                 y=train_Y,
                 method='kknn',
                 metric='Mean_ROC', 
                 maximize=TRUE, #FALSE,
                 trControl=tr_ctrl)

out_c50 <- train(x=train_X,
                 y=train_Y,
                 method='C5.0Tree',
                 metric='Mean_ROC', 
                 maximize=TRUE, #FALSE,
                 trControl=tr_ctrl)

stopCluster(cl)

models <- list(rf=out_rf,rrf=out_rrf,en=out_en,svmrbf=out_svmrbf,svmlinear=out_svmlinear,
               gbm=out_gbm,pls=out_pls,knn=out_knn,c50=out_c50)
resample_models <- resamples(models)
summary(resample_models,metric=c('Kappa','Mean_Balanced_Accuracy'))

saveRDS(list(models=models,train_X=train_X,train_Y=train_Y,test_X=test_X,test_Y=test_Y),file.path(out_path,'ml_region_overlap.rds'))

pdf(file.path(out_path,'ml_region_overlap.pdf'),height=6,width=10)
bwplot(resample_models,metric=c('Kappa','Mean_Balanced_Accuracy'))
dev.off()

rm(list=ls()[grepl('out_',ls()) & ls() != 'out_path'])


# classification workflow on dataset performance
# multiclass: region
# assess classifiers on training data repeated CV

ncores <- 60
cl <- makeCluster(ncores,type='SOCK')
registerDoSNOW(cl)

OTU <- OTU_OVERLAP
META <- META_OVERLAP

sMETA <- META[rownames(OTU),] %>%
  dplyr::select(id=SampleID,surface=sample.Surface,region=region,FR=FR,IE=IE) %>%
  group_by(region) %>%
  filter(n() > 10) 

sOTU <- OTU[sMETA$id,]
FILT <- (colSums(sOTU>0) < 3) # so no empties in test set
sOTU <- sOTU[,!FILT]
sOTU <- sOTU[rowSums(sOTU) > 0,]
sMETA <- sMETA %>% filter(id %in% rownames(sOTU))
Y <- as.factor(sMETA$region)

qOTU <- apply(sOTU,2,qnormalize)
zOTU <- apply(sOTU,2,znormalize)

Xz <- zOTU[sMETA$id,]
Xq <- qOTU[sMETA$id,]

names(Y) <- sMETA$id

set.seed(78) 

train_Y <- Y[names(Y) %in% train_names]
test_Y <- Y[!(names(Y) %in% train_names)]

X <- Xz
train_X <- X[rownames(X) %in% train_names,]
test_X <- X[!(rownames(X) %in% train_names),]

tr_ctrl <- trainControl(method='repeatedCV', number=10, repeats=10, 
                        classProbs=TRUE,
                        summaryFunction=multiClassSummary,
                        sampling='up',
                        seeds=rf_seeds,
                        allowParallel=TRUE)

param_sweep <- expand.grid(mtry=floor(seq(ncol(train_X)^.25,ncol(train_X)^.75,length=5)))
out_rf_overlap <- train(x=train_X,
                        y=train_Y,
                        method='rf',
                        ntree=128,
                        tuneGrid=param_sweep,
                        metric='Mean_ROC', 
                        maximize=TRUE, #FALSE,
                        trControl=tr_ctrl)

OTU <- OTU_METAPHLAN
META <- META_METAPHLAN

sMETA <- META[rownames(OTU),] %>%
  dplyr::select(id=SampleID,surface=sample.Surface,region=region,FR=FR,IE=IE) %>%
  group_by(region) %>%
  filter(n() > 10) 

sOTU <- OTU[sMETA$id,]
FILT <- (colSums(sOTU>0) < 3) # so no empties in test set
sOTU <- sOTU[,!FILT]
sOTU <- sOTU[rowSums(sOTU) > 0,]
sMETA <- sMETA %>% filter(id %in% rownames(sOTU))
Y <- as.factor(sMETA$region)


qOTU <- apply(sOTU,2,qnormalize)
zOTU <- apply(sOTU,2,znormalize)

Xz <- zOTU[sMETA$id,]
Xq <- qOTU[sMETA$id,]

names(Y) <- sMETA$id

set.seed(78) 

train_Y <- Y[names(Y) %in% train_names]
test_Y <- Y[!(names(Y) %in% train_names)]

X <- Xz
train_X <- X[rownames(X) %in% train_names,]
test_X <- X[!(rownames(X) %in% train_names),]

tr_ctrl <- trainControl(method='repeatedCV', number=10, repeats=10, 
                        classProbs=TRUE,
                        summaryFunction=multiClassSummary,
                        sampling='up',
                        seeds=rf_seeds,
                        allowParallel=TRUE)

param_sweep <- expand.grid(mtry=floor(seq(ncol(train_X)^.25,ncol(train_X)^.75,length=5)))
out_rf_metaphlan <- train(x=train_X,
                          y=train_Y,
                          method='rf',
                          ntree=128,
                          tuneGrid=param_sweep,
                          metric='Mean_ROC', 
                          maximize=TRUE, #FALSE,
                          trControl=tr_ctrl)

OTU <- HUM
META <- META_HUM

sMETA <- META[rownames(OTU),] %>%
  dplyr::select(id=SampleID,surface=sample.Surface,region=region,FR=FR,IE=IE) %>%
  group_by(region) %>%
  filter(n() > 10) 

sOTU <- OTU[sMETA$id,]
FILT <- (colSums(sOTU>0) < 3) # so no empties in test set
sOTU <- sOTU[,!FILT]
sOTU <- sOTU[rowSums(sOTU) > 0,]
sMETA <- sMETA %>% filter(id %in% rownames(sOTU))
Y <- as.factor(sMETA$region)


qOTU <- apply(sOTU,2,qnormalize)
zOTU <- apply(sOTU,2,znormalize)

Xz <- zOTU[sMETA$id,]
Xq <- qOTU[sMETA$id,]

names(Y) <- sMETA$id

set.seed(78) 

train_Y <- Y[names(Y) %in% train_names]
test_Y <- Y[!(names(Y) %in% train_names)]

X <- Xz
train_X <- X[rownames(X) %in% train_names,]
test_X <- X[!(rownames(X) %in% train_names),]

tr_ctrl <- trainControl(method='repeatedCV', number=10, repeats=10, 
                        classProbs=TRUE,
                        summaryFunction=multiClassSummary,
                        sampling='up',
                        seeds=rf_seeds,
                        allowParallel=TRUE)

param_sweep <- expand.grid(mtry=floor(seq(ncol(train_X)^.25,ncol(train_X)^.75,length=5)))
out_rf_humann <- train(x=train_X,
                       y=train_Y,
                       method='rf',
                       ntree=128,
                       tuneGrid=param_sweep,
                       metric='Mean_ROC', 
                       maximize=TRUE, #FALSE,
                       trControl=tr_ctrl)

OTU <- MEG
META <- META_RAMEG

sMETA <- META[rownames(OTU),] %>%
  dplyr::select(id=SampleID,surface=sample.Surface,region=region,FR=FR,IE=IE) %>%
  group_by(region) %>%
  filter(n() > 10) 

sOTU <- OTU[sMETA$id,]
FILT <- (colSums(sOTU>0) < 3) # so no empties in test set
sOTU <- sOTU[,!FILT]
sOTU <- sOTU[rowSums(sOTU) > 0,]
sMETA <- sMETA %>% filter(id %in% rownames(sOTU))
Y <- as.factor(sMETA$region)


qOTU <- apply(sOTU,2,qnormalize)
zOTU <- apply(sOTU,2,znormalize)

Xz <- zOTU[sMETA$id,]
Xq <- qOTU[sMETA$id,]

names(Y) <- sMETA$id

set.seed(78) 

train_Y <- Y[names(Y) %in% train_names]
test_Y <- Y[!(names(Y) %in% train_names)]

X <- Xz
train_X <- X[rownames(X) %in% train_names,]
test_X <- X[!(rownames(X) %in% train_names),]

tr_ctrl <- trainControl(method='repeatedCV', number=10, repeats=10, 
                        classProbs=TRUE,
                        summaryFunction=multiClassSummary,
                        sampling='up',
                        seeds=rf_seeds,
                        allowParallel=TRUE)

param_sweep <- expand.grid(mtry=floor(seq(ncol(train_X)^.25,ncol(train_X)^.75,length=5)))
out_rf_megan <- train(x=train_X,
                      y=train_Y,
                      method='rf',
                      ntree=128,
                      tuneGrid=param_sweep,
                      metric='Mean_ROC', 
                      maximize=TRUE, #FALSE,
                      trControl=tr_ctrl)


stopCluster(cl)


models <- list(rf_overlap=out_rf_overlap,rf_metaphlan=out_rf_metaphlan,rf_megan=out_rf_megan,rf_humann=out_rf_humann)
resample_models <- resamples(models)
summary(resample_models,metric=c('Kappa','Mean_Balanced_Accuracy'))

saveRDS(list(models=models,train_X=train_X,train_Y=train_Y,test_X=test_X,test_Y=test_Y),file.path(out_path,'ml_region_datasets.rds'))

pdf(file.path(out_path,'ml_region_datasets.pdf'),height=6,width=10)
bwplot(resample_models,metric=c('Kappa','Mean_Balanced_Accuracy'))
dev.off()


rm(list=ls()[grepl('out_',ls()) & ls() != 'out_path'])

# assess generalization error via test set after training/validation cv fit
# multiclass: region

ncores <- 60
cl <- makeCluster(ncores,type='SOCK')
registerDoSNOW(cl)

OTU <- OTU_OVERLAP
META <- META_OVERLAP

sMETA <- META[rownames(OTU),] %>%
  dplyr::select(id=SampleID,surface=sample.Surface,region=region,FR=FR,IE=IE) %>%
  group_by(region) %>%
  filter(n() > 10) 

sOTU <- OTU[sMETA$id,]
FILT <- (colSums(sOTU>0) < 3) # so no empties in test set
sOTU <- sOTU[,!FILT]
sOTU <- sOTU[rowSums(sOTU) > 0,]
sMETA <- sMETA %>% filter(id %in% rownames(sOTU))
Y <- as.factor(sMETA$region)

qOTU <- apply(sOTU,2,qnormalize)
zOTU <- apply(sOTU,2,znormalize)

Xz <- zOTU[sMETA$id,]
Xq <- qOTU[sMETA$id,]

names(Y) <- sMETA$id

set.seed(78) 

train_Y <- Y[names(Y) %in% train_names]
test_Y <- Y[!(names(Y) %in% train_names)]

X <- Xz
train_X <- X[rownames(X) %in% train_names,]
test_X <- X[!(rownames(X) %in% train_names),]

tr_ctrl <- trainControl(method='repeatedCV', number=10, repeats=10, 
                        classProbs=TRUE,
                        sampling='up',
                        seeds=rf_seeds,
                        summaryFunction=multiClassSummary,
                        allowParallel=TRUE)

param_sweep <- expand.grid(mtry=floor(seq(ncol(train_X)^.25,ncol(train_X)^.75,length=5)))
out_rf <- train(x=train_X,
                y=train_Y,
                method='rf',
                ntree=128,
                tuneGrid=param_sweep,
                importance=TRUE,
                proximity=TRUE,
                metric='Mean_ROC', 
                maximize=TRUE, #FALSE,
                trControl=tr_ctrl)

stopCluster(cl)

out_rf

pred_p <- predict(out_rf,test_X,type='prob')
pred_Y <- factor(colnames(pred_p)[apply(pred_p,1,which.max)],levels=colnames(pred_p))
pred_out <- data.frame(obs=test_Y,
                       pred=pred_Y,
                       pred_p)
multiClassSummary(pred_out,lev=levels(pred_out$pred))
caret::confusionMatrix(pred_Y,test_Y)

ROC <- roc(response = test_Y, 
           predictor = pred_p[,1], 
           levels = rev(unique(test_Y)))
plot(ROC, col='red',lwd = 2)


fm <- out_rf$finalModel
pred_p <- fm$votes
pred_Y <- fm$predicted
pred_out <- data.frame(obs=fm$y,
                       pred=pred_Y,
                       pred_p)
multiClassSummary(pred_out,lev=levels(pred_out$pred))
caret::confusionMatrix(pred_Y,fm$y)

# plot features from rf for region for best performing classes
regions <- c('s_e','s_w_w_coast','w','w_coast')
top_taxa <- sapply(seq_along(regions), function(i) fm$xNames[order(varImp(fm)[,regions[i]],decreasing=TRUE)[1]])
names(top_taxa) <- regions
p1 <- data.frame(OTU_OVERLAP[,top_taxa],
                 Region=META_OVERLAP[rownames(OTU_OVERLAP),'region']) %>%
  group_by(Region) %>%
  filter(Region %in% regions) %>%
  ungroup() %>%
  gather(Top,Abundance,-Region) %>%
  group_by(Region) %>%
  dplyr::mutate(N=length(Region)/length(regions)) %>%
  ungroup() %>%
  dplyr::mutate(Bin=cut(Abundance,breaks=c(-1,seq(min(OTU_OVERLAP[OTU_OVERLAP!=0])/10,max(OTU_OVERLAP),length=10)),ordered_result=TRUE)) %>%
  group_by(Region,Top,Bin) %>%
  dplyr::mutate(Samples=length(Bin)/unique(N)) %>%
  ungroup() %>%
  dplyr::mutate(Pair=ifelse(Top == top_taxa[as.character(Region)],'A','B')) %>%
  ggplot() + geom_bar(aes(x=Bin,y=Samples,fill=Pair),colour='black',stat='identity') + facet_grid(Region ~ Top) +
  theme(legend.position='none', axis.text.x = element_text(angle=45,hjust=1)) + xlab('RPK Bin')

pdf(file.path(out_path,'rf_region_final_import_abund.pdf'),height=8,width=8)
p1
dev.off()

# make proximity plots
rf_prox <- cmdscale(1 - out_rf$finalModel$proximity)
p1 <- data.frame(rf_prox,Region=out_rf$finalModel$y) %>%
  ggplot(aes(x = X1, y = X2, colour = Region)) +
  geom_point(size=4,alpha=.7) +
  facet_wrap(~ Region) +
  geom_hline(yintercept=0) + geom_vline(xintercept=0) +
  labs(x = "Axis 1", y = "Axis 2")

pdf(file.path(out_path,'rf_region_final_prox.pdf'),height=10,width=10)
p1
dev.off()

saveRDS(list(model=out_rf,test_X=test_X,test_Y=test_Y),file.path(out_path,'rf_region_final.rds'))

pdf(file.path(out_path,'rf_region_final_features.pdf'),height=15,width=10)
plot(varImp(out_rf),top=25)
dev.off()

rm(list=ls()[grepl('out_',ls()) & ls() != 'out_path'])

# classification workflow
# metaphlan
# multiclass: city
# assess classifiers on training data repeated CV

ncores <- 60
cl <- makeCluster(ncores,type='SOCK')
registerDoSNOW(cl)

OTU <- OTU_METAPHLAN
META <- META_METAPHLAN

# keep cities with N>=10
sMETA <- META[rownames(OTU),] %>%
  dplyr::select(id=SampleID,surface=sample.Surface,region=region,FR=FR,IE=IE,city=city) %>%
  group_by(city) %>%
  filter(n() >= 10) 

sOTU <- OTU[sMETA$id,]
FILT <- (colSums(sOTU>0) < 3) # so no empties in test set
sOTU <- sOTU[,!FILT]
sOTU <- sOTU[rowSums(sOTU) > 0,]
sMETA <- sMETA %>% filter(id %in% rownames(sOTU))
Y <- as.factor(sMETA$city)

qOTU <- apply(sOTU,2,qnormalize)
zOTU <- apply(sOTU,2,znormalize)

Xz <- zOTU[sMETA$id,]
Xq <- qOTU[sMETA$id,]

names(Y) <- sMETA$id

set.seed(78) 

p <- .8
train_idx <- createDataPartition(Y,times=1,p=p)$Resample1
train_names <- names(Y)[train_idx]
train_Y <- Y[names(Y) %in% train_names]
test_Y <- Y[!(names(Y) %in% train_names)]

X <- Xz #Xz
train_X <- X[rownames(X) %in% train_names,]
test_X <- X[!(rownames(X) %in% train_names),]

# upsampling to overcome low represented class
tr_ctrl <- trainControl(method='repeatedCV', number=10, repeats=10, 
                        classProbs=TRUE,
                        summaryFunction=multiClassSummary,
                        sampling='up',
                        allowParallel=TRUE)

param_sweep <- expand.grid(mtry=floor(seq(ncol(train_X)^.25,ncol(train_X)^.75,length=5)))
out_rf <- train(x=train_X,
                y=train_Y,
                method='rf',
                ntree=128,
                tuneGrid=param_sweep,
                metric='Mean_ROC', 
                maximize=TRUE, #FALSE,
                trControl=tr_ctrl)

param_sweep <- expand.grid(alpha=c(.1,.5,.9),lambda=c(.0005,.005,.05))
out_en <- train(x=train_X,
                y=train_Y,
                method='glmnet',
                tuneGrid=param_sweep,
                maxit=1000000,
                metric='Mean_ROC', 
                maximize=TRUE, #FALSE,
                trControl=tr_ctrl)

param_sweep <- expand.grid(mtry=floor(seq(ncol(train_X)^.25,ncol(train_X)^.75,length=3)),coefReg=c(.1,.5,.9))
out_rrf <- train(x=train_X,
                 y=train_Y,
                 method='RRFglobal',
                 ntree=128,
                 tuneGrid=param_sweep,
                 metric='Mean_ROC', 
                 maximize=TRUE, #FALSE,
                 trControl=tr_ctrl)

out_svmlinear <- train(x=train_X,
                       y=train_Y,
                       method='svmLinear',
                       metric='Mean_ROC', 
                       maximize=TRUE, #FALSE,
                       trControl=tr_ctrl)

out_svmrbf <- train(x=train_X,
                    y=train_Y,
                    method='svmRadial',
                    metric='Mean_ROC', 
                    maximize=TRUE, #FALSE,
                    trControl=tr_ctrl)

out_gbm <- train(x=train_X,
                 y=train_Y,
                 method='gbm',
                 metric='Mean_ROC', 
                 maximize=TRUE, #FALSE,
                 trControl=tr_ctrl)

out_pls <- train(x=train_X,
                 y=train_Y,
                 method='pls',
                 metric='Mean_ROC', 
                 maximize=TRUE, #FALSE,
                 trControl=tr_ctrl)

out_knn <- train(x=train_X,
                 y=train_Y,
                 method='kknn',
                 metric='Mean_ROC', 
                 maximize=TRUE, #FALSE,
                 trControl=tr_ctrl)

out_c50 <- train(x=train_X,
                 y=train_Y,
                 method='C5.0Tree',
                 metric='Mean_ROC', 
                 maximize=TRUE, #FALSE,
                 trControl=tr_ctrl)

stopCluster(cl)


models <- list(rf=out_rf,rrf=out_rrf,en=out_en,svmrbf=out_svmrbf,svmlinear=out_svmlinear,
               gbm=out_gbm,pls=out_pls,knn=out_knn,c50=out_c50)
resample_models <- resamples(models)
summary(resample_models,metric=c('Kappa','Mean_Balanced_Accuracy'))

saveRDS(list(models=models,train_X=train_X,train_Y=train_Y,test_X=test_X,test_Y=test_Y),file.path(out_path,'ml_city_metaphlan.rds'))

pdf(file.path(out_path,'ml_city_metaphlan.pdf'),height=6,width=10)
bwplot(resample_models,metric=c('Kappa','Mean_Balanced_Accuracy'))
dev.off()

rm(list=ls()[grepl('out_',ls()) & ls() != 'out_path'])

# classification workflow
# overlap
# multiclass: city
# assess classifiers on training data repeated CV

ncores <- 60
cl <- makeCluster(ncores,type='SOCK')
registerDoSNOW(cl)

OTU <- OTU_OVERLAP
META <- META_OVERLAP

sMETA <- META[rownames(OTU),] %>%
  dplyr::select(id=SampleID,surface=sample.Surface,region=region,FR=FR,IE=IE,city=city) %>%
  group_by(city) %>%
  filter(n() >= 10) 

sOTU <- OTU[sMETA$id,]
FILT <- (colSums(sOTU>0) < 3) # so no empties in test set
sOTU <- sOTU[,!FILT]
sOTU <- sOTU[rowSums(sOTU) > 0,]
sMETA <- sMETA %>% filter(id %in% rownames(sOTU))
Y <- as.factor(sMETA$city)

qOTU <- apply(sOTU,2,qnormalize)
zOTU <- apply(sOTU,2,znormalize)

Xz <- zOTU[sMETA$id,]
Xq <- qOTU[sMETA$id,]

names(Y) <- sMETA$id

set.seed(78) 

p <- .8
train_idx <- createDataPartition(Y,times=1,p=p)$Resample1
train_names <- names(Y)[train_idx]
train_Y <- Y[names(Y) %in% train_names]
test_Y <- Y[!(names(Y) %in% train_names)]

X <- Xz #Xz
train_X <- X[rownames(X) %in% train_names,]
test_X <- X[!(rownames(X) %in% train_names),]

tr_ctrl <- trainControl(method='repeatedCV', number=10, repeats=10, 
                        classProbs=TRUE,
                        summaryFunction=multiClassSummary,
                        sampling='up',
                        allowParallel=TRUE)

param_sweep <- expand.grid(mtry=floor(seq(ncol(train_X)^.25,ncol(train_X)^.75,length=5)))
out_rf <- train(x=train_X,
                y=train_Y,
                method='rf',
                ntree=128,
                tuneGrid=param_sweep,
                metric='Mean_ROC', 
                maximize=TRUE, #FALSE,
                trControl=tr_ctrl)

#tr_ctrl$seeds <- out_rf$control$seeds
rf_seeds <- out_rf$control$seeds

param_sweep <- expand.grid(alpha=c(.1,.5,.9),lambda=c(.0005,.005,.05))
out_en <- train(x=train_X,
                y=train_Y,
                method='glmnet',
                tuneGrid=param_sweep,
                maxit=1000000,
                metric='Mean_ROC', 
                maximize=TRUE, #FALSE,
                trControl=tr_ctrl)

param_sweep <- expand.grid(mtry=floor(seq(ncol(train_X)^.25,ncol(train_X)^.75,length=3)),coefReg=c(.1,.5,.9))
out_rrf <- train(x=train_X,
                 y=train_Y,
                 method='RRFglobal',
                 ntree=128,
                 tuneGrid=param_sweep,
                 metric='Mean_ROC', 
                 maximize=TRUE, #FALSE,
                 trControl=tr_ctrl)

out_svmlinear <- train(x=train_X,
                       y=train_Y,
                       method='svmLinear',
                       metric='Mean_ROC', 
                       maximize=TRUE, #FALSE,
                       trControl=tr_ctrl)

out_svmrbf <- train(x=train_X,
                    y=train_Y,
                    method='svmRadial',
                    metric='Mean_ROC', 
                    maximize=TRUE, #FALSE,
                    trControl=tr_ctrl)

out_gbm <- train(x=train_X,
                 y=train_Y,
                 method='gbm',
                 metric='Mean_ROC', 
                 maximize=TRUE, #FALSE,
                 trControl=tr_ctrl)

out_pls <- train(x=train_X,
                 y=train_Y,
                 method='pls',
                 metric='Mean_ROC', 
                 maximize=TRUE, #FALSE,
                 trControl=tr_ctrl)

out_knn <- train(x=train_X,
                 y=train_Y,
                 method='kknn',
                 metric='Mean_ROC', 
                 maximize=TRUE, #FALSE,
                 trControl=tr_ctrl)

out_c50 <- train(x=train_X,
                 y=train_Y,
                 method='C5.0Tree',
                 metric='Mean_ROC', 
                 maximize=TRUE, #FALSE,
                 trControl=tr_ctrl)

stopCluster(cl)


models <- list(rf=out_rf,rrf=out_rrf,en=out_en,svmrbf=out_svmrbf,svmlinear=out_svmlinear,
               gbm=out_gbm,pls=out_pls,knn=out_knn,c50=out_c50)
resample_models <- resamples(models)
summary(resample_models,metric=c('Kappa','Mean_Balanced_Accuracy'))

saveRDS(list(models=models,train_X=train_X,train_Y=train_Y,test_X=test_X,test_Y=test_Y),file.path(out_path,'ml_city_overlap.rds'))

pdf(file.path(out_path,'ml_city_overlap.pdf'),height=6,width=10)
bwplot(resample_models,metric=c('Kappa','Mean_Balanced_Accuracy'))
dev.off()


rm(list=ls()[grepl('out_',ls()) & ls() != 'out_path'])

# classification workflow for dataset performance
# multiclass: city
# assess classifiers on training data repeated CV

ncores <- 60
cl <- makeCluster(ncores,type='SOCK')
registerDoSNOW(cl)

OTU <- OTU_OVERLAP
META <- META_OVERLAP

sMETA <- META[rownames(OTU),] %>%
  dplyr::select(id=SampleID,surface=sample.Surface,region=region,FR=FR,IE=IE,city=city) %>%
  group_by(city) %>%
  filter(n() >= 10) 

sOTU <- OTU[sMETA$id,]
FILT <- (colSums(sOTU>0) < 3) # so no empties in test set
sOTU <- sOTU[,!FILT]
sOTU <- sOTU[rowSums(sOTU) > 0,]
sMETA <- sMETA %>% filter(id %in% rownames(sOTU))
Y <- as.factor(sMETA$city)

qOTU <- apply(sOTU,2,qnormalize)
zOTU <- apply(sOTU,2,znormalize)

Xz <- zOTU[sMETA$id,]
Xq <- qOTU[sMETA$id,]

names(Y) <- sMETA$id

set.seed(78) 

p <- .8
train_idx <- createDataPartition(Y,times=1,p=p)$Resample1
train_names <- names(Y)[train_idx]
train_Y <- Y[names(Y) %in% train_names]
test_Y <- Y[!(names(Y) %in% train_names)]

X <- Xz #Xz
train_X <- X[rownames(X) %in% train_names,]
test_X <- X[!(rownames(X) %in% train_names),]

tr_ctrl <- trainControl(method='repeatedCV', number=10, repeats=10, 
                        classProbs=TRUE,
                        seeds=rf_seeds,
                        summaryFunction=multiClassSummary,
                        sampling='up',
                        allowParallel=TRUE)

param_sweep <- expand.grid(mtry=floor(seq(ncol(train_X)^.25,ncol(train_X)^.75,length=5)))
out_rf_overlap <- train(x=train_X,
                        y=train_Y,
                        method='rf',
                        ntree=128,
                        tuneGrid=param_sweep,
                        metric='Mean_ROC', 
                        maximize=TRUE, #FALSE,
                        trControl=tr_ctrl)

OTU <- OTU_METAPHLAN
META <- META_METAPHLAN

sMETA <- META[rownames(OTU),] %>%
  dplyr::select(id=SampleID,surface=sample.Surface,region=region,FR=FR,IE=IE,city=city) %>%
  group_by(city) %>%
  filter(n() >= 10) 

sOTU <- OTU[sMETA$id,]
FILT <- (colSums(sOTU>0) < 3) # so no empties in test set
sOTU <- sOTU[,!FILT]
sOTU <- sOTU[rowSums(sOTU) > 0,]
sMETA <- sMETA %>% filter(id %in% rownames(sOTU))
Y <- as.factor(sMETA$city)

qOTU <- apply(sOTU,2,qnormalize)
zOTU <- apply(sOTU,2,znormalize)

Xz <- zOTU[sMETA$id,]
Xq <- qOTU[sMETA$id,]

names(Y) <- sMETA$id

set.seed(78) 

p <- .8
train_idx <- createDataPartition(Y,times=1,p=p)$Resample1
train_names <- names(Y)[train_idx]
train_Y <- Y[names(Y) %in% train_names]
test_Y <- Y[!(names(Y) %in% train_names)]

X <- Xz #Xz
train_X <- X[rownames(X) %in% train_names,]
test_X <- X[!(rownames(X) %in% train_names),]

tr_ctrl <- trainControl(method='repeatedCV', number=10, repeats=10, 
                        classProbs=TRUE,
                        seeds=rf_seeds,
                        summaryFunction=multiClassSummary,
                        sampling='up',
                        allowParallel=TRUE)

param_sweep <- expand.grid(mtry=floor(seq(ncol(train_X)^.25,ncol(train_X)^.75,length=5)))
out_rf_metaphlan <- train(x=train_X,
                          y=train_Y,
                          method='rf',
                          ntree=128,
                          tuneGrid=param_sweep,
                          metric='Mean_ROC', 
                          maximize=TRUE, #FALSE,
                          trControl=tr_ctrl)


OTU <- HUM
META <- META_HUM

sMETA <- META[rownames(OTU),] %>%
  dplyr::select(id=SampleID,surface=sample.Surface,region=region,FR=FR,IE=IE,city=city) %>%
  group_by(city) %>%
  filter(n() >= 10) 

sOTU <- OTU[sMETA$id,]
FILT <- (colSums(sOTU>0) < 3) # so no empties in test set
sOTU <- sOTU[,!FILT]
sOTU <- sOTU[rowSums(sOTU) > 0,]
sMETA <- sMETA %>% filter(id %in% rownames(sOTU))
Y <- as.factor(sMETA$city)

qOTU <- apply(sOTU,2,qnormalize)
zOTU <- apply(sOTU,2,znormalize)

Xz <- zOTU[sMETA$id,]
Xq <- qOTU[sMETA$id,]

names(Y) <- sMETA$id

set.seed(78) 

p <- .8
train_idx <- createDataPartition(Y,times=1,p=p)$Resample1
train_names <- names(Y)[train_idx]
train_Y <- Y[names(Y) %in% train_names]
test_Y <- Y[!(names(Y) %in% train_names)]

X <- Xz #Xz
train_X <- X[rownames(X) %in% train_names,]
test_X <- X[!(rownames(X) %in% train_names),]

tr_ctrl <- trainControl(method='repeatedCV', number=10, repeats=10, 
                        classProbs=TRUE,
                        seeds=rf_seeds,
                        summaryFunction=multiClassSummary,
                        sampling='up',
                        allowParallel=TRUE)

param_sweep <- expand.grid(mtry=floor(seq(ncol(train_X)^.25,ncol(train_X)^.75,length=5)))
out_rf_humann <- train(x=train_X,
                       y=train_Y,
                       method='rf',
                       ntree=128,
                       tuneGrid=param_sweep,
                       metric='Mean_ROC', 
                       maximize=TRUE, #FALSE,
                       trControl=tr_ctrl)


OTU <- MEG
META <- META_RAMEG

sMETA <- META[rownames(OTU),] %>%
  dplyr::select(id=SampleID,surface=sample.Surface,region=region,FR=FR,IE=IE,city=city) %>%
  group_by(city) %>%
  filter(n() >= 10) 

sOTU <- OTU[sMETA$id,]
FILT <- (colSums(sOTU>0) < 3) # so no empties in test set
sOTU <- sOTU[,!FILT]
sOTU <- sOTU[rowSums(sOTU) > 0,]
sMETA <- sMETA %>% filter(id %in% rownames(sOTU))
Y <- as.factor(sMETA$city)

qOTU <- apply(sOTU,2,qnormalize)
zOTU <- apply(sOTU,2,znormalize)

Xz <- zOTU[sMETA$id,]
Xq <- qOTU[sMETA$id,]

names(Y) <- sMETA$id

set.seed(78) 

p <- .8
train_idx <- createDataPartition(Y,times=1,p=p)$Resample1
train_names <- names(Y)[train_idx]
train_Y <- Y[names(Y) %in% train_names]
test_Y <- Y[!(names(Y) %in% train_names)]

X <- Xz #Xz
train_X <- X[rownames(X) %in% train_names,]
test_X <- X[!(rownames(X) %in% train_names),]

tr_ctrl <- trainControl(method='repeatedCV', number=10, repeats=10, 
                        classProbs=TRUE,
                        seeds=rf_seeds,
                        summaryFunction=multiClassSummary,
                        sampling='up',
                        allowParallel=TRUE)

param_sweep <- expand.grid(mtry=floor(seq(ncol(train_X)^.25,ncol(train_X)^.75,length=5)))
out_rf_megan <- train(x=train_X,
                      y=train_Y,
                      method='rf',
                      ntree=128,
                      tuneGrid=param_sweep,
                      metric='Mean_ROC', 
                      maximize=TRUE, #FALSE,
                      trControl=tr_ctrl)


stopCluster(cl)


models <- list(rf_overlap=out_rf_overlap,rf_metaphlan=out_rf_metaphlan,rf_megan=out_rf_megan,rf_humann=out_rf_humann)
resample_models <- resamples(models)
summary(resample_models,metric=c('Kappa','Mean_Balanced_Accuracy'))

saveRDS(list(models=models,train_X=train_X,train_Y=train_Y,test_X=test_X,test_Y=test_Y),file.path(out_path,'ml_city_datasets.rds'))

pdf(file.path(out_path,'ml_city_datasets.pdf'),height=6,width=10)
bwplot(resample_models,metric=c('Kappa','Mean_Balanced_Accuracy'))
dev.off()


rm(list=ls()[grepl('out_',ls()) & ls() != 'out_path'])

# generalization error fitting training/validation then using test set
# overlap
# multiclass: city
# assess classifiers on training data repeated CV

ncores <- 60
cl <- makeCluster(ncores,type='SOCK')
registerDoSNOW(cl)

OTU <- OTU_OVERLAP
META <- META_OVERLAP

sMETA <- META[rownames(OTU),] %>%
  dplyr::select(id=SampleID,surface=sample.Surface,region=region,FR=FR,IE=IE,city=city) %>%
  group_by(city) %>%
  filter(n() >= 10) 

sOTU <- OTU[sMETA$id,]
FILT <- (colSums(sOTU>0) < 3) # so no empties in test set
sOTU <- sOTU[,!FILT]
sOTU <- sOTU[rowSums(sOTU) > 0,]
sMETA <- sMETA %>% filter(id %in% rownames(sOTU))
Y <- as.factor(sMETA$city)

qOTU <- apply(sOTU,2,qnormalize)
zOTU <- apply(sOTU,2,znormalize)

Xz <- zOTU[sMETA$id,]
Xq <- qOTU[sMETA$id,]

names(Y) <- sMETA$id

set.seed(78) 

p <- .8
train_idx <- createDataPartition(Y,times=1,p=p)$Resample1
train_names <- names(Y)[train_idx]
train_Y <- Y[names(Y) %in% train_names]
test_Y <- Y[!(names(Y) %in% train_names)]

X <- Xz #Xz
train_X <- X[rownames(X) %in% train_names,]
test_X <- X[!(rownames(X) %in% train_names),]

tr_ctrl <- trainControl(method='repeatedCV', number=10, repeats=10, 
                        #method='LOOCV',
                        classProbs=TRUE,
                        seeds=rf_seeds,
                        summaryFunction=multiClassSummary,
                        sampling='up',
                        allowParallel=TRUE)

param_sweep <- expand.grid(mtry=floor(seq(ncol(train_X)^.25,ncol(train_X)^.75,length=5)))
out_rf <- train(x=train_X,
                y=train_Y,
                method='rf',
                ntree=128,
                tuneGrid=param_sweep,
                importance=TRUE,
                proximity=TRUE,
                metric='Mean_ROC', 
                maximize=TRUE, #FALSE,
                trControl=tr_ctrl)

stopCluster(cl)

out_rf

pred_p <- predict(out_rf,test_X,type='prob')
pred_Y <- factor(colnames(pred_p)[apply(pred_p,1,which.max)],levels=colnames(pred_p))
pred_out <- data.frame(obs=test_Y,
                       pred=pred_Y,
                       pred_p)
multiClassSummary(pred_out,lev=levels(pred_out$pred))
caret::confusionMatrix(pred_Y,test_Y)

fm <- out_rf$finalModel
pred_p <- fm$votes
pred_Y <- fm$predicted
pred_out <- data.frame(obs=fm$y,
                       pred=pred_Y,
                       pred_p)
multiClassSummary(pred_out,lev=levels(pred_out$pred))
caret::confusionMatrix(pred_Y,fm$y)

# assess RF importance of features
cities <- as.character(unique(train_Y))
top_taxa <- sapply(seq_along(cities), function(i) fm$xNames[order(varImp(fm)[,cities[i]],decreasing=TRUE)[1]])
names(top_taxa) <- cities
p1 <- data.frame(OTU_OVERLAP[,top_taxa],
                 City=META_OVERLAP[rownames(OTU_OVERLAP),'city']) %>%
  group_by(City) %>%
  filter(City %in% cities) %>%
  ungroup() %>%
  gather(Top,Abundance,-City) %>%
  group_by(City) %>%
  dplyr::mutate(N=length(City)/length(cities)) %>%
  ungroup() %>%
  dplyr::mutate(Bin=cut(Abundance,breaks=c(-1,seq(min(OTU_OVERLAP[OTU_OVERLAP!=0])/10,max(OTU_OVERLAP),length=10)),ordered_result=TRUE)) %>%
  group_by(City,Top,Bin) %>%
  dplyr::summarise(Samples=length(Bin)/unique(N)) %>%
  ungroup() %>%
  dplyr::mutate(Pair=ifelse(Top == top_taxa[as.character(City)],'A','B')) %>%
  ggplot() + geom_bar(aes(x=Bin,y=Samples,fill=Pair),colour='black',stat='identity') + facet_grid(City ~ Top) +
  theme(legend.position='none', axis.text.x = element_text(angle=45,hjust=1)) + xlab('RPK Bin')

pdf(file.path(out_path,'rf_city_final_import_abund.pdf'),height=15,width=15)
p1
dev.off()

# make proximity plots
rf_prox <- cmdscale(1 - out_rf$finalModel$proximity)
p1 <- data.frame(rf_prox,City=out_rf$finalModel$y) %>%
  ggplot(aes(x = X1, y = X2, colour = City)) +
  geom_point(size=4,alpha=.7) +
  facet_wrap(~ City) +
  geom_hline(yintercept=0) + geom_vline(xintercept=0) +
  labs(x = "Axis 1", y = "Axis 2")

pdf(file.path(out_path,'rf_city_final_prox.pdf'),height=10,width=10)
p1
dev.off()

saveRDS(list(model=out_rf,test_X=test_X,test_Y=test_Y),file.path(out_path,'rf_city_final.rds'))

pdf(file.path(out_path,'rf_city_final_features.pdf'),height=15,width=10)
plot(varImp(out_rf),top=25)
dev.off()


rm(list=ls()[grepl('out_',ls()) & ls() != 'out_path'])

# classification: train/val/test
# overlap
# 2-class: front/rear
# assess classifiers on training data repeated CV

ncores <- 60
cl <- makeCluster(ncores,type='SOCK')
registerDoSNOW(cl)

OTU <- OTU_OVERLAP
META <- META_OVERLAP

sMETA <- META[rownames(OTU),] %>%
  dplyr::select(id=SampleID,surface=sample.Surface,region=region,FR=FR,IE=IE) 

sOTU <- OTU[sMETA$id,]
FILT <- (colSums(sOTU>0) < 3) # so no empties in test set
sOTU <- sOTU[,!FILT]
sOTU <- sOTU[rowSums(sOTU) > 0,]
sMETA <- sMETA %>% filter(id %in% rownames(sOTU))
Y <- as.factor(sMETA$FR)

qOTU <- apply(sOTU,2,qnormalize)
zOTU <- apply(sOTU,2,znormalize)

Xz <- zOTU[sMETA$id,]
Xq <- qOTU[sMETA$id,]

names(Y) <- sMETA$id

set.seed(78) 

p <- .8
train_idx <- createDataPartition(Y,times=1,p=p)$Resample1
train_names <- names(Y)[train_idx]
train_Y <- Y[names(Y) %in% train_names]
test_Y <- Y[!(names(Y) %in% train_names)]

X <- Xz
train_X <- X[rownames(X) %in% train_names,]
test_X <- X[!(rownames(X) %in% train_names),]

tr_ctrl <- trainControl(method='repeatedCV', number=10, repeats=10, 
                        classProbs=TRUE,
                        sampling='down',
                        seeds=rf_seeds,
                        summaryFunction=twoClassSummary,
                        allowParallel=TRUE)

param_sweep <- expand.grid(mtry=floor(seq(ncol(train_X)^.25,ncol(train_X)^.75,length=5)))
out_rf <- train(x=train_X,
                y=train_Y,
                method='rf',
                ntree=128,
                tuneGrid=param_sweep,
                importance=TRUE,
                proximity=TRUE,
                metric='ROC', 
                maximize=TRUE, #FALSE,
                trControl=tr_ctrl)

stopCluster(cl)



out_rf

pred_p <- predict(out_rf,test_X,type='prob')
pred_Y <- factor(colnames(pred_p)[apply(pred_p,1,which.max)],levels=colnames(pred_p))
pred_out <- data.frame(obs=test_Y,
                       pred=pred_Y,
                       pred_p)
twoClassSummary(pred_out,lev=levels(pred_out$pred))
caret::confusionMatrix(pred_Y,test_Y)

ROC <- roc(response = test_Y, 
           predictor = pred_p[,1], 
           levels = rev(unique(test_Y)))
plot(ROC, col='red',lwd = 2)


fm <- out_rf$finalModel
pred_p <- fm$votes
pred_Y <- fm$predicted
pred_out <- data.frame(obs=fm$y,
                       pred=pred_Y,
                       pred_p)
twoClassSummary(pred_out,lev=levels(pred_out$pred))
caret::confusionMatrix(pred_Y,fm$y)

saveRDS(list(model=out_rf,test_X=test_X,test_Y=test_Y),file.path(out_path,'rf_FR_final.rds'))

# classification: train/val/test
# overlap
# 2-class: interior/exterior
# assess classifiers on training data repeated CV

ncores <- 60
cl <- makeCluster(ncores,type='SOCK')
registerDoSNOW(cl)

OTU <- OTU_OVERLAP
META <- META_OVERLAP

sMETA <- META[rownames(OTU),] %>%
  dplyr::select(id=SampleID,surface=sample.Surface,region=region,FR=FR,IE=IE)

sOTU <- OTU[sMETA$id,]
FILT <- (colSums(sOTU>0) < 3) # so no empties in test set
sOTU <- sOTU[,!FILT]
sOTU <- sOTU[rowSums(sOTU) > 0,]
sMETA <- sMETA %>% filter(id %in% rownames(sOTU))
Y <- as.factor(sMETA$IE)

qOTU <- apply(sOTU,2,qnormalize)
zOTU <- apply(sOTU,2,znormalize)

Xz <- zOTU[sMETA$id,]
Xq <- qOTU[sMETA$id,]

names(Y) <- sMETA$id

set.seed(78) 

p <- .8
train_idx <- createDataPartition(Y,times=1,p=p)$Resample1
train_names <- names(Y)[train_idx]
train_Y <- Y[names(Y) %in% train_names]
test_Y <- Y[!(names(Y) %in% train_names)]

X <- Xz
train_X <- X[rownames(X) %in% train_names,]
test_X <- X[!(rownames(X) %in% train_names),]

tr_ctrl <- trainControl(method='repeatedCV', number=10, repeats=10, 
                        classProbs=TRUE,
                        sampling='down',
                        seeds=rf_seeds,
                        summaryFunction=twoClassSummary,
                        allowParallel=TRUE)

param_sweep <- expand.grid(mtry=floor(seq(ncol(train_X)^.25,ncol(train_X)^.75,length=5)))
out_rf <- train(x=train_X,
                y=train_Y,
                method='rf',
                ntree=128,
                tuneGrid=param_sweep,
                importance=TRUE,
                proximity=TRUE,
                metric='ROC', 
                maximize=TRUE, #FALSE,
                trControl=tr_ctrl)

stopCluster(cl)


out_rf

pred_p <- predict(out_rf,test_X,type='prob')
pred_Y <- factor(colnames(pred_p)[apply(pred_p,1,which.max)],levels=colnames(pred_p))
pred_out <- data.frame(obs=test_Y,
                       pred=pred_Y,
                       pred_p)
twoClassSummary(pred_out,lev=levels(pred_out$pred))
caret::confusionMatrix(pred_Y,test_Y)

ROC <- roc(response = test_Y, 
           predictor = pred_p[,1], 
           levels = rev(unique(test_Y)))
plot(ROC, col='red',lwd = 2)


fm <- out_rf$finalModel
pred_p <- fm$votes
pred_Y <- fm$predicted
pred_out <- data.frame(obs=fm$y,
                       pred=pred_Y,
                       pred_p)
twoClassSummary(pred_out,lev=levels(pred_out$pred))
caret::confusionMatrix(pred_Y,fm$y)

saveRDS(list(model=out_rf,test_X=test_X,test_Y=test_Y),file.path(out_path,'rf_IE_final.rds'))

### perform differential abundance analysis on human data
### classes to make comparison are based on best performing classes
### during classification
### using deseq2 pipeline 

# humann
# surface

PS <- phyloseq(otu_table(HUM,taxa_are_rows = FALSE),
               sample_data(META_HUM),
               tax_table(matrix(colnames(HUM),ncol=1,dimnames=list(colnames(HUM),'Pathway'))))
surfaces <- c('Stethoscope','RearLights_controlPanel','RearBench_seats')
PS1 <- subset_samples(PS, sample.Surface %in% surfaces)
PS1 <- filter_taxa(PS1,function(x) sum(x) > 0, TRUE)

set.seed(54)
comb <- expand.grid(surfaces,surfaces,stringsAsFactors = FALSE)
comb <- unique(t(apply(comb[comb[,1] != comb[,2],],1,sort)))
RES <- vector(mode='list',length=nrow(comb))
for (i in 1:nrow(comb)){
  PS2 <- subset_samples(PS1, sample.Surface %in% comb[i,])
  PS2 <- filter_taxa(PS2,function(x) sum(x) > 0, TRUE)
  diagdds <- phyloseq_to_deseq2(PS2, ~ sample.Surface)
  diagdds <- DESeq(diagdds, test='Wald', fitType='parametric')
  res <- results(diagdds, cooksCutoff = FALSE)
  RES[[i]] <- res
}

N <- sum(sapply(RES,nrow))
alpha <- .01 # alpha for BH FDR correction
for (i in seq_along(RES)) RES[[i]]$padj2 <- p.adjust(RES[[i]]$pvalue,method='BH',n=N)
RES3 <- sapply(seq_along(RES),function(i) RES[[i]][which(RES[[i]]$padj2 < alpha),])

saveRDS(list(res=RES,sig=RES3),file.path(out_path,'deseq2_humann.rds'))

sink(file.path(out_path,'deseq2_humann.txt'))
for (i in seq_along(RES)){
  comp <- sort(unlist(comb[i,]),decreasing=TRUE)
  cat(comp[1],'vs',comp[2],'\n\n')
  
  R <- round(as.data.frame(RES[[i]]),5)
  R_names <- rownames(R)
  cat('PW',paste0(colnames(R)),'\n',sep='\t')
  for(j in 1:nrow(R)){
    cat(c(R_names[j],paste0(R[j,])),'\n',sep='\t')
  }
  cat('\n')
  
}
sink()

pdf(file.path(out_path,'volc_deseq.pdf'),width=20,height=10)
for(i in seq_along(RES)) print(plot_volcano(RES[[i]],HUM,HUMPW,1,.01,7,.95))
for(i in seq_along(RES)) print(plot_volcano(RES[[i]],HUM,HUMPW,2,.01,7,.95))
dev.off()

pdf(file.path(out_path,'hm_deseq.pdf'),width=20,height=20)
plot_hm(RES3,comb,HUM,HUMPW,7,lasCol=1,cexCol=1.5)
dev.off()

# humann
# region

PS <- phyloseq(otu_table(HUM,taxa_are_rows = FALSE),
               sample_data(META_HUM),
               tax_table(matrix(colnames(HUM),ncol=1,dimnames=list(colnames(HUM),'Pathway'))))
regions <- c('s_e','s_w_w_coast','w','w_coast')
PS1 <- subset_samples(PS, region %in% regions)
PS1 <- filter_taxa(PS1,function(x) sum(x) > 0, TRUE)

set.seed(54)
comb <- expand.grid(regions,regions,stringsAsFactors = FALSE)
comb <- unique(t(apply(comb[comb[,1] != comb[,2],],1,sort)))
RES <- vector(mode='list',length=nrow(comb))
for (i in 1:nrow(comb)){
  PS2 <- subset_samples(PS1, region %in% comb[i,])
  PS2 <- filter_taxa(PS2,function(x) sum(x) > 0, TRUE)
  diagdds <- phyloseq_to_deseq2(PS2, ~ region)
  diagdds <- DESeq(diagdds, test='Wald', fitType='parametric')
  res <- results(diagdds, cooksCutoff = FALSE)
  RES[[i]] <- res
}

N <- sum(sapply(RES,nrow))
alpha <- .01
for (i in seq_along(RES)) RES[[i]]$padj2 <- p.adjust(RES[[i]]$pvalue,method='BH',n=N)
RES3 <- sapply(seq_along(RES),function(i) RES[[i]][which(RES[[i]]$padj2 < alpha),])

saveRDS(list(res=RES,sig=RES3),file.path(out_path,'deseq2_region_humann.rds'))

sink(file.path(out_path,'deseq2_region_humann.txt'))
for (i in seq_along(RES)){
  comp <- sort(unlist(comb[i,]),decreasing=TRUE)
  cat(comp[1],'vs',comp[2],'\n\n')
  
  R <- round(as.data.frame(RES[[i]]),5)
  R_names <- rownames(R)
  cat('PW',paste0(colnames(R)),'\n',sep='\t')
  for(j in 1:nrow(R)){
    cat(c(R_names[j],paste0(R[j,])),'\n',sep='\t')
  }
  cat('\n')
  
}
sink()

pdf(file.path(out_path,'volc_region_deseq.pdf'),width=20,height=10)
for(i in seq_along(RES)) print(plot_volcano(RES[[i]],HUM,HUMPW,1,.01,7,.95))
for(i in seq_along(RES)) print(plot_volcano(RES[[i]],HUM,HUMPW,2,.01,7,.95))
dev.off()

pdf(file.path(out_path,'hm_region_deseq.pdf'),width=20,height=20)
plot_hm(RES3,comb,HUM,HUMPW,7,lasCol=1,cexCol=1.5)
dev.off()

# humann
# city

PS <- phyloseq(otu_table(HUM,taxa_are_rows = FALSE),
               sample_data(META_HUM),
               tax_table(matrix(colnames(HUM),ncol=1,dimnames=list(colnames(HUM),'Pathway'))))
cities <- ##############################################
PS1 <- subset_samples(PS, city %in% cities)
PS1 <- filter_taxa(PS1,function(x) sum(x) > 0, TRUE)

set.seed(54)
comb <- expand.grid(cities,cities,stringsAsFactors = FALSE)
comb <- unique(t(apply(comb[comb[,1] != comb[,2],],1,sort)))
RES <- vector(mode='list',length=nrow(comb))
for (i in 1:nrow(comb)){
  PS2 <- subset_samples(PS1, city %in% comb[i,])
  PS2 <- filter_taxa(PS2,function(x) sum(x) > 0, TRUE)
  diagdds <- phyloseq_to_deseq2(PS2, ~ city)
  diagdds <- DESeq(diagdds, test='Wald', fitType='parametric')
  res <- results(diagdds, cooksCutoff = FALSE)
  RES[[i]] <- res
}

N <- sum(sapply(RES,nrow))
alpha <- .01
for (i in seq_along(RES)) RES[[i]]$padj2 <- p.adjust(RES[[i]]$pvalue,method='BH',n=N)
RES3 <- sapply(seq_along(RES),function(i) RES[[i]][which(RES[[i]]$padj2 < alpha),])

saveRDS(list(res=RES,sig=RES3),file.path(out_path,'deseq2_city_humann.rds'))

sink(file.path(out_path,'deseq2_city_humann.txt'))
for (i in seq_along(RES)){
  comp <- sort(unlist(comb[i,]),decreasing=TRUE)
  cat(comp[1],'vs',comp[2],'\n\n')
  
  R <- round(as.data.frame(RES[[i]]),5)
  R_names <- rownames(R)
  cat('PW',paste0(colnames(R)),'\n',sep='\t')
  for(j in 1:nrow(R)){
    cat(c(R_names[j],paste0(R[j,])),'\n',sep='\t')
  }
  cat('\n')
  
}
sink()

pdf(file.path(out_path,'volc_city_deseq.pdf'),width=20,height=10)
for(i in seq_along(RES)) print(plot_volcano(RES[[i]],HUM,HUMPW,1,.01,7,.95))
for(i in seq_along(RES)) print(plot_volcano(RES[[i]],HUM,HUMPW,2,.01,7,.95))
dev.off()

pdf(file.path(out_path,'hm_city_deseq.pdf'),width=10,height=10)
plot_hm(RES3,comb,HUM,HUMPW,7,lasCol=2,cexCol=1)
dev.off()