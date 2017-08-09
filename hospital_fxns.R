qnormalize <- function(x) qnorm(rank(x,ties.method='average')/(length(x)+1))
znormalize <- function(x) (x-mean(x))/(max(x)-min(x))

gm_mean <- function(x, na.rm=TRUE) exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))

anobrayss <- function(x,y){
  x <- x[rowSums(x)>0,]
  BRAYMETA <- y[rownames(x),]
  BRAYMETA <- BRAYMETA[!is.na(BRAYMETA$SampleID),]
  SURFACE_NAMES <- names(which(table(BRAYMETA$sample.Surface) > 25))
  BRAYMETA <- BRAYMETA %>% filter(sample.Surface %in% SURFACE_NAMES) %>% as.data.frame()
  rownames(BRAYMETA) <- BRAYMETA$SampleID
  BRAYOTU <- x[BRAYMETA$SampleID,]
  BRAYOTU <- BRAYOTU[,colSums(BRAYOTU)>0]
  BRAYMETA <- BRAYMETA[rownames(BRAYOTU),]
  SS <- BRAYMETA$sample.Surface
  BRAY <- vegdist(BRAYOTU,'bray')
  
  out <- anosim(BRAY,SS)
  print(plot(out))
  out
}
anobrayfr <- function(x,y){
  x <- x[rowSums(x)>0,]
  BRAYMETA <- y[rownames(x),]
  BRAYMETA <- BRAYMETA[!is.na(BRAYMETA$SampleID),]
  rownames(BRAYMETA) <- BRAYMETA$SampleID
  BRAYOTU <- x[BRAYMETA$SampleID,]
  BRAYOTU <- BRAYOTU[,colSums(BRAYOTU)>0]
  BRAYMETA <- BRAYMETA[rownames(BRAYOTU),]
  FR <- ifelse(BRAYMETA$FR=='Front',1,0)
  BRAY <- vegdist(BRAYOTU,'bray')
  
  out <- anosim(BRAY,FR)
  print(plot(out))
  out
}
anobrayie <- function(x,y){
  x <- x[rowSums(x)>0,]
  BRAYMETA <- y[rownames(x),]
  BRAYMETA <- BRAYMETA[!is.na(BRAYMETA$SampleID),]
  rownames(BRAYMETA) <- BRAYMETA$SampleID
  BRAYOTU <- x[BRAYMETA$SampleID,]
  BRAYOTU <- BRAYOTU[,colSums(BRAYOTU)>0]
  BRAYMETA <- BRAYMETA[rownames(BRAYOTU),]
  IE <- ifelse(BRAYMETA$IE=='Interior',1,0)
  BRAY <- vegdist(BRAYOTU,'bray')
  
  out <- anosim(BRAY,IE)
  print(plot(out))
  out
}
anobrayreg <- function(x,y){
  x <- x[rowSums(x)>0,]
  BRAYMETA <- y[rownames(x),]
  BRAYMETA <- BRAYMETA[!is.na(BRAYMETA$SampleID),]
  rownames(BRAYMETA) <- BRAYMETA$SampleID
  BRAYOTU <- x[BRAYMETA$SampleID,]
  BRAYOTU <- BRAYOTU[,colSums(BRAYOTU)>0]
  BRAYMETA <- BRAYMETA[rownames(BRAYOTU),]
  REGION <- BRAYMETA$region
  BRAY <- vegdist(BRAYOTU,'bray')
  
  out <- anosim(BRAY,REGION)
  print(plot(out))
  out
}

gm_mean <- function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

plot_volcano <- function(res,abund,pw,superclass,alpha=.01,n_group=7,lab_percentile=.95){
  
  name <- gsub('^.*\\(MAP\\)\\: ','',res@elementMetadata$description[2])
  name <- gsub('^sample\\.','',name)
  
  res$padj2 <- ifelse(is.na(res$padj2),1,res$padj2)
  df <- data.frame(fc=res$log2FoldChange,pval=res$padj2,class=rownames(res)) %>%
    left_join(pw,by='class') %>%
    left_join(data.frame(abundance=colSums(abund),class=colnames(abund)),by='class')%>%
    dplyr::mutate(pval2=-log10(pval))
  
  if (superclass==2){
    top_sc <- pw %>%
      select(superclass2) %>%
      filter(!is.na(superclass2)) %>%
      group_by(superclass2) %>%
      dplyr::summarise(N=n()) %>%
      arrange(desc(N)) %>%
      select(superclass2) %>%
      unlist() 
    
    top_sc <- top_sc[1:min(n_group,length(top_sc))]
    df$superclass_clean <- ifelse(df$superclass2 %in% top_sc,top_sc,'Other')
    
  }else if(superclass==1){
    top_sc <- pw %>%
      select(superclass1) %>%
      filter(!is.na(superclass1)) %>%
      dplyr::mutate(superclass1=ifelse(superclass1=='Degradation/Utilization/Assimilatio','Degradation/Utilization/Assimilation',superclass1)) %>%
      group_by(superclass1) %>%
      dplyr::summarise(N=n()) %>%
      arrange(desc(N)) %>%
      select(superclass1) %>%
      unlist() 
    
    top_sc <- top_sc[1:min(n_group,length(top_sc))]
    df$superclass_clean <- ifelse(df$superclass1 %in% top_sc,top_sc,'Other')
    
  }
  
  df$superclass_clean <- factor(df$superclass_clean,levels=sort(unique(df$superclass_clean)))
  df$superclass_clean <- relevel(df$superclass_clean,ref='Other')
  
  prop <- prop.table(table(df$superclass_clean[df$pval < alpha]))
  df$prop <- prop[df$superclass_clean]
  
  roof <- ceiling(-log10(min(df$pval)))
  roof <- max(ifelse(roof %% 2 == 0, roof, roof + 1),10)
  
  df %>%
    ggplot() + 
    geom_hline(yintercept=-log10(alpha),colour='red',size=2) +
    geom_point(data = . %>% filter(pval >= alpha),
               aes(fc,pval2),
               colour='gray',size=2,alpha=.7) + 
    geom_point(data = . %>% filter(pval < alpha),
               aes(fc,pval2,colour=superclass_clean,size=prop),
               alpha=.7) +
    scale_color_brewer(type='qual',palette=2,drop=FALSE) +
    geom_text_repel(data = . %>% filter(abs(fc) > quantile(abs(fc),lab_percentile,na.rm=TRUE) & pval < alpha), 
                    aes(fc,pval2,label=class,colour=superclass_clean),
                    size=3,segment.size=0.5,fontface='bold') +
    scale_size_continuous(breaks = seq(0,1,by=.2),labels = seq(0,1,by=.2),range=c(1,10),limits=c(0,1)) +
    scale_y_continuous(breaks=seq(0,roof,by=2),labels=10^(-seq(0,roof,by=2)),limits=c(0,roof)) + 
    xlim(-max(abs(df$fc),na.rm=TRUE),max(abs(df$fc),na.rm=TRUE)) +
    labs(x='log2 Fold Change',y='p-value',colour='Superclass',size='Proportion') +
    ggtitle(name) +
    theme(axis.text = element_text(size=14),
          axis.title = element_text(size=18),
          title = element_text(size=22,face='bold'))
}

plot_hm <- function(res,comb,abund,pw,n_group=7,...){

  top_sc <- pw %>%
    select(superclass2) %>%
    filter(!is.na(superclass2)) %>%
    group_by(superclass2) %>%
    dplyr::summarise(N=n()) %>%
    arrange(desc(N)) %>%
    select(superclass2) %>%
    unlist() 
  
  top_sc <- top_sc[1:min(n_group,length(top_sc))]
  superclass_clean2 <- ifelse(pw$superclass2 %in% top_sc,top_sc,'Other')
  names(superclass_clean2) <- pw$class

  top_sc <- pw %>%
    select(superclass1) %>%
    filter(!is.na(superclass1)) %>%
    dplyr::mutate(superclass1=ifelse(superclass1=='Degradation/Utilization/Assimilatio','Degradation/Utilization/Assimilation',superclass1)) %>%
    group_by(superclass1) %>%
    dplyr::summarise(N=n()) %>%
    arrange(desc(N)) %>%
    select(superclass1) %>%
    unlist() 
  
  top_sc <- top_sc[1:min(n_group,length(top_sc))]
  superclass_clean1 <- ifelse(pw$superclass1 %in% top_sc,top_sc,'Other')
  names(superclass_clean1) <- pw$class
  
  sig_pw_names <- unique(unlist(sapply(res,rownames)))
  sig_pw <- do.call('cbind',lapply(RES,function(x) x[sig_pw_names,'log2FoldChange']))
  dimnames(sig_pw) <- list(sig_pw_names,apply(comb,1,paste0,collapse='_vs_'))
  
  bar2 <- superclass_clean2[rownames(sig_pw)]
  bar1 <- superclass_clean1[rownames(sig_pw)]
  
  bar2[is.na(bar2)] <- 'Other'
  bar1[is.na(bar1)] <- 'Other'
  
  #rownames(sig_pw) <- bar1
  roword <- order(bar1)
  sig_pw <- sig_pw[roword,]
  sig_pw_names <- sig_pw_names[roword]
  bar2 <- bar2[roword]
  bar1 <- bar1[roword]
  
  rowlabs <- rownames(sig_pw)
  collabs <- colnames(sig_pw)
  
  bar2 <- factor(bar2,levels=sort(unique(bar2)))
  if ('Other' %in% bar2) bar2 <- relevel(bar2,ref='Other')
  bar1 <- factor(bar1,levels=sort(unique(bar1)))
  if ('Other' %in% bar1) bar1 <- relevel(bar1,ref='Other')
  
  rowcolidx2 <- as.integer(bar2)
  rowcolidx1 <- as.integer(bar1)
  rowcols2 <- rainbow(length(unique(rowcolidx2)))[rowcolidx2]
  rowcols1 <- rainbow(length(unique(rowcolidx1)))[rowcolidx1]
  
  colleft <- colorRampPalette(c('green','darkgreen'))(25)
  colcent <- colorRampPalette(c('darkgreen','black','darkred'))(100)
  colright <- colorRampPalette(c('darkred','red'))(25)
  colors <- c(colleft,colcent,colright)

  sig_pw <- as.matrix(scale(sig_pw,FALSE,FALSE))
  
  heatmap3::heatmap3(sig_pw,
                     
                     #labRow='',
                     #labCol='',
                     
                     margins=c(15,25),
                     
                     cexRow=.8,
                     #cexCol=1.5,
                     
                     #lasCol=1,
                     
                     Colv=NA,
                     Rowv=NA,
                     
                     RowSideColors=cbind(Superclass1=rowcols1,Superclass2=rowcols2),
                     ColSideWidth=1.5,
                     
                     col=colors,
                     balanceColor=TRUE,
                     
                     scale='none',
                     
                     ...,
                     
  )

}



plot_volcano_amr <- function(res,abund,tax,alpha=.01,n_group=7,lab_percentile=.95){
  
  tax <- data.frame(class=rownames(tax),tax,stringsAsFactors=FALSE)
  name <- gsub('^.*\\(MAP\\)\\: ','',res@elementMetadata$description[2])
  name <- gsub('^sample\\.','',name)
  
  res$padj2 <- ifelse(is.na(res$padj2),1,res$padj2)
  df <- data.frame(fc=res$log2FoldChange,pval=res$padj2,class=rownames(res)) %>%
    left_join(tax,by='class') %>%
    left_join(data.frame(abundance=colSums(abund),class=colnames(abund)),by='class')%>%
    dplyr::mutate(pval2=-log10(pval))
  
  top_sc <- tax %>%
    select(Genus) %>%
    filter(!is.na(Genus)) %>%
    group_by(Genus) %>%
    dplyr::summarise(N=n()) %>%
    arrange(desc(N)) %>%
    select(Genus) %>%
    unlist() 
  
  top_sc <- top_sc[1:min(n_group,length(top_sc))]
  df$superclass_clean <- ifelse(df$Genus %in% top_sc,top_sc,'Other')

  
  df$superclass_clean <- factor(df$superclass_clean,levels=sort(unique(df$superclass_clean)))
  df$superclass_clean <- relevel(df$superclass_clean,ref='Other')
  
  prop <- prop.table(table(df$superclass_clean[df$pval < alpha]))
  df$prop <- prop[df$superclass_clean]
  
  roof <- ceiling(-log10(min(df$pval)))
  roof <- max(ifelse(roof %% 2 == 0, roof, roof + 1),10)
  
  df %>%
    ggplot() + 
    geom_hline(yintercept=-log10(alpha),colour='red',size=2) +
    geom_point(data = . %>% filter(pval >= alpha),
               aes(fc,pval2),
               colour='gray',size=2,alpha=.7) + 
    geom_point(data = . %>% filter(pval < alpha),
               aes(fc,pval2,colour=superclass_clean,size=prop),
               alpha=.7) +
    scale_color_brewer(type='qual',palette=2,drop=FALSE) +
    geom_text_repel(data = . %>% filter(abs(fc) > quantile(abs(fc),lab_percentile,na.rm=TRUE) & pval < alpha), 
                    aes(fc,pval2,label=class,colour=superclass_clean),
                    size=3,segment.size=0.5,fontface='bold') +
    scale_size_continuous(breaks = seq(0,1,by=.2),labels = seq(0,1,by=.2),range=c(1,10),limits=c(0,1)) +
    scale_y_continuous(breaks=seq(0,roof,by=2),labels=10^(-seq(0,roof,by=2)),limits=c(0,roof)) + 
    xlim(-max(abs(df$fc),na.rm=TRUE),max(abs(df$fc),na.rm=TRUE)) +
    labs(x='log2 Fold Change',y='p-value',colour='Genus',size='Proportion') +
    ggtitle(name) +
    theme(axis.text = element_text(size=14),
          axis.title = element_text(size=18),
          title = element_text(size=22,face='bold'))
}

plot_hm_amr <- function(res,comb,abund,tax,n_group=7,...){
  
  tax <- data.frame(class=rownames(tax),tax,stringsAsFactors=FALSE)
  
  top_sc <- tax %>%
    select(Genus) %>%
    filter(!is.na(Genus)) %>%
    group_by(Genus) %>%
    dplyr::summarise(N=n()) %>%
    arrange(desc(N)) %>%
    select(Genus) %>%
    unlist() 
  
  top_sc <- top_sc[1:min(n_group,length(top_sc))]
  superclass_clean <- ifelse(tax$Genus %in% top_sc,top_sc,'Other')
  names(superclass_clean) <- tax$class

  sig_pw_names <- unique(unlist(sapply(res,rownames)))
  sig_pw <- do.call('cbind',lapply(RES,function(x) x[sig_pw_names,'log2FoldChange']))
  dimnames(sig_pw) <- list(sig_pw_names,apply(comb,1,paste0,collapse='_vs_'))
  
  bar <- superclass_clean[rownames(sig_pw)]
  
  bar[is.na(bar)] <- 'Other'
  
  roword <- order(bar)
  sig_pw <- sig_pw[roword,]
  sig_pw_names <- sig_pw_names[roword]
  bar <- bar[roword]
  
  rowlabs <- rownames(sig_pw)
  collabs <- colnames(sig_pw)
  
  bar <- factor(bar,levels=sort(unique(bar)))
  if ('Other' %in% bar) bar <- relevel(bar,ref='Other')
  
  rowcolidx <- as.integer(bar)
  rowcols <- rainbow(length(unique(rowcolidx)))[rowcolidx]
  
  colleft <- colorRampPalette(c('green','darkgreen'))(25)
  colcent <- colorRampPalette(c('darkgreen','black','darkred'))(100)
  colright <- colorRampPalette(c('darkred','red'))(25)
  colors <- c(colleft,colcent,colright)
  
  sig_pw <- as.matrix(scale(sig_pw,FALSE,FALSE))
  
  heatmap3::heatmap3(sig_pw,
                     
                     #labRow='',
                     #labCol='',
                     
                     margins=c(15,25),
                     
                     cexRow=.8,
                     #cexCol=1.5,
                     
                     #lasCol=1,
                     
                     Colv=NA,
                     Rowv=NA,
                     
                     RowSideColors=cbind(Genus=rowcols),
                     ColSideWidth=1.5,
                     
                     col=colors,
                     balanceColor=TRUE,
                     
                     scale='none',
                     
                     ...,
                     
  )
  
}


