setwd('~/Dropbox/Code/R/TxnSiteAnalysis/')
library(reshape2)
library(ggplot2)

otown <- read.delim('alldata_crl.txt')

otown <- dcast(otown, '... ~ variable')

numericcols <- c('area','avgExon','avgIntron','cytoGAPDH','cytoRNA','nucArea','nucGAPDH','nucRNA',
                 'numCyclin','numRnaPerTxnSiteAvgExon','numTxnSites','totalGAPDH','totalRNA','volume')

for (col in numericcols){
  otown[[col]] <- as.numeric(otown[[col]])
}

otown$id <- NA
otown$id[otown$date==120501] <- 'Senescent'
otown$id[otown$date==120801] <- 'Quiescent'
otown$id[otown$date==120809] <- 'Cycling'

ggplot(otown, 
       aes(x=volume,y=totalRNA,color=cellCycleStage)) + 
  geom_point() + facet_wrap(~gene,scales = "free_y")+geom_rug()

subset(otown,gene=='EEF2')[,c('totalRNA','date','cellType')]

# 4 genes that corr with volume

ggplot(subset(otown,gene %in% c("EEF2","LMNA","BABAM1","ZNF444")), 
       aes(x=volume,y=totalRNA,color=gene)) + 
  geom_point() + facet_wrap(~gene, scales = 'free_y') + expand_limits(y=0) + 
  xlab('Volume (femtoliter)') + ylab('Total RNA')

# quiescent, cycling, senescent

ggplot(subset(otown,id %in% c('Quiescent')), 
       aes(x=volume,y=totalRNA,color=id)) + 
  geom_point() + facet_wrap(~gene,scales = "free_y") + expand_limits(y=0)

ggplot(subset(otown,id %in% c('Quiescent','Cycling')), 
       aes(x=volume,y=totalRNA,color=id)) + 
  geom_point() + facet_wrap(~gene,scales = "free_y") + expand_limits(y=0)

ggplot(subset(otown,id %in% c('Quiescent','Cycling','Senescent')), 
       aes(x=volume,y=totalRNA,color=id)) + 
  geom_point() + facet_wrap(~gene,scales = "free_y")+ expand_limits(y=0) + 
  xlab('Volume (femtoliter)') + ylab('Total RNA')

ggplot(subset(otown,gene %in% "EEF2"), 
       aes(x=volume,y=totalRNA,color=id)) + 
  geom_point() + facet_wrap(~gene,scales = "free_y")+geom_rug()

# CORR

data <- otown[!is.na(otown$volume),]
newData <- matrix(nrow=length(unique(otown$gene)),ncol=5)
colnames(newData) <- c('gene','corr','freqNeg','freqPos','abund')
j <- 0
for(i in 1:length(unique(otown$gene))) {
  geneiter <- as.character(unique(otown$gene))[i]
  j <- j + 1
  tmp <- subset(otown, gene %in% geneiter)
  
  ### CORR
  
  #geneData <- tmp[,c('volume','cytoRNA')]
  fm <- lm(as.numeric(cytoRNA) ~ as.numeric(volume), data=tmp)
  rs <- summary(fm)[c("r.squared", "adj.r.squared")]$r.squared
  
  ### FREQ
  
  gone <- subset(tmp, cellCycleStage %in% 'G1')
  freqNeg <- mean(as.numeric(gone$numTxnSites[gone$numTxnSites<3]))/2
  
  gtwo <- subset(tmp, cellCycleStage %in% 'G2')
  freqPos <- mean(as.numeric(gtwo$numTxnSites[gtwo$numTxnSites<5]))/4
  
  ### ABUND
  
  abund <- mean(as.numeric(tmp$cytoRNA))
  
  newData[j,] <- c(geneiter, sqrt(rs), freqNeg, freqPos, abund)
}

newData <- as.data.frame(newData)
plot(as.vector(newData$freqNeg),as.vector(newData$corr))
identify(as.vector(newData$freqNeg),as.vector(newData$corr),labels = as.vector(newData$gene))
#plot(as.vector(newData$freq),as.vector(newData$corr),xlim=c(-.5,1.5),ylim=c(-.5,1.5))
#plot(as.vector(newData$abund),as.vector(newData$corr))

fm <- lm(as.numeric(as.character(freqNeg)) ~ as.numeric(as.character(freqPos)), data=newData)
fm
plot(as.vector(newData$freqPos),as.vector(newData$freqNeg))
identify(as.vector(newData$freqPos),as.vector(newData$freqNeg),newData$gene)
abline(fm)

tmp <- subset(data,gene %in% 'ACTN4')
fm <- lm(as.numeric(cytoRNA) ~ as.numeric(volume), data=tmp)
rs <- summary(fm)[c("r.squared", "adj.r.squared")]$r.squared
sqrt(rs)