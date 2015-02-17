source('~/Dropbox/densitypaper/densitypaperdataanalysis/multiplot.R')
library(reshape2)
library(ggplot2)

data <- read.delim('~/Dropbox/densitypaper/densitypaperdataanalysis/ExtractedData/CRL.txt',header=T,stringsAsFactors=F)


numericcols <- c('volume','area','totalGAPDH','cytoGAPDH','nucGAPDH','nucArea',
                 'nucIntensityTotal','nucIntensityAvg','numCyclin','totalRNA','cytoRNA',
                 'nucRNA','numTxnSites','avgIntron','avgExon','numRnaPerTxnSiteAvgExon');
                 
  
for (col in numericcols){
  data[[col]] <- as.numeric(data[[col]])
}

data$volume <- data$volume/1000

seedFromFile <- read.table('~/Dropbox/densitypaper/densitypaperdataanalysis/densityfigures/seed1')
.Random.seed <- t(seedFromFile)
gapdhDat <- data[runif(dim(data)[1])<.08,]

gapdhDat$gene <- 'GAPDH'
gapdhDat$cytoRNA <- gapdhDat$cytoGAPDH



# High abundance genes

highFit <- lm(cytoRNA ~ volume, data=subset(data,gene %in% 'EEF2'))
medFit <- lm(cytoRNA ~ volume, data=subset(data,gene %in% 'LMNA'))
lowFit <- lm(cytoRNA ~ volume, data=subset(data,gene %in% 'TBCB'))

selectDat <- subset(data,gene %in% c('EEF2','LMNA','TBCB'))

clrs <- c('royalblue','lightseagreen','mediumpurple2')

pdf(file='~/Dropbox/densitypaper/densitypaperdataanalysis/densityfigures/Figure1/MultipleGenes/high_low_LMNA_EEF2_TBCB.pdf',height=2.5,width=2.7)
ggplot(subset(selectDat,volume<5),
       aes(x=volume,y=cytoRNA,color=gene)) +
  geom_point(size=1.3) + expand_limits(y=0,x=0) +
  theme_bw() + theme_classic() + theme(legend.position = "none") + 
  xlab('Volume (picoliter)') + ylab('Cytoplasmic RNA') +
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6')) +
  scale_colour_manual(values=c('EEF2'=clrs[2],
                               'LMNA'=clrs[1],
                               'TBCB'=clrs[3])) +
  geom_abline(slope=coef(highFit)[[2]],intercept=coef(highFit)[[1]],col=clrs[2],linetype='dashed') +
  geom_abline(slope=coef(medFit)[[2]],intercept=coef(medFit)[[1]],col=clrs[1],linetype='dashed') +
  geom_abline(slope=coef(lowFit)[[2]],intercept=coef(lowFit)[[1]],col=clrs[3],linetype='dashed')
  
dev.off()