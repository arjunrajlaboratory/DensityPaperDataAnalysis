library(reshape2)
library(ggplot2)

data <- read.delim('~/Dropbox/densitypaper/ExtractedData/CRL.txt',header=T,stringsAsFactors=T)

numericcols <- c('volume','area','totalGAPDH','cytoGAPDH','nucGAPDH','nucArea',
                 'nucIntensityTotal','nucIntensityAvg','numCyclin','totalRNA','cytoRNA',
                 'nucRNA','numTxnSites','avgIntron','avgExon','numRnaPerTxnSiteAvgExon');


for (col in numericcols){
  data[[col]] <- as.numeric(data[[col]])
}

data$volume <- data$volume/1000

data <- subset(data,nucArea > 100)

data$stage <- 'G2'
data$stage[data$numCyclin<231] <- 'S'
data$stage[data$numCyclin<21] <- 'G1'

pdf('~/Dropbox/densitypaper/SupplementaryFigures/NuclearSizeVsVolume/NuclearSizeVsVolume_CellCycle.pdf',width=4,height=3)
ggplot(data,aes(x=volume,y=nucArea,col=stage)) +
  geom_point(size=1) + expand_limits(y=0,x=0) +
  theme_bw() + theme_classic() + theme(legend.position = "none") + 
  xlab('Volume (picoliter)') + ylab('Nuclear Area (um^2)') +
  theme(axis.title=element_text(size="10"), axis.text=element_text(size='10'))
dev.off()

pdf('~/Dropbox/densitypaper/SupplementaryFigures/NuclearSizeVsVolume/NuclearSize_Stage.pdf',width=3.3,height=3)
ggplot(data,aes(x=nucArea,fill=factor(stage))) +
  geom_density(alpha=.25) +
  xlab('Nuclear area') + ylab('Density') +
  theme_bw() + theme_classic() + 
  theme(axis.title=element_text(size="12"), axis.text=element_text(size='12'),legend.position='none')
dev.off()