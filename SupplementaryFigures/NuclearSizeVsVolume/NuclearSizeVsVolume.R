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

pdf('~/Dropbox/densitypaper/SupplementaryFigures/NuclearSizeVsVolume/NuclearSizeVsVolume.pdf',width=5,height=3.5)
ggplot(data,aes(x=volume,y=nucArea)) +
  geom_point(size=1) + expand_limits(y=0,x=0) +
  theme_bw() + theme_classic() + theme(legend.position = "none") + 
  xlab('Volume (picoliter)') + ylab('Nuclear Area (um^2)') +
  theme(axis.title=element_text(size="10"), axis.text=element_text(size='10'))
dev.off()
