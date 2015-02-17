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

data$stage <- 'G2'
data$stage[data$numCyclin<231] <- 'S'
data$stage[data$numCyclin<21] <- 'G1'

pdf(file='~/Dropbox/densitypaper/SupplementaryFigures/G1_G2_Vols/VolumeDistStages.pdf',height=3,width=3.3)
ggplot(subset(data,stage %in% c('G1','G2')),aes(x=volume,fill=factor(stage))) +
  geom_density(alpha=.25) +
  xlab('Volume (picoliter)') + ylab('Density') +
  theme_bw() + theme_classic() + 
  theme(axis.title=element_text(size="12"), axis.text=element_text(size='12'),legend.position='none')
dev.off()