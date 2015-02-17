library(reshape2)
library(ggplot2)

data <- read.delim('~/Dropbox/densitypaper/ExtractedData/rRNA.txt',header=T,stringsAsFactors=F)

numericcols <- c('volume','area','totalGAPDH','cytoGAPDH','nucGAPDH','nucArea',
                 'dapiIntensityTotal','dapiIntensityAvg','numCyclin',
                 'rrnaIntensityTotal','rrnaIntensityAvg','itsIntensityTotal','itsIntensityAvg');


for (col in numericcols){
  data[[col]] <- as.numeric(data[[col]])
}

data$volume <- data$volume/1000

tmp <- subset(data,date==unique(data$date)[2])

pdf('~/Dropbox/densitypaper/SupplementaryFigures/rRNA/rRNA.pdf',height=3,width=3.3)
ggplot(tmp,aes(x=volume,y=rrnaIntensityTotal)) +
  geom_point(size=1.3) +
  expand_limits(x=0,y=0) +
  xlab('Volume (picoliter)') + ylab('rRNA Intensity') +
  theme_classic() + 
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6'))
dev.off()

pdf('~/Dropbox/densitypaper/SupplementaryFigures/rRNA/rRNA_ITS.pdf',height=3,width=3.3)
ggplot(tmp,aes(x=volume,y=itsIntensityTotal)) +
  geom_point(size=1.3) +
  expand_limits(x=0,y=0) +
  xlab('Volume (picoliter)') + ylab('rRNA ITS Intensity') +
  theme_classic() + 
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6'))
dev.off()