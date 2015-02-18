library(ggplot2)

data <- read.delim('~/Dropbox/densitypaper/densitypaperdataanalysis/ExtractedData/CRL.txt',
                   header=T,stringsAsFactors=F)

data$stage <- 'G1'
data$stage[data$numCyclin>20] <- 'S'
data$stage[data$numCyclin>230] <- 'G2'

tmp <- subset(data,stage %in% c('G1','G2'))

pdf('~/Dropbox/densitypaper/densitypaperdataanalysis/SupplementaryFigures/CyclinHistone/VolumeDist.pdf',width=3,height=3)
ggplot(tmp,aes(x=volume/1000,fill=stage)) + geom_density(alpha=.25) +
  theme_classic() +
  theme(legend.position='') +
  theme(axis.title=element_text(size="12"), axis.text=element_text(size='12')) +
  xlab('Volume (picoliters)') + ylab('Fraction of cells')
dev.off()