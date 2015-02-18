library(ggplot2)
library(plyr)

data <- read.delim('~/Dropbox/densitypaper/densitypaperdataanalysis/ExtractedData/triptolide.txt',header=T,stringsAsFactors=F)

data <- data[,1:7]
data <- unique(data)

freqDat <- ddply(data,.(gene,conc),summarise,meanNumTxnSites=mean(numTxnSites),seTxnSites=sd(numTxnSites)/sqrt(length(numTxnSites)))

tmp <- subset(freqDat,gene %in% 'MYC')

pdf('~/Dropbox/densitypaper/densitypaperdataanalysis/SupplementaryFigures/PolII/MYCfreq.pdf',width=2,height=4)
ggplot(tmp,aes(x=conc,y=meanNumTxnSites,fill=conc)) + geom_bar(stat='identity') + expand_limits(y=2) +
  geom_errorbar(aes(ymin=meanNumTxnSites-seTxnSites,ymax=meanNumTxnSites+seTxnSites),width=.2) +
  theme_classic() + theme(legend.position='') +
  theme(axis.title=element_text(size="12"), axis.text=element_text(size='12')) +
  xlab('Triptolide Conc.') +
  ylab('Mean number of transcription sites per cell') +
  scale_fill_manual(values = c('deepskyblue4','deepskyblue1'))
dev.off()

t.test(subset(data,gene %in% 'MYC' & conc %in% 'Ctrl')$numTxnSites,subset(data,gene %in% 'MYC' & conc %in% '100nM')$numTxnSites)

#############################

tmp <- subset(freqDat,gene %in% 'UBC')

freqDatUBC <- ddply(subset(data,gene %in% 'UBC' & date %in% '141105'),.(conc),summarise,meanNumTxnSites=mean(numTxnSites),seTxnSites=sd(numTxnSites)/sqrt(length(numTxnSites)))

pdf('~/Dropbox/densitypaper/densitypaperdataanalysis/SupplementaryFigures/PolII/UBCfreq.pdf',width=2,height=4)
ggplot(subset(tmp,conc %in% c('Ctrl','100nM')),aes(x=conc,y=meanNumTxnSites,fill=conc)) + geom_bar(stat='identity') + expand_limits(y=2) +
  geom_errorbar(aes(ymin=meanNumTxnSites-seTxnSites,ymax=meanNumTxnSites+seTxnSites),width=.2) +
  theme_classic() + theme(legend.position='') +
  theme(axis.title=element_text(size="12"), axis.text=element_text(size='12')) +
  xlab('Triptolide Conc.') +
  ylab('Mean number of transcription sites per cell') +
  scale_fill_manual(values = c('darkslateblue','darkorchid3'))
dev.off()

t.test(subset(data,gene %in% 'UBC' & conc %in% 'Ctrl')$numTxnSites,subset(data,gene %in% 'UBC' & conc %in% '100nM')$numTxnSites)

