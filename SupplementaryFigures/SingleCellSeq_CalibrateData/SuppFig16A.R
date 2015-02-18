library(plyr)
library(ggplot2)

concData <- read.delim('~/Dropbox/densitypaper/densitypaperdataanalysis/SingleCellSeqData/ercc_concentrations.txt',header=T,stringsAsFactors=F)
colnames(concData) <- c('ID','gene_id','subgroup','concMix1','concMix2','foldChange','log2')
erccReads <- read.delim('~/Dropbox/densitypaper/densitypaperdataanalysis/SingleCellSeqData/erccDataFpkm.txt',header=T,stringsAsFactors=F)

erccSummary <- ddply(erccReads,.(gene_id),summarise,
                     meanCountPerLength = mean(counts/exLength),
                     meanFpkm = mean(fpkm))

erccSummary <- merge(erccSummary,concData)

erccReadsVsConcentration <- ggplot(erccSummary,aes(x=log10(concMix1),y=log10(meanCountPerLength))) + geom_point(size=1) + 
  theme_classic() +
  xlab('ERCC reference concentration (attomoles/microliter)') +
  ylab('Mean reads per length') +
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6'))



erccFpkmVsConcentration <- ggplot(erccSummary,aes(x=concMix1,y=meanFpkm)) + geom_point(size=1) + 
  scale_x_log10() + scale_y_log10() +
  theme_classic() +
  geom_hline(yintercept=10) +
  xlab('ERCC reference concentration (attomoles/microliter)') +
  ylab('Mean FPKM') +
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6'))

erccFpkmVsConcentration

pdf(file='~/Dropbox/densitypaper/densitypaperdataanalysis/SupplementaryFigures/SingleCellSeq_CalibrateData/erccFpkmVsConcentration.pdf',height=3,width=3.3)
erccFpkmVsConcentration
dev.off()
