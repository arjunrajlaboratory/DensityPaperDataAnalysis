# Finding genes that correlate inversely with volume
library(plyr)
library(ggplot2)

data <- read.table('~/Dropbox/densitypaper/SingleCellSeqData/seqData_volume_and_counts_allGenes.txt',header=T,stringsAsFactors=F)
nmDat <- read.delim('~/Dropbox/densitypaper/densityfigures/Figure5/Nm_all_genes.txt',header=T,stringsAsFactors=F)

meanTable <- ddply(data,.(gene_id),summarise,
                   meanFpkm = mean(fpkm),
                   maxFpkm = max(fpkm),
                   medianFpkm = median(fpkm))

corTable <- ddply(subset(data,gene_id %in% subset(meanTable,medianFpkm>0)$gene_id),.(gene_id),summarise,
                  r = cor(volume,actualCount),
                  meanFpkm = mean(fpkm))
colnames(corTable) <- c('gene','r')
tmp <- merge(nmDat,corTable)
tmp$colorCol <- 'notOfInterest'
tmp$colorCol[(tmp$Nm_seq<0.2 & tmp$r<0.1 & tmp$r>-0.1)] <- 'lowNoiseZeroR'

pdf('~/Dropbox/densitypaper/SupplementaryFigures/SingleCellSeq_zero_anticorr/rnearzero_criteria.pdf',width=3.3,height=3)
ggplot(tmp,aes(x=Nm_seq,y=r,col=colorCol)) + geom_point(size=.5) + scale_x_log10() +
  theme_classic() +
  xlab('Noise measure') +
  ylab('Correlation coefficient, RNA count with volume') +
  geom_hline(yintercept=0) +
  theme(legend.position='none') +
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6'))
dev.off()

pdf('~/Dropbox/densitypaper/SupplementaryFigures/SingleCellSeq_zero_anticorr/R_hist.pdf',width=3.3,height=3)
ggplot(corTable,aes(x=r)) + geom_density(fill='pink') + theme_classic() +
  xlab('Correlation coefficient (RNA count and volume)') +
  ylab('') +
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6'))
dev.off()

lowNoiseHighR <- subset(tmp,Nm_seq<0.2 & r<0.1 & r>-0.1)
lowNoiseHighR <- lowNoiseHighR[order(lowNoiseHighR$Nm_seq),]

geneList <- lowNoiseHighR$gene[c(1:3)]

pdf('~/Dropbox/densitypaper/SupplementaryFigures/SingleCellSeq_zero_anticorr/zerocorr_BSG.pdf',width=2.2,height=2.2)
ggplot(subset(data,gene_id %in% geneList[1]),aes(x=volume/1000,y=actualCount)) + geom_point(size=1,col='steelblue3') + 
  expand_limits(x=0,y=0) + ggtitle(geneList[1]) +
  xlab('Volume (picoliter)') +
  ylab('Inferred RNA count') +
  theme_classic() +
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6'))
dev.off()

pdf('~/Dropbox/densitypaper/SupplementaryFigures/SingleCellSeq_zero_anticorr/zerocorr_HNRNPL.pdf',width=2.2,height=2.2)
ggplot(subset(data,gene_id %in% geneList[2]),aes(x=volume/1000,y=actualCount)) + geom_point(size=1,col='steelblue3') + 
  expand_limits(x=0,y=0) + ggtitle(geneList[2]) +
  xlab('Volume (picoliter)') +
  ylab('Inferred RNA count') +
  theme_classic() +
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6'))
dev.off()

pdf('~/Dropbox/densitypaper/SupplementaryFigures/SingleCellSeq_zero_anticorr/zerocorr_LY6E.pdf',width=2.2,height=2.2)
ggplot(subset(data,gene_id %in% geneList[3]),aes(x=volume/1000,y=actualCount)) + geom_point(size=1,col='steelblue3') + 
  expand_limits(x=0,y=0) + ggtitle(geneList[3]) +
  xlab('Volume (picoliter)') +
  ylab('Inferred RNA count') +
  theme_classic() +
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6'))
dev.off()