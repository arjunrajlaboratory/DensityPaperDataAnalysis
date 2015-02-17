# Compare Nm and CV between FISH and seq
library(plyr)
library(ggplot2)

#seqNmData <- read.table('~/Dropbox/densitypaper/densityfigures/Figure5/Nm_all_genes.txt',header=T,stringsAsFactors=F)
#seqData <- read.delim('~/Dropbox/densitypaper/SingleCellSeqData/genomicDataFpkm.txt',header=T,stringsAsFactors=F)
#fishData <- read.delim('~/Dropbox/densitypaper/ExtractedData/CRL.txt',header=T,stringsAsFactors=F)
#fishNmData <- read.delim('~/Dropbox/densitypaper/densityfigures/Figure4/CV_Nm/CV_Nm_Table.txt',header=T,stringsAsFactors=F)

#seqDataTrunc <- subset(seqData,gene_id %in% unique(fishData$gene))

#cvSeq <- ddply(seqDataTrunc,.(gene_id),summarise,cvSeq=sd(fpkm)/mean(fpkm))
#colnames(cvSeq) <- c('gene','cvSeq')
#cvFish <- ddply(fishData,.(gene),summarise,cvFish=sd(totalRNA)/mean(totalRNA))
#cvData <- merge(cvSeq,cvFish)

cvData <- read.delim('~/Dropbox/densitypaper/densityfigures/Figure5/CV_FISH_Seq_compare.txt',header=T,stringsAsFactors=F)

#seqDat <- read.delim('~/Dropbox/densitypaper/SingleCellSeqData/genomicDataFpkm.txt',header=T,stringsAsFactors=F)
seqDat <- read.table('~/Dropbox/densitypaper/SingleCellSeqData/seqData_volume_and_counts_fishGenes.txt',header=T,stringsAsFactors=F)
seqDat <- subset(seqDat,gene_id %in% cvData$gene)
#volumeDat <- read.delim('~/Dropbox/densitypaper/SingleCellSeqData/fpkmRatio.txt',header=T,stringsAsFactors=T)
#seqDat <- merge(seqDat,volumeDat)
meanSeqDat <- ddply(seqDat,.(gene_id),summarise,meanFpkm = mean(fpkm))
genesHighFpkm <- subset(meanSeqDat,meanFpkm>10)$gene_id
lowNoiseCvDat <- subset(cvData,gene %in% genesHighFpkm)

cvData <- subset(cvData,gene %in% genesHighFpkm)


p1 <- ggplot(cvData,aes(x=CV2_fish,y=CV2_seq_fpkm)) + 
  geom_point(size=1.3) + 
  geom_errorbar(aes(ymin=lower_seq_fpkm,ymax=upper_seq_fpkm)) +
  geom_errorbarh(aes(xmin=lower_fish,xmax=upper_fish)) +
  scale_x_log10() + scale_y_log10() +
  expand_limits(x=c(.5,10),y=c(.5,10)) +
  theme_classic() +
  xlab('Count CV^2 (FISH)') + ylab('FPKM CV^2 (Seq)') +
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6'))

cor(log10(cvData$CV2_fish),log10(cvData$CV2_seq_fpkm))

p2 <- ggplot(cvData,aes(x=CV2_fish,y=CV2_seq_count)) + 
  geom_point(size=1.3) + 
  geom_errorbar(aes(ymin=lower_seq_count,ymax=upper_seq_count)) +
  geom_errorbarh(aes(xmin=lower_fish,xmax=upper_fish)) +
  scale_x_log10() + scale_y_log10() +
  expand_limits(x=c(.5,10),y=c(.5,10)) +
  theme_classic() +
  xlab('Count CV^2 (FISH)') + ylab('Count CV^2 (Seq)') +
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6'))

cor(log10(cvData$CV2_fish),log10(cvData$CV2_seq_count))
cor(log10(lowNoiseCvDat$CV2_fish),log10(lowNoiseCvDat$CV2_seq_count))

nmData2 <- read.delim('~/Dropbox/densitypaper/densityfigures/Figure5/Nm_FISH_Seq_compare.txt',header=T,stringsAsFactors=F)
lowNoiseNmDat <- subset(nmData2,gene %in% genesHighFpkm)

p3 <- ggplot(lowNoiseNmDat,aes(x=Nm_fish,y=Nm_seq)) + 
  geom_point(size=1.3) + 
  geom_errorbar(aes(ymin=lower_seq,ymax=upper_seq)) +
  geom_errorbarh(aes(xmin=lower_fish,xmax=upper_fish)) +
  scale_x_log10() + scale_y_log10() +
  expand_limits(x=c(.01,10),y=c(.01,10)) +
  theme_classic() +
  xlab('Nm (FISH)') + ylab('Nm (Seq)') +
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6'))

cor(log(nmData2$Nm_fish),log(nmData2$Nm_seq))
cor(log(lowNoiseNmDat$Nm_fish),log(lowNoiseNmDat$Nm_seq))

pdf(file='~/Dropbox/densitypaper/SupplementaryFigures/SingleCellSeq_CV_Nm_Compare/CV_fpkm.pdf',height=2,width=2.2)
p1
dev.off()

pdf(file='~/Dropbox/densitypaper/SupplementaryFigures/SingleCellSeq_CV_Nm_Compare/CV_count.pdf',height=2,width=2.2)
p2
dev.off()

pdf(file='~/Dropbox/densitypaper/SupplementaryFigures/SingleCellSeq_CV_Nm_Compare/Nm.pdf',height=2,width=2.2)
p3
dev.off()