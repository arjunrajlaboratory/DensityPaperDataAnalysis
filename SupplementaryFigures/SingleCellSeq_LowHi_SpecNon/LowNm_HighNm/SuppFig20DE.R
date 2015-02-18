library(Hmisc)
library(dendroextras)
library(reshape2)
library(plyr)
library(stringr)
library(RColorBrewer)
library(gplots)
library(ggplot2)
library(data.table)
library(msarc)
library(GO.db)

# Pick out genes for GO analysis

data <- read.delim('~/Dropbox/densitypaper/densitypaperdataanalysis/densityfigures/Figure6/Nm_all_genes.txt',header=T,stringsAsFactors=F)
fpkmDat <- read.table('~/Dropbox/densitypaper/densitypaperdataanalysis/SingleCellSeqData/seqData_volume_and_counts_allGenes.txt',header=T,stringsAsFactors=F)
volumeDat <- read.delim('~/Dropbox/densitypaper/densitypaperdataanalysis/SingleCellSeqData/fpkmRatio.txt',header=T,stringsAsFactors=F)
fpkmDat <- merge(fpkmDat,volumeDat)
fpkmDat <- subset(fpkmDat,fpkmRatio>30)
meanFpkm <- ddply(fpkmDat,.(gene_id),summarise,meanFpkm = mean(fpkm))
colnames(meanFpkm) <- c('gene','meanFpkm')
data <- merge(data,meanFpkm)
highFpkmData <- subset(data,meanFpkm>10)

fit <- lm(data=highFpkmData,log10(Nm_seq) ~ log10(meanFpkm))
fitSlope <- coef(fit)[[2]]
fitInt <- coef(fit)[[1]]

pdf('~/Dropbox/densitypaper/densitypaperdataanalysis/SupplementaryFigures/SingleCellSeq_LowHi_SpecNon/LowNm_HighNm/UntransformedNm.pdf',height=2,width=2)
ggplot(highFpkmData,aes(x=log10(meanFpkm),y=log10(Nm_seq))) + geom_point(size=.3) + 
  geom_abline(slope=coef(fit)[[2]],intercept=coef(fit)[[1]],col='red') +
  theme_classic() +
  xlab('Log 10 (mean FPKM)') + ylab('Log 10 (Nm)') +
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6'))
dev.off()

tmp <- highFpkmData
tmp$log10Nm <- log10(tmp$Nm_seq)
tmp$log10Fpkm <- log10(tmp$meanFpkm)


tmp$yToAdd <- -(fitSlope*tmp$log10Fpkm+fitInt)
tmp$transformedLog10Nm <- tmp$log10Nm + tmp$yToAdd

pdf('~/Dropbox/densitypaper/densitypaperdataanalysis/SupplementaryFigures/SingleCellSeq_LowHi_SpecNon/LowNm_HighNm/TransformedNm.pdf',height=2,width=2)
ggplot(tmp,aes(x=log10Fpkm,y=transformedLog10Nm)) + geom_point(size=.3) + 
  geom_hline(yintercept=0.5) + geom_hline(yintercept=-0.5) + geom_hline(yintercept=0,col='red') +
  theme_classic() +
  xlab('Log 10 (mean FPKM)') + ylab('Log 10 (Transformed Nm)') +
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6'))
dev.off()

highNmGenes <- subset(tmp,transformedLog10Nm>0.5)$gene
write.table(highNmGenes,'~/Dropbox/densitypaper/densitypaperdataanalysis/SupplementaryFigures/SingleCellSeq_LowHi_SpecNon/LowNm_HighNm/high_Nm_after_transform.txt',quote=F,sep='\t',row.names=F,col.names=F)
lowNmGenes <- subset(tmp,transformedLog10Nm<(-0.5))$gene
write.table(lowNmGenes,'~/Dropbox/densitypaper/densitypaperdataanalysis/SupplementaryFigures/SingleCellSeq_LowHi_SpecNon/LowNm_HighNm/low_Nm_after_transform.txt',quote=F,sep='\t',row.names=F,col.names=F)
