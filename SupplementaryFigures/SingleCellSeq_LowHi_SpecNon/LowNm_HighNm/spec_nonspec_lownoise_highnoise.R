library(ggplot2)
library(plyr)
library('scales')

data <- read.delim('~/Dropbox/densitypaper/densityfigures/Figure5/Nm_all_genes.txt',header=T,stringsAsFactors=F)

bulkData <- read.delim('~/Dropbox/densitypaper/BulkSeqData/CRL_A549_Bulk_Fpkm.txt',header=T,stringsAsFactors=F)

#fpkmDat <- read.delim('~/Dropbox/densitypaper/SingleCellSeqData/genomicDataFpkm.txt',header=T,stringsAsFactors=F)
#volumeDat <- read.delim('~/Dropbox/densitypaper/SingleCellSeqData/fpkmRatio.txt',header=T,stringsAsFactors=F)
#fpkmDat <- merge(fpkmDat,volumeDat)
#fpkmDat <- subset(fpkmDat,fpkmRatio>30)
fpkmDat <- read.table('~/Dropbox/densitypaper/SingleCellSeqData/seqData_volume_and_counts_allGenes.txt',header=T,stringsAsFactors=F)
meanFpkmDat <- ddply(fpkmDat,.(gene_id),summarise,meanFpkm=mean(fpkm))
genesHighFpkm <- subset(meanFpkmDat,meanFpkm>10)$gene_id

bulkData$group <- 'NotConsidered'

bulkData$group[bulkData$fpkmCrl>5*bulkData$fpkmA549 & bulkData$fpkmCrl>5] <- 'crlSpecific'
bulkData$group[bulkData$fpkmA549>5*bulkData$fpkmCrl & bulkData$fpkmA549>5] <- 'a549Specific'
bulkData$group[bulkData$fpkmA549<2*bulkData$fpkmCrl & bulkData$fpkmCrl<2*bulkData$fpkmA549 & bulkData$fpkmCrl>5 & bulkData$fpkmA549>5] <- 'equalExpression'
colnames(bulkData) <- c('gene','fpkmCrl','fpkmA549','group')
bulkData <- subset(bulkData,gene %in% genesHighFpkm)

specificGenes <- subset(bulkData,group %in% 'crlSpecific')$gene
nonspecificGenes <- subset(bulkData,group %in% 'equalExpression')$gene

lowNoise <- read.delim('~/Dropbox/densitypaper/SupplementaryFigures/SingleCellSeq_LowHi_SpecNon/LowNm_HighNm/low_Nm_after_transform.txt',header=F,stringsAsFactors=F)
colnames(lowNoise) <- 'gene'
lowNoiseGenes <- lowNoise$gene

highNoise <- read.delim('~/Dropbox/densitypaper/SupplementaryFigures/SingleCellSeq_LowHi_SpecNon/LowNm_HighNm/high_Nm_after_transform.txt',header=F,stringsAsFactors=F)
colnames(highNoise) <- 'gene'
highNoiseGenes <- highNoise$gene


sum(highNoiseGenes %in% specificGenes)/length(specificGenes)
sum(highNoiseGenes %in% nonspecificGenes)/length(nonspecificGenes)
sum(lowNoiseGenes %in% specificGenes)/length(specificGenes)
sum(lowNoiseGenes %in% nonspecificGenes)/length(nonspecificGenes)

sum(highNoiseGenes %in% specificGenes)/length(highNoiseGenes)
sum(highNoiseGenes %in% nonspecificGenes)/length(highNoiseGenes)
sum(lowNoiseGenes %in% specificGenes)/length(lowNoiseGenes)
sum(lowNoiseGenes %in% nonspecificGenes)/length(lowNoiseGenes)

highNoiseWithGroups <- merge(highNoise,bulkData[,c('gene','group')])
highNoiseWithGroups$noise <- 'high'
lowNoiseWithGroups <- merge(lowNoise,bulkData[,c('gene','group')])
lowNoiseWithGroups$noise <- 'low'

dat <- rbind(highNoiseWithGroups,lowNoiseWithGroups)
dat$group[dat$group %in% c('a549Specific','NotConsidered')] <- 'other'

#pdf('~/Dropbox/densitypaper/densityfigures/Figure5/GO_Stuff/Spec_Nonspec_LowNoise_HighNoise.pdf',height=3,width=1.5)
ggplot(dat,aes(x=noise,fill=group)) + geom_bar() + 
  scale_fill_manual(values=c('midnightblue','mediumvioletred','mediumturquoise')) +
  theme_classic() + theme(legend.position='none') +
  xlab('Noise Level') + 
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6'))
#dev.off()

meanFpkmDat$specNon <- 'other'
meanFpkmDat$specNon[meanFpkmDat$gene_id %in% specificGenes] <- 'specific'
meanFpkmDat$specNon[meanFpkmDat$gene_id %in% nonspecificGenes] <- 'equalExpression'

meanFpkmDat$noise <- 'other'
meanFpkmDat$noise[meanFpkmDat$gene_id %in% lowNoiseGenes] <- 'low'
meanFpkmDat$noise[meanFpkmDat$gene_id %in% highNoiseGenes] <- 'high'

t.test(subset(meanFpkmDat,specNon %in% 'specific')$meanFpkm,subset(meanFpkmDat,specNon %in% 'equalExpression')$meanFpkm)
t.test(subset(meanFpkmDat,noise %in% 'low')$meanFpkm,subset(meanFpkmDat,noise %in% 'high')$meanFpkm)

ggplot(subset(meanFpkmDat,noise %in% c('low','high')),aes(x=noise,y=meanFpkm)) + geom_boxplot() + ylim(c(0,100))
ggplot(subset(meanFpkmDat,noise %in% c('low','high')),aes(x=noise,y=meanFpkm)) + geom_boxplot()

ggplot(subset(meanFpkmDat,noise %in% c('low','high')),aes(x=meanFpkm,fill=noise)) + geom_density(alpha=.25) + xlim(c(0,100))

mean(subset(meanFpkmDat,noise %in% 'other')$meanFpkm)
mean(subset(meanFpkmDat,noise %in% 'low')$meanFpkm)
mean(subset(meanFpkmDat,noise %in% 'high')$meanFpkm)
