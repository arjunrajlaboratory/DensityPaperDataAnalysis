library(ggplot2)
library(plyr)

data <- read.delim('~/Dropbox/densitypaper/densityfigures/Figure5/Nm_all_genes.txt',header=T,stringsAsFactors=F)

bulkData <- read.delim('~/Dropbox/densitypaper/BulkSeqData/CRL_A549_Bulk_Fpkm.txt',header=T,stringsAsFactors=F)

fpkmDat <- read.delim('~/Dropbox/densitypaper/SingleCellSeqData/genomicDataFpkm.txt',header=T,stringsAsFactors=F)
meanFpkmDat <- ddply(fpkmDat,.(gene_id),summarise,meanFpkm=mean(fpkm))
genesHighFpkm <- subset(meanFpkmDat,meanFpkm>10)$gene_id

bulkData$group <- 'NotConsidered'

bulkData$group[bulkData$fpkmCrl>5*bulkData$fpkmA549 & bulkData$fpkmCrl>5] <- 'crlSpecific'
bulkData$group[bulkData$fpkmA549>5*bulkData$fpkmCrl & bulkData$fpkmA549>5] <- 'a549Specific'
bulkData$group[bulkData$fpkmA549<2*bulkData$fpkmCrl & bulkData$fpkmCrl<2*bulkData$fpkmA549 & bulkData$fpkmCrl>5 & bulkData$fpkmA549>5] <- 'equalExpression'
colnames(bulkData) <- c('gene','fpkmCrl','fpkmA549','group')
bulkData <- subset(bulkData,gene %in% genesHighFpkm)

data <- merge(data,bulkData[,c(1,4)])
data <- subset(data,gene %in% genesHighFpkm)

ggplot(bulkData,aes(x=fpkmCrl,y=fpkmA549,color=group)) + geom_point() + 
  theme_classic() + theme(legend.position='none') +
  scale_x_log10() + scale_y_log10() +
  scale_color_manual(values=c('gray80','blueviolet','darkblue','gray80')) +
  xlab('CRL Bulk FPKM') + ylab('A549 Bulk FPKM')

dataSubset <- subset(data,!(group %in% 'NotConsidered' | group %in% 'a549Specific'))

dataSubset$meanCountInverse <- 1/dataSubset$meanCount

ggplot(dataSubset,aes(x=meanCount,y=Nm_seq,color=group)) + geom_point()
ggplot(dataSubset,aes(x=meanCountInverse,y=Nm_seq,color=group)) + geom_point()

fit <- lm(Nm_seq ~ meanCountInverse + factor(group), data=dataSubset)
summary(fit)
anova(fit)

nmBoxplot <- ggplot(dataSubset,aes(x=group,y=Nm_seq,fill=group)) + geom_boxplot(outlier.size=.5,weight=.5) +
  theme_classic() + theme(legend.position='none') +
  ylab('Noise measure') + xlab('') +
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6'))

countBoxplot <- ggplot(dataSubset,aes(x=group,y=meanCount,fill=group)) + geom_boxplot(outlier.size=.5,weight=.5) +
  theme_classic() + theme(legend.position='none') +
  ylab('Noise measure') + xlab('') +
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6'))

t.test(subset(dataSubset,group %in% 'equalExpression')$Nm_seq,subset(dataSubset,group %in% 'crlSpecific')$Nm_seq)

t.test(subset(dataSubset,group %in% 'equalExpression')$meanCount,subset(dataSubset,group %in% 'crlSpecific')$meanCount)

nmDensityPlot <- ggplot(dataSubset,aes(x=Nm_seq,fill=group)) + geom_density(alpha=.25) + xlim(c(0,10)) +
  theme_classic() + theme(legend.position='none') +
  xlab('Noise measure') + ylab('Density') +
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6'))

countDensityPlot <- ggplot(dataSubset,aes(x=meanCount,fill=group)) + geom_density(alpha=.25) + xlim(c(0,1000)) +
  theme_classic() + theme(legend.position='none') +
  xlab('Mean mRNA count') + ylab('Density') +
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6'))

pdf(file='~/Dropbox/densitypaper/densityfigures/Figure5/Spec_NonSpec_Genes/nmBoxplot.pdf',height=2,width=2.2)
nmBoxplot
dev.off()

pdf(file='~/Dropbox/densitypaper/densityfigures/Figure5/Spec_NonSpec_Genes/countBoxplot.pdf',height=2,width=2.2)
countBoxplot
dev.off()

pdf(file='~/Dropbox/densitypaper/densityfigures/Figure5/Spec_NonSpec_Genes/nmDensityPlot.pdf',height=2,width=2.2)
nmDensityPlot
dev.off()

pdf(file='~/Dropbox/densitypaper/densityfigures/Figure5/Spec_NonSpec_Genes/countDensityPlot.pdf',height=2,width=2.2)
countDensityPlot
dev.off()