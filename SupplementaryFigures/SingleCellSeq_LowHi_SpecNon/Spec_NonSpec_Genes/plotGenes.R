library(ggplot2)
library(plyr)

data <- read.delim('~/Dropbox/densitypaper/densityfigures/Figure5/Nm_all_genes.txt',header=T,stringsAsFactors=F)

bulkData <- read.delim('~/Dropbox/densitypaper/BulkSeqData/CRL_A549_Bulk_Fpkm.txt',header=T,stringsAsFactors=F)

#fpkmDat <- read.table('~/Dropbox/densitypaper/SingleCellSeqData/seqData_volume_and_counts_allGenes.txt',header=T,stringsAsFactors=F)
#meanFpkmDat <- ddply(fpkmDat,.(gene_id),summarise,meanFpkm=mean(fpkm))
#highAbundGenes <- subset(meanFpkmDat,meanFpkm>10)$gene_id

#bulkData <- subset(bulkData,gene %in% highAbundGenes)

bulkData$group <- 'NotConsidered'

bulkData$group[bulkData$fpkmCrl>5*bulkData$fpkmA549 & bulkData$fpkmCrl>5] <- 'crlSpecific'
bulkData$group[bulkData$fpkmA549>5*bulkData$fpkmCrl & bulkData$fpkmA549>5] <- 'a549Specific'
bulkData$group[bulkData$fpkmA549<2*bulkData$fpkmCrl & bulkData$fpkmCrl<2*bulkData$fpkmA549 & bulkData$fpkmCrl>5 & bulkData$fpkmA549>5] <- 'equalExpression'
colnames(bulkData) <- c('gene','fpkmCrl','fpkmA549','group')

data <- merge(data,bulkData[,c(1,4)])

gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}

png(file='~/Dropbox/densitypaper/SupplementaryFigures/SingleCellSeq_LowHi_SpecNon/Spec_NonSpec_Genes/plotGenes.png',width=500,height=550)
ggplot(bulkData,aes(x=fpkmCrl,y=fpkmA549,color=group)) + geom_point(size=0.4) + 
  theme_classic() + theme(legend.position='none') +
  scale_x_log10() + scale_y_log10() +
  scale_color_manual(values=c('gray80',gg_color_hue(4)[1],gg_color_hue(4)[3],'gray80')) +
  xlab('CRL Bulk FPKM') + ylab('A549 Bulk FPKM') +
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6'))
dev.off()
