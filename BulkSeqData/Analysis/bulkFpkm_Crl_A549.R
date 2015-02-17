# calculate FPKM (again) for bulk A549 and CRL, using same method as I did for single-cell seq

library(plyr)
library(ggplot2)

crlCountData <- read.delim('~/Dropbox/densitypaper/BulkSeqData/CRL_total1_HTSeq_out_Exon_newGtf.txt',header=F,stringsAsFactors=F)
colnames(crlCountData) <- c('gene_id','count')
a549CountData <- read.delim('~/Dropbox/densitypaper/BulkSeqData/A549_total1_HTSeq_out_Exon.txt',header=F,stringsAsFactors=F)
colnames(a549CountData) <- c('gene_id','count')
lengthData <- read.delim('~/Dropbox/densitypaper/SingleCellSeqData/exon_lengths.txt',header=T,stringsAsFactors=F)
colnames(lengthData) <- c('gene_id','exLength')

crlCountData <- merge(crlCountData,lengthData)
a549CountData <- merge(a549CountData,lengthData)

crlCountData$totalReads <- sum(crlCountData$count)
a549CountData$totalReads <- sum(a549CountData$count)

crlCountData$fpkmCrl <- crlCountData$count/((crlCountData$exLength/1000)*(crlCountData$totalReads/1000000))
a549CountData$fpkmA549 <- a549CountData$count/((a549CountData$exLength/1000)*(a549CountData$totalReads/1000000))

totalFpkm <- merge(crlCountData[,c('gene_id','fpkmCrl')],a549CountData[,c('gene_id','fpkmA549')])

write.table(totalFpkm,'~/Dropbox/densitypaper/BulkSeqData/CRL_A549_Bulk_Fpkm.txt',quote=F,sep='\t',row.names=F,col.names=T)

ggplot(totalFpkm,aes(x=fpkmCrl,y=fpkmA549)) + 
  geom_point() +
  scale_y_log10() + scale_x_log10() +
  geom_abline(slope=1,intercept=0,color='red')

highFpkm <- subset(totalFpkm, fpkmCrl > 5 & fpkmA549 > 5)
outlierData <- subset(highFpkm,fpkmA549 > 10*fpkmCrl | fpkmCrl > 10*fpkmA549)
ggplot(outlierData,aes(x=fpkmCrl,y=fpkmA549)) + 
  geom_point() +
  scale_y_log10() + scale_x_log10() +
  geom_abline(slope=1,intercept=0,color='red')
