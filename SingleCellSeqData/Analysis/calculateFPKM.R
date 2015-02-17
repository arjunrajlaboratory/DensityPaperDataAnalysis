library(plyr)
library(ggplot2)
library(data.table)

countData <- read.delim('~/Dropbox/densitypaper/SingleCellSeqData/singlecellseq_star_HTSeqCounts.tsv',header=T,stringsAsFactors=F)
lengthData <- read.delim('~/Dropbox/densitypaper/SingleCellSeqData/exon_lengths.txt',header=T,stringsAsFactors=F)
colnames(lengthData) <- c('gene_id','exLength')

data <- merge(countData,lengthData)
data$type <- NA
data$type[grep('Ctrl*',data$sampleID)] <- 'Ctrl'
data$type[grep('Fixed*',data$sampleID)] <- 'Fixed'
data$type[grep('Live*',data$sampleID)] <- 'Live'

data$rnaType <- 'Genomic'
data$rnaType[grep('ERCC-00*',data$gene_id)] <- 'ERCC'

totalReadTable <- ddply(data,.(sampleID),summarize,totalReads = sum(counts))
totalData <- merge(data,totalReadTable)

erccData <- data[grep('ERCC-00*',data$gene_id),]
erccReadTable <- ddply(erccData,.(sampleID),summarize,sumCountsErcc=sum(counts))
erccData <- merge(erccData,erccReadTable)

genomicData <- subset(data,!(gene_id %in% unique(erccData$gene_id)) )
genomicReadTable <- ddply(genomicData,.(sampleID),summarize,sumCountsGenomic=sum(counts))
genomicData <- merge(genomicData,genomicReadTable)

totalData$fpkm <- totalData$counts/((totalData$exLength/1000)*(totalData$totalReads/1000000))
erccData$fpkm <- erccData$counts/((erccData$exLength/1000)*(erccData$sumCountsErcc/1000000))
genomicData$fpkm <- genomicData$counts/((genomicData$exLength/1000)*(genomicData$sumCountsGenomic/1000000))

write.table(totalData,'~/Dropbox/densitypaper/SingleCellSeqData/totalDataFpkm.txt',quote=F,sep='\t',row.names=F,col.names=T)
write.table(genomicData,'~/Dropbox/densitypaper/SingleCellSeqData/genomicDataFpkm.txt',quote=F,sep='\t',row.names=F,col.names=T)
write.table(erccData,'~/Dropbox/densitypaper/SingleCellSeqData/erccDataFpkm.txt',quote=F,sep='\t',row.names=F,col.names=T)