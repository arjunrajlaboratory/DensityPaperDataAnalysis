library(plyr)
library(ggplot2)
library(data.table)
library(reshape)

data <- read.delim('~/Dropbox/densitypaper/SingleCellSeqData/totalDataFpkm.txt',header=T,stringsAsFactors=F)

data <- subset(data,totalReads>1000000)

erccData <- read.delim('~/Dropbox/singleCellSeqData/singlecellseq/analysis/erccDataFpkm.txt',header=T,stringsAsFactors=F)
erccCounts <- unique(erccData[,c('sampleID','sumCountsErcc')])
genomicData <- read.delim('~/Dropbox/singleCellSeqData/singlecellseq/analysis/genomicDataFpkm.txt',header=T,stringsAsFactors=F)
genomicCounts <- unique(genomicData[,c('sampleID','sumCountsGenomic')])

countsTable <- merge(erccCounts,genomicCounts)
countsTable$countRatio <- countsTable$sumCountsGenomic/countsTable$sumCountsErcc

fpkmTable <- data.table(data)
totalFpkmTable <- fpkmTable[ , list(totalFpkm = sum(fpkm)), by = c('sampleID','rnaType','type')]
totalFpkmTable <- data.frame(totalFpkmTable)

genomicFpkmTable <- subset(totalFpkmTable,rnaType %in% 'Genomic')
colnames(genomicFpkmTable) <- c('sampleID','rnaType','type','genomicFpkm')

erccFpkmTable <- subset(totalFpkmTable,rnaType %in% 'ERCC')
colnames(erccFpkmTable) <- c('sampleID','rnaType','type','erccFpkm')

totalFpkmTable <- merge(genomicFpkmTable[,c(1,3,4)],erccFpkmTable[,c(1,3,4)])

totalFpkmTable$fpkmRatio <- totalFpkmTable$genomicFpkm/totalFpkmTable$erccFpkm

totalFpkmTableToPrint <- totalFpkmTable[,c('sampleID','fpkmRatio')]
write.table(totalFpkmTableToPrint,'~/Dropbox/densitypaper/SingleCellSeqData/fpkmRatio.txt',quote=F,sep='\t',row.names=F,col.names=T)