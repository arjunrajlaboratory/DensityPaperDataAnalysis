library(plyr)
library(ggplot2)
library(data.table)

data <- read.delim('~/Dropbox/densitypaper/SingleCellSeqData/refSeq_2014-06-17_refGene_ERCC.gtf',header=F,stringsAsFactors=F)
colnames(data) <- c(
  'chr',
  'source',
  'feature',
  'start',
  'end',
  'score',
  'strand',
  'frame',
  'group'
)

dataEx <- subset(data,feature %in% 'exon')

tosplit <- dataEx$group

groupVec <- gsub('; $', '', tosplit)
groupList <- strsplit(groupVec, split='; ')
groupTmpTable <- do.call(rbind, groupList)
dataEx$transcript_id <- gsub('transcript_id ', '', groupTmpTable[,2])
dataEx$gene_name <- gsub('gene_id', '', groupTmpTable[,1])

trim <- function (x) gsub("^\\s+|\\s+$", "", x)

dataEx$gene_name <- trim(dataEx$gene_name)

calcTotalExonLength <- function(start,end) {
  #subDat <- data[,c('start','end')]
  subDat <- cbind(start,end)
  subDat <- data.frame(subDat)
  colnames(subDat) <- c('start','end')
  
  orderDat <- subDat[order(subDat$start),]
  orderDat <- unique(orderDat)
  
  geneLength <- 0
  for (i in 1:nrow(orderDat)) {
    if (i == 1) {
      geneLength <- orderDat[1,2] - orderDat[1,1]
    } else if (orderDat[i,1] >= orderDat[i-1,2]) { # this exon does not overlap previous at all
      geneLength <- geneLength + (orderDat[i,2] - orderDat[i,1])
    } else if (orderDat[i,2] > orderDat[i-1, 2]) { # this exon starts in the middle of previous, but extends beyond it
      geneLength <- geneLength + (orderDat[i,2] - orderDat[i-1,2])
    }
  }
  
  return(geneLength)
}

dataExTable <- data.table(dataEx)
exonLengthTable <- dataExTable[ , list(exontotal = calcTotalExonLength(start,end)), by = 'gene_name']
exonLengthTable <- as.data.frame(exonLengthTable)

write.table(exonLengthTable,'~/Dropbox/densitypaper/SingleCellSeqData/exon_lengths.txt',quote=F,sep='\t',row.names=F,col.names=T)
