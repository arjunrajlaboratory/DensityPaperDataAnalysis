setwd('Documents/ForeskinData/ForeskinData121206/')
annoTable <- read.delim('refseq_annotations_121206.gtf',header=F,stringsAsFactors=F)
colnames(annoTable) <- c(
  'seqname',
  'source',
  'feature',
  'start',
  'end',
  'score',
  'strand',
  'frame',
  'group'
)

groupVec <- annoTable$group
groupVec <- gsub('; $', '', groupVec)
groupList <- strsplit(groupVec, split='; ')
#all(unlist(lapply(groupList, length)) == 2)
groupTmpTable <- do.call(rbind, groupList)
annoTable$transcript_id <- gsub('transcript_id ', '', groupTmpTable[,2])
annoTable$gtf_idx <- 1:nrow(annoTable)

myGeneTranscriptPair <- read.delim('refflat_annotations_121209.txt',stringsAsFactors = F, header=T)
myGeneTranscriptPair <- myGeneTranscriptPair[,1:2]
colnames(myGeneTranscriptPair) <- c('gene_name','transcript_id')
myGeneTranscriptPair <- myGeneTranscriptPair[!duplicated(myGeneTranscriptPair$transcript_id),]
annoTable <- merge(myGeneTranscriptPair,annoTable)
annoTable$group <- paste('gene_id ', annoTable$gene_name, '; transcript_id ', annoTable$transcript_id, '; gene_name ', annoTable$gene_name, ';', sep='')
annoTable <- annoTable[order(annoTable$gtf_idx),]
outGtfTable <- annoTable[,c('seqname','source','feature','start','end','score','strand','frame','group')]
write.table(outGtfTable,'refseq_annotations_121206_for_cufflinks.gtf',quote=F,sep='\t',row.names=F,col.names=F)

#x <- c('a', 'b', 'c'); x %in% c('a', 'c')
#x[c(1,1,3, 2, 1)]
#order(order(y))
#y[order(y)]



#######

exonTable <- annoTable[annoTable$feature=="exon",]

bedTable <- data.frame(
  chrom=annoTable$seqname, 
  chromStart=annoTable$start,
  chromEnd=annoTable$end,
  name=annoTable$transcript_id,
  score=annoTable$score,
  strand=annoTable$strand,
  stringsAsFactors=F
)
#write.table(bedTable, 'my_first_bed_file_121206.bed', quote=F, sep='\t', row.names=F, col.names=F)

#bedTableShort = bedTable[1:50,]
#write.table(bedTableShort, 'my_first_bed_file_short_121206.bed', quote=F, sep='\t', row.names=F, col.names=F)
