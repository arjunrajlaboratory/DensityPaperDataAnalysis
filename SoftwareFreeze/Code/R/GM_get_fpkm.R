setwd('~/Documents/GM_Seq_tmp/')

data <- read.delim('feature_quantifications_GM12878_Seq_transcriptsonly',header=F,stringsAsFactors=F)
colnames(data) <- c(
  'seqname',
  'strand',
  'type',
  'location',
  'fpkm',
  'max',
  'Ucount',
  'NUcount',
  'length'
)

table1 <- data[,c(1,5)]

coltosplit <- table1[,1]
newcols <- strsplit(coltosplit,'\\({1}')
newcolsTable <- do.call(rbind, newcols)
table1$ref_id <- newcolsTable[,1]

setwd('~/Downloads/')

vegaTable <- read.delim('vega_annotations_hg19_130107.gtf',header=T,stringsAsFactors=F)
vegaNamesTable <- vegaTable[,c('name','name2')]
colnames(vegaNamesTable) <- c('ref_id','gene_name')

refseqTable <- read.delim('refseq_annotations_hg19_130107.gtf',header=T,stringsAsFactors=F)
refseqNamesTable <- refseqTable[,c('name','name2')]
colnames(refseqNamesTable) <- c('refseq_id','gene_name')

ucscTable <- read.delim('ucsc_to_refseq_hg19_130110.gtf',header=T,stringsAsFactors=F)
colnames(ucscTable) <- c('ucsc_id','refseq_id')
#ucscTable <- ucscTable[!duplicated(ucscTable$refseq_id),]

ucscRefseqTable <- merge(refseqNamesTable,ucscTable)
ucscTable <- ucscRefseqTable[,c('ucsc_id','gene_name')]

colnames(refseqNamesTable) <- c('ref_id','gene_name')
colnames(ucscTable) <- c('ref_id','gene_name')
catTable <- rbind(refseqNamesTable,ucscTable,vegaNamesTable)

outTable <- merge(table1,catTable)
finalTable <- outTable[,c(4,3)]

setwd('~/Documents/GM_Seq_tmp/')
write.table(finalTable,'GM_name_counts.txt',quote=F,sep='\t',row.names=F,col.names=T)
