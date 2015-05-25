setwd('/Volumes/Untitled/')
data_noDrug <- read.delim('NoDrug/feature_quantifications_NoDrug_transcriptsonly',header=F,stringsAsFactors=F)
colnames(data_noDrug) <- c(
  'seqname',
  'strand',
  'type',
  'location',
  'min',
  'max',
  'Ucount',
  'NUcount',
  'length'
)

data_1uM1Day <- read.delim('1uM1Day/feature_quantifications_1uM1Day_transcriptsonly',header=F,stringsAsFactors=F)
colnames(data_1uM1Day) <- c(
  'seqname',
  'strand',
  'type',
  'location',
  'min',
  'max',
  'Ucount',
  'NUcount',
  'length'
)

data_1uM5Day <- read.delim('1uM5Day/feature_quantifications_1uM5Day_transcriptsonly',header=F,stringsAsFactors=F)
colnames(data_1uM5Day) <- c(
  'seqname',
  'strand',
  'type',
  'location',
  'min',
  'max',
  'Ucount',
  'NUcount',
  'length'
)

data_50uM1Day <- read.delim('50uM1Day/feature_quantifications_50uM1Day_transcriptsonly',header=F,stringsAsFactors=F)
colnames(data_50uM1Day) <- c(
  'seqname',
  'strand',
  'type',
  'location',
  'min',
  'max',
  'Ucount',
  'NUcount',
  'length'
)

data_50uM5Day <- read.delim('50uM5Day/feature_quantifications_50uM5Day_transcriptsonly',header=F,stringsAsFactors=F)
colnames(data_50uM5Day) <- c(
  'seqname',
  'strand',
  'type',
  'location',
  'min',
  'max',
  'Ucount',
  'NUcount',
  'length'
)

input_to_DESeq <- data_noDrug[,c('seqname','Ucount')]
input_to_DESeq[,3] <- data_1uM1Day[,c('Ucount')]
input_to_DESeq[,4] <- data_1uM5Day[,c('Ucount')]
input_to_DESeq[,5] <- data_50uM1Day[,c('Ucount')]
input_to_DESeq[,6] <- data_50uM5Day[,c('Ucount')]
colnames(input_to_DESeq) <- c(
  'seqname',
  'NoDrug',
  '1uM1Day',
  '1uM5Day',
  '50uM1Day',
  '50uM5Day'
)

write.table(input_to_DESeq,'counts_five_conditions.txt',quote=F,sep='\t',row.names=F,col.names=T)

source("http://www.bioconductor.org/biocLite.R")
biocLite("DESeq")
library( "DESeq" )

DESeq_table <- read.delim('counts_five_conditions.txt',header=T,row.names=1)
condition <- factor(c('NoDrug','1uM1Day','1uM5Day','50uM1Day','50uM5Day'))
cds <- newCountDataSet(DESeq_table,condition)
cds <- estimateSizeFactors( cds )
sizeFactors( cds )
head( counts( cds, normalized=TRUE ) )
cds = estimateDispersions( cds , method='blind',sharingMode='fit-only')
str( fitInfo(cds) )
plotDispEsts( cds )
head( fData(cds) )

res = nbinomTest( cds, "NoDrug", "50uM1Day" )
plotMA(res)
hist(res$pval, breaks=100, col="skyblue", border="slateblue", main="")
resSig = res[ res$padj < 0.1, ]
head( resSig[ order(resSig$pval), ] )
most_sig <- resSig[ order(resSig$pval), ]
most_up <- resSig[ order( -resSig$foldChange, -resSig$baseMean ), ]
most_down <- resSig[ order( resSig$foldChange, -resSig$baseMean ), ]
write.table(most_sig,'noDrug_50uM1Day_mostsig.txt',quote=F,sep='\t',row.names=F,col.names=T)

table1 <- read.delim('noDrug_1uM1Day_mostsig.txt',header=T,stringsAsFactors=F)
coltosplit <- table1[,1]
newcols <- strsplit(coltosplit,'\\({1}')
newcolsTable <- do.call(rbind, newcols)
table1$ref_id <- newcolsTable[,1]

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

x <- read.delim('NoDrug/feature_quantifications_NoDrug_transcriptsonly',header=F,stringsAsFactors=F)
colnames(x) <- c(
  'seqname',
  'strand',
  'type',
  'location',
  'min',
  'max',
  'Ucount',
  'NUcount',
  'length'
)
coltosplit <- x[,1]
newcols <- strsplit(coltosplit,'\\({1}')
newcolsTable <- do.call(rbind, newcols)
x$ref_id <- newcolsTable[,1]
xtable <- merge(x,catTable)