library(ggplot2)
library(plyr)
library(reshape)

genomicData <- read.delim('~/Dropbox/densitypaper/densitypaperdataanalysis/SingleCellSeqData/genomicDataFpkm.txt',header=T,stringsAsFactors=F)
genomicDataSub <- genomicData[,c('sampleID','type','sumCountsGenomic')]
genomicDataSub <- unique(genomicDataSub)
erccData <- read.delim('~/Dropbox/densitypaper/densitypaperdataanalysis/SingleCellSeqData/erccDataFpkm.txt',header=T,stringsAsFactors=F)
erccDataSub <- erccData[,c('sampleID','type','sumCountsErcc')]
erccDataSub <- unique(erccDataSub)

data <- merge(genomicDataSub,erccDataSub)
data <- subset(data,type %in% 'Live' | type %in% 'Ctrl')
data$category <- 'greaterThan30'
data$category[data$sumCountsGenomic/data$sumCountsErcc<30] <- 'lessThan30'

ggplot(data,aes(x=sumCountsGenomic,y=sumCountsErcc,color=type)) + geom_point()
ggplot(subset(data,type %in% 'Live'),aes(x=sumCountsGenomic,y=sumCountsErcc,color=type)) + geom_point()



erccGenomicCounts <- ggplot(subset(data,type %in% 'Live'),aes(x=sumCountsGenomic,y=sumCountsErcc,color=category)) + 
  geom_point(size=1.1) +
  theme_classic() + theme(legend.position='none') +
  xlab('Summed genomic counts') + ylab('Summed ERCC counts') +
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6'))


hist(data$sumCountsGenomic/data$sumCountsErcc)

volumeData <- read.delim('~/Dropbox/densitypaper/densitypaperdataanalysis/SingleCellSeqData/fpkmRatio.txt',header=T,stringsAsFactors=F)
genomicWithVolume <- merge(genomicData,volumeData)
genomicWithVolume$category <- 'greaterThan30'
genomicWithVolume$category[genomicWithVolume$fpkmRatio<30] <- 'lessThan30'
gapdhData <- subset(genomicWithVolume,gene_id %in% 'GAPDH' & type %in% 'Live')


gapdhWithFunkyCells <- ggplot(gapdhData,aes(x=fpkmRatio,y=fpkm*fpkmRatio,color=category)) + 
  geom_point(size=1.1) +
  theme_classic() + theme(legend.position='none') +
  xlab('ERCC/Genomic count ratio ("volume")') + ylab('GAPDH FPKM * count ratio ("count")') +
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6'))
  






# Just use FPKM - not sure if we can presume to measure count and volume accurately for the funky cells
genomicWithVolumeLive <- subset(genomicWithVolume,type %in% 'Live')

goodCell <- 'Live-PL2-F06'
badCell <- 'Live-PL2-G07'

countAndVolMean <- ddply(genomicWithVolumeLive,.(gene_id,category),summarize,
                         meanFpkm=mean(fpkm),
                         meanCount=mean(counts))

goodCellsMean <- subset(countAndVolMean,category %in% 'greaterThan30')
goodCellsMean <- goodCellsMean[,c(1,3,4)]
colnames(goodCellsMean) <- c('gene','meanFpkmGood','meanCountGood')
weirdCellsMean <- subset(countAndVolMean,category %in% 'lessThan30')
weirdCellsMean <- weirdCellsMean[,c(1,3,4)]
colnames(weirdCellsMean) <- c('gene','meanFpkmWeird','meanCountWeird')

meanData <- merge(goodCellsMean,weirdCellsMean)

weirdVsNormal <- ggplot(meanData,aes(x=meanFpkmGood,y=meanFpkmWeird)) + 
  geom_point(size=.4) + scale_x_log10() + scale_y_log10() +
  theme_classic() + theme(legend.position='none') +
  expand_limits(x=0,y=0) +
  xlab('Mean FPKM "high ratio" cells') + ylab('Mean FPKM "low ratio" cells') +
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6'))

# compare a single cell to the "good" mean

goodSingleCell <- subset(genomicWithVolumeLive,sampleID %in% goodCell)
goodSingleCell <- goodSingleCell[,c('gene_id','counts','fpkm')]
colnames(goodSingleCell) <- c('gene','singleCountGood','singleFpkmGood')
weirdSingleCell <- subset(genomicWithVolumeLive,sampleID %in% badCell)
weirdSingleCell <- weirdSingleCell[,c('gene_id','counts','fpkm')]
colnames(weirdSingleCell) <- c('gene','singleCountWeird','singleFpkmWeird')

compData <- merge(meanData,goodSingleCell)
compData <- merge(compData,weirdSingleCell)

ggplot(compData,aes(x=meanFpkmGood,y=singleFpkmGood)) + geom_point() + scale_x_log10() + scale_y_log10()
ggplot(compData,aes(x=meanFpkmGood,y=singleFpkmWeird)) + geom_point() + scale_x_log10() + scale_y_log10()





pdf(file='~/Dropbox/densitypaper/densitypaperdataanalysis/SupplementaryFigures/SingleCellSeq_ClassifyWeirdCells/gapdhWithFunkyCells.pdf',height=2,width=2.2)
gapdhWithFunkyCells
dev.off()

pdf(file='~/Dropbox/densitypaper/densitypaperdataanalysis/SupplementaryFigures/SingleCellSeq_ClassifyWeirdCells/erccGenomicCounts.pdf',height=2,width=2.2)
erccGenomicCounts
dev.off()