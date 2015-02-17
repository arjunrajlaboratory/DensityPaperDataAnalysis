# 
library(plyr)
library(ggplot2)

fishData <- read.delim('~/Dropbox/densitypaper/ExtractedData/CRL.txt',header=T,stringsAsFactors=F)
seqData <- read.delim('~/Dropbox/densitypaper/SingleCellSeqData/genomicDataFpkm.txt',header=T,stringsAsFactors=F)
seqData <- subset(seqData,gene_id %in% unique(fishData$gene))
seqData <- subset(seqData,type %in% 'Live')

volData <- read.delim('~/Dropbox/densitypaper/SingleCellSeqData/fpkmRatio.txt',header=T,stringsAsFactors=F)

seqData <- merge(seqData,volData)
seqData <- subset(seqData,fpkmRatio>30)

#pdf('~/Dropbox/Lab Meeting 140829/gapdh1.pdf',height=2,width=2.2)
ggplot(subset(seqData,gene_id %in% 'GAPDH'),aes(x=fpkmRatio,y=fpkm*fpkmRatio)) + geom_point(size=1) + theme_classic() + 
  xlab('\"volume\" = genomic/ercc count ratio') +
  ylab('\"count\" = fpkm * ratio') +
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6'))
#dev.off()


# normalize seq "volume" so that its mean is the same as GAPDH, measured by FISH
# check that normalized "volume" and total GAPDH have similar CDFs

seqData$gapdhEquiv <- seqData$fpkmRatio*mean(fishData$totalGAPDH)/mean(seqData$fpkmRatio)

seqTotalRna <- data.frame(unique(seqData$gapdhEquiv))
seqTotalRna$type <- 'seq'
colnames(seqTotalRna) <- c('totalRNA','type')

fishTotalRna <- data.frame(unique(fishData$totalGAPDH))
fishTotalRna$type <- 'fish'
colnames(fishTotalRna) <- c('totalRNA','type')

bigData <- rbind(seqTotalRna,fishTotalRna)

cdfPlot <- ggplot(bigData,aes(x=totalRNA,color=type)) + stat_ecdf() + 
  theme_classic() + theme(legend.position='none') +
  xlab('Total RNA per cell equivalent') +
  ylab('Cumulative density') +
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6'))

ggplot(bigData,aes(x=totalRNA,fill=type)) + geom_density(alpha=.25)
ggplot(bigData,aes(x=type,y=totalRNA)) + geom_boxplot()


# Get equation to find real volume, given "volume" (FPKM ratio) - use PCA

fishDataForPCA <- fishData[,c('totalGAPDH','volume')]

covMat <- cov(fishDataForPCA)

eigenValues <- eigen(covMat)$values
eigenVectors <- eigen(covMat)$vectors

slopeFromEigen <- eigenVectors[2,1]/eigenVectors[1,1]

slopeForNonCentered <- slopeFromEigen
interceptForNonCentered <- -slopeFromEigen*mean(fishDataForPCA$totalGAPDH)+mean(fishDataForPCA$volume)

volumeVsRnaPlot <- ggplot(fishDataForPCA,aes(x=totalGAPDH,y=volume)) + 
  geom_point(size=.5) +
  theme_classic() +
  expand_limits(x=0,y=0) +
  geom_abline(slope=slopeForNonCentered,intercept=interceptForNonCentered,color='red') +
  ylab('Volume (picoliter)') +
  xlab('GAPDH mRNA count (FISH)') +
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6'))

seqData$volume <- seqData$gapdhEquiv*slopeForNonCentered+interceptForNonCentered

ggplot(seqData,aes(x=volume,y=gapdhEquiv)) + geom_point()




# Now turn fpkm*fpkmRatio into actual counts, also using PCA for fit

meanSeqDat <- ddply(seqData,.(gene_id),summarise,
                    logMeanCountEquiv = log(mean(fpkm*fpkmRatio)),
                    log10MeanFpkm = log10(mean(fpkm)),
                    meanFpkm = mean(fpkm),
                    seFpkm = sd(fpkm)/sqrt(length(fpkm)))
colnames(meanSeqDat) <- c('gene','logMeanCountEquiv','log10MeanFpkm','meanFpkm','seFpkm')

meanFishDat <- ddply(fishData,.(gene),summarise,
                     logMeanFishCount = log(mean(totalRNA)),
                     log10MeanFishCount = log10(mean(totalRNA)),
                     meanFishCount = mean(totalRNA),
                     seCount = sd(totalRNA)/sqrt(length(totalRNA)))

mergedMeanDat <- merge(meanSeqDat,meanFishDat)

fishVsSeqForPCA <- mergedMeanDat[,c('logMeanCountEquiv','logMeanFishCount')]

sl <- coef(lm(data=mergedMeanDat,log10MeanFpkm ~ log10MeanFishCount))[[2]]
intc <- coef(lm(data=mergedMeanDat,log10MeanFpkm ~ log10MeanFishCount))[[1]]


pcaForFig <- mergedMeanDat[,c('log10MeanFishCount','log10MeanFpkm')]
covMat <- cov(pcaForFig)

eigenValues <- eigen(covMat)$values
eigenVectors <- eigen(covMat)$vectors

slopeFromEigen <- eigenVectors[2,1]/eigenVectors[1,1]

slopeFig <- slopeFromEigen
interceptFig <- -slopeFromEigen*mean(pcaForFig$log10MeanFishCount)+mean(pcaForFig$log10MeanFpkm)


fishVsFpkm <- ggplot(mergedMeanDat,aes(y=meanFpkm,x=meanFishCount)) + geom_point(size=1) + 
  scale_x_log10() + scale_y_log10() +
  theme_classic() +
  #geom_abline(slope=sl,intercept=intc) +
  geom_errorbar(aes(ymin=meanFpkm-seFpkm,ymax=meanFpkm+seFpkm)) +
  geom_errorbarh(aes(xmin=meanFishCount-seCount,xmax=meanFishCount+seCount)) +
  #expand_limits(x=1,y=1) +
  ylab('Mean FPKM') +
  xlab('Mean FISH Count') +
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6'))

fishVsFpkm2 <- ggplot(mergedMeanDat,aes(y=log10MeanFpkm,x=log10MeanFishCount)) + geom_point(size=1) + 
  theme_classic() +
  geom_abline(slope=slopeFig,intercept=interceptFig) +
  expand_limits(x=1,y=1) +
  ylab('Log Mean FPKM') +
  xlab('Log Mean FISH Count') +
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6'))



ggplot(fishVsSeqForPCA,aes(x=logMeanCountEquiv,y=logMeanFishCount)) + geom_point()

ggplot(fishVsSeqForPCA,aes(x=logMeanCountEquiv,y=logMeanFishCount,label=mergedMeanDat$gene)) + geom_text()

covMat <- cov(fishVsSeqForPCA)

eigenValues <- eigen(covMat)$values
eigenVectors <- eigen(covMat)$vectors

slopeFromEigen <- eigenVectors[2,1]/eigenVectors[1,1]

slopeForNonCenteredFishVsSeq <- slopeFromEigen
interceptForNonCenteredFishVsSeq <- -slopeFromEigen*mean(fishVsSeqForPCA$logMeanCountEquiv)+mean(fishVsSeqForPCA$logMeanFishCount)

ggplot(fishVsSeqForPCA,aes(x=logMeanCountEquiv,y=logMeanFishCount)) + geom_point() +
  geom_abline(slope=slopeForNonCenteredFishVsSeq,intercept=interceptForNonCenteredFishVsSeq,color='blue')

ggplot(fishVsSeqForPCA,aes(x=logMeanCountEquiv,y=logMeanFishCount,label=mergedMeanDat$gene)) + geom_text() + geom_point(color='red') +
  geom_abline(slope=slopeForNonCenteredFishVsSeq,intercept=interceptForNonCenteredFishVsSeq,color='blue')

seqData$actualCount <- exp(slopeForNonCenteredFishVsSeq*log(seqData$fpkm*seqData$fpkmRatio)+interceptForNonCenteredFishVsSeq)

#write.table(seqData,'~/Dropbox/densitypaper/SingleCellSeqData/seqData_volume_and_counts_fishGenes.txt')

summarizeSeq <- ddply(seqData,.(gene_id),summarise,
                      meanActualCount=mean(actualCount),
                      meanFpkm = mean(fpkm))
colnames(summarizeSeq) <- c('gene','meanActualCount','meanFpkm')
summarizeFish <- ddply(fishData,.(gene),summarise,
                       meanTotalRna = mean(totalRNA))

compDat <- merge(summarizeSeq,summarizeFish)

ggplot(compDat,aes(x=meanTotalRna,y=meanActualCount)) + geom_point() + scale_x_log10() + scale_y_log10()
ggplot(compDat,aes(x=meanTotalRna,y=meanFpkm)) + geom_point() + scale_x_log10() + scale_y_log10()
ggplot(compDat,aes(x=meanActualCount,y=meanFpkm)) + geom_point() + scale_x_log10() + scale_y_log10()





# do this for all genes, not just FISHed genes
seqDataTotal <- read.delim('~/Dropbox/densitypaper/SingleCellSeqData/genomicDataFpkm.txt',header=T,stringsAsFactors=F)
seqDataTotal <- subset(seqDataTotal,type %in% 'Live')
seqDataTotal <- merge(seqDataTotal,volData)
seqDataTotal <- subset(seqDataTotal,fpkmRatio>30)

seqDataTotal$gapdhEquiv <- seqDataTotal$fpkmRatio*mean(fishData$totalGAPDH)/mean(seqDataTotal$fpkmRatio)

seqDataTotal$volume <- seqDataTotal$gapdhEquiv*slopeForNonCentered+interceptForNonCentered
seqDataTotal$actualCount <- exp(slopeForNonCenteredFishVsSeq*log(seqDataTotal$fpkm*seqDataTotal$fpkmRatio)+interceptForNonCenteredFishVsSeq)

#write.table(seqDataTotal,'~/Dropbox/densitypaper/SingleCellSeqData/seqData_volume_and_counts_allGenes.txt')

summarizeSeq <- ddply(seqDataTotal,.(gene_id),summarise,
                      meanActualCount=mean(actualCount),
                      meanFpkm = mean(fpkm))

ggplot(summarizeSeq,aes(x=meanActualCount,y=meanFpkm)) + geom_point() + scale_x_log10() + scale_y_log10()




#pdf(file='~/Dropbox/densitypaper/SupplementaryFigures/SingleCellSeq_CalibrateData/cdfPlot.pdf',height=3,width=3.3)
cdfPlot
#dev.off()

#pdf(file='~/Dropbox/densitypaper/SupplementaryFigures/SingleCellSeq_CalibrateData/volumeVsRnaPlot.pdf',height=3,width=3.3)
volumeVsRnaPlot
#dev.off()

pdf(file='~/Dropbox/densitypaper/SupplementaryFigures/SingleCellSeq_CalibrateData/fishVsFpkmPlot.pdf',height=3,width=3.3)
fishVsFpkm
dev.off()

pdf(file='~/Dropbox/densitypaper/SupplementaryFigures/SingleCellSeq_CalibrateData/fishVsFpkmPlot2.pdf',height=3,width=3.3)
fishVsFpkm2
dev.off()