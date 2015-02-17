# PCA vs LM for FISH GAPDH vs. Volume

library(plyr)
library(ggplot2)

fishData <- read.delim('~/Dropbox/densitypaper/ExtractedData/CRL.txt',header=T,stringsAsFactors=F)

fishDataForPCA <- fishData[,c('totalGAPDH','volume')]

covMat <- cov(fishDataForPCA)

eigenValues <- eigen(covMat)$values
eigenVectors <- eigen(covMat)$vectors

slopeFromEigen <- eigenVectors[2,1]/eigenVectors[1,1]

slopeForNonCentered <- slopeFromEigen
interceptForNonCentered <- -slopeFromEigen*mean(fishDataForPCA$totalGAPDH)+mean(fishDataForPCA$volume)

fitGapdhVsVolume <- lm(totalGAPDH ~ volume,data=fishDataForPCA)
slopeGapdhVsVolume <- coef(fitGapdhVsVolume)[[2]]
interceptGapdhVsVolume <- coef(fitGapdhVsVolume)[[1]]

fitVolumeVsGapdh <- lm(volume ~ totalGAPDH,data=fishDataForPCA)
slopeVolumeVsGapdh <- coef(fitVolumeVsGapdh)[[2]]
interceptVolumeVsGapdh <- coef(fitVolumeVsGapdh)[[1]]

volumeVsRnaPlot <- ggplot(fishDataForPCA,aes(x=totalGAPDH,y=volume)) + 
  geom_point(size=.5) +
  theme_classic() +
  expand_limits(x=0,y=0) +
  #geom_abline(slope=slopeForNonCentered,intercept=interceptForNonCentered,color='red') +
  #geom_abline(slope=1/slopeGapdhVsVolume,intercept=-interceptGapdhVsVolume/slopeGapdhVsVolume,color='blue') +
  geom_abline(slope=slopeVolumeVsGapdh,intercept=interceptVolumeVsGapdh,color='green') +
  ylab('Volume (picoliter)') +
  xlab('GAPDH mRNA count (FISH)') +
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6'))

ggplot(fishDataForPCA,aes(x=volume,y=totalGAPDH)) + 
  geom_point(size=.5) +
  theme_classic() +
  expand_limits(x=0,y=0) +
  geom_abline(slope=1/slopeForNonCentered,intercept=-interceptForNonCentered/slopeForNonCentered,color='red') +
  geom_abline(slope=slopeGapdhVsVolume,intercept=interceptGapdhVsVolume,color='blue') +
  #geom_abline(slope=slopeVolumeVsGapdh,intercept=interceptVolumeVsGapdh,color='green') +
  ylab('Volume (picoliter)') +
  xlab('GAPDH mRNA count (FISH)') +
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6'))


eef2Dat <- subset(fishData,gene %in% 'EEF2')

eef2Fit <- lm(cytoRNA ~ volume, data=eef2Dat)
eef2Slope <- coef(eef2Fit)[[2]]
eef2Int <- coef(eef2Fit)[[1]]

ggplot(eef2Dat,aes(x=volume,y=cytoRNA))