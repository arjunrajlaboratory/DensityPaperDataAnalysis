library(ggplot2)
library(plyr)

data <- read.delim('~/Dropbox/densitypaper/densitypaperdataanalysis/ExtractedData/triptolide.txt',header=T,stringsAsFactors=F)
data <- subset(data,numTxnSites>0)

tmp <- subset(data,date %in% '141105')

ggplot(tmp, aes(x=intensity,fill=conc)) + geom_density(alpha=0.5)

intensityDat <- ddply(tmp,.(conc,dataNum,objNum),summarise,
                      meanInt = mean(intensity),
                      maxInt = max(intensity),
                      numTxnSites = max(numTxnSites))



tmp <- subset(data,gene %in% 'MYC')

pdf('~/Dropbox/densitypaper/densitypaperdataanalysis/densityfigures/Figure4/Triptolide/MYC.pdf',height=2.5,width=2)
ggplot(tmp, aes(x=intensity/1000,fill=conc)) + geom_density(alpha=0.5) + 
  scale_fill_manual(values = c('deepskyblue4','deepskyblue1')) +
  theme_classic() + theme(legend.position='none') +
  xlab('Intensity of each \'ON\' transcription site') + ylab('Density') +
  geom_vline(xintercept=10) +
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6'))
dev.off()

thresh <- 10000
ks.test(subset(tmp,conc %in% 'Ctrl' & intensity>thresh)$intensity,subset(tmp,conc %in% '100nM' & intensity>thresh)$intensity)
ks.test(subset(tmp,conc %in% 'Ctrl')$intensity,subset(tmp,conc %in% '100nM')$intensity)

ctrl <- subset(tmp,conc %in% 'Ctrl')
kd <- subset(tmp,conc %in% '100nM')

numBelowThreshCtrl <- length(subset(ctrl,intensity<thresh)$intensity)
numAboveThreshCtrl <- length(subset(ctrl, intensity>thresh)$intensity)
percentAboveThreshCtrl <- numAboveThreshCtrl/numBelowThreshCtrl

numBelowThreshKd <- length(subset(kd, intensity<thresh)$intensity)
numAboveThreshKd <- length(subset(kd, intensity>thresh)$intensity)
percentAboveThreshKd <- numAboveThreshKd/numBelowThreshKd

allIntensities <- tmp$intensity
numCtrlIntens <- dim(ctrl)[1]
numKdIntens <- dim(kd)[1]

percentAboveKd <- NA
percentAboveCtrl <- NA

for(i in 1:10000) {
  sampleNums <- sample(length(allIntensities),numCtrlIntens,replace=F)
  ctrlSample <- allIntensities[sampleNums]
  kdSample <- allIntensities[setdiff(1:length(allIntensities),sampleNums)]
  
  numBelowThreshKd <- sum(kdSample<thresh)
  numAboveThreshKd <- sum(kdSample>thresh)
  percentAboveKd[i] <- numAboveThreshKd/numBelowThreshKd
  
  numBelowThreshCtrl <- sum(ctrlSample<thresh)
  numAboveThreshCtrl <- sum(ctrlSample>thresh)
  percentAboveCtrl[i] <- numAboveThreshCtrl/numBelowThreshCtrl
}

pValKdMyc <- sum(percentAboveKd<percentAboveThreshKd)/sum(percentAboveKd)
pValCtrlMyc <- sum(percentAboveCtrl>percentAboveThreshCtrl)/sum(percentAboveCtrl)

#############################

tmp <- subset(data,gene %in% 'UBC' & conc %in% c('Ctrl','100nM'))

pdf('~/Dropbox/densitypaper/densitypaperdataanalysis/densityfigures/Figure4/Triptolide/UBC.pdf',height=2.5,width=2)
ggplot(tmp, aes(x=intensity/1000,fill=conc)) + geom_density(alpha=0.5) + 
  scale_fill_manual(values = c('darkslateblue','darkorchid1')) +
  theme_classic() + theme(legend.position='none') +
  xlab('Intensity of each \'ON\' transcription site') + ylab('Density') +
  geom_vline(xintercept=5) +
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6'))
dev.off()

thresh <- 5000
ks.test(subset(tmp,conc %in% 'Ctrl' & intensity>thresh)$intensity,subset(tmp,conc %in% '100nM' & intensity>thresh)$intensity)
ks.test(subset(tmp,conc %in% 'Ctrl')$intensity,subset(tmp,conc %in% '100nM')$intensity)

ctrl <- subset(tmp,conc %in% 'Ctrl')
kd <- subset(tmp,conc %in% '100nM')

numBelowThreshCtrl <- length(subset(ctrl,intensity<thresh)$intensity)
numAboveThreshCtrl <- length(subset(ctrl, intensity>thresh)$intensity)
percentAboveThreshCtrl <- numAboveThreshCtrl/numBelowThreshCtrl

numBelowThreshKd <- length(subset(kd, intensity<thresh)$intensity)
numAboveThreshKd <- length(subset(kd, intensity>thresh)$intensity)
percentAboveThreshKd <- numAboveThreshKd/numBelowThreshKd

allIntensities <- tmp$intensity
numCtrlIntens <- dim(ctrl)[1]
numKdIntens <- dim(kd)[1]

percentAboveKd <- NA
percentAboveCtrl <- NA

for(i in 1:10000) {
  sampleNums <- sample(length(allIntensities),numCtrlIntens,replace=F)
  ctrlSample <- allIntensities[sampleNums]
  kdSample <- allIntensities[setdiff(1:length(allIntensities),sampleNums)]
  
  numBelowThreshKd <- sum(kdSample<thresh)
  numAboveThreshKd <- sum(kdSample>thresh)
  percentAboveKd[i] <- numAboveThreshKd/numBelowThreshKd
  
  numBelowThreshCtrl <- sum(ctrlSample<thresh)
  numAboveThreshCtrl <- sum(ctrlSample>thresh)
  percentAboveCtrl[i] <- numAboveThreshCtrl/numBelowThreshCtrl
}

pValKdUbc <- sum(percentAboveKd<percentAboveThreshKd)/sum(percentAboveKd)
pValCtrlUbc <- sum(percentAboveCtrl>percentAboveThreshCtrl)/sum(percentAboveCtrl)

pValCtrlMyc
pValKdMyc