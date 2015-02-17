source('~/Dropbox/densitypaper/densitypaperdataanalysis/multiplot.R')
library(ggplot2)
library(plyr)

data <- read.delim('~/Dropbox/densitypaper/densitypaperdataanalysis/ExtractedData/si.txt',header=T,stringsAsFactors=F)

numDat <- unique(data[,c('date','type','dataNum','objNum','numRNA')])
numSummarize <- ddply(numDat,.(type),summarize,
                      meanRNA = mean(numRNA),
                      sterrRNA = sd(numRNA)/sqrt(length(numRNA)))

cols <- c('deepskyblue3','deepskyblue4')


p1 <- ggplot(numSummarize,aes(x=type,y=meanRNA,fill=type)) +
  geom_bar(stat='identity') +
  geom_errorbar(aes(ymin=meanRNA-sterrRNA,ymax=meanRNA+sterrRNA),width=.25) +
  expand_limits(x=0,y=c(0,max(numSummarize$meanRNA)*(4/3))) +
  theme_classic() +
  scale_fill_manual(values=cols) +
  xlab('') + ylab('LMNA mRNA') +
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6'),legend.position='') +
  theme(axis.text.y = element_text(angle = 90, hjust = 1))

numSitesDat <- unique(data[,c('date','type','dataNum','objNum','numTxnSites')])
numSitesSummarize <- ddply(numSitesDat,.(type),summarize,
                           meanNum = mean(numTxnSites),
                           sterrNum = sd(numTxnSites)/sqrt(length(numTxnSites)))

p2 <- ggplot(numSitesSummarize,aes(x=type,y=meanNum,fill=type)) +
  geom_bar(stat='identity') +
  geom_errorbar(aes(ymin=meanNum-sterrNum,ymax=meanNum+sterrNum),width=.25) +
  expand_limits(x=0,y=c(0,max(numSitesSummarize$meanNum)*(4/3))) +
  theme_classic() +
  scale_fill_manual(values=cols) +
  xlab('') + ylab('Txn sites') +
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6'),legend.position='') +
  theme(axis.text.y = element_text(angle = 90, hjust = 1))


intensityDat <- ddply(subset(data,numTxnSites>0),.(date,type,dataNum,objNum),summarize,
                            intPerCell = mean(intensity))
intensitySummarize <- ddply(intensityDat,.(type),summarize,
                            meanI = mean(intPerCell),
                            sterrI = sd(intPerCell)/sqrt(length(intPerCell)))


p3 <- ggplot(intensitySummarize,aes(x=type,y=meanI,fill=type)) +
  geom_bar(stat='identity') +
  geom_errorbar(aes(ymin=meanI-sterrI,ymax=meanI+sterrI),width=.25) +
  expand_limits(x=0,y=c(0,max(intensitySummarize$meanI)*(4/3))) +
  theme_classic() +
  scale_fill_manual(values=cols) +
  xlab('') + ylab('Txn site intensity') +
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6'),legend.position='') +
  theme(axis.text.y = element_text(angle = 90, hjust = 1))

pdf('~/Dropbox/densitypaper/densitypaperdataanalysis/densityfigures/Figure2/Knockdown/si.pdf',width=3.9,height=1.75)
multiplot(p1,p2,p3,cols=3)
dev.off()