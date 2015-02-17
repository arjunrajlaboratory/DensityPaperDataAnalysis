library(ggplot2)
library(plyr)

data <- read.delim('~/Dropbox/densitypaper/densitypaperdataanalysis/ExtractedData/FreqIntens.txt',header=T,stringsAsFactors=F)

numericcols <- c('numTxnSites','numCyclin','xVar','intensity');

for (col in numericcols){
  data[[col]] <- as.numeric(data[[col]])
}

data$stage <- 'G1'
data$stage[data$numCyclin>20] <- 'S'
data$stage[data$numCyclin>230] <- 'G2'

data$stage[data$gene %in% 'TUSC3'] <- 'G1'
data$stage[data$numCyclin>40 & data$gene %in% 'TUSC3'] <- 'S'
data$stage[data$numCyclin>200 & data$gene %in% 'TUSC3'] <- 'G2'

data$stage[data$date %in% '141202'] <- 'G1'
data$stage[data$numCyclin>30 & data$date %in% '141202'] <- 'S'
data$stage[data$numCyclin>130 & data$date %in% '141202'] <- 'G2'

data$stageFactor <- 2
data$stageFactor[data$stage %in% 'S'] <- 4
data$stageFactor[data$stage %in% 'G2'] <- 4

data$stageFactor[(data$stage %in% 'S' & data$gene %in% 'TUSC3')] <- 2

data$normalizedTxnSites <- data$numTxnSites/data$stageFactor

intDat <- subset(data, (gene %in% 'EEF2' & date %in% '140114') | (gene %in% 'UBC' & date %in% '131220') |
                   (gene %in% 'TUSC3') | (gene %in% 'MYC' & date %in% '131220') )

colE <- c('springgreen4','springgreen3','springgreen1')
colM <- c('deepskyblue4','deepskyblue3','deepskyblue1')
colU <- c('darkorchid4','darkorchid3','darkorchid1')
colR <- c('orangered4','orangered3','orangered2')
cols <- c('springgreen4','deepskyblue4','orangered4','darkorchid4')


# Intensity by stage

ibspcData <- ddply(intDat,.(gene,dataNum,objNum,date,stage), function(df) {
  intensity <- df$intensity
  stage <- df$stage[1]
  if(stage %in% 'G1' & length(intensity)>2) {        ## if G1 cells have >2 txn sites, need to get rid of dimmest
    tmp <- intensity[order(-intensity)[1:2]]
    meanInt <- sum(tmp)/2
  }
  else {
    meanInt <- mean(intensity)
  }
  meanInt
})

colnames(ibspcData) <- c('gene','dataNum','objNum','date','stage','intPerCell')

intByStagePerCell <- ddply(ibspcData,.(gene,stage),summarize,
                           meanInt = mean(intPerCell),
                           seInt = sd(intPerCell)/sqrt(length(intPerCell)))

intByStagePerCell$stageHack[intByStagePerCell$stage %in% 'G1'] <- 'first'
intByStagePerCell$stageHack[intByStagePerCell$stage %in% 'S'] <- 'second'
intByStagePerCell$stageHack[intByStagePerCell$stage %in% 'G2'] <- 'third'
dodge <- position_dodge(.6)

pdf('~/Dropbox/densitypaper/densitypaperdataanalysis/densityfigures/Figure4/TxnIntensity/intensbystage_EEF2.pdf',width=1,height=2)
ggplot(subset(intByStagePerCell,gene %in% 'EEF2'),aes(paste(gene,stageHack,sep='.'), meanInt, fill=paste(gene,stageHack,sep='.'))) +
  geom_bar(stat='identity',position=dodge,width=.9) +
  geom_errorbar(aes(ymax = meanInt+seInt, ymin = meanInt-seInt,), width = .2, position = dodge) +
  xlab('Gene') + ylab('Mean transcription site intensity') +
  scale_fill_manual(values=colE) +
  theme_bw() + theme_classic() + theme(legend.position='') +
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6'))
dev.off()

pdf('~/Dropbox/densitypaper/densitypaperdataanalysis/densityfigures/Figure4/TxnIntensity/intensbystage_MYC.pdf',width=1,height=2)
ggplot(subset(intByStagePerCell,gene %in% 'MYC'),aes(paste(gene,stageHack,sep='.'), meanInt, fill=paste(gene,stageHack,sep='.'))) +
  geom_bar(stat='identity',position=dodge,width=.9) +
  geom_errorbar(aes(ymax = meanInt+seInt, ymin = meanInt-seInt,), width = .2, position = dodge) +
  xlab('Gene') + ylab('Mean transcription site intensity') +
  scale_fill_manual(values=colM) +
  theme_bw() + theme_classic() + theme(legend.position='') +
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6'))
dev.off()

pdf('~/Dropbox/densitypaper/densitypaperdataanalysis/densityfigures/Figure4/TxnIntensity/intensbystage_UBC.pdf',width=1,height=2)
ggplot(subset(intByStagePerCell,gene %in% 'UBC'),aes(paste(gene,stageHack,sep='.'), meanInt, fill=paste(gene,stageHack,sep='.'))) +
  geom_bar(stat='identity',position=dodge,width=.9) +
  geom_errorbar(aes(ymax = meanInt+seInt, ymin = meanInt-seInt,), width = .2, position = dodge) +
  xlab('Gene') + ylab('Mean transcription site intensity') +
  scale_fill_manual(values=colU) +
  theme_bw() + theme_classic() + theme(legend.position='') +
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6'))
dev.off()

pdf('~/Dropbox/densitypaper/densitypaperdataanalysis/densityfigures/Figure4/TxnIntensity/intensbystage_TUSC3.pdf',width=1,height=2)
ggplot(subset(intByStagePerCell,gene %in% 'TUSC3'),aes(paste(gene,stageHack,sep='.'), meanInt, fill=paste(gene,stageHack,sep='.'))) +
  geom_bar(stat='identity',position=dodge,width=.9) +
  geom_errorbar(aes(ymax = meanInt+seInt, ymin = meanInt-seInt,), width = .2, position = dodge) +
  xlab('Gene') + ylab('Mean transcription site intensity') +
  scale_fill_manual(values=colR) +
  theme_bw() + theme_classic() + theme(legend.position='') +
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6'))
dev.off()

# Intensity by volume

ibvpcData <- ddply(intDat,.(gene,date,dataNum,objNum,stage,xVar), function(df) {
  intensity <- df$intensity
  stage <- df$stage[1]
  if(stage %in% 'G1' & length(intensity)>2) {        ## if G1 cells have >2 txn sites, need to get rid of dimmest
    tmp <- intensity[order(-intensity)[1:2]]
    meanInt <- sum(tmp)/2
  }
  else {
    meanInt <- mean(intensity)
  }
  meanInt
})

colnames(ibvpcData) <- c('gene','date','dataNum','objNum','stage','xVar','intPerCell')

qTable <- ddply(ibvpcData,'gene',summarize,
                q1=quantile(xVar,.25),
                q2=quantile(xVar,.5),
                q3=quantile(xVar,.75))

ibvpcData <- merge(ibvpcData,qTable)
ibvpcData$quartile <- 'q4'
ibvpcData$quartile[ibvpcData$xVar < ibvpcData$q3] <- 'q3'
ibvpcData$quartile[ibvpcData$xVar < ibvpcData$q2] <- 'q2'
ibvpcData$quartile[ibvpcData$xVar < ibvpcData$q1] <- 'q1'

intByVolPerCell <- ddply(ibvpcData,.(quartile,gene),summarize,
                         meanInt = mean(intPerCell),
                         seInt = sd(intPerCell)/sqrt(length(intPerCell)),
                         meanVol = mean(xVar),
                         seVol = sd(xVar)/sqrt(length(xVar)))

pdf('~/Dropbox/densitypaper/densitypaperdataanalysis/densityfigures/Figure4/TxnIntensity/intensbyvol.pdf',width=2.25,height=4)
ggplot(intByVolPerCell, aes(meanVol, meanInt, col=gene)) +
  geom_point(size=1.3) +
  geom_errorbar(aes(ymax = meanInt+seInt, ymin = meanInt-seInt,), width = 75) +
  geom_errorbarh(aes(xmax = meanVol+seVol, xmin = meanVol-seVol,), height = 50) + 
  facet_wrap(~gene,scales='free') +
  expand_limits(y=c(0,1500)) +
  xlab('GAPDH mRNA count') + ylab('Mean transcription site intensity') +
  theme_classic() + 
  theme(legend.position = "none",strip.text = element_blank(),strip.background=element_blank()) + 
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6')) +
  scale_colour_manual(values=cols)
dev.off()