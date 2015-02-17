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

freqDat <- subset(data, (gene %in% c('EEF2','UBC','TUSC3')) | (gene %in% 'MYC' & date %in% '131220') )

colE <- c('springgreen4','springgreen3','springgreen1')
colM <- c('deepskyblue4','deepskyblue3','deepskyblue1')
colU <- c('darkorchid4','darkorchid3','darkorchid1')
colR <- c('orangered4','orangered3','orangered2')
cols <- c('springgreen4','deepskyblue4','orangered4','darkorchid4')


## Frequency by stage

fbsData <- unique(freqDat[,c('dataNum','objNum','date','gene','normalizedTxnSites','stage')])
fbsData$normalizedTxnSites[fbsData$normalizedTxnSites>1] <- 1                ### Change >2 spots to 2!!

frequencyByStage <- ddply(fbsData,.(stage,gene),summarize,
                          meanFreq = mean(normalizedTxnSites),
                          seFreq = sd(normalizedTxnSites)/sqrt(length(normalizedTxnSites)))

frequencyByStage$stageHack[frequencyByStage$stage %in% 'G1'] <- 'first'
frequencyByStage$stageHack[frequencyByStage$stage %in% 'S'] <- 'second'
frequencyByStage$stageHack[frequencyByStage$stage %in% 'G2'] <- 'third'

dodge <- position_dodge(.6)

pdf('~/Dropbox/densitypaper/densitypaperdataanalysis/densityfigures/Figure5/freqbystage.pdf',width=2.5,height=2)
ggplot(frequencyByStage,aes(gene, meanFreq, fill=paste(gene,stageHack,sep='.'))) +
  geom_bar(stat='identity',position=dodge,width=.55) +
  geom_errorbar(aes(ymax = meanFreq+seFreq, ymin = meanFreq-seFreq,), width = .2, position = dodge) +
  xlab('Gene') + ylab('Transcription sites per allele') +
  scale_fill_manual(values=c(colE,colM,colR,colU)) +
  theme_bw() + theme_classic() + theme(legend.position='') +
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6')) +
  expand_limits(y=c(0,1)) +
  geom_hline(yintercept=frequencyByStage[1:4,3]/2)
dev.off()


# Frequency by volume

fbvData <- unique(freqDat[,c('dataNum','objNum','date','gene','normalizedTxnSites','xVar')])
fbvData$normalizedTxnSites[fbvData$normalizedTxnSites>1] <- 1                ### Change >2 spots to 2!!

qTable <- ddply(fbvData,'gene',summarize,
                q1=quantile(xVar,.25),
                q2=quantile(xVar,.5),
                q3=quantile(xVar,.75))

fbvData <- merge(fbvData,qTable)
fbvData$quartile <- 'q4'
fbvData$quartile[fbvData$xVar < fbvData$q3] <- 'q3'
fbvData$quartile[fbvData$xVar < fbvData$q2] <- 'q2'
fbvData$quartile[fbvData$xVar < fbvData$q1] <- 'q1'

frequencyByVolume <- ddply(fbvData,.(gene,quartile),summarize,
                           meanFreq = mean(normalizedTxnSites),
                           seFreq = sd(normalizedTxnSites)/sqrt(length(normalizedTxnSites)),
                           meanVol = mean(xVar),
                           seVol = sd(xVar)/sqrt(length(xVar)))

pdf('~/Dropbox/densitypaper/densitypaperdataanalysis/densityfigures/Figure5/freqbyvol.pdf',width=2,height=4)
ggplot(frequencyByVolume, aes(meanVol, meanFreq, col=gene)) +
  geom_point(size=1.3) +
  geom_errorbar(aes(ymax = meanFreq+seFreq, ymin = meanFreq-seFreq), width=75) +
  geom_errorbarh(aes(xmax = meanVol+seVol, xmin = meanVol-seVol),height=.02) +
  facet_wrap(~gene,scales='free_x') +
  expand_limits(y=c(0,1)) +
  xlab('GAPDH mRNA count') + ylab('Transcription sites per allele') +
  theme_classic() + 
  theme(legend.position = "none",strip.text = element_blank(),strip.background=element_blank()) + 
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6')) +
  scale_colour_manual(values=cols)
dev.off()
