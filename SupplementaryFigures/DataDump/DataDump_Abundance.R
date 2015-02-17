library(data.table)
library(ggplot2)
library(plyr)

crlDat <- read.delim('~/Dropbox/densitypaper/ExtractedData/CRL.txt',header=T,stringsAsFactors=F)
crlDat$type <- 'cycling'

ssDat <- read.delim('~/Dropbox/densitypaper/ExtractedData/SerumStarved.txt',header=T,stringsAsFactors=F)
ssDat$type <- 'quiescent'

senDat <- read.delim('~/Dropbox/densitypaper/ExtractedData/Senescent.txt',header=T,stringsAsFactors=F)
#senDat <- subset(senDat,!(date %in% 140128))
senDat$type <- 'senescent'

a549Dat <- read.delim('~/Dropbox/densitypaper/ExtractedData/A549.txt',header=T,stringsAsFactors=F)
a549Dat$type <- 'a549'

numericcols <- c('volume','area','totalGAPDH','cytoGAPDH','nucGAPDH','nucArea',
                 'nucIntensityTotal','nucIntensityAvg','numCyclin','totalRNA','cytoRNA',
                 'nucRNA','numTxnSites','avgIntron','avgExon','numRnaPerTxnSiteAvgExon');


for (col in numericcols){
  crlDat[[col]] <- as.numeric(crlDat[[col]])
  ssDat[[col]] <- as.numeric(ssDat[[col]])
  senDat[[col]] <- as.numeric(senDat[[col]])
  a549Dat[[col]] <- as.numeric(a549Dat[[col]])
}

crlDat$volume <- crlDat$volume/1000
ssDat$volume <- ssDat$volume/1000
senDat$volume <- senDat$volume/1000
a549Dat$volume <- a549Dat$volume/1000

pdf('~/Dropbox/densitypaper/SupplementaryFigures/DataDump/CRL_cycling.pdf',width=7.5,height=4)
ggplot(crlDat, 
       aes(x=volume,y=cytoRNA)) + 
  geom_point(size=1) + facet_wrap(~gene, scales = 'free_y') + expand_limits(x=0,y=0) + 
  stat_smooth(method="lm", se=FALSE, fullrange = TRUE,col='red') +
  theme_classic() +
  theme(axis.title=element_text(size="10"), axis.text=element_text(size='10')) +
  xlab('Volume (picoliter)') + ylab('mRNA Abundance') +
  #theme(strip.background = element_blank())
  theme(strip.background = element_rect(fill='white'),strip.text = element_text(face='italic',size=10))
dev.off()

pdf('~/Dropbox/densitypaper/SupplementaryFigures/DataDump/CRL_quiescent.pdf',width=7.5,height=4)
ggplot(ssDat, 
       aes(x=volume,y=cytoRNA)) + 
  geom_point(size=1) + facet_wrap(~gene, scales = 'free_y') + expand_limits(x=0,y=0) + 
  stat_smooth(method="lm", se=FALSE, fullrange = TRUE,col='red') +
  theme_classic() +
  theme(axis.title=element_text(size="10"), axis.text=element_text(size='10')) +
  xlab('Volume (picoliter)') + ylab('mRNA Abundance') +
  theme(strip.background = element_rect(fill='white'),strip.text = element_text(face='italic',size=10))
dev.off()

pdf('~/Dropbox/densitypaper/SupplementaryFigures/DataDump/CRL_senescent.pdf',width=7.5,height=4)
ggplot(senDat, 
       aes(x=volume,y=cytoRNA)) + 
  geom_point(size=1.5) + facet_wrap(~gene, scales = 'free_y') + expand_limits(x=0,y=0) + 
  stat_smooth(method="lm", se=FALSE, fullrange = TRUE,col='red') +
  theme_classic() +
  theme(axis.title=element_text(size="10"), axis.text=element_text(size='10')) +
  xlab('Volume (picoliter)') + ylab('mRNA Abundance') +
  theme(strip.background = element_rect(fill='white'),strip.text = element_text(face='italic',size=10))
dev.off()

pdf('~/Dropbox/densitypaper/SupplementaryFigures/DataDump/A549_cycling.pdf',width=7.5,height=4)
ggplot(a549Dat, 
       aes(x=volume,y=cytoRNA)) + 
  geom_point(size=1) + facet_wrap(~gene, scales = 'free_y') + expand_limits(x=0,y=0) + 
  stat_smooth(method="lm", se=FALSE, fullrange = TRUE,col='red') +
  theme_classic() +
  theme(axis.title=element_text(size="10"), axis.text=element_text(size='10')) +
  xlab('Volume (picoliter)') + ylab('mRNA Abundance') +
  theme(strip.background = element_rect(fill='white'),strip.text = element_text(face='italic',size=10))
dev.off()