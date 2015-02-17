library(reshape2)
library(ggplot2)
library(plyr)

crl <- read.delim('~/Dropbox/ExtractedData_131216/alldata_CRL.txt',header=T,stringsAsFactors=F)
sen <- read.delim('~/Dropbox/ExtractedData_131216/alldata_Senescent.txt',header=T,stringsAsFactors=F)
ss <- read.delim('~/Dropbox/ExtractedData_131216/alldata_SS.txt',header=T,stringsAsFactors=F)

data <- rbind(crl,sen,ss)

data$volume <- data$volume/1000

numericcols <- c('volume','area','totalGAPDH','cytoGAPDH','nucGAPDH','nucArea',
                 'nucIntensityTotal','nucIntensityAvg','numCyclin','totalRNA','cytoRNA',
                 'nucRNA','numTxnSites','avgIntron','avgExon','numRnaPerTxnSiteAvgExon');


for (col in numericcols){
  data[[col]] <- as.numeric(data[[col]])
}

slopeIntTable <- ddply(data,.(cellType,gene),function(x) {
  tmp <- summary(lm(cytoRNA~volume,x))$coefficients; 
  out <- c(tmp[1,1:2], tmp[2,1:2]);
  names(out) <- paste(rep(c("inter", "slope"), each=ncol(tmp)/2), 
                      c('estimate','stderr'), sep=".");
  out
})

meanVolTable <- ddply(data,.(cellType),summarize,meanVol = mean(volume))

slopeIntTable <- merge(slopeIntTable,meanVolTable)

#summary(lm(cytoRNA~volume,subset(data,gene %in% 'EEF2' & cellType %in% 'CRL2097')))$coefficients

crlTable <- subset(slopeIntTable,cellType %in% 'CRL2097')
ssTable <- subset(slopeIntTable,cellType %in% 'CRL2097_SerumStarved')

colnames(crlTable)[3:6] <- paste(colnames(crlTable)[3:6],'CRL2097',sep='.')
crlTable <- crlTable[,c(2:6)]
colnames(ssTable)[3:6] <- paste(colnames(ssTable)[3:6],'SS',sep='.')
ssTable <- ssTable[,c(2:6)]

crlSsDat <- merge(crlTable,ssTable)

genesInBoth <- intersect(unique(crl$gene),unique(ss$gene))
volCRL <- mean(crl$volume)
volSS <- mean(ss$volume)

# scatter plot vol-indep
ggplot(subset(crlSsDat,gene %in% genesInBoth),aes(x=inter.estimate.CRL2097,y=inter.estimate.SS)) +
  geom_point() +
  geom_abline(slope=1,intercept=0) +
  scale_x_log10() +
  scale_y_log10() +
  geom_errorbar(aes(ymin=(inter.estimate.SS-inter.stderr.SS),
                    ymax=(inter.estimate.SS+inter.stderr.SS))) +
  geom_errorbarh(aes(xmin=(inter.estimate.CRL2097-inter.stderr.CRL2097),
                     xmax=(inter.estimate.CRL2097+inter.stderr.CRL2097))) +
  expand_limits(x=c(1,1000),y=c(1,1000)) +
  theme_classic()

# scatter plot vol-dep
ggplot(subset(crlSsDat,gene %in% genesInBoth),
       aes(x=slope.estimate.CRL2097*volCRL,y=slope.estimate.SS*volSS)) +
  geom_point() +
  geom_abline(slope=1,intercept=0) +
  scale_x_log10() +
  scale_y_log10() +
  geom_errorbar(aes(ymin=(slope.estimate.SS-slope.stderr.SS)*volSS,
                    ymax=(slope.estimate.SS+slope.stderr.SS)*volSS)) +
  geom_errorbarh(aes(xmin=(slope.estimate.CRL2097-slope.stderr.CRL2097)*volCRL,
                    xmax=(slope.estimate.CRL2097+slope.stderr.CRL2097)*volCRL)) +
  expand_limits(x=c(3e3,3e6),y=c(3e3,3e6)) +
  theme_classic()

# bar plot vol-indep
ggplot(subset(slopeIntTable,cellType %in% c('CRL2097','CRL2097_SerumStarved') & gene %in% genesInBoth),
       aes(x=gene,y=inter.estimate,fill=cellType)) +
  geom_bar(stat='identity',position='dodge') +  
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(legend.position='none')

# bar plot vol-dep
ggplot(subset(slopeIntTable,cellType %in% c('CRL2097','CRL2097_SerumStarved') & gene %in% genesInBoth),
       aes(x=gene,y=slope.estimate*meanVol,fill=cellType)) +
  geom_bar(stat='identity',position='dodge') +  
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(legend.position='none')
