library(reshape2)
library(ggplot2)
library(plyr)

crl <- read.delim('~/Dropbox/densitypaper/densitypaperdataanalysis/ExtractedData/CRL.txt',header=T,stringsAsFactors=F)
ss <- read.delim('~/Dropbox/densitypaper/densitypaperdataanalysis/ExtractedData/SerumStarved.txt',header=T,stringsAsFactors=F)

data <- rbind(crl,ss)

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

meanVolTable <- ddply(data,.(cellType),summarise,meanVol = mean(volume))

slopeIntTable <- merge(slopeIntTable,meanVolTable)


crlTable <- subset(slopeIntTable,cellType %in% 'CRL2097')
ssTable <- subset(slopeIntTable,cellType %in% 'CRL2097_SerumStarved')

colnames(crlTable)[3:6] <- paste(colnames(crlTable)[3:6],'CRL2097',sep='.')
crlTable <- crlTable[,c(2:6)]
colnames(ssTable)[3:6] <- paste(colnames(ssTable)[3:6],'SS',sep='.')
ssTable <- ssTable[,c(2:6)]

crlSsDat <- merge(crlTable,ssTable)

genesInBoth <- intersect(unique(crl$gene),unique(ss$gene))
volCRL <- mean(crl$volume)/1000
volSS <- mean(ss$volume)/1000

crlSsDat$inter.min.SS <- crlSsDat$inter.estimate.SS-crlSsDat$inter.stderr.SS
crlSsDat$inter.min.SS[crlSsDat$inter.min.SS<0] <- crlSsDat$inter.estimate.SS[crlSsDat$inter.min.SS<0]
crlSsDat$inter.max.SS <- crlSsDat$inter.estimate.SS+crlSsDat$inter.stderr.SS

crlSsDat$inter.min.CRL2097 <- crlSsDat$inter.estimate.CRL2097-crlSsDat$inter.stderr.CRL2097
crlSsDat$inter.min.CRL2097[crlSsDat$inter.min.CRL2097<0] <- 
  crlSsDat$inter.estimate.CRL2097[crlSsDat$inter.min.CRL2097<0]
crlSsDat$inter.max.CRL2097 <- crlSsDat$inter.estimate.CRL2097+crlSsDat$inter.stderr.CRL2097

crlSsDat$slope.min.SS <- crlSsDat$slope.estimate.SS-crlSsDat$slope.stderr.SS
crlSsDat$slope.min.SS[crlSsDat$slope.min.SS<0] <- crlSsDat$slope.estimate.SS[crlSsDat$slope.min.SS<0]
crlSsDat$slope.max.SS <- crlSsDat$slope.estimate.SS+crlSsDat$slope.stderr.SS

crlSsDat$slope.min.CRL2097 <- crlSsDat$slope.estimate.CRL2097-crlSsDat$slope.stderr.CRL2097
crlSsDat$slope.min.CRL2097[crlSsDat$slope.min.CRL2097<0] <- 
  crlSsDat$slope.estimate.CRL2097[crlSsDat$slope.min.CRL2097<0]
crlSsDat$slope.max.CRL2097 <- crlSsDat$slope.estimate.CRL2097+crlSsDat$slope.stderr.CRL2097

# scatter plot vol-indep
pdf('~/Dropbox/densitypaper/densitypaperdataanalysis/SupplementaryFigures/VolumeDepIndepTxn/Vol_Indep_Cyc_Qui.pdf',width=3.3,height=3)
ggplot(subset(crlSsDat,gene %in% genesInBoth),aes(x=inter.estimate.CRL2097,y=inter.estimate.SS)) +
  geom_point(size=1.5) +
  geom_abline(slope=1,intercept=0) +
  scale_x_log10() +
  scale_y_log10() +
  geom_errorbar(aes(ymin=inter.min.SS,ymax=inter.max.SS),width=.05) +
  geom_errorbarh(aes(xmin=inter.min.CRL2097,xmax=inter.max.CRL2097),height=.05) +
  expand_limits(x=c(1,1000),y=c(1,1000)) +
  theme_classic() +
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6')) +
  xlab('Volume-independent Expression (Cycling)') + ylab('Volume-independent Expression (Senescent)')
dev.off()

# scatter plot vol-dep
pdf('~/Dropbox/densitypaper/densitypaperdataanalysis/SupplementaryFigures/VolumeDepIndepTxn/Vol_Dep_Cyc_Qui.pdf',width=3.3,height=3)
ggplot(subset(crlSsDat,gene %in% genesInBoth),
       aes(x=slope.estimate.CRL2097*volCRL,y=slope.estimate.SS*volSS)) +
  geom_point(size=1.5) +
  geom_abline(slope=1,intercept=0) +
  scale_x_log10() +
  scale_y_log10() +
  geom_errorbar(aes(ymin=slope.min.SS*volSS,ymax=slope.max.SS*volSS),width=.05) +
  geom_errorbarh(aes(xmin=slope.min.CRL2097*volCRL,xmax=slope.max.CRL2097*volCRL),height=.05) +
  expand_limits(x=c(10,1500),y=c(10,1500)) +
  theme_classic() +
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6')) +
  xlab('Volume-dependent Expression (Cycling)') + ylab('Volume-dependent Expression (Senescent)')
dev.off()

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
