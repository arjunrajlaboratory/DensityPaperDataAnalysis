data <- read.delim('~/Documents/Data/131110_CRL_MeOH_EEF2/TxnSiteIntensityManual.txt',
                   header=F,stringsAsFactors=F)

colnames(data) <- c('datNum','objNum','numInt','avgInt','int1','int2','int3','int4','volume','cyclin')
#colnames(dataU) <- colnames(dataE)

#dataE$gene <- 'EEF2'
#dataU$gene <- 'UBC'

#data <- rbind(dataE,dataU)

data <- data[!is.nan(data$avgInt),]


g1 <- data[data$cyclin<20,]
g1 <- g1[g1$numInt<3,]
g1$stage <- 'G1'
g2 <- data[data$cyclin>230,]
g2$stage <- 'G2'
s <- data[(data$cyclin>20 & data$cyclin<230),]
s$stage <- 'S'

outDat <- rbind(g1,g2,s)

ggplot(subset(outDat, avgInt > 0),aes(x=volume,y=avgInt)) +
  geom_point() + geom_smooth(method = 'lm')

lin.out <- lm(avgInt ~ volume + factor(stage), subset(outDat, avgInt > 0))
