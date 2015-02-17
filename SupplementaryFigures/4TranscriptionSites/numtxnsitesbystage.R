library(ggplot2)
library(plyr)

data <- read.delim('~/Dropbox/densitypaper/ExtractedData/FreqIntens.txt',header=T,stringsAsFactors=F)

data$stage <- 'G1'
data$stage[data$numCyclin>20] <- 'S'
data$stage[data$numCyclin>230] <- 'G2'

data <- unique(data[,c('gene','date','fixative','dataNum','objNum','numTxnSites','stage')])
reducedDataNumSites <- subset(data,gene %in% 'EEF2')

reducedData <- subset(data,intensity>750)
reducedDataNumSites <- ddply(reducedData,.(datNum,objNum,stage),summarize,numTxnSites=length(intensity))
ggplot(reducedDataNumSites,aes(x=numTxnSites,fill=stage)) +
  geom_bar(position='dodge')

ggplot(data,aes(x=xVar,y=intensity)) +
  geom_point()
#ggplot(d=reducedDataNumSites,aes(x=numTxnSites,fill=stage)) +
#  geom_bar(position='dodge')

## THIS IS AWESOME. FREQUENCY TABLE ##
tab <- as.data.frame(prop.table(table(reducedDataNumSites$numTxnSites, reducedDataNumSites$stage), 2))
ggplot(tab, aes(y = Freq,  x = Var1, fill = Var2)) + 
  geom_bar(stat = "identity", position = position_dodge())

#tmp <- ddply(data,.(datNum,objNum,stage),summarize,numTxnSites=mean(numInt))

ggplot(dataNumSites,aes(x=numTxnSites,fill=stage)) +
  geom_bar(position='dodge')
