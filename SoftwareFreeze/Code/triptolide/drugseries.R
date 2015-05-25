library(ggplot2)
library(plyr)

dataCtrl <- read.delim('~/Documents/Data/141105_Triptolide/141105_Triptolide_100nM/141105_Triptolide_100nM.txt',header=T,stringsAsFactors=F)
data100 <- read.delim('~/Documents/Data/141105_Triptolide/141105_Triptolide_Ctrl/141105_Triptolide_Ctrl.txt',header=T,stringsAsFactors=F)
data10 <- read.delim('~/Documents/Data/141105_Triptolide/141105_Triptolide_10nM/141105_Triptolide_10nM.txt',header=T,stringsAsFactors=F)
data1 <- read.delim('~/Documents/Data/141105_Triptolide/141105_Triptolide_1nM/141105_Triptolide_1nM.txt',header=T,stringsAsFactors=F)
data10um <- read.delim('~/Documents/Data/141105_Triptolide/141105_Triptolide_10uM/141105_Triptolide_10uM.txt',header=T,stringsAsFactors=F)
data1um <- read.delim('~/Documents/Data/141105_Triptolide/141105_Triptolide_1uM/141105_Triptolide_1uM.txt',header=T,stringsAsFactors=F)

dat <- rbind(dataCtrl,data100,data10,data1,data10um,data1um)
dat <- rbind(dataCtrl,data100,data10,data1)

ggplot(dat,aes(x=intensity,fill=drugConc)) + geom_density(alpha=.25)

meanDat <- ddply(dat,.(drugConc),summarise,
                 meanInt = mean(intensity))

ggplot(dat,aes(y=intensity,x=drugConc)) + geom_boxplot()

intensityDat <- ddply(dat,.(drugConc,dataNum,objNum),summarise,
                      meanInt = mean(intensity),
                      maxInt = max(intensity),
                      numTxnSites = max(numTxnSites))

ggplot(intensityDat,aes(x=maxInt,fill=drugConc)) + geom_density(alpha=.25)

ggplot(intensityDat,aes(y=maxInt,x=drugConc)) + geom_boxplot()

dodge = position_dodge(0.2)
ggplot(intensityDat,aes(x=numTxnSites,fill=drugConc)) + geom_bar(position=dodge,width=.2)

ggplot(intensityDat,aes(x=numTxnSites,fill=drugConc)) + geom_density(alpha=.25)
