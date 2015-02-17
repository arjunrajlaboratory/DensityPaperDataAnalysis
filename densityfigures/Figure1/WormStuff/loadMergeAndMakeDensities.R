library("ggplot2")
library("plyr")
library(reshape2)
data1 <- read.csv('~/Dropbox/densitypaper/WormData/140616_wormdensity_analysis_CB502_fixed140610/totalextraction.csv',header=T,stringsAsFactors=F)
data2 <- read.csv('~/Dropbox/densitypaper/WormData/140616_wormdensity_analysis_CB502_fixed140610/WormIdentities.csv',header=T,stringsAsFactors=F)
combinedData1 <- merge(data1, data2)
data3 <- read.csv('~/Dropbox/densitypaper/WormData/140616_wormdensity_analysis_N2_fixed140609/WormIdentities.csv',header=T,stringsAsFactors=F)
data4 <- read.csv('~/Dropbox/densitypaper/WormData/140616_wormdensity_analysis_N2_fixed140609/totalextraction.csv',header=T,stringsAsFactors=F)
combinedData2 <- merge(data3, data4)
data5 <- read.csv('~/Dropbox/densitypaper/WormData/140620_wormdensity_analysis_CB502_fixed140610/WormIdentities.csv',header=T,stringsAsFactors=F)
data6 <- read.csv('~/Dropbox/densitypaper/WormData/140620_wormdensity_analysis_CB502_fixed140610/totalextraction.csv',header=T,stringsAsFactors=F)
combinedData3 <- merge(data5, data6)
data7 <- read.csv('~/Dropbox/densitypaper/WormData/140626_wormdensity_analysis_N2_fixed140623/WormIdentities.csv',header=T,stringsAsFactors=F)
data8 <- read.csv('~/Dropbox/densitypaper/WormData/140626_wormdensity_analysis_N2_fixed140623/totalextraction.csv',header=T,stringsAsFactors=F)
combinedData4 <- merge(data7, data8)
data9 <- read.csv('~/Dropbox/densitypaper/WormData/140630_wormdensity_analysis_N2_fixed140623/WormIdentities.csv',header=T,stringsAsFactors=F)
data10 <- read.csv('~/Dropbox/densitypaper/WormData/140630_wormdensity_analysis_N2_fixed140623/totalextraction.csv',header=T,stringsAsFactors=F)
combinedData5 <- merge(data9, data10)
data11 <- read.csv('~/Dropbox/densitypaper/WormData/140702_wormdensity_analysis_CB502_fixed140625/WormIdentities.csv',header=T,stringsAsFactors=F)
data12 <- read.csv('~/Dropbox/densitypaper/WormData/140702_wormdensity_analysis_CB502_fixed140625/totalextraction.csv',header=T,stringsAsFactors=F)
combinedData6 <- merge(data11, data12)
combinedData <- rbind(combinedData1, combinedData2, combinedData3, combinedData4, combinedData5, combinedData6)
maxZdiffTable <- ddply(combinedData,.(objArrayNum,objNum,WormNum,Strain,Fixed),summarize,
                       maxZdiff = max(tmrZdiff,alexaZdiff))
combinedData <- merge(combinedData,maxZdiffTable)
combinedData$Volume <- combinedData$area*0.125*0.125*0.3*combinedData$maxZdiff
combinedData$VolDensityAlexa <- combinedData$alexaspots/combinedData$Volume
combinedData$VolDensityTmr <- combinedData$tmrspots/combinedData$Volume
combinedData$PerNucDensityAlexa <- combinedData$alexaspots/combinedData$numnuclei
combinedData$PerNucDensityTmr <- combinedData$tmrspots/combinedData$numnuclei
combinedData$nucleiPerVol <- combinedData$numnuclei/combinedData$Volume
combinedData$nadVolPerCell <- combinedData$Volume/combinedData$numnuclei