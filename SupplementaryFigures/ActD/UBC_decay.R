data <- read.delim('~/Dropbox/densitypaper/ExtractedData/ActD_SerumStarved.txt',header=T,stringsAsFactors=F)

data <- subset(data,gene %in% 'UBC')

meanDat <- ddply(data,.(timePt),summarise,meanRNA = mean(totalRNA),seRNA = sd(totalRNA)/sqrt(length(totalRNA)))

m0 = meanDat$meanRNA[1]



xvals = seq(from=0,to=6,by=.01)
yvals = m0*exp(-(1/4.2)*xvals)

modelDat <- data.frame(xvals,yvals)

pdf('~/Dropbox/densitypaper/SupplementaryFigures/ActD/UBC_decay.pdf',height=3,width=3.2)
ggplot(meanDat,aes(x=timePt,y=meanRNA)) + geom_point(size=1.3) + 
  geom_errorbar(aes(ymin=meanRNA-seRNA,ymax=meanRNA+seRNA),width=.1) +
  geom_line(data=modelDat,aes(x=xvals,y=yvals),col='cadetblue') +
  theme_classic() +
  xlab('Time after Actinomycin D addition (hours)') +
  ylab('Mean UBC mRNA count') +
  theme(axis.title=element_text(size="12"), axis.text=element_text(size='12'))
dev.off()