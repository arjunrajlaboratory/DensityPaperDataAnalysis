library(boot)
library(plyr)

crl <- read.delim('~/Dropbox/ExtractedData_131216/alldata_CRL.txt',header=T,stringsAsFactors=F)
a549 <- read.delim('~/Dropbox/ExtractedData_131216/alldata_A549.txt',header=T,stringsAsFactors=F)
quies <- read.delim('~/Dropbox/ExtractedData_131216/alldata_SS.txt',header=T,stringsAsFactors=F)
sen <- read.delim('~/Dropbox/ExtractedData_131216/alldata_Senescent.txt',header=T,stringsAsFactors=F)

data <- rbind(crl,a549,quies,sen)

# Only calculate for G1, and only take genes that have >20 cells in G1
data <- data[data$numCyclin<20,]
data <- ddply(data,.(gene,cellType),subset,length(cytoRNA)>20)

coefs <- ddply(data, .(gene,cellType), function(df) {
  m <- lm(cytoRNA ~ volume, data=df)
  data.frame(intercept = coef(m)[1], slope = coef(m)[2])
})

data <- merge(data,coefs)

boot.nm <- function(mydata,i){
  cytoRNA <- mydata$cytoRNA[i]
  slope <- mydata$slope[i]
  intercept <- mydata$intercept[i]
  volume <- mydata$volume[i]
  boot.nm <- (sd(cytoRNA)/mean(cytoRNA))^2 - 
    slope*mean(volume)/(intercept+slope*mean(volume)) * cov(cytoRNA,volume)/(mean(cytoRNA)*mean(volume))
}

boots <- ddply(data, .(gene,cellType), function(df) {
  results <- boot(data=df, statistic = boot.nm, 1000)
  out <- c(results$t0[1],boot.ci(results,type='bca')$bca[4],boot.ci(results,type='bca')$bca[5]);
  names(out) <- c('est','lower','upper');
  out
})

crl <- subset(boots,cellType %in% 'CRL2097')[,c('gene','est','lower','upper')]
sen <- subset(boots,cellType %in% 'CRL2097_Senescent')[,c('gene','est','lower','upper')]
quies <- subset(boots,cellType %in% 'CRL2097_SerumStarved')[,c('gene','est','lower','upper')]
a549 <- subset(boots,cellType %in% 'A549')[,c('gene','est','lower','upper')]

colnames(crl) <- c(colnames(crl)[1],paste(colnames(crl)[2:4],'crl',sep='.'))
colnames(sen) <- c(colnames(sen)[1],paste(colnames(sen)[2:4],'sen',sep='.'))
colnames(quies) <- c(colnames(quies)[1],paste(colnames(quies)[2:4],'quies',sep='.'))
colnames(a549) <- c(colnames(a549)[1],paste(colnames(a549)[2:4],'a549',sep='.'))

dat1 <- merge(crl,sen)
dat2 <- merge(crl,quies)
dat3 <- merge(crl,a549)

### Now plot

cols = c('purple4','darkorchid3','magenta4')

#pdf('~/Dropbox/densitypaper/densitiyfigures/Figure4/cyc_vs_a549.pdf',width=2.5,height=2.2)
p1 <- ggplot(dat3, aes(x=est.crl,y=est.a549,label=gene)) +
  geom_abline(slope=1,intercept=0,col='gray78') +
  geom_point(size=1.5,col=cols[1]) +
  scale_y_log10() + scale_x_log10() +
  geom_errorbarh(aes(xmin=lower.crl,xmax=upper.crl),height=.1,col=cols[1]) +
  geom_errorbar(aes(ymin=lower.a549,ymax=upper.a549),width=.1,col=cols[1]) +
  expand_limits(x=c(.001,10),y=c(.001,20)) +
  xlab('Noise measure (CRL)') + ylab('Noise Measure (A549)') +
  theme_classic() + 
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6'))
#dev.off()

#pdf('~/Dropbox/densitypaper/densitiyfigures/Figure4/cyc_vs_sen.pdf',width=2.5,height=2.2)
p2 <- ggplot(dat1, aes(x=est.crl,y=est.sen)) +
  geom_abline(slope=1,intercept=0,col='gray78') +
  geom_point(size=1.5,col=cols[2]) +
  geom_errorbarh(aes(xmin=lower.crl,xmax=upper.crl),height=.1,col=cols[2]) +
  geom_errorbar(aes(ymin=lower.sen,ymax=upper.sen),width=.1,col=cols[2]) +
  expand_limits(x=c(.001,10),y=c(.001,20)) +
  scale_y_log10() + scale_x_log10() +
  xlab('Noise measure (Cycling)') + ylab('Noise Measure (Senescent)') +
  theme_classic() + 
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6'))
#dev.off()

#pdf('~/Dropbox/densitypaper/densitiyfigures/Figure4/cyc_vs_quies.pdf',width=2.5,height=2.2)
p3 <- ggplot(dat2, aes(x=est.crl,y=est.quies)) +
  geom_abline(slope=1,intercept=0,col='gray78') +
  geom_point(size=1.5,col=cols[3]) +
  geom_errorbarh(aes(xmin=lower.crl,xmax=upper.crl),height=.1,col=cols[3]) +
  geom_errorbar(aes(ymin=lower.quies,ymax=upper.quies),width=.1,col=cols[3]) +
  expand_limits(x=c(.001,10),y=c(.001,20)) +
  scale_y_log10() + scale_x_log10() +
  xlab('Noise measure (Cycling)') + ylab('Noise Measure (Quiescent)') +
  theme_classic() + 
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6'))
#dev.off()

#pdf('~/Dropbox/densitypaper/densitiyfigures/Figure4/celltypes.pdf',width=8,height=2.2)
multiplot(p1,p3,p2,cols=3)
#dev.off()