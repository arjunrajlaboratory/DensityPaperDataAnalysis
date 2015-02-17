library(boot)
library(plyr)
library(ggplot2)

data <- read.delim('~/Dropbox/ExtractedData_131216/alldata_CRL.txt',header=T,stringsAsFactors=F)

data$stage <- 'G1'
data$stage[data$numCyclin>20] <- 'SG2'

data <- subset(data,stage %in% c('G1','SG2'))

data <- ddply(data,.(gene,stage),subset,length(cytoRNA)>20)

coefs <- ddply(data, .(gene,stage), function(df) {
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

boots <- ddply(data, .(gene,stage), function(df) {
  results <- boot(data=df, statistic = boot.nm, 5000)
  out <- c(results$t0[1],boot.ci(results,type='bca')$bca[4],boot.ci(results,type='bca')$bca[5]);
  names(out) <- c('est','lower','upper');
  out
})

g1 <- subset(boots,stage %in% 'G1')[,c('gene','est','lower','upper')]
g2 <- subset(boots,stage %in% 'SG2')[,c('gene','est','lower','upper')]

colnames(g1) <- c(colnames(g1)[1],paste(colnames(g1)[2:4],'g1',sep='.'))
colnames(g2) <- c(colnames(g2)[1],paste(colnames(g2)[2:4],'g2',sep='.'))

dat <- merge(g1,g2)

# Now plot

#pdf('~/Dropbox/densitypaper/densitiyfigures/Figure4/g1_g2.pdf',width=2.5,height=2.2)
ggplot(dat, aes(x=est.g1,y=est.g2)) +
  geom_point(size=1.5) +
  geom_errorbarh(aes(xmin=lower.g1,xmax=upper.g1),height=.1) +
  geom_errorbar(aes(ymin=lower.g2,ymax=upper.g2),width=.1) +
  geom_abline(slope=1,intercept=0) +
  expand_limits(x=c(.001,10),y=c(.001,10)) +
  scale_y_log10() + scale_x_log10() +
  xlab('Noise measure (G1)') + ylab('Noise Measure (G2)') +
  theme_classic() + 
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6'))
#dev.off()