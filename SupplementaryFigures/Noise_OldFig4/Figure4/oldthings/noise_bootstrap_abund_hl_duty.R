library(boot)
library(plyr)
library(ggplot2)

data <- read.delim('~/Dropbox/ExtractedData_131216/alldata_CRL.txt',header=T,stringsAsFactors=F)

# Only calculate for G1, and only take genes that have >20 cells in G1
data <- data[data$numCyclin<20 & data$numTxnSites<3,]
data <- ddply(data,.(gene),subset,length(cytoRNA)>20)

hl <- read.delim('~/Dropbox/Foreskin/half_life_with_gene_names.txt',header=F,stringsAsFactors=F)

colnames(hl) <- c('nm','gene','hl')
hl <- hl[,c('gene','hl')]
hl <- hl[hl$hl!='N.D.',]
hl$hl <- gsub(' ', '', hl$hl)
hl$hl <- gsub('>24', '24', hl$hl)
hl$hl <- as.numeric(hl$hl)

dutyNoisefloorAbund <- ddply(data,.(gene),summarize,
                   dutycycle = mean(numTxnSites)/2,
                   noisefloor = 1/mean(cytoRNA)+.02,
                   meanRNA = mean(cytoRNA))

coefs <- ddply(data, .(gene), function(df) {
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

boots <- ddply(data, .(gene), function(df) {
  results <- boot(data=df, statistic = boot.nm, 10000)
  out <- c(results$t0[1],boot.ci(results,type='bca')$bca[4],boot.ci(results,type='bca')$bca[5]);
  names(out) <- c('est','lower','upper');
  out
})

# Merge hl, duty cycle, noise, noise floor

outData <- merge(dutyNoisefloorAbund, hl)
outData <- merge(outData,boots)
outData <- unique(outData)

### Now make the plots...

cols = c('navy','turquoise4','dodgerblue3')

#pdf('~/Dropbox/densitypaper/densitiyfigures/Figure4/noise_vs_mean.pdf',width=2.5,height=2.2)
p1 <- ggplot(outData,aes(x=meanRNA,y=est,label=gene)) +
  #geom_text() +
  geom_point(size=1.5,col=cols[1]) +
  geom_errorbar(aes(ymin=lower,ymax=upper),width=.10,col=cols[1]) +
  scale_y_log10() + scale_x_log10() +
  geom_line(data=outData, aes(x=meanRNA, y=noisefloor), colour="gray78") +
  geom_line(data=outData, aes(x=meanRNA, y=1/meanRNA), col = 'gray78') +
  xlab('Mean RNA') + ylab('Noise Measure') +
  theme_classic() + 
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6')) +
  theme(legend.position='')
#dev.off()

#pdf('~/Dropbox/densitypaper/densitiyfigures/Figure4/noise_vs_prob.pdf',width=2.5,height=2.2)
p2 <- ggplot(outData,aes(x=dutycycle,y=est)) +
  geom_point(size=1.5,col=cols[2]) +
  geom_errorbar(aes(ymin=lower,ymax=upper),width=.045,col=cols[2]) +
  scale_y_log10() +
  expand_limits(x=c(0,1)) +
  xlab('Transcriptional ON Time') + ylab('Noise Measure') +
  theme_classic() + 
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6')) +
  theme(legend.position='')
#dev.off()

#pdf('~/Dropbox/densitypaper/densitiyfigures/Figure4/noise_vs_hl.pdf',width=2.5,height=2.2)
p3 <- ggplot(outData,aes(x=hl,y=est)) +
  geom_point(size=1.5,col=cols[3]) +
  geom_errorbar(aes(ymin=lower,ymax=upper),width=1,col=cols[3]) +
  scale_y_log10() +
  xlab('Halflife') + ylab('Noise Measure') +
  theme_classic() + 
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6')) +
  theme(legend.position='')
#dev.off()

#pdf('~/Dropbox/densitypaper/densitiyfigures/Figure4/abund_hl_duty.pdf',width=8,height=2.2)
multiplot(p1,p2,p3,cols=3)
#dev.off()

ggplot(outData,aes(x=dutycycle,y=hl,col=est,label=round(meanRNA))) +
  geom_text() +
  scale_color_gradient(low='red',high='blue',trans='log') +
  theme_classic()

summary(lm(est ~ dutycycle + hl + meanRNA, data=outData))
summary(lm(est ~ hl, data=outData))
summary(lm(est ~ dutycycle, data=outData))
summary(lm(log(est) ~ log(meanRNA), data=outData))
