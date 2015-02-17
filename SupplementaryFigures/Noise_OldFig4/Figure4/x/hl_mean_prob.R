library(plyr)
library(ggplot2)

data <- read.delim('~/Dropbox/densitypaper/densitiyfigures/Figure4/Table1_NoiseMeasureCycling.txt',
                   header=T,stringsAsFactors=F)

data <- data[data$Gene != 'DNMT1',]

# Noise and abundance

pdf('~/Dropbox/densitypaper/densitiyfigures/Figure4/noise_vs_mean.pdf',width=2.5,height=2.2)
ggplot(data,aes(x=MeanRna,y=NoiseMeasure)) +
  geom_point(size=1.5) +
  geom_errorbar(aes(ymin=Lower,ymax=Upper),width=.10) +
  scale_y_log10() + scale_x_log10() +
  geom_line(data=data, aes(x=MeanRna, y=NoiseFloor), colour="red") +
  expand_limits(y=c(0.015,10)) +
  xlab('Mean RNA') + ylab('Noise Measure') +
  theme_classic() + 
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6'))
dev.off()

pdf('~/Dropbox/densitypaper/densitiyfigures/Figure4/noise_vs_prob.pdf',width=2.5,height=2.2)
ggplot(data,aes(x=TranscriptionProbability,y=NoiseMeasure)) +
  geom_point(size=1.5) +
  geom_errorbar(aes(ymin=Lower,ymax=Upper),width=.080) +
  scale_y_log10() +
  expand_limits(y=c(0.015,10)) +
  xlab('Transcriptional ON Time') + ylab('Noise Measure') +
  theme_classic() + 
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6'))
dev.off()

pdf('~/Dropbox/densitypaper/densitiyfigures/Figure4/noise_vs_hl.pdf',width=2.5,height=2.2)
ggplot(data,aes(x=Halflife,y=NoiseMeasure)) +
  geom_point(size=1.5) +
  geom_errorbar(aes(ymin=Lower,ymax=Upper),width=.5) +
  scale_y_log10() +
  expand_limits(y=c(0.015,10),x=c(0,1)) +
  xlab('mRNA Half-life') + ylab('Noise Measure') +
  theme_classic() + 
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6'))
dev.off()