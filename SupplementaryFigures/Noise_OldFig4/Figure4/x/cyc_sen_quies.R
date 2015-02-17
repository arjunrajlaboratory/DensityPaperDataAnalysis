cycling <- read.delim('~/Dropbox/densitypaper/densitiyfigures/Figure4/Table1_NoiseMeasureCycling.txt',
                      header = T, stringsAsFactors = F)
sen <- read.delim('~/Dropbox/densitypaper/densitiyfigures/Figure4/Noise_Senescent.txt',
                  header = T, stringsAsFactors = F)
quies <- read.delim('~/Dropbox/densitypaper/densitiyfigures/Figure4/Noise_Quiescent.txt',
                    header = T, stringsAsFactors = F)
colnames(sen) <- c(colnames(sen)[1],paste(colnames(sen)[2:6],'sen',sep='.'))
colnames(quies) <- c(colnames(quies)[1],paste(colnames(quies)[2:6],'quies',sep='.'))

dat1 <- merge(cycling,sen)
dat2 <- merge(cycling,quies)

pdf('~/Dropbox/densitypaper/densitiyfigures/Figure4/cyc_vs_sen.pdf',width=2.5,height=2.2)
ggplot(dat1,aes(x=NoiseMeasure,y=NoiseMeasure.sen)) +
  geom_point(size=1.3) +
  geom_errorbarh(aes(xmin=Lower,xmax=Upper),height=.01) +
  geom_errorbar(aes(ymin=Lower.sen,ymax=Upper.sen),width=.01) +
  xlab('Noise Measure (Cycling)') + ylab('Noise Measure (Senescent)') +
  theme_classic() + 
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6'))
dev.off()

pdf('~/Dropbox/densitypaper/densitiyfigures/Figure4/cyc_vs_quies.pdf',width=2.5,height=2.2)
ggplot(dat2,aes(x=NoiseMeasure,y=NoiseMeasure.quies)) +
  geom_point(size=1.3) +
  geom_errorbarh(aes(xmin=Lower,xmax=Upper),height=.01) +
  geom_errorbar(aes(ymin=Lower.quies,ymax=Upper.quies),width=.01) +
  xlab('Noise Measure (Cycling)') + ylab('Noise Measure (Quiescent)') +
  theme_classic() + 
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6'))
dev.off()