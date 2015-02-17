library(ggplot2)

data <- read.delim('~/Dropbox/densitypaper/SupplementaryFigures/VolumeControls/VolumeComparison.txt',
                   header=T,stringsAsFactors=F)

pdf('~/Dropbox/densitypaper/SupplementaryFigures/VolumeControls/EEF2_Vol.pdf',width=3,height=2.7)
ggplot(data,aes(x=GAPDH_volume/1000,y=EEF2_volume/1000)) +
  geom_point(size=1.3) +
  expand_limits(x=0,y=c(0,5)) +
  xlab('Volume Using GAPDH (pL)') + ylab('Volume Using EEF2 (pL)') +
  theme_classic() +
  geom_abline(slope=1,intercept=0) +
  theme(axis.title=element_text(size="10"), axis.text=element_text(size='10'))
dev.off()

pdf('~/Dropbox/densitypaper/SupplementaryFigures/VolumeControls/Half_GAPDH_Vol.pdf',width=3,height=2.7)
ggplot(data,aes(x=GAPDH_volume/1000,y=Half_GAPDH/1000)) +
  geom_point(size=1.3) +
  expand_limits(x=0,y=0) +
  xlab('Volume Using GAPDH (pL)') + ylab('Volume Using Half GAPDH (pL)') +
  theme_classic() +
  geom_abline(slope=1,intercept=0) +
  theme(axis.title=element_text(size="10"), axis.text=element_text(size='10'))
dev.off()