library(ggplot2)
library(plyr)

cycling <- read.delim('~/Dropbox/densitypaper/ExtractedData/CRL.txt',header=T,stringsAsFactors=F)
ss <- read.delim('~/Dropbox/densitypaper/ExtractedData/SerumStarved.txt',header=T,stringsAsFactors=F)

densityTableCycling <- ddply(cycling,.(gene),summarize,
                             densCycling = mean(cytoRNA/volume),
                             seCycling = sd(cytoRNA/volume)/sqrt(length(cytoRNA)))
densityTableSS <- ddply(ss,.(gene),summarize,
                        densSS = mean(cytoRNA/volume),
                        seSS = sd(cytoRNA/volume)/sqrt(length(cytoRNA)))

densityTable <- merge(densityTableCycling,densityTableSS)

pdf('~/Dropbox/densitypaper/SupplementaryFigures/DensityCyclingQuiescent/Density.pdf',width=3.3,height=3)
ggplot(densityTable,aes(x=densCycling,y=densSS)) +
  geom_point(size=1.3) +
  #expand_limits(x=c(0,1.4),y=c(0,1.4)) +
  scale_x_log10() +
  scale_y_log10() +
  geom_errorbar(aes(ymin=densSS-seSS,ymax=densSS+seSS)) +
  geom_errorbarh(aes(xmin=densCycling-seCycling,xmax=densCycling+seCycling)) +
  theme_classic() +
  geom_abline(slope=1,intercept=0) +
  xlab('Cytoplasmic RNA concentration (cycling)') + ylab('Cytoplasmic RNA concentration (quiescent)') +
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6'))
dev.off()