library(ggplot2)
library(plyr)

cycling <- read.delim('~/Dropbox/densitypaper/ExtractedData/CRL.txt',header=T,stringsAsFactors=F)
ss <- read.delim('~/Dropbox/densitypaper/ExtractedData/SerumStarved.txt',header=T,stringsAsFactors=F)

corTableCycling <- ddply(cycling,.(gene),summarize,
                  r2cycling = cor(cytoRNA,volume)^2)
corTableSS <- ddply(ss,.(gene),summarize,
                    r2ss = cor(cytoRNA,volume)^2)

corTable <- merge(corTableCycling,corTableSS)

ggplot(corTable,aes(x=r2cycling,y=r2ss)) +
  geom_point() +
  expand_limits(x=0,y=0) +
  theme_classic() +
  geom_abline(slope=1,intercept=0)

densityTableCycling <- ddply(cycling,.(gene),summarize,
                             densCycling = mean(cytoRNA/volume),
                             seCycling = sd(cytoRNA/volume)/sqrt(length(cytoRNA)))
densityTableSS <- ddply(ss,.(gene),summarize,
                             densSS = mean(cytoRNA/volume),
                        seSS = sd(cytoRNA/volume)/sqrt(length(cytoRNA)))

densityTable <- merge(densityTableCycling,densityTableSS)

ggplot(densityTable,aes(x=densCycling,y=densSS)) +
  geom_point() +
  expand_limits(x=0,y=0) +
  geom_errorbar(aes(ymin=densSS-seSS,ymax=densSS+seSS)) +
  geom_errorbarh(aes(xmin=densCycling-seCycling,xmax=densCycling+seCycling)) +
  theme_classic() +
  geom_abline(slope=1,intercept=0)

cvTable <- ddply(cycling,.(gene),summarize,
                 cv = sd(cytoRNA)/mean(cytoRNA),
                 mean = mean(cytoRNA))

ggplot(cvTable, aes(x=mean,y=cv)) +
  geom_point() +
  expand_limits(x=0,y=0) +
  theme_classic()
