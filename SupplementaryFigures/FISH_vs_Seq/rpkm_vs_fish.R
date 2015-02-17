library(plyr)
library(ggplot2)

fishData <- read.delim('~/Dropbox/densitypaper/ExtractedData/CRL.txt',header=T,stringsAsFactors=F)
seqData <- read.delim('~/Dropbox/NucSeqData/htseq_CRL/CRL_rpkm.txt',header=T,stringsAsFactors=F)

seqReduced <- ddply(seqData,.(gene,total1_ex_rpkm),summarize,exonRPKM=sum(total1_ex_rpkm,total2_ex_rpkm)/2)

fishReduced <- ddply(fishData,.(gene),summarize,meanRNA=mean(cytoRNA))

tmp <- merge(fishReduced,seqReduced)

pdf('~/Dropbox/densitypaper/SupplementaryFigures/FISH_vs_Seq/rpkm_vs_fish.pdf',width=4.3,height=4)
ggplot(tmp,aes(x=meanRNA,y=exonRPKM)) +
  geom_point(size=1.5) +
  scale_x_log10() + scale_y_log10() +
  theme_classic() +
  xlab('Mean RNA Count (RNA-FISH)') + ylab('RPKM') +
  theme(axis.title=element_text(size="12"), axis.text=element_text(size='12'))
dev.off()

summary(lm(exonRPKM ~ meanRNA, data=tmp))
