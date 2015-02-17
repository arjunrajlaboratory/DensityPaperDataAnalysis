library(ggplot2)

data <- read.delim('~/Dropbox/ExtractedData_131216/alldata_CRL.txt',header=T,stringsAsFactors=F)

lmna <- subset(data,gene %in% 'LMNA')

pdf('~/Dropbox/densitypaper/SupplementaryFigures/LMNA/LMNA.pdf',height=3,width=3.3)
ggplot(lmna,aes(x=volume/1000,y=cytoRNA)) +
  geom_point(col='palevioletred3',size=1.3) +
  expand_limits(x=0,y=0) +
  xlab('Volume (picoliter)') + ylab('LMNA mRNA') +
  theme_classic() + 
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6'))
dev.off()