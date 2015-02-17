library(ggplot2)

data <- read.delim('~/Dropbox/ExtractedData_131216/alldata_CRL.txt',header=T,stringsAsFactors=T)

cols <- c('green4','salmon','green4')

tmp <- subset(data,gene %in% 'USF2' & volume < 5000)

p1 <- ggplot(tmp,aes(x=volume,y=cytoRNA)) +
  geom_point(size=1.3,col=cols[1]) +
  expand_limits(x=0,y=0) +
  xlab('Volume (picoliter)') + ylab('RNA') +
  theme_classic() + 
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6'))

tmp <- subset(data,gene %in% 'RND3' & volume < 5000)

p2 <- ggplot(tmp,aes(x=volume,y=cytoRNA)) +
  geom_point(size=1.3,col=cols[2]) +
  expand_limits(x=0,y=0) +
  xlab('Volume (picoliter)') + ylab('RNA') +
  theme_classic() + 
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6'))

tmp <- subset(data,gene %in% 'EEF2' & volume < 5000)

p3 <- ggplot(tmp,aes(x=volume,y=cytoRNA)) +
  geom_point(size=1.3,col=cols[3]) +
  expand_limits(x=0,y=0) +
  xlab('Volume (picoliter)') + ylab('RNA') +
  theme_classic() + 
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6'))

pdf('~/Dropbox/densitypaper/densitiyfigures/Figure4/3genes.pdf',width=5,height=1.5)
multiplot(p1,p2,p3,cols=3)
dev.off()
