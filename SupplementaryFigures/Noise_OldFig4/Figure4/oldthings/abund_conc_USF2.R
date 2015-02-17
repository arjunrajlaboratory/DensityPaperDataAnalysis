library(ggplot2)
library(gridExtra)

data <- read.delim('~/Dropbox/ExtractedData_131216/alldata_CRL.txt',header=T,stringsAsFactors=F)

numericcols <- c('volume','area','totalGAPDH','cytoGAPDH','nucGAPDH','nucArea',
                 'nucIntensityTotal','nucIntensityAvg','numCyclin','totalRNA','cytoRNA',
                 'nucRNA','numTxnSites','avgIntron','avgExon','numRnaPerTxnSiteAvgExon');


for (col in numericcols){
  data[[col]] <- as.numeric(data[[col]])
}

data$volume <- data$volume/1000

usf2 <- subset(data,gene %in% 'USF2')

fit <- lm(cytoRNA ~ volume, data = usf2)
slope <- fit$coef[[2]]
intercept <- fit$coeff[[1]]

cols = c('slateblue4','turquoise4')

#### This is all abundance

datPlot <- ggplot(usf2, aes(x=volume,y=cytoRNA))
histPlot <- ggplot(usf2, aes(x=cytoRNA))

scatter <- datPlot + 
  geom_abline(slope=slope,intercept=intercept,col='gray87') +
  geom_point(size=1.2,col=cols[1]) + 
  theme_bw() + theme_classic() + theme(legend.position = "none") + 
  expand_limits(x=c(0,6),y=c(0,125)) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  xlab('Volume (picoliter)') + ylab('USF2 mRNA Abundance') +
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6'))

hist_right <- histPlot+geom_histogram(binwidth=4,fill=cols[1])+coord_flip()+theme_classic() +
  theme(legend.position = "none") +
  expand_limits(x=c(0,125),y=0) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6'))

pdf('~/Dropbox/densitypaper/densityfigures/Figure4/abund_hist.pdf',width=3.5,height=2)
grid.arrange(scatter, hist_right, ncol=2, nrow=1, widths=c(3, 1.5), heights=c(1.7, 3.5))
dev.off()

####




### This is all concentration

datPlot <- ggplot(usf2, aes(x=volume,y=cytoRNA/volume))
histPlot <- ggplot(usf2, aes(x=cytoRNA/volume))

scatter <- datPlot + 
  geom_abline(slope=0,intercept=slope,col='gray87') +
  geom_point(size=1.2,col=cols[2]) +  
  theme_bw() + theme_classic() + theme(legend.position = "none") + 
  expand_limits(x=c(0,6),y=c(0,60)) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  xlab('Volume (picoliter)') + ylab('USF2 mRNA Concentration') + 
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6'))

hist_right <- histPlot+geom_histogram(binwidth=2,fill=cols[2])+coord_flip()+theme_classic() +
  theme(legend.position = "none") +
  expand_limits(x=c(0,60),y=0) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6'))

pdf('~/Dropbox/densitypaper/densityfigures/Figure4/concentration_hist.pdf',width=3.5,height=2)
grid.arrange(scatter, hist_right, ncol=2, nrow=1, widths=c(3, 1.5), heights=c(1.7, 3.5))
dev.off()
####