library(ggplot2)
library(gridExtra)

data <- read.delim('~/Dropbox/densitypaper/densitypaperdataanalysis/ExtractedData/CRL.txt',header=T,stringsAsFactors=F)

numericcols <- c('volume','area','totalGAPDH','cytoGAPDH','nucGAPDH','nucArea',
                 'nucIntensityTotal','nucIntensityAvg','numCyclin','totalRNA','cytoRNA',
                 'nucRNA','numTxnSites','avgIntron','avgExon','numRnaPerTxnSiteAvgExon');


for (col in numericcols){
  data[[col]] <- as.numeric(data[[col]])
}

data$volume <- data$volume/1000

usf2 <- subset(data,gene %in% 'TBCB')

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
  expand_limits(x=c(0,7),y=c(0,230)) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  xlab('Volume (picoliter)') + ylab('TBCB mRNA Abundance') +
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6'))

hist_right <- histPlot+geom_histogram(binwidth=8,fill=cols[1])+coord_flip()+theme_classic() +
  theme(legend.position = "none") +
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6'))

pdf('~/Dropbox/densitypaper/densitypaperdataanalysis/densityfigures/Figure6/Abund_Conc/abund_hist_TBCB.pdf',width=3.5,height=2)
grid.arrange(scatter, hist_right, ncol=2, nrow=1, widths=c(3, 1.5), heights=c(1.7, 3.5))
dev.off()