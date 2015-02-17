source('~/Dropbox/densitypaper/densitypaperdataanalysis/multiplot.R')
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

data$cellCycleStage <- 'G1'
data$cellCycleStage[data$numCyclin>20] <- 'S'
data$cellCycleStage[data$numCyclin>230] <- 'G2'

seedFromFile <- read.table('~/Dropbox/densitypaper/densitypaperdataanalysis/densityfigures/seed_fig1_2')
.Random.seed <- t(seedFromFile)
gapdhDat <- data[runif(dim(data)[1])<.15,]

gapdhDat$cellCycleStage <- factor(gapdhDat$cellCycleStage, levels = c('G1','S','G2'))

fit <- lm(cytoGAPDH~volume,data=gapdhDat)
meanvol <- mean(data$volume)

p <- ggplot(gapdhDat, 
            aes(x=volume,y=cytoGAPDH,col=factor(cellCycleStage)))
q <- ggplot(gapdhDat, aes(x=volume,fill=cellCycleStage))
r <- ggplot(gapdhDat, aes(x=cytoGAPDH,fill=cellCycleStage))

cols <- c('darkorchid3','skyblue3','hotpink')

hist_top <- q+geom_histogram(binwidth=.2)+theme_classic() +
  theme(legend.position = "none") +
  scale_fill_manual(values=cols) +
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6'))

empty <- ggplot()+geom_point(aes(1,1), colour="white")+
  theme(axis.ticks=element_blank(), 
       panel.background=element_blank(), 
       axis.text.x=element_blank(), axis.text.y=element_blank(),           
       axis.title.x=element_blank(), axis.title.y=element_blank())

scatter <- p + geom_point(size=1.3,alpha=1) + 
  geom_abline(intercept=fit$coef[[1]],slope=fit$coef[[2]],col='gray50',linetype='dashed') +
  geom_vline(xintercept=meanvol,col='gray50',linetype='dashed') +
  scale_color_manual(values=cols) + 
  theme_bw() + theme_classic() + theme(legend.position = "none") + 
  expand_limits(x=c(0,8),y=c(0,7000)) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  xlab('Volume') + ylab('mRNA') +
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6'))

hist_right <- r+geom_histogram(binwidth=200)+coord_flip()+theme_classic() +
  theme(legend.position = "none") +
  scale_fill_manual(values=cols) +
  theme(axis.title=element_text(size="6"), axis.text=element_text(size='6'))

pdf('~/Dropbox/densitypaper/densitypaperdataanalysis/densityfigures//Figure1/OneGeneDetail/MarginalHist_Data_nolines.pdf',width=4,height=3.5)
grid.arrange(hist_top, empty, scatter, hist_right, ncol=2, nrow=2, widths=c(3, 1.4), heights=c(1.7, 3.5))
dev.off()

