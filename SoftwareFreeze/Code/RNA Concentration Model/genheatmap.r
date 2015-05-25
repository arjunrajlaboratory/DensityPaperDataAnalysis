setwd('~/Dropbox/Code/RNA Concentration Model/')
data10 <- read.delim('abund10.txt',header=F,stringsAsFactors=F)
data10 <- as.matrix(data10)

x <- 1:ncol(data10)
y <- 1:nrow(data10)

library("RColorBrewer", lib.loc="/Library/Frameworks/R.framework/Versions/2.15/Resources/library")
colors <- brewer.pal(4, "BuPu")
pal <- colorRampPalette(colors)

filled.contour(x,y,t(data10),col=pal(18))
points(c(42.86),c(60),col=pal(18)[[6]],pch=19,cex=2)

####
data100 <- read.delim('abund100.txt',header=F,stringsAsFactors=F)
data100 <- as.matrix(data100)

x <- .1*1:ncol(data100)
y <- (log10(24)/6)*1:nrow(data100)

library("RColorBrewer", lib.loc="/Library/Frameworks/R.framework/Versions/2.15/Resources/library")
colors <- brewer.pal(4, "BuPu")
pal <- colorRampPalette(colors)

myGenes100 <- myGenesTotal[myGenesTotal$avg_abund>50,]
myGenes100 <- myGenes100[myGenes100$avg_abund<500,]

library(calibrate)
filled.contour(x,y,t(data100),col=pal(18),
               xlab='intron spot frequency',
               ylab='log(halflife)',
               main='Expresison level = 100')
points(myGenes100$neg-.05,log10(myGenes100$halflife),
       pch=19,cex=2,col=as.matrix(pal(18))[round(myGenes100$corr_gapdh*18),],)
points(myGenes100$neg-.05,log10(myGenes100$halflife),cex=2)
text(c())

###
library(fields)
library(akima)

s <- interp(data10)