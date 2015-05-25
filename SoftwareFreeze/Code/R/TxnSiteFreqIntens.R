setwd('~/Documents/Data/131220_CRL_MeOH_UBC/')
library(vioplot)
library(reshape2)
library(ggplot2)

data <- read.delim('TxnSiteIntensityManual.txt',header=F,stringsAsFactors=T)

#CHECK ME
colnames(data) <- c('datNum','objNum','numInt','avgInt','int1','int2','int3','int4','volume','cyclin')

sortData <- data[order(data$volume),]

dimSort <- dim(sortData)[1]
d <- ceiling(dimSort/4)
mid <- round(dimSort/2)

dat1 <- sortData[1:d,]
dat2 <- sortData[(d+1):mid,]
dat3 <- sortData[(mid+1):(dimSort-d),]
dat4 <- sortData[(dimSort-d+1):dimSort,]

g1 <- data[data$cyclin<20,]
g1 <- g1[g1$numInt<3,]
g1$stage <- 'G1'
g2 <- data[data$cyclin>230,]
g2$stage <- 'G2'

sortG1 <- g1[order(g1$volume),]
dimSortG1 <- dim(sortG1)[1]
dg1 <- ceiling(dimSortG1/4)
midg1 <- round(dimSortG1/2)

dat1_g1 <- sortG1[1:dg1,]
dat1_g1$numInt <- dat1_g1$numInt/2
dat2_g1 <- sortG1[(dg1+1):midg1,]
dat2_g1$numInt <- dat2_g1$numInt/2
dat3_g1 <- sortG1[(midg1+1):(dimSortG1-dg1),]
dat3_g1$numInt <- dat3_g1$numInt/2
dat4_g1 <- sortG1[(dimSortG1-dg1+1):dimSortG1,]
dat4_g1$numInt <- dat4_g1$numInt/2

### Small/large: transcription intensity

boxplot(dat1$avgInt,dat4$avgInt, main="Txn intensity, small and large cells", 
        xlab="GAPDH", ylab="Avg txn site intensity", 
        names = c('Smallest 1/4 cells','Largest 1/4 cells'))



df <- data.frame(
  vols = c(mean(dat1$volume),mean(dat4$volume)),
  ints = c(mean(dat1$avgInt),mean(dat4$avgInt)),
  seV = c(sd(dat1$volume)/sqrt(length(dat1$volume)),sd(dat4$volume)/sqrt(length(dat4$volume))),
  seI = c(sd(dat1$avgInt)/sqrt(length(dat1$avgInt)),sd(dat4$avgInt)/sqrt(length(dat4$avgInt)))
)

df <- data.frame(
  vols = c(mean(dat1$volume),mean(dat2$volume),mean(dat3$volume),mean(dat4$volume)),
  ints = c(mean(dat1$avgInt),mean(dat2$avgInt),mean(dat3$avgInt),mean(dat4$avgInt)),
  seV = c(sd(dat1$volume)/sqrt(length(dat1$volume)),
          sd(dat2$volume)/sqrt(length(dat2$volume)),
          sd(dat3$volume)/sqrt(length(dat3$volume)),
          sd(dat4$volume)/sqrt(length(dat4$volume))),
  seI = c(sd(dat1$avgInt)/sqrt(length(dat1$avgInt)),
          sd(dat2$avgInt)/sqrt(length(dat2$avgInt)),
          sd(dat3$avgInt)/sqrt(length(dat3$avgInt)),
          sd(dat4$avgInt)/sqrt(length(dat4$avgInt)))
)

p1 <- ggplot(df, aes(vols, ints)) +
  geom_point(size=5) +
  geom_errorbarh(aes(xmax = vols + seV, xmin = vols - seV), height = max(df$ints)/40) +
  geom_errorbar(aes(ymax = ints + seI, ymin = ints - seI), width = max(df$vols)/40) +
  xlab('GAPDH mRNA') + ylab('Txn Site Intensity') +
  theme(axis.title=element_text(size="20"), axis.text=element_text(size='15'))

t.test(dat1$avgInt,dat4$avgInt)

### Small/large: transcription frequency

boxplot(dat1$numInt,dat4$numInt, main="Txn frequency, small and large cells", 
        xlab="GAPDH", ylab="Num txn sites per cell", 
        names = c('Smallest 1/4 cells','Largest 1/4 cells'))

df2 <- data.frame(
  vols = c(mean(dat1$volume),mean(dat4$volume)),
  nums = c(mean(dat1$numInt),mean(dat4$numInt)),
  seV = c(sd(dat1$volume)/sqrt(length(dat1$volume)),sd(dat4$volume)/sqrt(length(dat4$volume))),
  seN = c(sd(dat1$numInt)/sqrt(length(dat1$numInt)),sd(dat4$numInt)/sqrt(length(dat4$numInt)))
)

df2 <- data.frame(
  vols = c(mean(dat1$volume),mean(dat2$volume),mean(dat3$volume),mean(dat4$volume)),
  nums = c(mean(dat1$numInt),mean(dat2$numInt),mean(dat3$numInt),mean(dat4$numInt)),
  seV = c(sd(dat1$volume)/sqrt(length(dat1$volume)),
          sd(dat2$volume)/sqrt(length(dat2$volume)),
          sd(dat3$volume)/sqrt(length(dat3$volume)),
          sd(dat4$volume)/sqrt(length(dat4$volume))),
  seN = c(sd(dat1$numInt)/sqrt(length(dat1$numInt)),
          sd(dat2$numInt)/sqrt(length(dat2$numInt)),
          sd(dat3$numInt)/sqrt(length(dat3$numInt)),
          sd(dat4$numInt)/sqrt(length(dat4$numInt)))
)

p <- ggplot(df2, aes(vols, nums))
p + geom_point(size=5) + expand_limits(y=0) +
  geom_errorbarh(aes(xmax = vols + seV, xmin = vols - seV), height = max(df2$nums)/40) +
  geom_errorbar(aes(ymax = nums + seN, ymin = nums - seN,), width = max(df2$vols)/40) +
  xlab('GAPDH mRNA') + ylab('Txn Sites Per Cell') +
  theme(axis.title=element_text(size="20"), axis.text=element_text(size='15'))

# Do this again, only with cells in G1:

dfg1 <- data.frame(
  vols = c(mean(dat1_g1$volume),mean(dat2_g1$volume),mean(dat3_g1$volume),mean(dat4_g1$volume)),
  nums = c(mean(dat1_g1$numInt),mean(dat2_g1$numInt),mean(dat3_g1$numInt),mean(dat4_g1$numInt)),
  seV = c(sd(dat1_g1$volume)/sqrt(length(dat1_g1$volume)),
          sd(dat2_g1$volume)/sqrt(length(dat2_g1$volume)),
          sd(dat3_g1$volume)/sqrt(length(dat3_g1$volume)),
          sd(dat4_g1$volume)/sqrt(length(dat4_g1$volume))),
  seN = c(sd(dat1_g1$numInt)/sqrt(length(dat1_g1$numInt)),
          sd(dat2_g1$numInt)/sqrt(length(dat2_g1$numInt)),
          sd(dat3_g1$numInt)/sqrt(length(dat3_g1$numInt)),
          sd(dat4_g1$numInt)/sqrt(length(dat4_g1$numInt)))
)

p2 <- ggplot(dfg1, aes(vols, nums)) +
  geom_point(size=5) + expand_limits(y=0) +
  geom_errorbarh(aes(xmax = vols + seV, xmin = vols - seV), height = max(dfg1$nums)/40) +
  geom_errorbar(aes(ymax = nums + seN, ymin = nums - seN,), width = max(dfg1$vols)/40) +
  xlab('GAPDH mRNA') + ylab('Txn Sites Per Allele') +
  theme(axis.title=element_text(size="20"), axis.text=element_text(size='15'))

t.test(dat1$numInt,dat4$numInt)

### G1/G2: transcription intensity

dat <- rbind(g1,g2)

p <- ggplot(dat,aes(factor(stage),avgInt))
p + geom_boxplot(size=.75) + 
  xlab('Cell Cycle Stage') + ylab('Txn Site Intensity') +
  theme(axis.title=element_text(size="20"), axis.text=element_text(size='15'))

df4 <- data.frame(
  vols = c('G1','G2'),
  nums = c(mean(g1$avgInt),mean(g2$avgInt)),
  seN = c(sd((g1$avgInt))/sqrt(length(g1$avgInt)),sd((g2$avgInt))/sqrt(length(g2$avgInt)))
)

p3 <- ggplot(df4, aes(vols, nums)) +
  geom_point(size=5) +
  geom_errorbar(aes(ymax = nums + seN, ymin = nums - seN,), width = .2) +
  xlab('Cell Cycle Stage') + ylab('Txn Site Intensity') + expand_limits(y=0) +
  theme(axis.title=element_text(size="20"), axis.text=element_text(size='15'))

boxplot(g1$avgInt,g2$avgInt, main="Txn intensity, G1 and G2", 
        xlab="GAPDH", ylab="Avg txn site intensity", 
        names = c('G1','G2'))

t.test(g1$avgInt,g2$avgInt)

### G1/G2: transcription frequency

p <- ggplot(dat,aes(factor(stage),numInt))
p + geom_boxplot(size=.75) + 
  xlab('Cell Cycle Stage') + ylab('Txn Sites Per Cell') +
  theme(axis.title=element_text(size="20"), axis.text=element_text(size='15'))

boxplot(g1$numInt/2,g2$numInt/4, main="Txn frequency per allele, G1 and G2", 
        ylab="Num txn sites per allele", 
        names = c('G1','G2'))

df3 <- data.frame(
  vols = c('G1','G2'),
  nums = c(mean(g1$numInt)/2,mean(g2$numInt)/4),
  seN = c(sd((g1$numInt)/2)/sqrt(length(g1$numInt)),sd((g2$numInt)/4)/sqrt(length(g2$numInt)))
)

p4 <- ggplot(df3, aes(vols, nums)) +
  geom_point(size=5) + 
  geom_errorbar(aes(ymax = nums + seN, ymin = nums - seN,), width = .2) +
  xlab('Cell Cycle Stage') + ylab('Txn Sites Per Allele') + expand_limits(y=0) +
  theme(axis.title=element_text(size="20"), axis.text=element_text(size='15'))

t.test(g1$numInt/2,g2$numInt/4)

### G1/G2: volume

boxplot(g1$volume,g2$volume, main="Volume, G1 and G2", 
        ylab="Volume", 
        names = c('G1','G2'))

t.test(g1$volume,g2$volume)



### Arrange data by individual transcription site

meltDat <- melt.data.frame(data=data,measure.vars=c('int1','int2','int3','int4'))
tmp <- meltDat[meltDat$value>0,]
plot(tmp$volume,tmp$value)
cor(tmp$volume,tmp$value)