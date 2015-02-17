library(boot)
library(plyr)

data <- read.delim('~/Dropbox/ExtractedData_131216/alldata_CRL.txt',header=T,stringsAsFactors=F)

# Only calculate for G1, and only take genes that have >20 cells in G1
data <- data[data$numCyclin<20,]
data <- ddply(data,.(gene),subset,length(cytoRNA)>20)

coefs <- ddply(data, .(gene), function(df) {
  m <- lm(cytoRNA ~ volume, data=df)
  data.frame(intercept = coef(m)[1], slope = coef(m)[2])
})

data <- merge(data,coefs)

boot.nm <- function(mydata,i){
  cytoRNA <- mydata$cytoRNA[i]
  slope <- mydata$slope[i]
  intercept <- mydata$intercept[i]
  volume <- mydata$volume[i]
  boot.nm <- (sd(cytoRNA)/mean(cytoRNA))^2 - 
    slope*mean(volume)/(intercept+slope*mean(volume)) * cov(cytoRNA,volume)/(mean(cytoRNA)*mean(volume))
}

boots <- ddply(data, .(gene), function(df) {
  results <- boot(data=df, statistic = boot.nm, 10000)
  out <- c(results$t0[1],boot.ci(results,type='bca')$bca[4],boot.ci(results,type='bca')$bca[5]);
  names(out) <- c('est','lower','upper');
  out
})