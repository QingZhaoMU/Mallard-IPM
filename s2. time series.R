rm(list=ls())
library(rstan)

path <- 'd:/RESEARCH/A. Waterfowl/1. Mallard/4.ipm/'

#========================
# load posterior samples
#========================
bpop <- read.csv(paste(c(path, 'data/data/bpop.csv'), collapse=''))
se <- read.csv(paste(c(path, 'data/data/bpop_se.csv'), collapse=''))
area <- bpop$minor
bpop <- as.matrix(bpop[,-1]) / 2
se <- as.matrix(se[,-1]) / 2
low <- bpop - 1.96 * se
upp <- bpop + 1.96 * se
low[which(low < 0)] <- 0
year <- as.numeric(substr(colnames(bpop),3,6))
year2 <- year[-length(year)]

load(paste(c(path, 'results/2. estimation/ipm additive.RData'), collapse=''))
ext <- extract(fit)
#save(ext, file=paste(c(path, 'ext additive.RData'), collapse=''))

D <- apply(ext$D, 2:3, quantile, probs=c(.5, .025, .975))
gamma <- apply(ext$gamma, 2:3, quantile, probs=c(.5, .025, .975))
omega <- apply(ext$omega, 2:3, quantile, probs=c(.5, .025, .975))

order <- c(1,2,3,4,5,6,7,8,9,16,17,13,14,15,10,11,12)
bpop <- bpop[order,]
se <- se[order,]
low <- low[order,]
upp <- upp[order,]
D <- D[,order,]
gamma <- gamma[,order,]
omega <- omega[,order,]
area <- area[order]
bcr <- c('NIF','NIF','BTP','BTP','BTP','PPR','PPR','PPR','PPR','PPR','PPR','BAP','BAP','BAP','BSS','BSS','BSS')

ymax <- c(300, 160, 140, 240, 360, 600, 700, 800, 460, 900, 700, 280, 440, 400, 200, 240, 160)

pdf(file=paste(c(path, 'bpop.pdf'), collapse=''), width=12, height=8)

par(mfrow=c(3,6))
par(mar=c(1,2,4,1))
par(oma=c(4,4,0,0))

plot(1, type='n', xlab='', ylab='', axes=F)
legend('center', col=c('navy','tomato'), pch=c(NA,19), lty=c(1,NA), 
       legend=c('IPM','Survey'), bty='n', cex=1.6, lwd=2)

for (i in 1:length(area)) {
  plot(D[1,i,] ~ year, type='n', ylim=c(0,ymax[i]), 
       xlab='', ylab='', main='', axes=F)
  polygon(x=c(year,rev(year)), y=c(D[2,i,],rev(D[3,i,])), col='lightskyblue', border=NA)
  points(bpop[i,] ~ year, pch=19, cex=.6, col='tomato')
  for (t in 1:length(year)) {
    lines(x=c(year[t],year[t]), y=c(low[i,t],upp[i,t]), col='tomato', lwd=.6)
  }
  lines(D[1,i,] ~ year, lwd=3, col='navy')
  title(main=paste(c(area[i], ' (', bcr[i], ')'), collapse=' '), line=.6)
  box()
  axis(1, at=seq(1980,2010,10))
  axis(2, at=seq(0,ymax[i],length.out=3), las=2)
}
title(xlab='Year', outer=T, line=2, cex.lab=2)
title(ylab='Female Population Density (N / 100 km2)', outer=T, line=1.2, cex.lab=2)

dev.off()

pdf(file=paste(c(path, 'recruitment-survival.pdf'), collapse=''), width=12, height=8)

par(mfrow=c(3,6))
par(mar=c(1,2,4,1))
par(oma=c(4,4,0,0))

plot(1, type='n', xlab='', ylab='', axes=F)
legend('center', col=c('navy','firebrick'), lty=1, lwd=2, cex=1.6, 
       legend=c('Recruitment','Survival'), bty='n')

for (i in 1:length(area)) {
  plot(gamma[1,i,] ~ year2, type='n', ylim=c(0,1.2), 
       xlab='', ylab='', main='', axes=F)
  polygon(x=c(year2,rev(year2)), y=c(gamma[2,i,],rev(gamma[3,i,])), 
          col='lightskyblue', border=NA)
  polygon(x=c(year2,rev(year2)), y=c(omega[2,i,],rev(omega[3,i,])), 
          col=rgb(255,0,0,100,maxColorValue=255), border=NA)
  lines(gamma[1,i,] ~ year2, lwd=2, col='navy')
  lines(omega[1,i,] ~ year2, lwd=2, col='firebrick')
  title(main=paste(c(area[i], ' (', bcr[i], ')'), collapse=' '), line=.6)
  box()
  axis(1, at=seq(1980,2010,10))
  axis(2, at=seq(0,1.2,length.out=3), las=2)
}
title(xlab='Year', outer=T, line=2, cex.lab=2)
title(ylab='Survival Rate / per capita Recruitment', outer=T, line=1.2, cex.lab=2)

dev.off()


