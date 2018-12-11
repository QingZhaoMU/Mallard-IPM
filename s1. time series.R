rm(list=ls())
library(rstan)
library(maptools)

path <- 'd:/RESEARCH/A. Waterfowl/1. Mallard/4.ipm/'

#========================
# load posterior samples
#========================
probs <- c(.5, .1, .9)

year <- 1972:2013
year_sim <- 2013:2100

load(paste(c(path, 'results/2. estimation/ext additive.RData'), collapse=''))
est <- apply(ext$D, 2:3, quantile, probs=probs)

load(paste(c(path, 'results/3. forecast/forecast_', 4, '_', -0.1, '.RData'), collapse=''))
sim1 <- apply(D_sim, 2:3, quantile, probs=probs)

load(paste(c(path, 'results/3. forecast/forecast_', 2, '_', 0.1, '.RData'), collapse=''))
sim2 <- apply(D_sim, 2:3, quantile, probs=probs)

#==========
# graphing
#==========
minor <- c(11,12,21,22,23,31,41,51,61,71,72,73,121,122,124,131,132)
bcr <- c(1, 1, 2, 2, 2, 4, 4, 4, 4, 3, 3, 3, 5, 5, 5, 4, 4)
bcr.name <- c('NIF', 'BTP', 'BSS', 'PPR', 'BAP')
order <- c(1,2,3,4,5,6,7,8,9,16,17,13,14,15,10,11,12)

minor <- minor[order]
bcr <- bcr[order]
narea <- length(minor)

est <- est[,order,]
sim1 <- sim1[,order,]
sim2 <- sim2[,order,]

ymax <- c(500,140,60,160,200,500,600,600,400,900,600,160,300,200,140,160,100)

pdf(file=paste(c(path, 'forecast_time series.pdf'), collapse=''), width=12, height=8)

par(mfrow=c(3,6))
par(mar=c(1,2,4,1))
par(oma=c(4,4,0,0))

plot(1, type='n', axes=F, xlab='', ylab='')
legend('left', bty='n', lty=1, lwd=2, cex=1.2, 
       col=c('grey12','goldenrod','navy','firebrick'), 
       legend=c('Historic\n','Historic Mean\n',
                'Precip. -10%\nTemp. +4 deg.\n',
                'Precip. +10%\nTemp. +2 deg.\n'))

for (i in 1:narea) {
  plot(1, type='n', xlim=range(c(year,year_sim)), ylim=c(0,ymax[i]), axes=F, 
       xlab='', ylab='', main='')
  polygon(x=c(year,rev(year)), y=c(est[2,i,],rev(est[3,i,])), 
          col='grey66', border=NA)
  polygon(x=c(year_sim,rev(year_sim)), y=c(sim1[2,i,],rev(sim1[3,i,])), 
          col='skyblue', border=NA)
  polygon(x=c(year_sim,rev(year_sim)), y=c(sim2[2,i,],rev(sim2[3,i,])), 
          col=rgb(255,0,0,66,maxColorValue=255), border=NA)
  lines(est[1,i,] ~ year, lwd=2, col='grey12')
  lines(sim1[1,i,] ~ year_sim, lwd=2, col='navy')
  lines(sim2[1,i,] ~ year_sim, lwd=2, col='firebrick')
  abline(h=mean(ext$D[,order,][,i,]), col='goldenrod', lwd=2)
  axis(1, at=seq(1980,2100,40))
  axis(2, at=seq(0,ymax[i],length.out=3), las=2)
  box()
  title(main=paste(c(minor[i], ' (', bcr.name[bcr[i]], ')'), collapse=''), line=.6)
}
title(xlab='Year', outer=T, line=2, cex.lab=2)
title(ylab='Female Population Density (N / 100 km2)', outer=T, line=1.2, cex.lab=2)

dev.off()


