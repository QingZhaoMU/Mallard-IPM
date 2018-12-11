rm(list=ls())
library(rstan)

path <- 'd:/RESEARCH/A. Waterfowl/1. Mallard/4.ipm/'

#========================
# load posterior samples
#========================
bpop <- read.csv(paste(c(path, 'data/data/bpop.csv'), collapse=''))
se <- read.csv(paste(c(path, 'data/data/bpop_se.csv'), collapse=''))
area <- bpop$minor
bcr <- c(1,1,2,2,2,4,4,4,4,3,3,3,5,5,5,4,4)
bpop <- as.matrix(bpop[,-1]) / 2
se <- as.matrix(se[,-1]) / 2
narea <- dim(bpop)[1]
nyear <- dim(bpop)[2]

load(paste(c(path, 'results/2. estimation/ext additive.RData'), collapse=''))
nmcmc <- dim(ext$D)[1]

dsim <- dobs <- array(, dim=c(nmcmc, narea, nyear-1))
for (i in 1:narea) {
  for (t in 1:(nyear-1)) {
    exp <- ext$D[,i,t]
    sim <- rnorm(nmcmc, exp, se[i,t])
    sim[which(sim < 0)] <- NA
    dsim[,i,t] <- (sim - exp) ^ 2 / var(exp)
    dobs[,i,t] <- (bpop[i,t] - exp) ^ 2 / var(exp)
  }
}

dobs[which(is.na(dsim))] <- NA

dsim_area <- apply(dsim, 1:2, sum, na.rm=T)
dobs_area <- apply(dobs, 1:2, sum, na.rm=T)

cols <- c('tomato','goldenrod','lightpink','darkgreen','steelblue')
ppp <- length(which(dsim_area > dobs_area)) / nmcmc / narea

order <- c(1,2,3,4,5,10,11,12,6,7,8,9,16,17,13,14,15)
bcr <- bcr[order]
dobs_area <- dobs_area[,order]
dsim_area <- dsim_area[,order]

pdf(file=paste(c(path, 'ppp.pdf'), collapse=''), width=6, height=6)
par(mar=c(4,4,4,1))

plot(1, xlim=c(0,360), ylim=c(0,360), type='n', axes=F, 
     xlab='', ylab='', main='')
for (i in 1:narea) {
  points(dobs_area[,i], dsim_area[,i], col=cols[bcr[i]])
}
axis(1, at=seq(0, 360, length.out=3))
axis(2, at=seq(0, 360, length.out=3), las=1)
box()
abline(0, 1, col='grey36', lwd=2)
title(xlab='Discrepancy of Observations', cex.lab=1.6, line=2.6)
title(ylab='Discrepancy of Simulations', cex.lab=1.6, line=2.6)
title(main='Posterior Predictive Check', cex.main=2)
text(x=300, y=0, cex=1.6, 
     labels=paste(c('P = ', round(ppp, digits=3)), collapse=''))

dev.off()


