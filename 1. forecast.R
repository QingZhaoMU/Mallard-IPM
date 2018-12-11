rm(list=ls())
library(rstan)
library(boot)

path <- 'd:/RESEARCH/A. Waterfowl/1. Mallard/4.ipm/'

temp_inc <- 2 # change in temperature degree by 2100
prec_inc <- .1 # change in precipitation proportion by 2100
nyear_sim <- 88 # number of years from 2013 to 2100

#========================
# load posterior samples
#========================
load(paste(c(path, 'results/2. estimation/ext additive.RData'), collapse=''))
nmcmc <- dim(ext$D)[1]

#==============
# read in data
#==============
temp <- read.csv(paste(c(path, 'data/data/temp.csv'), collapse=''))
prec <- read.csv(paste(c(path, 'data/data/prec.csv'), collapse=''))
area <- temp$minor

temp <- as.matrix(temp[,-1])
prec <- as.matrix(prec[,-1])
narea <- dim(temp)[1]
nyear <- dim(temp)[2]

#============
# simulation
#============
# simulate climate change
temp_sim <- prec_sim <- array(, dim=c(nmcmc, narea, nyear_sim))
for (k in 1:nmcmc) {
  for (i in 1:narea) {
    temp_sim[k,i,] <- sample(temp[i,], nyear_sim, replace=T) + seq(0,temp_inc,length.out=nyear_sim)
    prec_sim[k,i,] <- sample(prec[i,], nyear_sim, replace=T) + seq(0,mean(prec[i,])*prec_inc,length.out=nyear_sim)
  } # i
} # k

#par(mfrow=c(3,6))
#par(mar=c(1,1,2,1))
#for (i in 1:narea) {
#  plot(temp[i,] ~ c(1:nyear), type='l', axes=F, 
#       xlim=c(1, nyear+nyear_sim), ylim=range(c(temp[i,], temp_sim[k,i,])))
#  box()
#  lines(temp_sim[k,i,] ~ c((nyear+1):(nyear+nyear_sim)), col=2)
#  title(main=area[i])
#}

#par(mfrow=c(3,6))
#par(mar=c(1,1,2,1))
#for (i in 1:narea) {
#  plot(prec[i,] ~ c(1:nyear), type='l', axes=F, 
#       xlim=c(1, nyear+nyear_sim), ylim=range(c(prec[i,], prec_sim[k,i,])))
#  box()
#  lines(prec_sim[k,i,] ~ c((nyear+1):(nyear+nyear_sim)), col=2)
#  title(main=area[i])
#}

for (i in 1:narea) {
  temp_sim[,i,] <- (temp_sim[,i,] - mean(temp[i,])) / sd(temp[i,])
  prec_sim[,i,] <- (log(prec_sim[,i,]) - mean(log(prec[i,]))) / sd(log(prec[i,]))
}

# simulate population dynamics
logD_mean <- apply(log(ext$D), 1:2, mean)
logD_sd <- apply(log(ext$D), 1:2, sd)
D_max <- apply(ext$D, 2, max)

D_sim <- array(, dim=c(nmcmc, narea, nyear_sim))

for (i in 1:narea) {
  D_sim[,i,1] <- ext$D[,i,dim(ext$D)[3]]

  for (t in 2:nyear_sim) {
    gamma <- exp(ext$gamma_alpha_area[,i] + 
                 ext$gamma_bpop_area[,i] * (log(D_sim[,i,t-1]) - logD_mean[,i]) / logD_sd[,i] + 
                 ext$gamma_temp_area[,i] * temp_sim[,i,t] + 
                 ext$gamma_prec_area[,i] * prec_sim[,i,t] + 
                 rnorm(nmcmc, 0, ext$gamma_sigma))
    R <- rnorm(nmcmc, D_sim[,i,t-1]*gamma, sqrt(D_sim[,i,t-1]*gamma))
    R[which(R < 0)] <- 1

    omega <- inv.logit(ext$omega_alpha_area[,i] + 
                       ext$omega_bpop_area[,i] * (log(D_sim[,i,t-1]) - logD_mean[,i]) / logD_sd[,i] + 
                       ext$omega_temp_area[,i] * temp_sim[,i,t] + 
                       ext$omega_prec_area[,i] * prec_sim[,i,t] + 
                       rnorm(nmcmc, 0, ext$omega_sigma))
    S <- rnorm(nmcmc, D_sim[,i,t-1]*omega, sqrt(D_sim[,i,t-1]*omega*(1-omega)))
    S[which(S < 0)] <- 1

    D <- R + S
    D[which(D > D_max[i]*2)] <- D_max[i]*2
    D_sim[,i,t] <- D
  } # t
} # i

probs <- c(.5, .01, .9)
D_sim_qt <- apply(D_sim, 2:3, quantile, probs=probs)
D_est_qt <- apply(ext$D, 2:3, quantile, probs=probs)

year <- 1972:2013
year_sim <- 2013:2100

minor <- c(11,12,21,22,23,31,41,51,61,71,72,73,121,122,124,131,132)
bcr <- c(1, 1, 2, 2, 2, 4, 4, 4, 4, 3, 3, 3, 5, 5, 5, 4, 4)
area.name <- c('Northwestern Interior Forest','Boreal Taiga Plains',
               'Boreal Softwood Shield','Prairie Potholes','Badlands and Prairies')
bcr.name <- c('NIF', 'BTP', 'BSS', 'PPR', 'BAP')
order <- c(1,2,3,4,5,6,7,8,9,16,17,13,14,15,10,11,12)

minor <- minor[order]
bcr <- bcr[order]

D_sim_qt <- D_sim_qt[,order,]
D_est_qt <- D_est_qt[,order,]
ymax <- c(500,140,60,160,200,500,600,600,400,900,600,160,300,200,140,160,100)

pdf(file=paste(c(path, 'forecast_', temp_inc, '_', prec_inc, '_temp.pdf'), collapse=''), width=8, height=6)

par(mfrow=c(3,6))
par(mar=c(2,2,2,1))
par(oma=c(2,2,0,0))

plot(1, type='n', axes=F, xlab='', ylab='')
legend('left', bty='n', col=c('grey12','goldenrod','navy'), lty=1, lwd=2, cex=1.2, 
       legend=c('Historic','Historic Mean','Forecast'))

for (i in 1:narea) {
  plot(1, type='n', xlim=range(c(year,year_sim)), ylim=c(0,ymax[i]), axes=F, 
       xlab='', ylab='', main='')
  polygon(x=c(year,rev(year)), y=c(D_est_qt[2,i,],rev(D_est_qt[3,i,])), col='grey66', border=NA)
  polygon(x=c(year_sim,rev(year_sim)), y=c(D_sim_qt[2,i,],rev(D_sim_qt[3,i,])), col='lightskyblue', border=NA)
  lines(D_est_qt[1,i,] ~ year, lwd=2, col='grey12')
  lines(D_sim_qt[1,i,] ~ year_sim, lwd=2, col='navy')
  abline(h=mean(ext$D[,order,][,i,]), col='goldenrod', lwd=2)
  if (i > 11) {
    axis(1, at=seq(1980,2100,40))
  } else {
    axis(1, at=seq(1980,2100,40), labels=rep('',4))
  }
  axis(2, at=seq(0,ymax[i],length.out=3), las=2)
  box()
  title(main=paste(c(minor[i], ' (', bcr.name[bcr[i]], ')'), collapse=''))
}

dev.off()

#=====================
# save the simulation
#=====================
save(D_sim, file=paste(c(path, 'forecast_', temp_inc, '_', prec_inc, '.RData'), collapse=''))


