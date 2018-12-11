rm(list=ls())
library(rstan)
library(maptools)

path <- 'd:/RESEARCH/A. Waterfowl/1. Mallard/4.ipm/'

#==========================================
# load estimated and simulated populations
#==========================================
nscale <- 100

load(paste(c(path, 'results/2. estimation/ext additive.RData'), collapse=''))
his <- apply(ext$D, 2, mean)

load(paste(c(path, 'results/3. forecast/forecast_', 4, '_', -0.1, '.RData'), collapse=''))
sim1_all <- D_sim[,,dim(D_sim)[3]]

load(paste(c(path, 'results/3. forecast/forecast_', 2, '_', 0.1, '.RData'), collapse=''))
sim2_all <- D_sim[,,dim(D_sim)[3]]

#==================
# calculate scales
#==================
# density
sim1 <- apply(sim1_all, 2, mean)
sim2 <- apply(sim2_all, 2, mean)

his_scale <- round((log(his) - min(log(c(his, sim1, sim2)))) / 
             (max(log(c(his, sim1, sim2))) - min(log(c(his, sim1, sim2)))) * 
             (nscale - 1)) + 1
sim1_scale <- round((log(sim1) - min(log(c(his, sim1, sim2)))) / 
              (max(log(c(his, sim1, sim2))) - min(log(c(his, sim1, sim2)))) * 
              (nscale - 1)) + 1
sim2_scale <- round((log(sim2) - min(log(c(his, sim1, sim2)))) / 
              (max(log(c(his, sim1, sim2))) - min(log(c(his, sim1, sim2)))) * 
              (nscale - 1)) + 1

# change
c1 <- (sim1 - his) / his
c2 <- (sim2 - his) / his

c1p <- ifelse(c1 > 0, c1, NA)
c1n <- ifelse(c1 < 0, c1, NA)
c2p <- ifelse(c2 > 0, c2, NA)
c2n <- ifelse(c2 < 0, c2, NA)

c1p_scale <- round(c1p / max(c(c1p,c2p), na.rm=T) * (nscale/2 - 1)) + 1
c1n_scale <- round((c1n - min(c(c1n,c2n), na.rm=T)) / abs(min(c(c1n,c2n), na.rm=T)) * (nscale/2 - 1)) + 1
c1_scale <- numeric(length(c1))
for (i in 1:length(c1)) {
  if (is.na(c1p_scale[i])) {
    c1_scale[i] <- c1n_scale[i]
  } else {
    c1_scale[i] <- c1p_scale[i] + nscale/2
  }
}

c2p_scale <- round(c2p / max(c(c1p,c2p), na.rm=T) * (nscale/2 - 1)) + 1
c2n_scale <- round((c2n - min(c(c1n,c2n), na.rm=T)) / abs(min(c(c1n,c2n), na.rm=T)) * (nscale/2 - 1)) + 1
c2_scale <- numeric(length(c2))
for (i in 1:length(c2)) {
  if (is.na(c2p_scale[i])) {
    c2_scale[i] <- c2n_scale[i]
  } else {
    c2_scale[i] <- c2p_scale[i] + nscale/2
  }
}

# uncertainty
u1 <- apply(sim1_all, 2, sd) / his
u2 <- apply(sim2_all, 2, sd) / his

u1_scale <- round((u1 - min(c(u1,u2))) / (max(c(u1,u2)) - min(c(u1,u2))) * (nscale - 1)) + 1
u2_scale <- round((u2 - min(c(u1,u2))) / (max(c(u1,u2)) - min(c(u1,u2))) * (nscale - 1)) + 1

#===========
# draw maps
#===========
shape <- readShapeSpatial(paste(c(path, 'data/shape/minor.shp'), collapse=''))

minor <- sort(shape@data$MINOR)
bcr <- c('NIF', 'NIF', 'BTP', 'BTP', 'BTP', 'PPR', 'PPR', 'PPR', 'PPR', 
         'BSS', 'BSS', 'BSS', 'BAP', 'BAP', 'BAP', 'PPR', 'PPR')
info <- data.frame(MINOR=minor, bcr, his_scale, sim1_scale, sim2_scale, c1_scale, c2_scale, u1_scale, u2_scale)
shape@data <- merge(shape@data, info, all=T, sort=F)

dcol <- rev(heat.colors(nscale+10))[-(1:10)]
ccol <- c(colorRampPalette(c('firebrick','white'))(nscale/2+5)[1:(nscale/2)], 
          colorRampPalette(c('white','navy'))(nscale/2+5)[-(1:5)])
ucol <- colorRampPalette(c('white','darkgreen'))(nscale+10)[-(1:10)]

pdf(file=paste(c(path, 'forecast_by area.pdf'), collapse=''), width=6, height=6)

par(mfrow=c(3,3))
par(mar=c(0,0,0,0))
par(oma=c(0,2,2,0))

plot(shape, col=dcol[shape@data$his_scale], border='grey36', lwd=.6)

plot(shape, col=dcol[shape@data$sim1_scale], border='grey36', lwd=.6)

plot(shape, col=dcol[shape@data$sim2_scale], border='grey36', lwd=.6)

plot(shape, border='white')

plot(shape, col=ccol[shape@data$c1_scale], border='grey36', lwd=.6)

plot(shape, col=ccol[shape@data$c2_scale], border='grey36', lwd=.6)

plot(shape, border='white')

plot(shape, col=ucol[shape@data$u1_scale], border='grey36', lwd=.6)

plot(shape, col=ucol[shape@data$u2_scale], border='grey36', lwd=.6)

par(mfrow=c(3,1))
par(mar=c(6,0,6,0))

hist(x=1:1e6, col=dcol, border=dcol, breaks=nscale, axes=F, xlab='', main='')
hist(x=1:1e6, col=ccol, border=ccol, breaks=nscale, axes=F, xlab='', main='')
hist(x=1:1e6, col=ucol, border=ucol, breaks=nscale, axes=F, xlab='', main='')

dev.off()

range(c(his, sim1, sim2))
range(c(c1, c2))
range(c(u1, u2))


