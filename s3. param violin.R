rm(list=ls())
library(boot)
library(ggplot2)

path <- 'd:/RESEARCH/A. Waterfowl/1. Mallard/4.ipm/'

#========================
# load posterior samples
#========================
load(paste(c(path, 'results/2. estimation/ext additive.RData'), collapse=''))
nmcmc <- dim(ext$D)[1]

#==============
# violin polot
#==============
bpop <- read.csv(paste(c(path, 'data/data/bpop.csv'), collapse=''))
order <- c(1,2,3,4,5,10,11,12,6,7,8,9,16,17,13,14,15)
area <- bpop$minor[order]
narea <- length(area)
area_order <- paste(letters[1:length(area)], area, sep='.')
cols <- c(rep('tomato',2),rep('goldenrod',3),rep('lightpink',3),rep('darkgreen',6),rep('steelblue',3))

data_summary <- function(x) {
   m <- mean(x)
   ymin <- as.numeric(quantile(x, .025))
   ymax <- as.numeric(quantile(x, .975))
   return(c(y=m,ymin=ymin,ymax=ymax))
}

gamma_mean <- exp(ext$gamma_alpha_area[,order])
gamma_mean <- data.frame(value=as.vector(gamma_mean), area=as.factor(rep(area_order, each=nmcmc)))

gamma_bpop <- ext$gamma_bpop_area[,order]
gamma_bpop <- data.frame(value=as.vector(gamma_bpop), area=as.factor(rep(area_order, each=nmcmc)))

gamma_prec <- ext$gamma_prec_area[,order]
gamma_prec <- data.frame(value=as.vector(gamma_prec), area=as.factor(rep(area_order, each=nmcmc)))

gamma_temp <- ext$gamma_temp_area[,order]
gamma_temp <- data.frame(value=as.vector(gamma_temp), area=as.factor(rep(area_order, each=nmcmc)))

omega_mean <- inv.logit(ext$omega_alpha_area[,order])
omega_mean <- data.frame(value=as.vector(omega_mean), area=as.factor(rep(area_order, each=nmcmc)))

omega_bpop <- ext$omega_bpop_area[,order]
omega_bpop <- data.frame(value=as.vector(omega_bpop), area=as.factor(rep(area_order, each=nmcmc)))

omega_prec <- ext$omega_prec_area[,order]
omega_prec <- data.frame(value=as.vector(omega_prec), area=as.factor(rep(area_order, each=nmcmc)))

omega_temp <- ext$omega_temp_area[,order]
omega_temp <- data.frame(value=as.vector(omega_temp), area=as.factor(rep(area_order, each=nmcmc)))

gamma_bpop_p <- numeric(narea)
for (i in 1:narea) {
  tt <- gamma_bpop[which(gamma_bpop$area == area_order[i]), 1]
  pp <- length(which(tt < 0)) / nmcmc
  gamma_bpop_p[i] <- ifelse(pp < .5, pp, 1-pp)
}
cbind(area_order, gamma_bpop_p)

gamma_prec_p <- numeric(narea)
for (i in 1:narea) {
  tt <- gamma_prec[which(gamma_prec$area == area_order[i]), 1]
  pp <- length(which(tt < 0)) / nmcmc
  gamma_prec_p[i] <- ifelse(pp < .5, pp, 1-pp)
}
cbind(area_order, gamma_prec_p)
range(gamma_prec_p[1:8])

gamma_temp_p <- numeric(narea)
for (i in 1:narea) {
  tt <- gamma_temp[which(gamma_temp$area == area_order[i]), 1]
  pp <- length(which(tt < 0)) / nmcmc
  gamma_temp_p[i] <- ifelse(pp < .5, pp, 1-pp)
}
cbind(area_order, gamma_temp_p)
range(gamma_temp_p[6:17])
range(gamma_temp_p[1:5])

omega_bpop_p <- numeric(narea)
for (i in 1:narea) {
  tt <- omega_bpop[which(omega_bpop$area == area_order[i]), 1]
  pp <- length(which(tt < 0)) / nmcmc
  omega_bpop_p[i] <- ifelse(pp < .5, pp, 1-pp)
}
cbind(area_order, omega_bpop_p)

omega_prec_p <- numeric(narea)
for (i in 1:narea) {
  tt <- omega_prec[which(omega_prec$area == area_order[i]), 1]
  pp <- length(which(tt < 0)) / nmcmc
  omega_prec_p[i] <- ifelse(pp < .5, pp, 1-pp)
}
cbind(area_order, omega_prec_p)
range(omega_prec_p)

omega_temp_p <- numeric(narea)
for (i in 1:narea) {
  tt <- omega_temp[which(omega_temp$area == area_order[i]), 1]
  pp <- length(which(tt < 0)) / nmcmc
  omega_temp_p[i] <- ifelse(pp < .5, pp, 1-pp)
}
cbind(area_order, omega_temp_p)
range(omega_temp_p)

pdf(file=paste(c(path, 'violin.pdf'), collapse=''), width=5, height=5)

ggplot(gamma_mean, aes(x=area, y=value, color=area, fill=area)) + 
  geom_violin(trim=T) + 
  scale_color_manual(values=rep('grey36',17)) + 
  scale_fill_manual(values=cols) + 
  theme(legend.position="none") + 
  labs(title='Mean', x='Area', y='Posterior Value') + 
  stat_summary(fun.data=data_summary, col='grey6', size=.1) + 
  ylim(c(.3,.5)) + 
  scale_x_discrete(labels=area)

ggplot(gamma_bpop, aes(x=area, y=value, color=area, fill=area)) + 
  geom_violin(trim=T) + 
  scale_color_manual(values=rep('grey36',17)) + 
  scale_fill_manual(values=cols) + 
  theme(legend.position="none") + 
  labs(title='Density Dependence', x='Area', y='Posterior Value') + 
  stat_summary(fun.data=data_summary, col='grey6', size=.1) + 
  ylim(c(-.5,.5)) + 
  scale_x_discrete(labels=area)

ggplot(gamma_prec, aes(x=area, y=value, color=area, fill=area)) + 
  geom_violin(trim=T) + 
  scale_color_manual(values=rep('grey36',17)) + 
  scale_fill_manual(values=cols) + 
  theme(legend.position="none") + 
  labs(title='Effect of Precipitation', x='Area', y='Posterior Value') + 
  stat_summary(fun.data=data_summary, col='grey6', size=.1) + 
  ylim(c(-.5,.5)) + 
  scale_x_discrete(labels=area)

ggplot(gamma_temp, aes(x=area, y=value, color=area, fill=area)) + 
  geom_violin(trim=T) + 
  scale_color_manual(values=rep('grey36',17)) + 
  scale_fill_manual(values=cols) + 
  theme(legend.position="none") + 
  labs(title='Effect of Temperature', x='Area', y='Posterior Value') + 
  stat_summary(fun.data=data_summary, col='grey6', size=.1) + 
  ylim(c(-.5,.5)) + 
  scale_x_discrete(labels=area)

ggplot(omega_mean, aes(x=area, y=value, color=area, fill=area)) + 
  geom_violin(trim=T) + 
  scale_color_manual(values=rep('grey36',17)) + 
  scale_fill_manual(values=cols) + 
  theme(legend.position="none") + 
  labs(title='Mean', x='Area', y='Posterior Value') + 
  stat_summary(fun.data=data_summary, col='grey6', size=.1) + 
  ylim(c(.5,.7)) + 
  scale_x_discrete(labels=area)

ggplot(omega_bpop, aes(x=area, y=value, color=area, fill=area)) + 
  geom_violin(trim=T) + 
  scale_color_manual(values=rep('grey36',17)) + 
  scale_fill_manual(values=cols) + 
  theme(legend.position="none") + 
  labs(title='Density Dependence', x='Area', y='Posterior Value') + 
  stat_summary(fun.data=data_summary, col='grey6', size=.1) + 
  ylim(c(-.5,.5)) + 
  scale_x_discrete(labels=area)

ggplot(omega_prec, aes(x=area, y=value, color=area, fill=area)) + 
  geom_violin(trim=T) + 
  scale_color_manual(values=rep('grey36',17)) + 
  scale_fill_manual(values=cols) + 
  theme(legend.position="none") + 
  labs(title='Effect of Precipitation', x='Area', y='Posterior Value') + 
  stat_summary(fun.data=data_summary, col='grey6', size=.1) + 
  ylim(c(-.5,.5)) + 
  scale_x_discrete(labels=area)

ggplot(omega_temp, aes(x=area, y=value, color=area, fill=area)) + 
  geom_violin(trim=T) + 
  scale_color_manual(values=rep('grey36',17)) + 
  scale_fill_manual(values=cols) + 
  theme(legend.position="none") + 
  labs(title='Effect of Temperature', x='Area', y='Posterior Value') + 
  stat_summary(fun.data=data_summary, col='grey6', size=.1) + 
  ylim(c(-.5,.5)) + 
  scale_x_discrete(labels=area)

plot(1, type='n', xlab='', ylab='', axes=F)
legend('left', col=c('tomato','goldenrod','lightpink','darkgreen','steelblue'), 
       pch=15, cex=1.6, bty='n', 
       legend=c('Northwestern Interior Forest',
                'Boreal Taiga Plains','Boreal Softwood Shield',
                'Prairie Potholes','Bandlands and Prairies'))

dev.off()













