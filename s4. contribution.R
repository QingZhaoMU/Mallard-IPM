rm(list=ls())
library(boot)
library(hier.part)

path <- 'd:/RESEARCH/A. Waterfowl/1. Mallard/4.ipm/'

#========================
# load posterior samples
#========================
bpop <- read.csv(paste(c(path, 'data/data/bpop.csv'), collapse=''))
area <- bpop$minor
year <- as.numeric(substr(colnames(bpop[,-1]),3,6))
narea <- length(area)
nyear <- length(year)

load(paste(c(path, 'results/2. estimation/ext additive.RData'), collapse=''))
nmcmc <- dim(ext$D)[1]

order <- c(1,2,3,4,5,10,11,12,6,7,8,9,16,17,13,14,15)
area <- area[order]

D_est <- apply(ext$D[,order,], 2:3, median)
R_est <- apply(ext$R[,order,], 2:3, median)
S_est <- apply(ext$S[,order,], 2:3, median)

lg <- log(D_est[,-1] / D_est[,-nyear])
recruit <- apply(log(ext$gamma), 2:3, median)
survival <- apply(logit(ext$omega[,,-nyear]), 2:3, median)

D <- G <- R <- S <- numeric(narea)
for (i in 1:narea) {
  D[i] <- paste(c(
    format(round(median(D_est[i,]), digits=1), nsmall=1), 
    ' [', 
    format(round(min(D_est[i,]), digits=1), nsmall=1), 
    ', ', 
    format(round(max(D_est[i,]), digits=1), nsmall=1), 
    ']'), 
    collapse='')

  G[i] <- paste(c(
    format(round(median(exp(lg[i,])), digits=3), nsmall=3), 
    ' [', 
    format(round(min(exp(lg[i,])), digits=3), nsmall=3), 
    ', ', 
    format(round(max(exp(lg[i,])), digits=3), nsmall=3), 
    ']'), 
    collapse='')

  R[i] <- paste(c(
    format(round(median(exp(recruit[i,])), digits=3), nsmall=3), 
    ' [', 
    format(round(min(exp(recruit[i,])), digits=3), nsmall=3), 
    ', ', 
    format(round(max(exp(recruit[i,])), digits=3), nsmall=3), 
    ']'), 
    collapse='')

  S[i] <- paste(c(
    format(round(median(inv.logit(survival[i,])), digits=3), nsmall=3), 
    ' [', 
    format(round(min(inv.logit(survival[i,])), digits=3), nsmall=3), 
    ', ', 
    format(round(max(inv.logit(survival[i,])), digits=3), nsmall=3), 
    ']'), 
    collapse='')
}

ntest <- nmcmc
rcont_mat <- scont_mat <- matrix(, narea, ntest)

for (i in 1:narea) {
  for (k in 1:ntest) {
    x <- data.frame(ext$gamma[k,i,], ext$omega[k,i,-nyear])
    names(x) <- c('recruitment', 'survival')
    y <- log(ext$D[k,i,-1] / ext$D[k,i,-nyear])
    hpart <- hier.part(y=y, xcan=x, family="gaussian", gof="RMSPE", barplot=F)$I.perc
    rcont_mat[i,k] <- hpart[which(rownames(hpart)=='recruitment'),1]
    scont_mat[i,k] <- hpart[which(rownames(hpart)=='survival'),1]

    print(paste(i, k, sep='_'))
  } # k
} # i

rcont_qt <- apply(rcont_mat, 1, quantile, probs=c(.5, .025, .975))
scont_qt <- apply(scont_mat, 1, quantile, probs=c(.5, .025, .975))
rcont <- scont <- character(narea)
for (i in 1:narea) {
  rcont[i] <- paste(c(
    format(round(rcont_qt[1,i],digits=1), nsmall=1), '% ', 
    '[', format(round(rcont_qt[2,i],digits=1), nsmall=1), '%, ', 
    format(round(rcont_qt[3,i],digits=1), nsmall=1), '%]'), collapse='')
  scont[i] <- paste(c(
    format(round(scont_qt[1,i],digits=1), nsmall=1), '% ', 
    '[', format(round(scont_qt[2,i],digits=1), nsmall=1), '%, ', 
    format(round(scont_qt[3,i],digits=1), nsmall=1), '%]'), collapse='')
}

bcr <- c(rep('Northwestern Interior Forest', 2), rep('Boreal Taiga Plains', 3), rep('Boreal Softwood Shield', 3), 
         rep('Prairie Potholes', 6), rep('Bandlands and Prairies', 3))

out2 <- data.frame(
  area=area, 
  bcr=bcr, 
  D=D, 
  G=G, 
  R=R, 
  RCont=rcont, 
  S=S, 
  SCont=scont)
out2

write.csv(out2, file=paste(c(path, 'contribution.csv'), collapse=''), row.names=F)

mean(exp(recruit))
range(exp(recruit))
mean(rcont_mat)
range(rcont_mat)

mean(inv.logit(survival))
range(inv.logit(survival))
mean(scont_mat)
range(scont_mat)

rsratio <- rcont_mat / scont_mat
rsratio_med <- apply(rsratio, 1, median)
range(rsratio_med)


