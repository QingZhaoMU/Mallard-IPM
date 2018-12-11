rm(list=ls())
library(boot)
library(rstan)

path <- 'c:/Zhao/1. Mallard/4.ipm/'

#==============
# read in data
#==============
# bpop
bpop <- read.csv(paste(c(path, 'data/data/bpop.csv'), collapse=''))
bpop_se <- read.csv(paste(c(path, 'data/data/bpop_se.csv'), collapse=''))

area <- bpop$minor
year <- substr(names(bpop)[-1], 3, 6)
narea <- length(area)
nyear <- length(year)

bcr <- c(1, 1, 2, 2, 2, 4, 4, 4, 4, 3, 3, 3, 5, 5, 5, 4, 4)
nbcr <- max(bcr)

bpop_obs <- as.matrix(bpop[,-1]) / 2
bpop_obs_sigma <- as.matrix(bpop_se[,-1]) / 2

# climate
temp <- read.csv(paste(c(path, 'data/data/temp.csv'), collapse=''))
prec <- read.csv(paste(c(path, 'data/data/prec.csv'), collapse=''))

temp <- as.matrix(temp[,-1])
temp0 <- (temp[,1] - mean(temp[,1])) / sd(temp[,1])
for (i in 1:narea) {
  temp[i,] <- (temp[i,] - mean(temp[i,])) / sd(temp[i,])
}

prec <- log(as.matrix(prec[,-1]))
prec0 <- (prec[,1] - mean(prec[,1])) / sd(prec[,1])
for (i in 1:narea) {
  prec[i,] <- (prec[i,] - mean(prec[i,])) / sd(prec[i,])
}

# band-recovery
load(paste(c(path, 'data/data/marr.RData'), collapse=''))
marr <- marr$af0

#==============
# define model
#==============
sink(paste(c(path, 'ipm additive.stan'), collapse=''))
cat("

functions {
  real getprod(vector X, int beg, int end) {
    real p;     
    p <- 1;
    for (j in beg:end) p <- p * X[j];
    return p;
  }

  real getsum(vector X, int beg, int end) {
    real p;     
    p <- 0;
    for (j in beg:end) p <- p + X[j];
    return p;
  }
}

data {
  # basic
  int narea; # number of reference areas
  int nyear; # number of banding years
  int nbcr; # nubmer of ecoregions
  int bcr[narea]; # ecoregion for each reference area

  # bpop
  matrix[narea, nyear] bpop_obs; # estimated breeding population density
  matrix[narea, nyear] bpop_obs_sigma; # SE of estimated breeding population density

  # climate
  matrix[narea, nyear] temp; # standardized spring temperature
  matrix[narea, nyear] prec; # standardized winter precipitation
  vector[narea] temp0; # spring temperature in the first year
  vector[narea] prec0; # winter precipitation in the first year

  # band-recovery
  int marr[narea, nyear, nyear+1]; # band-recovery data, standard banding
} # data

parameters {
  # initial population density
  real D0_alpha; # intercept for initial population density
  real D0_temp; # temperature effect on initial population density
  real D0_prec; # precipitation effect on initial population density
  real<lower=0> D0_sigma; # SD of process error in initial population density
  vector[narea] log_D0; # log density of the first year

  # reproduction
  real gamma_alpha_mu; # grand mean of reproduction intercept
  real<lower=0> gamma_alpha_sigma_bcr; # SD of reproduction intercept for ecoregions
  real<lower=0> gamma_alpha_sigma_area; # SD of reproduction intercept for reference areas
  real gamma_alpha_bcr[nbcr]; # mean of reproduction intercept for ecoregions
  real gamma_alpha_area[narea]; # reproduction intercept for reference areas

  real gamma_bpop_mu; # grand mean of reproduction density dependence
  real<lower=0> gamma_bpop_sigma_bcr; # SD of reproduction density dependence for ecoregions
  real<lower=0> gamma_bpop_sigma_area; # SD of reproduction density dependence for reference areas
  real gamma_bpop_bcr[nbcr]; # mean of reproduction density dependence for ecoregions
  real gamma_bpop_area[narea]; # reproduction density dependence for reference areas

  real gamma_temp_mu; # grand mean of reproduction temperature effect 
  real<lower=0> gamma_temp_sigma_bcr; # SD of reproduction temperature effect for ecoregions
  real<lower=0> gamma_temp_sigma_area; # SD of reproduction temperature effect for reference areas
  real gamma_temp_bcr[nbcr]; # mean of reproduction temperature effect for ecoregions
  real gamma_temp_area[narea]; # reproduction temperature effect for reference areas

  real gamma_prec_mu; # grand mean of reproduction precipitation effect 
  real<lower=0> gamma_prec_sigma_bcr; # SD of reproduction precipitation effect for ecoregions
  real<lower=0> gamma_prec_sigma_area; # SD of reproduction precipitation effect for reference areas
  real gamma_prec_bcr[nbcr]; # mean of reproduction precipitation effect for ecoregions
  real gamma_prec_area[narea]; # reproduction precipitation effect for reference areas

  real<lower=0> gamma_sigma; # SD of process error for reproduction rate
  real log_gamma[narea, nyear-1]; # reproduction rate
  real<lower=0> R[narea, nyear-1]; # reproduction

  # survival
  real omega_alpha_mu; # grand mean of survival intercept
  real<lower=0> omega_alpha_sigma_bcr; # SD of survival intercept for ecoregions
  real<lower=0> omega_alpha_sigma_area; # SD of survival intercept for reference areas
  real omega_alpha_bcr[nbcr]; # mean of survival intercept for ecoregions
  real omega_alpha_area[narea]; # survival intercept for reference areas

  real omega_bpop_mu; # grand mean of survival density dependence
  real<lower=0> omega_bpop_sigma_bcr; # SD of survival density dependence for ecoregions
  real<lower=0> omega_bpop_sigma_area; # SD of survival density dependence for reference areas
  real omega_bpop_bcr[nbcr]; # mean of survival density dependence for ecoregions
  real omega_bpop_area[narea]; # survival density dependence for reference areas

  real omega_temp_mu; # grand mean of survival temperature effect 
  real<lower=0> omega_temp_sigma_bcr; # SD of survival temperature effect for ecoregions
  real<lower=0> omega_temp_sigma_area; # SD of survival temperature effect for reference areas
  real omega_temp_bcr[nbcr]; # mean of survival temperature effect for ecoregions
  real omega_temp_area[narea]; # survival temperature effect for reference areas

  real omega_prec_mu; # grand mean of survival precipitation effect 
  real<lower=0> omega_prec_sigma_bcr; # SD of survival precipitation effect for ecoregions
  real<lower=0> omega_prec_sigma_area; # SD of survival precipitation effect for reference areas
  real omega_prec_bcr[nbcr]; # mean of survival precipitation effect for ecoregions
  real omega_prec_area[narea]; # survival precipitation effect for reference areas

  real<lower=0> omega_sigma; # SD of process error for survival rate
  vector[nyear-1] logit_omega[narea]; # survival rate
  real<lower=0> S[narea, nyear-1]; # survival

  # recovery
  real logit_rec_mu; # mean of logit recovery rate
  real<lower=0> logit_rec_sigma; # SD of logit recovery rate
  vector[nyear] logit_rec; # logit recovery rate
} # parameters

transformed parameters {
  # population density
  vector[narea] log_D0_mu; # expected log initial population density
  matrix<lower=0>[narea, nyear] D; # population density in each site and year

  # reproduction
  real log_gamma_mu[narea, nyear-1]; # expected log reproduction rate
  real<lower=0> gamma[narea, nyear-1]; # reproduction rate

  # survival
  real logit_omega_mu[narea, nyear-1]; # expected logit survival rate
  vector<lower=0, upper=1>[nyear-1] omega[narea]; # survival rate

  # recovery
  vector<lower=0, upper=1>[nyear] rec; # recovery rate
  simplex[nyear+1] pr[narea, nyear]; # cell probability of band-recovery array 

  # population density
  for (i in 1:narea) {
    log_D0_mu[i] <- 
      D0_alpha + 
      D0_temp * temp0[i] + 
      D0_prec * prec0[i];
    D[i,1] <- exp(log_D0[i]);
    for (t in 2:nyear) {
      D[i,t] <- R[i,t-1] + S[i,t-1];
    } # t
  } # i

  # reproduction
  for (i in 1:narea) {
    for (t in 1:(nyear-1)) {
      log_gamma_mu[i,t] <- 
        gamma_alpha_area[i] + 
        gamma_bpop_area[i] * (log(D[i,t]) - mean(log(row(D,i)))) / sd(log(row(D,i))) + 
        gamma_temp_area[i] * temp[i,t+1] + 
        gamma_prec_area[i] * prec[i,t+1];
      gamma[i,t] <- exp(log_gamma[i,t]);
    } # t
  } # i

  # survival
  for (i in 1:narea) {
    for (t in 1:(nyear-1)) {
      logit_omega_mu[i,t] <- 
        omega_alpha_area[i] + 
        omega_bpop_area[i] * (log(D[i,t]) - mean(log(row(D,i)))) / sd(log(row(D,i))) + 
        omega_temp_area[i] * temp[i,t+1] + 
        omega_prec_area[i] * prec[i,t+1];
      omega[i,t] <- inv_logit(logit_omega[i,t]);
    } # t
  } # i

  # recovery
  for (t in 1:nyear) {
    rec[t] <- inv_logit(logit_rec[t]);
  } # t

  for (i in 1:narea) {
    for (t in 1:nyear) {
      pr[i,t,t] <- rec[t];

      for (s in (t+1):nyear) {
        pr[i,t,s] <- getprod(omega[i],t,s-1) * rec[s];
      } # s

      for (s in 1:(t-1)) {
        pr[i,t,s] <- 0;
      } # s
    } # t

    for (t in 1:nyear) {
      pr[i,t,nyear+1] <- 1 - getsum(pr[i,t],t,nyear);
    } # t
  } # i

} # transformed parameters

model {
  ### priors
  # population density
  D0_alpha ~ normal(0, 100);
  D0_temp ~ normal(0, 100); 
  D0_prec ~ normal(0, 100); 
  D0_sigma ~ cauchy(0, 5);

  # reproduction
  gamma_alpha_mu ~ normal(0, 100);
  gamma_alpha_sigma_bcr ~ cauchy(0, 5); 
  gamma_alpha_sigma_area ~ cauchy(0, 5); 
  for (j in 1:nbcr) {
    gamma_alpha_bcr[j] ~ normal(gamma_alpha_mu, gamma_alpha_sigma_bcr);
  } # j
  for (i in 1:narea) {
    gamma_alpha_area[i] ~ normal(gamma_alpha_bcr[bcr[i]], gamma_alpha_sigma_area);
  } # i

  gamma_bpop_mu ~ normal(0, 100);
  gamma_bpop_sigma_bcr ~ cauchy(0, 5); 
  gamma_bpop_sigma_area ~ cauchy(0, 5); 
  for (j in 1:nbcr) {
    gamma_bpop_bcr[j] ~ normal(gamma_bpop_mu, gamma_bpop_sigma_bcr);
  } # j
  for (i in 1:narea) {
    gamma_bpop_area[i] ~ normal(gamma_bpop_bcr[bcr[i]], gamma_bpop_sigma_area);
  } # i

  gamma_temp_mu ~ normal(0, 100);
  gamma_temp_sigma_bcr ~ cauchy(0, 5); 
  gamma_temp_sigma_area ~ cauchy(0, 5); 
  for (j in 1:nbcr) {
    gamma_temp_bcr[j] ~ normal(gamma_temp_mu, gamma_temp_sigma_bcr);
  } # j
  for (i in 1:narea) {
    gamma_temp_area[i] ~ normal(gamma_temp_bcr[bcr[i]], gamma_temp_sigma_area);
  } # i

  gamma_prec_mu ~ normal(0, 100);
  gamma_prec_sigma_bcr ~ cauchy(0, 5); 
  gamma_prec_sigma_area ~ cauchy(0, 5); 
  for (j in 1:nbcr) {
    gamma_prec_bcr[j] ~ normal(gamma_prec_mu, gamma_prec_sigma_bcr);
  } # j
  for (i in 1:narea) {
    gamma_prec_area[i] ~ normal(gamma_prec_bcr[bcr[i]], gamma_prec_sigma_area);
  } # i

  gamma_sigma ~ cauchy(0, 5);

  # survival
  omega_alpha_mu ~ normal(0, 100);
  omega_alpha_sigma_bcr ~ cauchy(0, 5); 
  omega_alpha_sigma_area ~ cauchy(0, 5); 
  for (j in 1:nbcr) {
    omega_alpha_bcr[j] ~ normal(omega_alpha_mu, omega_alpha_sigma_bcr);
  } # j
  for (i in 1:narea) {
    omega_alpha_area[i] ~ normal(omega_alpha_bcr[bcr[i]], omega_alpha_sigma_area);
  } # i

  omega_bpop_mu ~ normal(0, 100);
  omega_bpop_sigma_bcr ~ cauchy(0, 5); 
  omega_bpop_sigma_area ~ cauchy(0, 5); 
  for (j in 1:nbcr) {
    omega_bpop_bcr[j] ~ normal(omega_bpop_mu, omega_bpop_sigma_bcr);
  } # j
  for (i in 1:narea) {
    omega_bpop_area[i] ~ normal(omega_bpop_bcr[bcr[i]], omega_bpop_sigma_area);
  } # i

  omega_temp_mu ~ normal(0, 100);
  omega_temp_sigma_bcr ~ cauchy(0, 5); 
  omega_temp_sigma_area ~ cauchy(0, 5); 
  for (j in 1:nbcr) {
    omega_temp_bcr[j] ~ normal(omega_temp_mu, omega_temp_sigma_bcr);
  } # j
  for (i in 1:narea) {
    omega_temp_area[i] ~ normal(omega_temp_bcr[bcr[i]], omega_temp_sigma_area);
  } # i

  omega_prec_mu ~ normal(0, 100);
  omega_prec_sigma_bcr ~ cauchy(0, 5); 
  omega_prec_sigma_area ~ cauchy(0, 5); 
  for (j in 1:nbcr) {
    omega_prec_bcr[j] ~ normal(omega_prec_mu, omega_prec_sigma_bcr);
  } # j
  for (i in 1:narea) {
    omega_prec_area[i] ~ normal(omega_prec_bcr[bcr[i]], omega_prec_sigma_area);
  } # i

  omega_sigma ~ cauchy(0, 5);

  # recovery
  logit_rec_mu ~ normal(0,100);
  logit_rec_sigma ~ cauchy(0,5);

  ### likelihood
  # population density
  for (i in 1:narea) {
    log_D0[i] ~ normal(log_D0_mu[i], D0_sigma);
    for (t in 1:(nyear-1)) {
      R[i,t] ~ normal(D[i,t]*gamma[i,t], sqrt(D[i,t]*gamma[i,t])); 
      S[i,t] ~ normal(D[i,t]*omega[i,t], sqrt(D[i,t]*omega[i,t]*(1-omega[i,t])));
    } # t
  } # i

  # reproduction
  for (i in 1:narea) {
    for (t in 1:(nyear-1)) {
      log_gamma[i,t] ~ normal(log_gamma_mu[i,t], gamma_sigma);
    }
  }

  # survival
  for (i in 1:narea) {
    for (t in 1:(nyear-1)) {
      logit_omega[i,t] ~ normal(logit_omega_mu[i,t], omega_sigma);
    }
  }

  # recovery
  for (t in 1:nyear) {
    logit_rec[t] ~ normal(logit_rec_mu, logit_rec_sigma);
  }

  ### Observation
  # population density
  for (i in 1:narea) {
    for (t in 1:nyear) {
      bpop_obs[i,t] ~ normal(D[i,t], bpop_obs_sigma[i,t]);
    } # t
  } # k

  # band-recovery
  for (i in 1:narea) {
    for (t in 1:nyear) {
      marr[i,t] ~ multinomial(pr[i,t]);
    } # t
  } # i

} # model

", fill=TRUE)
sink()

#================
# setup for STAN
#================
data <- list(narea=narea, nyear=nyear, 
             nbcr=nbcr, bcr=bcr, 
             temp=temp, prec=prec, temp0=temp0, prec0=prec0, 
             marr=marr, 
             bpop_obs=bpop_obs, bpop_obs_sigma=bpop_obs_sigma)

init_fun <- function() {
              list(
                D0_alpha=0, D0_temp=0, D0_prec=0, D0_sigma=1, 
                gamma_alpha_mu=0, gamma_alpha_sigma_bcr=1, gamma_alpha_sigma_area=1, 
                gamma_alpha_bcr=rep(0,nbcr), gamma_alpha_area=rep(0,narea), 
                gamma_bpop_mu=0, gamma_bpop_sigma_bcr=1, gamma_bpop_sigma_area=1, 
                gamma_bpop_bcr=rep(0,nbcr), gamma_bpop_area=rep(0,narea), 
                gamma_temp_mu=0, gamma_temp_sigma_bcr=1, gamma_temp_sigma_area=1, 
                gamma_temp_bcr=rep(0,nbcr), gamma_temp_area=rep(0,narea), 
                gamma_prec_mu=0, gamma_prec_sigma_bcr=1, gamma_prec_sigma_area=1, 
                gamma_prec_bcr=rep(0,nbcr), gamma_prec_area=rep(0,narea), 
                gamma_sigma=1, 
                log_gamma=matrix(0, narea, nyear-1), 
                omega_alpha_mu=0, omega_alpha_sigma_bcr=1, omega_alpha_sigma_area=1, 
                omega_alpha_bcr=rep(0,nbcr), omega_alpha_area=rep(0,narea), 
                omega_bpop_mu=0, omega_bpop_sigma_bcr=1, omega_bpop_sigma_area=1, 
                omega_bpop_bcr=rep(0,nbcr), omega_bpop_area=rep(0,narea), 
                omega_temp_mu=0, omega_temp_sigma_bcr=1, omega_temp_sigma_area=1, 
                omega_temp_bcr=rep(0,nbcr), omega_temp_area=rep(0,narea), 
                omega_prec_mu=0, omega_prec_sigma_bcr=1, omega_prec_sigma_area=1, 
                omega_prec_bcr=rep(0,nbcr), omega_prec_area=rep(0,narea), 
                omega_sigma=1, 
                logit_omega=matrix(0, narea, nyear-1), 
                logit_rec_mu=0, logit_rec_sigma=1, 
                logit_rec=rep(0, nyear)
              ) # list
            } # function

#===========
# call STAN
#===========
rstan_options(auto_write=TRUE)
options(mc.cores=parallel::detectCores())

fit <- stan(file=paste(c(path, 'ipm additive.stan'), collapse=''), 
            data=data, init=init_fun, 
            par=c('D0_alpha', 'D0_temp', 'D0_prec', 'D0_sigma', 
                  'gamma_alpha_mu', 'gamma_alpha_sigma_bcr', 'gamma_alpha_sigma_area', 
                  'gamma_alpha_bcr', 'gamma_alpha_area', 
                  'gamma_bpop_mu', 'gamma_bpop_sigma_bcr', 'gamma_bpop_sigma_area', 
                  'gamma_bpop_bcr', 'gamma_bpop_area', 
                  'gamma_temp_mu', 'gamma_temp_sigma_bcr', 'gamma_temp_sigma_area', 
                  'gamma_temp_bcr', 'gamma_temp_area', 
                  'gamma_prec_mu', 'gamma_prec_sigma_bcr', 'gamma_prec_sigma_area', 
                  'gamma_prec_bcr', 'gamma_prec_area', 
                  'gamma_sigma', 
                  'omega_alpha_mu', 'omega_alpha_sigma_bcr', 'omega_alpha_sigma_area', 
                  'omega_alpha_bcr', 'omega_alpha_area', 
                  'omega_bpop_mu', 'omega_bpop_sigma_bcr', 'omega_bpop_sigma_area', 
                  'omega_bpop_bcr', 'omega_bpop_area', 
                  'omega_temp_mu', 'omega_temp_sigma_bcr', 'omega_temp_sigma_area', 
                  'omega_temp_bcr', 'omega_temp_area', 
                  'omega_prec_mu', 'omega_prec_sigma_bcr', 'omega_prec_sigma_area', 
                  'omega_prec_bcr', 'omega_prec_area', 
                  'omega_sigma', 
                  'logit_rec_mu', 'logit_rec_sigma', 
                  'rec', 
                  'gamma', 'omega', 
                  'D', 'R', 'S'), 
#            iter=200, warmup=100, thin=1, chains=1, 
            iter=3000, warmup=1000, thin=1, chains=5, 
            control=list(adapt_delta=0.80))

#==============
# save results
#==============
save(fit, file=paste(c(path, 'ipm additive.RData'), collapse=''))

#===============
# check results
#===============
print(fit, digits=3, 
      c('D0_alpha', 'D0_temp', 'D0_prec', 'D0_sigma', 
        'gamma_alpha_mu', 'gamma_alpha_sigma_bcr', 'gamma_alpha_sigma_area', 
        'gamma_alpha_bcr', 'gamma_alpha_area', 
        'gamma_bpop_mu', 'gamma_bpop_sigma_bcr', 'gamma_bpop_sigma_area', 
        'gamma_bpop_bcr', 'gamma_bpop_area', 
        'gamma_temp_mu', 'gamma_temp_sigma_bcr', 'gamma_temp_sigma_area', 
        'gamma_temp_bcr', 'gamma_temp_area', 
        'gamma_prec_mu', 'gamma_prec_sigma_bcr', 'gamma_prec_sigma_area', 
        'gamma_prec_bcr', 'gamma_prec_area', 
        'gamma_sigma',  
        'omega_alpha_mu', 'omega_alpha_sigma_bcr', 'omega_alpha_sigma_area', 
        'omega_alpha_bcr', 'omega_alpha_area', 
        'omega_bpop_mu', 'omega_bpop_sigma_bcr', 'omega_bpop_sigma_area', 
        'omega_bpop_bcr', 'omega_bpop_area', 
        'omega_temp_mu', 'omega_temp_sigma_bcr', 'omega_temp_sigma_area', 
        'omega_temp_bcr', 'omega_temp_area', 
        'omega_prec_mu', 'omega_prec_sigma_bcr', 'omega_prec_sigma_area', 
        'omega_prec_bcr', 'omega_prec_area', 
        'omega_sigma', 
        'logit_rec_mu', 'logit_rec_sigma'))

traceplot(fit, inc_warmup=T, 
          c('D0_alpha', 'D0_sigma', 
            'gamma_alpha_mu', 'gamma_sigma', 
            'omega_alpha_mu', 'omega_sigma', 
            'gamma_prec_sigma_area', 
            'omega_temp_sigma_area'))



