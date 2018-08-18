rm(list=ls())
library(boot)
library(rstan)

path <- 'd:/RESEARCH/A. Waterfowl/4.ipm/'

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

bpop_obs <- as.matrix(bpop[,-1])
bpop_obs_sigma <- as.matrix(bpop_se[,-1])

# climate
temp <- read.csv(paste(c(path, 'data/data/temp.csv'), collapse=''))

temp <- as.matrix(temp[,-1])
temp_mean <- apply(temp, 1, mean)
temp_mean_mat <- matrix(rep(temp_mean, each=nyear), nrow=narea, ncol=nyear, byrow=T)
temp_var <- temp - temp_mean_mat
temp <- temp_var / sd(temp_var)

# band-recovery
load(paste(c(path, 'data/data/marr.RData'), collapse=''))
marr <- marr$af0 + marr$am0

#==============
# define model
#==============
sink(paste(c(path, 'ipm.stan'), collapse=''))
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
  matrix[narea, nyear] temp; # spring temperature anomaly

  # band-recovery
  int marr[narea, nyear+1, nyear+2]; # band-recovery data, standard banding
} # data

parameters {
  # initial population density
  real D0_alpha; # intercept for initial population density
  real D0_beta_temp; # temperature effect on initial population density
  real<lower=0> D0_sigma; # SD of process error in initial population density
  vector[narea] log_D0; # log density of the first year

  # reproduction
  real gamma_alpha_mu[nbcr]; # mean of reproduction intercept for each ecoregion
  real<lower=0> gamma_alpha_sigma[nbcr]; # standard deviation of reproduction intercept for each ecoregion
  real gamma_alpha[narea]; # reproduction intercept for each reference area
  real gamma_beta_ddp; # density dependence on reproduction
  real gamma_beta_temp; # temperature effect on reproduction
  real gamma_beta_ddp_temp; # interaction of density and temperature on reproduction
  real<lower=0> gamma_sigma; # SD of process error for reproduction rate
  real log_gamma[narea, nyear-1]; # reproduction rate
  real<lower=0> R[narea, nyear-1]; # reproduction

  # survival
  real omega_alpha_mu[nbcr]; # mean of survival intercept for each ecoregion
  real<lower=0> omega_alpha_sigma[nbcr]; # standard deviation of survival intercept for each ecoregion
  real omega_alpha[narea]; # survival intercept for each reference area
  real omega_beta_ddp; # density dependence on survival
  real omega_beta_temp; # temperature effect on survival
  real omega_beta_ddp_temp; # interaction of density and temperature on survival
  real<lower=0> omega_sigma; # SD of process error for survival rate
  vector[nyear] logit_omega[narea]; # survival rate
  real<lower=0> S[narea, nyear-1]; # survival

  # recovery
  real logit_rec_mu; # mean of logit recovery rate
  real<lower=0> logit_rec_sigma; # SD of logit recovery rate
  vector[nyear+1] logit_rec; # logit recovery rate
} # parameters

transformed parameters {
  # population density
  vector[narea] log_D0_mu; # expected log initial population density
  matrix<lower=0>[narea, nyear] D; # population density in each site and year

  # reproduction
  real log_gamma_mu[narea, nyear-1]; # expected log reproduction rate
  real<lower=0> gamma[narea, nyear-1]; # reproduction rate

  # survival
  real logit_omega_mu[narea, nyear]; # expected logit survival rate
  vector<lower=0, upper=1>[nyear] omega[narea]; # survival rate

  # recovery
  vector<lower=0, upper=1>[nyear+1] rec; # recovery rate
  simplex[nyear+2] pr[narea, nyear+1]; # cell probability of band-recovery array 

  # population density
  for (i in 1:narea) {
    log_D0_mu[i] <- 
      D0_alpha + 
      D0_beta_temp * temp[i,1];
    D[i,1] <- exp(log_D0[i]);
    for (t in 2:nyear) {
      D[i,t] <- R[i,t-1] + S[i,t-1];
    } # t
  } # i

  # reproduction
  for (i in 1:narea) {
    for (t in 1:(nyear-1)) {
      log_gamma_mu[i,t] <- 
        gamma_alpha[i] + 
        gamma_beta_ddp * (log(D[i,t]) - mean(log(D[]))) / sd(log(D[])) + 
        gamma_beta_temp * temp[i,t] + 
        gamma_beta_ddp_temp * (log(D[i,t]) - mean(log(D[]))) / sd(log(D[])) * temp[i,t];
      gamma[i,t] <- exp(log_gamma[i,t]);
    } # t
  } # i

  # survival
  for (i in 1:narea) {
    for (t in 1:nyear) {
      logit_omega_mu[i,t] <- 
        omega_alpha[i] + 
        omega_beta_ddp * (log(D[i,t]) - mean(log(D[]))) / sd(log(D[])) + 
        omega_beta_temp * temp[i,t] + 
        omega_beta_ddp_temp * (log(D[i,t]) - mean(log(D[]))) / sd(log(D[])) * temp[i,t];
      omega[i,t] <- inv_logit(logit_omega[i,t]);
    } # t
  } # i

  # recovery
  for (t in 1:(nyear+1)) {
    rec[t] <- inv_logit(logit_rec[t]);
  } # t

  for (k in 1:narea) {
    # adult
    for (t in 1:(nyear+1)) {
      pr[k,t,t] <- rec[t];

      for (j in (t+1):(nyear+1)) {
        pr[k,t,j] <- getprod(omega[k],t,j-1) * rec[j];
      } # j

      for (j in 1:(t-1)) {
        pr[k,t,j] <- 0;
      } # j
    } # t

    for (t in 1:(nyear+1)) {
      pr[k,t,nyear+2] <- 1 - getsum(pr[k,t],t,nyear+1);
    } # t
  } # k

} # transformed parameters

model {
  ### priors
  # population density
  D0_alpha ~ normal(0, 100);
  D0_beta_temp ~ normal(0, 100); 
  D0_sigma ~ cauchy(0, 5);

  # reproduction
  for (j in 1:nbcr) {
    gamma_alpha_mu[j] ~ normal(0, 100);
    gamma_alpha_sigma[j] ~ cauchy(0, 5);
  } # j
  for (i in 1:narea) {
    gamma_alpha[i] ~ normal(gamma_alpha_mu[bcr[i]], gamma_alpha_sigma[bcr[i]]);
  } # i
  gamma_beta_ddp ~ normal(0, 100); 
  gamma_beta_temp ~ normal(0, 100); 
  gamma_beta_ddp_temp ~ normal(0, 100); 
  gamma_sigma ~ cauchy(0, 5);

  # survival
  for (j in 1:nbcr) {
    omega_alpha_mu[j] ~ normal(0, 100);
    omega_alpha_sigma[j] ~ cauchy(0, 5);
  } # j
  for (i in 1:narea) {
    omega_alpha[i] ~ normal(omega_alpha_mu[bcr[i]], omega_alpha_sigma[bcr[i]]);
  } # i
  omega_beta_ddp ~ normal(0, 100); 
  omega_beta_temp ~ normal(0, 100); 
  omega_beta_ddp_temp ~ normal(0, 100); 
  omega_sigma ~ cauchy(0, 5);

  # recovery
  logit_rec_mu ~ normal(0,100);
  logit_rec_sigma ~ cauchy(0,5);

  ### likelihood
  # population density
  for (i in 1:narea) {
    log_D0[i] ~ normal(log_D0_mu[i], D0_sigma);
    for (t in 2:nyear) {
      R[i,t-1] ~ normal(D[i,t-1]*gamma[i,t-1], sqrt(D[i,t-1]*gamma[i,t-1])); 
      S[i,t-1] ~ normal(D[i,t-1]*omega[i,t-1], sqrt(D[i,t-1]*omega[i,t-1]*(1-omega[i,t-1])));
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
    for (t in 1:nyear) {
      logit_omega[i,t] ~ normal(logit_omega_mu[i,t], omega_sigma);
    }
  }

  # recovery
  for (t in 1:(nyear+1)) {
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
  for (k in 1:narea) {
    for (t in 1:(nyear+1)) {
      marr[k,t] ~ multinomial(pr[k,t]);
    } # t
  } # k

} # model

", fill=TRUE)
sink()

#================
# setup for STAN
#================
data <- list(narea=narea, nyear=nyear, temp=temp, marr=marr, 
             bpop_obs=bpop_obs, bpop_obs_sigma=bpop_obs_sigma)

init_fun <- function() {
              list(
                D0_alpha=0, D0_beta_temp=0, D0_sigma=1, 
                gamma_alpha_mu=rep(0,nbcr), gamma_alpha_sigma=rep(1,nbcr), 
                gamma_alpha=rep(0,narea), 
                gamma_beta_ddp=0, gamma_beta_temp=0, gamma_beta_ddp_temp=0, 
                gamma_sigma=1, 
                log_gamma=matrix(0, narea, nyear-1), 
                omega_alpha_mu=rep(0,nbcr), omega_alpha_sigma=rep(1,nbcr), 
                omega_alpha=rep(0,narea), 
                omega_beta_ddp=0, omega_beta_temp=0, omega_beta_ddp_temp=0, 
                omega_sigma=1, 
                logit_omega=matrix(0, narea, nyear), 
                logit_rec_mu=0, logit_rec_sigma=1, 
                logit_rec=rep(0, nyear+1)
              ) # list
            } # function

#===========
# call STAN
#===========
rstan_options(auto_write=TRUE)
options(mc.cores=parallel::detectCores())

fit <- stan(file=paste(c(path, 'ipm.stan'), collapse=''), 
            data=data, init=init_fun, 
            par=c('D0_alpha', 'D0_beta_temp', 'D0_sigma', 
                  'gamma_alpha_mu', 'gamma_alpha_sigma', 'gamma_alpha', 
                  'gamma_beta_ddp', 'gamma_beta_temp', 'gamma_beta_ddp_temp',  
                  'gamma_sigma', 
                  'omega_alpha_mu', 'omega_alpha_sigma', 'omega_alpha', 
                  'omega_beta_ddp', 'omega_beta_temp', 'omega_beta_ddp_temp',  
                  'omega_sigma', 
                  'logit_rec_mu', 'logit_rec_sigma', 'rec', 
                  'D', 'R', 'S'), 
            iter=200, warmup=100, thin=1, chains=1)

#==============
# save results
#==============
save(fit, file=paste(c(path, 'ipm.RData'), collapse=''))

#===============
# check results
#===============
print(fit, digits=3, 
      c('D0_alpha', 'D0_beta_temp', 'D0_sigma', 
        'gamma_alpha_mu', 'gamma_alpha_sigma', 'gamma_alpha', 
        'gamma_beta_ddp', 'gamma_beta_temp', 'gamma_beta_ddp_temp',  
        'gamma_sigma',  
        'omega_alpha_mu', 'omega_alpha_sigma', 'omega_alpha', 
        'omega_beta_ddp', 'omega_beta_temp', 'omega_beta_ddp_temp',  
        'omega_sigma', 
        'logit_rec_mu', 'logit_rec_sigma'))

traceplot(fit, inc_warmup=T, 
          c('D0_alpha', 'D0_sigma', 
            'gamma_alpha', 'gamma_sigma', 
            'omega_alpha', 'omega_sigma'))



