####################################################################################################################################
### Fit a hidden Markov model to diving data from MALE Weddell seals.                                                            ###
### Used in the analysis presented in:                                                                                           ###
### "Sex-specific variation in the use of vertical habitat by a resident Antarctic top predator"                                 ###
### Theoni Photopoulou, Karine Heerah, Jennifer Pohle and Lars Boehme (2020)                                                     ###
####################################################################################################################################

### load required packages and scripts
library(boot)
library(Rcpp)
library(RcppArmadillo)
library(dplyr)
library(purrr)
library(lubridate)

# please insert the correct path of the forward_algorithm.cpp file 
# make sure that a C++ compiler (e.g. Rtools) is installed on your computer
sourceCpp("mllk_cov.cpp")

# load the functions to run the HMM
source("hmm7Rcpp.R")

# load the data
load("weddell_dive_data_2011.RData")
pdat <- weddell_dive_data_2011
head(pdat)

# there are about 500 dives without hunting (about 2.4% of dives)
table(pdat %>% filter(source=="dv") %>% select(hunt_avgdep) %>% is.na(.))

# separate out the data by sex
femaleREF <- pdat %>% filter(sex=="F") %>% distinct(ID); femaleREF
maleREF <- pdat %>% filter(sex=="M") %>% distinct(ID); maleREF
pdatF <- pdat[pdat$ID %in% femaleREF$ID,] %>% droplevels()
pdatM <- pdat[pdat$ID %in% maleREF$ID,] %>% droplevels()

# create list of observations (obslist) for model: one list entry for each individual
obslist <- split(x=pdatF, f=pdatF$ID)  # F only
map(obslist, .f=nrow)
class(obslist)
length(obslist)

# Functions that turn the parameters of the beta distribution (alpha and beta) 
# into mean and variance and vice versa.
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}
estBetaMuVar <- function(alpha, beta) {
  mu <- alpha / (alpha + beta)
  var <- (alpha * beta) / ( (alpha + beta)^2 * (alpha + beta + 1))
  return(params = list(mu = mu, var = var))
}




# I am labelling the states, but I cannot assume these correspond to functional behaviours. 
stateNames <- c("haulout", "surface", "shallow", "epipelagic", "pelagic", "benthic")
N <- length(stateNames)

### DEFINE START VALUES

# dive duration (seconds) - all N states have duration: N=6 values
dd.mu0 <- c(10148, 2277, 61, 447, 843, 1138)
dd.sigma0 <- c(12424, 2987, 36, 315, 257, 286)

# hunting depth (metres) - N-2 states have hunting depth: N-2=4 values
hd.mu0 <- c(9, 39, 194, 413)
hd.sigma0 <- c(3, 25, 88, 63)

# proportion of dive duration spent hunting (range 0-1) - only states 3 to N have hunting time: N-2=4 values
ph.mu0 <- c(0.51, 0.55, 0.49, 0.43)
ph.sigma0 <- c(0.03, 0.05, 0.03, 0.02)

# proportion of bathymetry reached (hunt_avgdep/bathymetry range 0-1) - only states 3 to N have prop_hunt_bathy: N-2=4 values
pb.pi0 <- c(0.01, 0.01, 0.01, 0.89)

# salinity at hunting depth - only states 3 to N have salinity: N-2=4 values
sal.mu0 <- c(34.18, 34.29, 34.44, 34.61)
sal.sigma0 <- c(0.18, 0.12, 0.12, 0.06)

# transition probability matrix
gamma0 <- diag(N)
gamma0[!gamma0] <- rep(0.5/(N-1),N^2-N)
diag(gamma0) <- rep(0.5,N) 

# the initial specification of the state transition probability matrix
delta0 <- rep(1/N,N) 

# covariate names: cosine of hour, sine of hour and standardised week of the year
cov.names <- c("coshr", "sinhr", "week_stnd")

# beta: coefficients of the covariate effects on the transition probability matrix
beta0 <- cbind(c(rep(-2,N*(N-1)), 
                 rep(rep(0,N*(N-1)),length(cov.names)+2))) # 3 covariates and 
                                                           # 1 interaction but with 2 terms (coshr*wk, sinhr*wk)

## run numerical maximum-likelihood estimation and return estimates
s <- Sys.time()
mod <- mle(obslist,
           dd.mu0,dd.sigma0,
           hd.mu0,hd.sigma0,
           ph.mu0,ph.sigma0,
           pb.pi0,
           sal.mu0,sal.sigma0,
           beta0,delta0,
           N,cov.names, return_hessian = TRUE)
mod
Sys.time()-s

# AIC = -2*mllk + 2*p
npar <- length(unlist(mod[c(1:6,9:13)]))
-2*(-mod$mllk) + 2*npar # ph.alpha and ph.beta aren't extra parameters that are estimated,
                        # they are derived from ph.mu and ph.sigma, so we still have 9 parameters

# maximum likelihood estimates
mod[1:11]

ct70_6st_hmm7Rcpp_F_hourwk <- mod

# source the Viterbi algoirthm for working out the most likely state sequence
source("viterbi_hmm7Rcpp_2int.R")
viterbi_list_ct70 <- viterbi(obslist,cov.names,mod,N)
length(viterbi_list_ct70)
names(viterbi_list_ct70) <- unique(pdatF$ID)

viterbi_vec_ct70 <- unlist(viterbi_list_ct70)
length(viterbi_vec_ct70)
dim(pdat)
pdatF$viterbi_state <- viterbi_vec_ct70

# look at the distribution of Viterbi-decoded states
table(pdatF$viterbi_state)
table(pdatF$ID, pdatF$viterbi_state)
# look at the individual variability in the state occupancy
indivF_tab <- round(table(pdatF$ID, pdatF$viterbi_state)/rowSums(table(pdatF$ID, pdatF$viterbi_state)), digits=3); indivF_tab
range(indivF_tab[,1])

# find average state occupancy
delta.avg <- as.vector(round(table(pdatF$viterbi_state)/length(pdatF$viterbi_state), digits=2)); delta.avg




# PLOT state-dependent densities

require(ggplot2)
require(RColorBrewer)

# define colours using a palette
nbStates <- length(stateNames)
brew.cols <- brewer.pal(nbStates+1, "Set2")




# DURATION state-dependent densities
x <- seq(min(pdatF$DURATION), max(pdatF$DURATION), length=1000)
mod[1:2]
mu1_dd <- mod$dd.mu[1]
sd1_dd <- mod$dd.sigma[1]
mu2_dd <- mod$dd.mu[2]
sd2_dd <- mod$dd.sigma[2]
mu3_dd <- mod$dd.mu[3]
sd3_dd <- mod$dd.sigma[3]
mu4_dd <- mod$dd.mu[4]
sd4_dd <- mod$dd.sigma[4]
mu5_dd <- mod$dd.mu[5]
sd5_dd <- mod$dd.sigma[5]
mu6_dd <- mod$dd.mu[6]
sd6_dd <- mod$dd.sigma[6]

d1_dd <- dgamma(x, shape = mu1_dd^2/sd1_dd^2, scale = sd1_dd^2/mu1_dd)*delta.avg[1]
d2_dd <- dgamma(x, shape = mu2_dd^2/sd2_dd^2, scale = sd2_dd^2/mu2_dd)*delta.avg[2]
d3_dd <- dgamma(x, shape = mu3_dd^2/sd3_dd^2, scale = sd3_dd^2/mu3_dd)*delta.avg[3]
d4_dd <- dgamma(x, shape = mu4_dd^2/sd4_dd^2, scale = sd4_dd^2/mu4_dd)*delta.avg[4]
d5_dd <- dgamma(x, shape = mu5_dd^2/sd5_dd^2, scale = sd5_dd^2/mu5_dd)*delta.avg[5]
d6_dd <- dgamma(x, shape = mu6_dd^2/sd6_dd^2, scale = sd6_dd^2/mu6_dd)*delta.avg[6]

dmarg_dd <- d1_dd + d2_dd + d3_dd + d4_dd + d5_dd + d6_dd

dd_dat <- pdatF$DURATION

quartz()
ggplot(data=data.frame(x=dd_dat), aes(x,..density..)) + 
  geom_histogram(boundary=0, binwidth=150, fill="grey90") + 
  ylim(0,0.0025) + theme_bw() + 
  geom_line(data=data.frame(x=x, d1_md=d1_dd), aes(x, d1_dd, colour="Haulout"), size=1.3) +
  geom_line(data=data.frame(x=x, d2_md=d2_dd), aes(x, d2_dd, colour="Surface"), size=1.3) +
  geom_line(data=data.frame(x=x, d3_md=d3_dd), aes(x, d3_dd, colour="Shallow dive"), size=1.3) +
  geom_line(data=data.frame(x=x, d4_md=d4_dd), aes(x, d4_dd, colour="Epipelagic dive"), size=1.3) +
  geom_line(data=data.frame(x=x, d5_md=d5_dd), aes(x, d5_dd, colour="Pelagic dive"), size=1.3) +
  geom_line(data=data.frame(x=x, d6_md=d6_dd), aes(x, d6_dd, colour="Benthic dive"), size=1.3) +
  geom_line(data=data.frame(x=x, dmarg_dd=dmarg_dd), aes(x, dmarg_dd, color="Marginal"), 
            size=.8, linetype="dashed") +
  
  scale_colour_manual(name="Densities", 
                      values = c("Haulout" = brew.cols[6], "Surface" = brew.cols[3], 
                                 "Shallow dive" = brew.cols[5], "Epipelagic dive" = brew.cols[1], 
                                 "Pelagic dive" = brew.cols[7], "Benthic dive" = brew.cols[4],
                                 "Marginal" = brew.cols[2]),
                      breaks=c("Haulout", "Surface", "Shallow dive", "Epipelagic dive", 
                               "Pelagic dive", "Benthic dive", "Marginal")) + 
  scale_linetype_manual(name="Densities", values=c("Haulout" = 3, "Surface" = 3, "Shallow dive"=1, 
                                                   "Epipelagic dive"=1, "Pelagic dive"=1, "Benthic dive"=1, 
                                                   "Marginal" = 2)) + 
  scale_x_continuous(breaks=c(0,300,600,900,1200,1800,3600), 
                     labels=c(0,5,10,15,20,30,60),
                     limits=c(0,3600)) + 
  xlab("Behaviour duration (min)") + ylab("Density") + ggtitle("State-Dependent Duration Densities") +
  theme(legend.position=c(.8, .7)) +
  theme(text=element_text(size=15))




# HUNTING DEPTH state-dependent densities
x <- seq(min(pdatF$hunt_avgdep), max(pdatF$hunt_avgdep), length=1000)
mod[3:4]
mu3_hd <- mod$hd.mu[1]
sd3_hd <- mod$hd.sigma[1]
mu4_hd <- mod$hd.mu[2]
sd4_hd <- mod$hd.sigma[2]
mu5_hd <- mod$hd.mu[3]
sd5_hd <- mod$hd.sigma[3]
mu6_hd <- mod$hd.mu[4]
sd6_hd <- mod$hd.sigma[4]

d3_hd <- dgamma(x, shape = mu3_hd^2/sd3_hd^2, scale = sd3_hd^2/mu3_hd)*(delta.avg[3]/sum(delta.avg[3:N]))
d4_hd <- dgamma(x, shape = mu4_hd^2/sd4_hd^2, scale = sd4_hd^2/mu4_hd)*(delta.avg[4]/sum(delta.avg[3:N]))
d5_hd <- dgamma(x, shape = mu5_hd^2/sd5_hd^2, scale = sd5_hd^2/mu5_hd)*(delta.avg[5]/sum(delta.avg[3:N]))
d6_hd <- dgamma(x, shape = mu6_hd^2/sd6_hd^2, scale = sd6_hd^2/mu6_hd)*(delta.avg[6]/sum(delta.avg[3:N]))

dmarg_hd <- d3_hd + d4_hd + d5_hd + d6_hd

hd_dat <- pdatF$hunt_avgdep

quartz()
ggplot(data=data.frame(x=hd_dat), aes(x,..density..)) + 
  geom_histogram(boundary=0, binwidth=10, fill="grey90") + 
  ylim(0,0.025) + theme_bw() + 
  geom_line(data=data.frame(x=x, d3_md=d3_hd), aes(x, d3_hd, colour="Shallow dive"), size=1.3) +
  geom_line(data=data.frame(x=x, d4_md=d4_hd), aes(x, d4_hd, colour="Epipelagic dive"), size=1.3) +
  geom_line(data=data.frame(x=x, d5_md=d5_hd), aes(x, d5_hd, colour="Pelagic dive"), size=1.3) +
  geom_line(data=data.frame(x=x, d6_md=d6_hd), aes(x, d6_hd, colour="Benthic dive"), size=1.3) +
  geom_line(data=data.frame(x=x, dmarg_hd=dmarg_hd), aes(x, dmarg_hd, color="Marginal"), 
            size=.8, linetype="dashed") +
  
  scale_colour_manual(name="Densities", 
                      values = c("Shallow dive" = brew.cols[5], "Epipelagic dive" = brew.cols[1], 
                                 "Pelagic dive" = brew.cols[7], "Benthic dive" = brew.cols[4],
                                 "Marginal" = brew.cols[2]),
                      breaks=c("Shallow dive", "Epipelagic dive", 
                               "Pelagic dive", "Benthic dive", "Marginal")) + 
  scale_linetype_manual(name="Densities", values=c("Shallow dive"=1, 
                                                   "Epipelagic dive"=1, "Pelagic dive"=1, "Benthic dive"=1, 
                                                   "Marginal" = 2)) + 
  scale_x_continuous(breaks=c(0,50,100,200,400,600,800), 
                     limits=c(0,800)) + 
  xlab("Hunting depth (m)") + ylab("Density") + ggtitle("State-Dependent Duration Densities") +
  theme(legend.position=c(.8, .7)) +
  theme(text=element_text(size=15))




# PROPORTION OF DIVE SPENT HUNTING state-dependent densities
x <- seq(min(pdatF$hunt_prop[pdatF$source=="dv"]),  max(pdatF$hunt_prop), length=1000) # avoid zero
mod[5:8]

a3_ph <- mod$ph.alpha[1]
b3_ph <- mod$ph.beta[1]
a4_ph <- mod$ph.alpha[2]
b4_ph <- mod$ph.beta[2]
a5_ph <- mod$ph.alpha[3]
b5_ph <- mod$ph.beta[3]
a6_ph <- mod$ph.alpha[4]
b6_ph <- mod$ph.beta[4]

d3_ph <- dbeta(x, shape1 = a3_ph, shape2 = b3_ph, ncp=0)*(delta.avg[3]/sum(delta.avg[3:N]))
d4_ph <- dbeta(x, shape1 = a4_ph, shape2 = b4_ph, ncp=0)*(delta.avg[4]/sum(delta.avg[3:N]))
d5_ph <- dbeta(x, shape1 = a5_ph, shape2 = b5_ph, ncp=0)*(delta.avg[5]/sum(delta.avg[3:N]))
d6_ph <- dbeta(x, shape1 = a6_ph, shape2 = b6_ph, ncp=0)*(delta.avg[6]/sum(delta.avg[3:N]))

dmarg_ph <- d3_ph + d4_ph + d5_ph + d6_ph

ph_dat <- pdatF %>% filter(source=="dv") %>% dplyr::select(hunt_prop) 
ph_dat <- ph_dat$hunt_prop

quartz()
ggplot(data=data.frame(x=ph_dat), aes(x,..density..)) + 
  geom_histogram(boundary=0, binwidth=0.07, fill="grey90") + 
  ylim(0,2.2) + theme_bw() + 
  geom_line(data=data.frame(x=x, d3_ph=d3_ph), aes(x, d3_ph, colour="Shallow dive"), size=1.3) +
  geom_line(data=data.frame(x=x, d4_ph=d4_ph), aes(x, d4_ph, colour="Epipelagic dive"), size=1.3) +
  geom_line(data=data.frame(x=x, d5_ph=d5_ph), aes(x, d5_ph, colour="Pelagic dive"), size=1.3) +
  geom_line(data=data.frame(x=x, d6_ph=d6_ph), aes(x, d6_ph, colour="Benthic dive"), size=1.3) +
  geom_line(data=data.frame(x=x, dmarg_ph=dmarg_ph), aes(x, dmarg_ph, color="Marginal"), size=.8, linetype="dashed") +
  
scale_colour_manual(name="Densities", values = c("Shallow dive" = brew.cols[5], "Epipelagic dive" = brew.cols[1],
                                                 "Pelagic dive" = brew.cols[7], "Benthic dive" = brew.cols[4], 
                                                 "Marginal" = brew.cols[2]),
                    breaks=c("Shallow dive", "Epipelagic dive", "Pelagic dive", "Benthic dive", "Marginal")) +
  scale_linetype_manual(name="Densities", values=c("Shallow dive"=1, "Epipelagic dive"=1, "Pelagic dive"=1, 
                                                   "Benthic dive"=1, "Marginal" = 2)) +
  
  xlab("Proportion of time spent hunting") + ylab("Density") + ggtitle("State-Dependent Densities") +
  theme(legend.position=c(.8, .7)) +
  theme(text=element_text(size=15))




# PROPORTION OF BATHYMERTY state-dependent densities
round(mod[7]$pb.pi, digits=2)




# SALINITY AT DEPTH state-dependent densities
x <- seq(min(pdatF$psal), max(pdatF$psal), length=1000)
mod[8:9]
mu3_ps <- mod$sal.mu[1]
sd3_ps <- mod$sal.sigma[1]
mu4_ps <- mod$sal.mu[2]
sd4_ps <- mod$sal.sigma[2]
mu5_ps <- mod$sal.mu[3]
sd5_ps <- mod$sal.sigma[3]
mu6_ps <- mod$sal.mu[4]
sd6_ps <- mod$sal.sigma[4]

d3_ps <- dnorm(x, mean = mu3_ps, sd = sd3_ps)*(delta.avg[3]/sum(delta.avg[3:N]))
d4_ps <- dnorm(x, mean = mu4_ps, sd = sd4_ps)*(delta.avg[4]/sum(delta.avg[3:N]))
d5_ps <- dnorm(x, mean = mu5_ps, sd = sd5_ps)*(delta.avg[5]/sum(delta.avg[3:N]))
d6_ps <- dnorm(x, mean = mu6_ps, sd = sd6_ps)*(delta.avg[6]/sum(delta.avg[3:N]))

dmarg_ps <- d3_ps + d4_ps + d5_ps + d6_ps

ps_dat <- pdatF$psal

quartz()
ggplot(data=data.frame(x=ps_dat), aes(x,..density..)) + 
  geom_histogram(boundary=0, binwidth=0.08, fill="grey90") + 
  ylim(0,3) + theme_bw() + 
  geom_line(data=data.frame(x=x, d3_md=d3_ps), aes(x, d3_ps, colour="Shallow dive"), size=1.3) +
  geom_line(data=data.frame(x=x, d4_md=d4_ps), aes(x, d4_ps, colour="Epipelagic dive"), size=1.3) +
  geom_line(data=data.frame(x=x, d5_md=d5_ps), aes(x, d5_ps, colour="Pelagic dive"), size=1.3) +
  geom_line(data=data.frame(x=x, d6_md=d6_ps), aes(x, d6_ps, colour="Benthic dive"), size=1.3) +
  geom_line(data=data.frame(x=x, dmarg_ps=dmarg_ps), aes(x, dmarg_ps, color="Marginal"), 
            size=.8, linetype="dashed") +
  
  scale_colour_manual(name="Densities", 
                      values = c("Shallow dive" = brew.cols[5], "Epipelagic dive" = brew.cols[1], 
                                 "Pelagic dive" = brew.cols[7], "Benthic dive" = brew.cols[4],
                                 "Marginal" = brew.cols[2]),
                      breaks=c("Shallow dive", "Epipelagic dive", 
                               "Pelagic dive", "Benthic dive", "Marginal")) + 
  scale_linetype_manual(name="Densities", values=c("Shallow dive"=1, 
                                                   "Epipelagic dive"=1, "Pelagic dive"=1, "Benthic dive"=1, 
                                                   "Marginal" = 2)) + 
  scale_x_continuous(breaks=c(33.6,33.8,34.0,34.2,34.4,34.6,34.8), 
                     limits=c(min(ps_dat),34.9)) + 
  xlab("Salinity") + ylab("Density") + ggtitle("State-Dependent Duration Densities") +
  theme(legend.position=c(.83, .7)) +
  theme(text=element_text(size=15))


