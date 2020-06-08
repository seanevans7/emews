####################################################################################################################################
### Functions used to fit a hidden Markov model to diving data from Weddell seals.                                               ###
### Used in the analysis presented in:                                                                                           ###
### "Sex-specific variation in the use of vertical habitat by a resident Antarctic top predator"                                 ###
### Theoni Photopoulou, Karine Heerah, Jennifer Pohle and Lars Boehme (2020)                                                     ###
####################################################################################################################################

## pn2pw
## function that converts 'natural' parameters (possibly constrained) to 'working' parameters 
## (all of which are real-valued) - this is only necessary since I use the unconstrained optimizer nlm() below 
pn2pw <- function(dd.mu,dd.sigma,
                  hd.mu,hd.sigma,
                  ph.mu,ph.sigma,
                  pb.pi,
                  sal.mu, sal.sigma,
                  beta,delta,N){
  
  tdd.mu <- log(dd.mu); # mean for dive duration (gamma distibuted)
  tdd.sigma <- log(dd.sigma);
  
  thd.mu <- log(hd.mu); # mean for hunting depth (gamma distibuted)
  thd.sigma <- log(hd.sigma);
  
  tph.mu <- log(ph.mu); # mean for time spent hunting (gamma distibuted)
  tph.sigma <- log(ph.sigma);      
  
  tpb.pi <- logit(pb.pi); # proportion of bathymetry reached (binomially distributed)
  
  tsal.mu <- sal.mu; # mean for salinity (normally distibuted)
  tsal.sigma <- sal.sigma;      
  
  tbeta <- as.vector(beta)
  
  teta <- log(delta[-1]/delta[1])
  
  parvect <- c(tdd.mu,tdd.sigma,                           # dive duration
               thd.mu, thd.sigma,                          # hunting depth
               tph.mu,tph.sigma,                           # proportion of dive time spent hunting
               tpb.pi,                                     # proportion of bathymetry
               tsal.mu, tsal.sigma,                        # salinity
               tbeta, teta)
  return(parvect)
}




## pw2pn
## function that performs the inverse transformation of the parameters
pw2pn <- function(parvect,N){ 
  
  dd.mu <- dd.sigma <- NULL;
  hd.mu <- hd.sigma <- NULL;
  ph.mu <- ph.sigma <- NULL;
  sal.mu <- sal.sigma <- NULL;
  pb.pi <- NULL;
  
  dd.mu <- exp(parvect[1:N]) # all N parameters
  dd.sigma <- exp(parvect[(N+1):(2*N)]) # all N parameters
  
  hd.mu <- exp(parvect[(2*N+1):(3*N-2)]) # N-2 parameters: (2*N+1):(2*N+(N-2))
  hd.sigma <- exp(parvect[(3*N-1):(4*N-4)]) # N-2 parameters: (2*N+((N-2+1))):(2*N+2*(N-2))
  
  ph.mu <- exp(parvect[(4*N-3):(5*N-6)]) # N-2 parameters: (2*N+(2*(N-2))+1):(2*N+(3*(N-2)))
  ph.sigma <- exp(parvect[(5*N-5):(6*N-8)]) # N-2 parameters: (2*N+(3*(N-2))+1):(2*N+(4*(N-2)))

#  ph.alpha <- ((1-log(ph.mu))/log(ph.sigma)-1/log(ph.mu))*log(ph.mu)^2  # this is the beta distribution's alpha parameter
#  ph.beta <- log(ph.alpha)*(1/log(ph.mu)-1)                   # this is the beta distribution's beta parameter
  
  ph.alpha <- ((1-ph.mu)/ph.sigma-1/ph.mu)*ph.mu^2  # this is the beta distribution's alpha parameter
  ph.beta <- ph.alpha*(1/ph.mu-1)                   # this is the beta distribution's beta parameter
  
  pb.pi <- inv.logit(parvect[(6*N-7):(7*N-10)]) # N-2 parameters: (2*N+(4*(N-2))+1):(2*N+(5*(N-2)))
  
  sal.mu <- parvect[(7*N-9):(8*N-12)] # N-2 parameters: (2*N+(5*(N-2))+1):(2*N+(6*(N-2)))
  sal.sigma <- parvect[(8*N-11):(9*N-14)] # N-2 parameters: (2*N+(6*(N-2))+1):(2*N+(7*(N-2)))
  
  # the beta vector includes: N each coeff, times 3 covars  
  # (3+1=4 coefficients, including intercept) plus 1 interaction term (with 2 elements) -> total 4+2=6 beta coefficients
  end_beta_ind <- ((9*N-14)+(3+(length(cov.names)))*N*(N-1))
  beta <- matrix(parvect[(9*N-13):end_beta_ind],nrow=N*(N-1)) 
  
  delta <- exp(c(0,parvect[(end_beta_ind+1):(end_beta_ind+N-1)])) # length N in total, but N-1 in parvect
  delta <- delta/sum(delta) 
  
  return(list(dd.mu=dd.mu,dd.sigma=dd.sigma,
              hd.mu=hd.mu,hd.sigma=hd.sigma,
              ph.mu=ph.mu,ph.sigma=ph.sigma,
              ph.alpha=ph.alpha, ph.beta=ph.beta,
              pb.pi=pb.pi,
              sal.mu=sal.mu,sal.sigma=sal.sigma,
              beta=beta,delta=delta))
}




## mllk
## function that computes minus the log-likelihood of the movement HMM
mllk <- function(parvect,obslist,N,cov.names){
  
  lpn <- pw2pn(parvect,N)
  beta.mat <- t(lpn$beta)
  mllkall <- 0
  K <- length(obslist)
  llk <- rep(NA, K) # one llk per segment
  lscale <- c()
  
  for (k in 1:K){ # loop through segments within each individual k
    
    cov.mat <- matrix(1, nrow=nrow(obslist[[k]]), ncol=length(cov.names)+3) # cov.names plus intercept (1) plus interaction (2)
    for (cov in 1:length(cov.names)){
      cov.mat[,cov+1] <- obslist[[k]][,cov.names[cov]]
      if(cov==length(cov.names)){
        cov.mat[,cov+2] <- obslist[[k]][,cov.names[cov-1]]*obslist[[k]][,cov.names[cov]] # column 3+2=5 holds the interaction of covariate 2 with covariate 3
        cov.mat[,cov+3] <- obslist[[k]][,cov.names[cov-2]]*obslist[[k]][,cov.names[cov]] # column 3+3=6 holds the interaction of covariate 1 with covariate 3
      }
    }
    
    ind.ho <- which(obslist[[k]]$source=="ho") # index for all haulout observations 
    ind.sf <- which(obslist[[k]]$source=="sf") # index for all surface observations
    ind.dv1 <- which(obslist[[k]]$source=="dv") # index for all diving observations
    ind.dv2 <- which(obslist[[k]]$source=="dv" & !is.na(obslist[[k]]$hunt_avgdep)) # index for diving observations without missing hunting depth
    
    n <- dim(obslist[[k]])[1] # number of obs within a segment
    allprobs <- matrix(rep(1, N*n), nrow=n) # one allprobs matrix per segment
    
    ## ## ## KNOWN STATES
    # in state 1 (haulout) only duration contributes to state density
    j <- 1
    
    # DURATION (gamma)
    dd.prob <- rep(0,n) # assume 0 probability if not in ho
    dd.prob[ind.ho] <- dgamma(obslist[[k]][ind.ho,"DURATION"], 
                              shape=lpn$dd.mu[j]^2/lpn$dd.sigma[j]^2, 
                              scale=lpn$dd.sigma[j]^2/lpn$dd.mu[j])
    allprobs[,j] <- dd.prob
    
    # in state 2 (surface) only duration contributes to state density
    j <- 2
    
    # DURATION (gamma)
    dd.prob <- rep(0,n) # assume 0 probability if not in sf
    dd.prob[ind.sf] <- dgamma(obslist[[k]][ind.sf,"DURATION"], 
                              shape=lpn$dd.mu[j]^2/lpn$dd.sigma[j]^2, 
                              scale=lpn$dd.sigma[j]^2/lpn$dd.mu[j])
    
    allprobs[,j] <- dd.prob
    ## ## ## KNOWN STATES
    
    for (j in 3:N){ # loop through dive states

      dd.prob <- hd.prob <- ph.prob <- pb.prob <- sal.prob <- rep(0,n) # note 0 (if not in a dive state then no depth etc)
      dd.prob[ind.dv1] <- hd.prob[ind.dv1] <- ph.prob[ind.dv1] <- pb.prob[ind.dv1] <- sal.prob[ind.dv1] <- 1 # 1 by default if in diving state
      
      # DURATION (gamma)
      dd.prob[ind.dv1] <- dgamma(obslist[[k]][ind.dv1,"DURATION"], 
                                 shape=lpn$dd.mu[j]^2/lpn$dd.sigma[j]^2, 
                                 scale=lpn$dd.sigma[j]^2/lpn$dd.mu[j])
      
      # depth of interest (gamma) - this either the average hunting depth OR the maxdep if there was no hunting
      hd.prob[ind.dv2] <- dgamma(obslist[[k]][ind.dv2,"hunt_avgdep"], 
                                 shape=lpn$hd.mu[j-2]^2/lpn$hd.sigma[j-2]^2, 
                                 scale=lpn$hd.sigma[j-2]^2/lpn$hd.mu[j-2])        
      
      # proportion time spent hunting (beta)
      ph.prob[ind.dv2] <- dbeta(obslist[[k]][ind.dv2,"hunt_prop"], # note j-2 (only states 3-N so 3 params)
                                shape1=lpn$ph.alpha[j-2],
                                shape2=lpn$ph.beta[j-2])
      
      # proportion bathymetry reached (beta)
      pb.prob[ind.dv2] <- obslist[[k]][ind.dv2,"benthic"]*lpn$pb.pi[j-2] + 
        obslist[[k]][ind.dv2,"notbenthic"]*(1-lpn$pb.pi[j-2])
      
      # salinity at hunting depth (normal)
      sal.prob[ind.dv2] <- dnorm(obslist[[k]][ind.dv2,"psal"], 
                                 mean=lpn$sal.mu[j-2], 
                                 sd=lpn$sal.sigma[j-2])
      
      allprobs[,j] <- dd.prob*hd.prob*ph.prob*pb.prob*sal.prob
      
    }
    # # Rcpp version
    # returns the negative log likelihood, so you don't need to take the negative below
    lscale[k] <- nLogLike_rcpp(nbStates=N, nbObs=n, allProbs=allprobs, 
                               beta=beta.mat, covs=cov.mat, delta=lpn$delta)  
    
    # # non-Rcpp version
    # lscale <- 0
    # 
    # for (i in 1:n){ # loop through all obs in a segment
    #   if (i==1){
    #     foo <- lpn$delta
    #     foo <- foo*allprobs[1,] # estimating initial distribution, comment out to assume fixed
    #   }
    #   gamma <- diag(N)
    #   
    #   gamma[!gamma] <- exp(lpn$beta[,1]+
    #                          lpn$beta[,2]*cov.mat[i,2]+ # coshr
    #                          lpn$beta[,3]*cov.mat[i,3]+ # sinhr
    #                          lpn$beta[,4]*cov.mat[i,4]+ # week_stnd
    #                          lpn$beta[,5]*cov.mat[i,5]+ # sinhr*week_stnd
    #                          lpn$beta[,6]*cov.mat[i,6]) # coshr*week_stnd
    #   gamma <- gamma/apply(gamma,1,sum)
    #   
    #   foo <- foo%*%gamma*allprobs[i,]
    #   sumfoo <- sum(foo); lscale <- lscale+log(sumfoo); foo <- foo/sumfoo # scaling
    # } # closes i loop
    # 
    # llk[k] <- lscale
  } # closes k loop (segment loop)
  
  llkallseg <- sum(lscale)
  
  return(llkallseg)
}



## mle
## function that runs the numerical minimization of 'mllk' 
## (i.e. searches for the Maximum Likelihood Estimates of the parameters of the objective function)
mle <- function(obslist,
                dd.mu0,dd.sigma0,
                hd.mu0,hd.sigma0,
                ph.mu0,ph.sigma0,
                pb.pi0,
                sal.mu0,sal.sigma0,
                beta0,delta0,
                N,cov.names,return_hessian){
  parvect <- pn2pw(dd.mu0,dd.sigma0,
                   hd.mu0,hd.sigma0,
                   ph.mu0,ph.sigma0,
                   pb.pi0,
                   sal.mu0,sal.sigma0,
                   beta0,delta0,N)
  mod <- nlm(mllk,parvect,obslist,N,cov.names,hessian=return_hessian,
             print.level=2,iterlim=20000,steptol=1e-8)
  pn <- pw2pn(mod$estimate,N)
  outlist <- list(dd.mu=pn$dd.mu,dd.sigma=pn$dd.sigma,
                   hd.mu=pn$hd.mu,hd.sigma=pn$hd.sigma,
                   ph.mu=pn$ph.mu,ph.sigma=pn$ph.sigma,
                   ph.alpha=pn$ph.alpha,ph.beta=pn$ph.beta,
                   pb.pi=pn$pb.pi,
                   sal.mu=pn$sal.mu,sal.sigma=pn$sal.sigma,
                   beta=pn$beta,delta=pn$delta,
                   mllk=mod$minimum,hessian=mod$hessian)
  return(outlist)
}










