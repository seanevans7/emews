####################################################################################################################################
### Function used to implement the Viterbi algorithm for calculating the most likely state sequence from a                       ###
### hidden Markov model to diving data from Weddell seals.                                                                       ### 
###                                                                                                                              ###
### Used in the analysis presented in:                                                                                           ###
### "Sex-specific variation in the use of vertical habitat by a resident Antarctic top predator"                                 ###
### Theoni Photopoulou, Karine Heerah, Jennifer Pohle and Lars Boehme (2020)                                                     ###
####################################################################################################################################

# Viterbi algorithm for hidden Markov model including two interactions
viterbi <- function(obslist,cov.names,mod,N){
  
  K <- length(obslist)
  iv.seg <- vector("list")
  beta.mat <- t(mod$beta)
  
  for (k in 1:K){ # loop through data from each individual
    
    cov.mat <- matrix(1, nrow=nrow(obslist[[k]]), ncol=length(cov.names)+3) # cov.names plus intercept (1) plus interactions (2)
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
    
    n <- dim(obslist[[k]])[1] # obs within a segment
    allprobs <- matrix(rep(1, N*n), nrow=n) # one allprobs matrix per segment
    
    ## ## ## KNOWN STATES
    # in state 1 (haulout) only duration contributes to state density
    j <- 1
    
    # DURATION (gamma)
    dd.prob <- rep(0,n) # assume 0 probability if not in ho
    dd.prob[ind.ho] <- dgamma(obslist[[k]][ind.ho,"DURATION"], 
                              shape=mod$dd.mu[j]^2/mod$dd.sigma[j]^2, 
                              scale=mod$dd.sigma[j]^2/mod$dd.mu[j])
    allprobs[,j] <- dd.prob
    
    # in state 2 (surface) only maxdep and duration contribute to state density
    j <- 2
    
    # DURATION (gamma)
    dd.prob <- rep(0,n) # assume 0 probability if not in sf
    dd.prob[ind.sf] <- dgamma(obslist[[k]][ind.sf,"DURATION"], 
                              shape=mod$dd.mu[j]^2/mod$dd.sigma[j]^2, 
                              scale=mod$dd.sigma[j]^2/mod$dd.mu[j])
    
    allprobs[,j] <- dd.prob
    ## ## ## KNOWN STATES
    
    for (j in 3:N){ # loop through dive states
      
      dd.prob <- hd.prob <- ph.prob <- pb.prob <- sal.prob <- rep(0,n) # note 0 (if not in a dive state then no depth etc)
      dd.prob[ind.dv1] <- hd.prob[ind.dv1] <- ph.prob[ind.dv1] <- pb.prob[ind.dv1] <- sal.prob[ind.dv1] <- 1 # 1 by default if in diving state
      
      # DURATION (gamma)
      dd.prob[ind.dv1] <- dgamma(obslist[[k]][ind.dv1,"DURATION"], 
                                 shape=mod$dd.mu[j]^2/mod$dd.sigma[j]^2, 
                                 scale=mod$dd.sigma[j]^2/mod$dd.mu[j])
      
      # depth of interest (gamma) - this either the average hunting depth OR the maxdep if there was no hunting
      hd.prob[ind.dv2] <- dgamma(obslist[[k]][ind.dv2,"hunt_avgdep"], 
                                 shape=mod$hd.mu[j-2]^2/mod$hd.sigma[j-2]^2, 
                                 scale=mod$hd.sigma[j-2]^2/mod$hd.mu[j-2])        
      
      # proportion time spent hunting (beta)
      ph.prob[ind.dv2] <- dbeta(obslist[[k]][ind.dv2,"hunt_prop"], # note j-2 (only states 3-N so 3 params)
                                shape1=mod$ph.alpha[j-2],
                                shape2=mod$ph.beta[j-2])
      
      # proportion bathymetry reached (beta)
      pb.prob[ind.dv2] <- obslist[[k]][ind.dv2,"benthic"]*mod$pb.pi[j-2] + 
        obslist[[k]][ind.dv2,"notbenthic"]*(1-mod$pb.pi[j-2])
      
      # salinity at hunting depth (normal)
      sal.prob[ind.dv2] <- dnorm(obslist[[k]][ind.dv2,"psal"], 
                                 mean=mod$sal.mu[j-2], 
                                 sd=mod$sal.sigma[j-2])
      
      allprobs[,j] <- dd.prob*hd.prob*ph.prob*pb.prob*sal.prob
      
    } # closes j loop
    
    xi <- matrix(0,n,N)
    foo <- mod$delta*allprobs[1,]
    xi[1,] <- foo/sum(foo)
    
    trMat <- moveHMM:::trMatrix_rcpp(N,beta.mat,cov.mat)
    
    for (i in 2:n){ # forward loop
      foo <- apply(xi[i-1,]*trMat[,,i],2,max)*allprobs[i,]
      xi[i,] <- foo/sum(foo)
    } # closes i loop (forward)
    
    iv <- numeric(n) # state sequence
    iv[n] <- which.max(xi[n,]) 
    
    for (h in (n-1):1){ # backward loop
      iv[h] <- which.max(trMat[,iv[h+1],h+1]*xi[h,])
    } # closes h loop (backward)
    
    iv.seg[[k]] <- iv
    
  } # closes k loop
  return(iv.seg)
} # end function

