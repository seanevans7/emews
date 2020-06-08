#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double nLogLike_rcpp(int nbStates, int nbObs, arma::mat allProbs , arma::mat beta, arma::mat covs, arma::rowvec delta)
  {

    //=======================================================//
    // 1. Computation of transition probability matrix trMat //
    //=======================================================//

    arma::cube trMat(nbStates,nbStates,nbObs);
    trMat.zeros();
    arma::mat rowSums(nbStates,nbObs);
    rowSums.zeros();

    arma::mat g(nbObs,nbStates*(nbStates-1));

    if(nbStates>1) {
        g = covs*beta;

        for(int k=0;k<nbObs;k++) {
            int cpt=0; // counter for diagonal elements
            for(int i=0;i<nbStates;i++) {
                for(int j=0;j<nbStates;j++) {
                    if(i==j) {
                        // if diagonal element, set to one and increment counter
                        trMat(i,j,k)=1;
                        cpt++;
                    }
                    else
                        trMat(i,j,k) = exp(g(k,i*nbStates+j-cpt));

                    // keep track of row sums, to normalize in the end
                    rowSums(i,k)=rowSums(i,k)+trMat(i,j,k);
                }
            }
        }

        // normalization
        for(int k=0;k<nbObs;k++)
            for(int i=0;i<nbStates;i++)
                for(int j=0;j<nbStates;j++)
                    trMat(i,j,k) = trMat(i,j,k)/rowSums(i,k);
    }


    //======================//
    // 3. Forward algorithm //
    //======================//

    arma::mat Gamma(nbStates,nbStates); // transition probability matrix
    double lscale = 0; // scaled log-likelihood
    //int k=0; // animal index
    arma::rowvec alpha(nbStates);

    for(unsigned int i=0;i<allProbs.n_rows;i++) {

            Gamma = trMat.slice(i);

        if(i==0) {  //k<aInd.size() && i==(unsigned)(aInd(k)-1)) {
            // if 'i' is the 'k'-th element of 'aInd', switch to the next animal
            //k++;
            alpha = delta % allProbs.row(i);
        } else {
            alpha = (alpha * Gamma) % allProbs.row(i);
        }

        lscale = lscale + log(sum(alpha));
        alpha = alpha/sum(alpha);
    }

    return -lscale;
}