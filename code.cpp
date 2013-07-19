#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat bigM(
		int nSite, arma::vec migRate, arma::vec fecund, arma::vec survJ, arma::vec survA
) {
	arma::mat M(nSite*2, nSite*2);
	M.zeros();
   for (int i=1; i <= nSite; ++i) {
     for (int j=1; j <= nSite; ++j) {
			 M((2*i-1)-1,(2*j)-1) = (fecund(j-1)*migRate(j-1)/(double(nSite)-1.0));
       M((2*i)-1,(2*j-1)-1) = (survJ(j-1)*migRate(j-1)/(double(nSite)-1.0));
       M((2*i)-1,(2*j)-1) = (survA(j-1)*migRate(j-1)/(double(nSite)-1.0)); 
     } 
	 	 M((2*i-1)-1,(2*i)-1) = fecund(i-1)*(1.0-migRate(i-1));
     M((2*i)-1,(2*i-1)-1) = survJ(i-1)*(1.0-migRate(i-1));
     M((2*i)-1,(2*i)-1) = survA(i-1)*(1.0-migRate(i-1));    
	 }
	return M;
}


