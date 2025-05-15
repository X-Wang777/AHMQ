//#include <RcppArmadillo.h>
#include <rgen.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(rgen)]]

using namespace Rcpp;

// [[Rcpp::export]]
double dina_logL(const arma::mat y,unsigned int N, unsigned int J, unsigned int K, const arma::vec s, const arma::vec g, const arma::mat alpha,const arma::mat Q)
{ double logL = 0.;
  double qq;
  arma::mat AQ = alpha * Q.t();
  for (unsigned int j = 0; j < J; j++) {
    
    arma::vec ETA(N,arma::fill::zeros);
    qq = arma::conv_to<double>::from(Q.row(j) * (Q.row(j)).t());
    ETA.elem(arma::find(AQ.col(j) == qq)).fill(1.0);
    for (unsigned int i = 0; i < N; i++) {
      double eta_ij = ETA(i);
      logL += log(((1-s(j))*eta_ij+g(j)*(1-eta_ij)))*y(i,j)+log((s[j]*eta_ij+(1-g[j])*(1-eta_ij)))*(1-y(i,j));
    }
    
  }
  return  logL;
}
// [[Rcpp::export]]
arma::vec dina_logL_vec(const arma::mat y,unsigned int N, unsigned int J, unsigned int K, int chain,
                        const arma::mat s, const arma::mat g, const arma::cube alpha,const arma::cube Q){
  arma::vec logL(chain);
  for(int i = 0; i < chain; i++){
    logL(i)=dina_logL(y,N, J, K, s.col(i),  g.col(i),alpha.slice(i), Q.slice(i));
  }
  return  logL;
}

// [[Rcpp::export]]
double bar_Dev_theta(const arma::mat y,unsigned int N, unsigned int J, unsigned int K, int chain,
                     const arma::mat s, const arma::mat g, const arma::cube alpha,const arma::cube Q){
  arma::vec logL(chain);
  for(int i = 0; i < chain; i++){
    logL(i)=dina_logL(y,N, J, K, s.col(i),  g.col(i),alpha.slice(i), Q.slice(i));
  }
  return  -2*mean(logL);
}
// [[Rcpp::export]]
double compute_DIC(const arma::mat y,unsigned int N, unsigned int J, unsigned int K, int chain,
                   const arma::mat s, const arma::mat g, const arma::cube alpha,const arma::cube Q){
  arma::vec logL(chain);
  for(int i = 0; i < chain; i++){
    logL(i)=dina_logL(y,N, J, K, s.col(i),  g.col(i),alpha.slice(i), Q.slice(i));
  }
  double bar_Dev_theta =  -2*mean(logL);
  double Dev_theta_hat = min(-2 * logL);
  double P_D = bar_Dev_theta - Dev_theta_hat;
  double DIC = Dev_theta_hat + 2 * P_D;
  return DIC;
}
