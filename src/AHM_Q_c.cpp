#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <rgen.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(rgen)]]

using namespace Rcpp;

// [[Rcpp::export]]
arma::vec inv_bijectionvector(unsigned int K,double CL){
  arma::vec alpha(K);
  for(unsigned int k=0;k<K;k++){
    double twopow = pow(2,K-k-1);
    alpha(k) = (twopow<=CL);
    CL = CL - twopow*alpha(k);
  }
  return alpha;
}

// [[Rcpp::export]]
double twoToten(arma::ivec x)
{
  int K=x.n_elem;
  double res=0;
  for(int i=0;i<K;i++){
    res+=x(K-i-1)*pow(2,i);
  }
  return res;
}
// [[Rcpp::export]]
arma::mat ff(arma::mat a){
  arma::mat temp;
  temp=a;temp.diag().zeros();
  return temp;
}
// [[Rcpp::export]]
double pYi(const arma::vec& ETA_i,const arma::vec& Y_i,const arma::vec& ss,
           const arma::vec& gs){
  arma::vec one_m_ss = 1. - ss;
  arma::vec one_m_gs = 1. - gs;
  arma::vec one_m_ETA_i = 1. - ETA_i;
  arma::vec one_m_Y_i = 1. - Y_i;
  
  arma::vec ps = Y_i%(one_m_ss%ETA_i + gs%one_m_ETA_i) + one_m_Y_i%(ss%ETA_i + one_m_gs%one_m_ETA_i);
  
  return arma::prod(ps);
}

// [[Rcpp::export]]
arma::mat Boolean(arma::mat A){
  arma::mat R;
  R=arma::conv_to<arma::mat>::from(A>0);
  return R;
}

// [[Rcpp::export]]
arma::vec Booleanvec(arma::vec A){
  arma::vec R;
  R=arma::conv_to<arma::vec>::from(A>0);
  return R;
}

// [[Rcpp::export]]
arma::mat Reachability(const arma::mat& StrucMat,unsigned int K){
  arma::mat R=StrucMat;
  R.diag().zeros();
  arma::mat Identity=arma::eye(K,K);
  arma::mat Rnext=R+Identity;
  while(accu(R==Rnext)<(K*K)){
    R=Rnext;
    Rnext=Boolean(R*(R+Identity));
  }
  R.diag().ones();
  return R;
}
// [[Rcpp::export]]
arma::mat ConnectMat(const arma::mat& R,unsigned int K){
  arma::mat Connect(K,K);Connect.eye();
  for(unsigned int i=0;i<(K-1);i++){
    for(unsigned int j=(i+1);j<K;j++){
      if((R(i,j)==1)||(R(j,i)==1)){
        Connect(i,j)=1;Connect(j,i)=1;
      }
    }
  }
  return Connect;
}

// [[Rcpp::export]]
arma::mat  Transitive(const arma::mat& G,unsigned int K){
  arma::mat G_red=G;
  arma::mat R_c=Reachability(G_red,K);
  arma::mat C_c=ConnectMat(R_c,K);
  for(unsigned int i=0;i<K;i++){
    for(unsigned int j=0;j<K;j++){
      if(G(i,j)==1){
        arma::mat Gtemp=G_red;Gtemp(i,j)=0;
        arma::mat Rtemp=Reachability(Gtemp,K);
        arma::mat Ctemp=ConnectMat(Rtemp,K);
        if(arma::all(arma::vectorise(Ctemp==C_c))){
          G_red=Gtemp;
        }
      }
    }
  }
  return(G_red);
}

// [[Rcpp::export]]
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}

// [[Rcpp::export]]
Rcpp::List AHM_sg(const arma::mat &Y, const arma::mat &Q,
                  const arma::mat &ALPHAS, const arma::vec &ss_old,
                  double as0, double bs0, double ag0, double bg0)
{
  
  unsigned int N = Y.n_rows;
  unsigned int J = Y.n_cols;
  
  arma::vec ETA;
  arma::vec ss_new(J);
  arma::vec gs_new(J);
  arma::mat AQ = ALPHAS * Q.t();
  double T, S, G, y_dot_eta, qq, ps, pg;
  double ug, us;
  
  for (unsigned int j = 0; j < J; j++) {
    us = R::runif(0, 1);
    ug = R::runif(0, 1);
    ETA = arma::zeros<arma::vec>(N);
    qq = arma::conv_to<double>::from(Q.row(j) * (Q.row(j)).t());
    ETA.elem(arma::find(AQ.col(j) == qq)).fill(1.0);
    
    y_dot_eta = arma::conv_to<double>::from((Y.col(j)).t() * ETA);
    T = sum(ETA);
    S = T - y_dot_eta;
    G = sum(Y.col(j)) - y_dot_eta;
    
    // sample s and g as linearly truncated bivariate beta
    
    // draw g conditoned upon s_t-1
    pg = R::pbeta(1.0 - ss_old(j), G + ag0, N - T - G + bg0, 1, 0);
    gs_new(j) = R::qbeta(ug * pg, G + ag0, N - T - G + bg0, 1, 0);
    // draw s conditoned upon g
    ps = R::pbeta(1.0 - gs_new(j), S + as0, T - S + bs0, 1, 0);
    ss_new(j) = R::qbeta(us * ps, S + as0, T - S + bs0, 1, 0);
  }
  return Rcpp::List::create(Rcpp::Named("ss_new") = ss_new,
                            Rcpp::Named("gs_new") = gs_new);
}

// [[Rcpp::export]]
arma::mat maipow(arma::mat a,arma::mat b){
  int m=a.n_rows;
  int n=a.n_cols;
  arma::mat re(m,n);
  for(int i=0; i<m ;i++){
    for(int j=0;j<n;j++){
      re(i,j) = pow(a(i,j),b(i,j));
    }
  }
  return re;
}
// [[Rcpp::export]]
arma::vec vecpow_1_to_vec(double a,arma::vec b){
  int n=b.n_elem;
  arma::vec re(n);
  for(int i=0;i<n;i++){
    re(i) = pow(a,b(i));
  }
  return re;
}
// [[Rcpp::export]]
arma::vec vecpow(arma::vec a,arma::vec b){
  int n=b.n_elem;
  arma::vec re(n);
  for(int i=0;i<n;i++){
    re(i) = pow(a(i),b(i));
  }
  return re;
}


// [[Rcpp::export]]
Rcpp::List add_path_sample1(arma::mat R,arma::mat G,int K)// sample an added edge
{
  arma::mat G_new=G;
  arma::mat R_double = R+R.t();
  //R_double.elem(arma::find(R_double > 0)).fill(1.0);
  arma::uvec index=arma::find(R_double == 0);
  int n_edge=index.n_elem;
  int change=1;
  int sample_location;
  double prob=0;
  if(n_edge==0){
    change = 0;
  }else{
    prob=1.0/n_edge;
    sample_location=rgen::rmultinomial(arma::ones(n_edge,1)/n_edge);
    G_new(index(sample_location))=1;
  }
  return Rcpp::List::create(Rcpp::Named("G_new") = G_new,Rcpp::Named("change") = change,Rcpp::Named("prob") = prob);
}


// [[Rcpp::export]]
Rcpp::List reduce_path_sample1(arma::mat G,int K)// sample an removed edge
{
  arma::uvec index=arma::find(G == 1);
  int n_edge=index.n_elem;
  arma::mat G_new = G;
  int sample_location;
  double prob=0.0;
  int change = 1;
  if(n_edge==0){
    change = 0;
  }else{
    sample_location=rgen::rmultinomial(arma::ones(n_edge,1)/n_edge);
    G_new(index(sample_location))=0;
    prob=1.0/n_edge;
  }
  return Rcpp::List::create(Rcpp::Named("G_new") = G_new,Rcpp::Named("change") = change,Rcpp::Named("prob") = prob);
}


// [[Rcpp::export]]
arma::vec Trans_10to2(unsigned int K,double CL){
  arma::vec alpha(K);
  for(unsigned int k=0;k<K;k++){
    double twopow = pow(2,K-k-1);
    alpha(k) = (twopow<=CL);
    CL = CL - twopow*alpha(k);
  }
  return alpha;
}

// [[Rcpp::export]]
arma::mat Trans_10to2_mat(unsigned int K,const arma::vec& CL) {
  unsigned int Col=CL.n_elem;
  arma::mat alpha(Col,K);
  for(unsigned int i=0;i<Col;i++){
    arma::colvec alphai(K);double cl=CL(i);
    for(unsigned int k=0;k<K;k++){
      double twopow = pow(2,K-k-1);
      alphai(k) = (twopow<=cl);
      cl = cl - twopow*alphai(k);
    }
    alpha.row(i)=alphai.t();
  }
  return alpha;
}


// [[Rcpp::export]]
double Accept_ratio_piprior(arma::mat y,arma::vec s,arma::vec g,arma::vec pi0,arma::vec delta0,arma::mat alpha,arma::mat alpha_new,
                            arma::mat Q,arma::mat Q_new,arma::mat G,arma::mat G_new,arma::vec alpha_possible_index,
                            arma::vec alpha_new_possible_index,arma::vec alpha_possible_code,arma::vec alpha_new_possible_code,
                            arma::mat alpha_possible,arma::mat alpha_new_possible,double pr_G,int action,double p1,double p2)
{ 
  int N=y.n_rows;
  int J=y.n_cols;
  int K=alpha.n_cols;
  int Cg = alpha_possible_index.n_elem;
  int Cg_new = alpha_new_possible_index.n_elem;
  arma::mat R_new = Reachability(G_new,K);
  arma::vec pi = arma::zeros(Cg,1);
  arma::vec pi_new = arma::zeros(Cg_new,1);
  arma::vec delta = arma::zeros(Cg,1);
  arma::vec delta_new = arma::zeros(Cg_new,1);
  arma::mat eta(J,Cg);
  arma::mat eta_new(J,Cg_new);
  // double Ratio;
  double Ratio1;
  double Ratio2;
  double Ratio3;
  arma::mat alpha_possible_all;
  arma::mat alpha_new_possible_all;
  arma::vec aa;
  
  alpha_possible_all = Trans_10to2_mat(K,alpha_possible_index);
  alpha_new_possible_all = Trans_10to2_mat(K,alpha_new_possible_index);
  for(int c=0;c<Cg;c++){
    pi(c) = sum(pi0.elem(arma::find(alpha_possible_code==alpha_possible_index(c))));
    delta(c) = sum(delta0.elem(arma::find(alpha_possible_code==alpha_possible_index(c))));
    for(int j=0;j<J;j++){
      aa = arma::cumprod(vecpow((alpha_possible_all.row(c)).t(),(Q.row(j)).t()));
      eta(j,c)=aa(K-1);
    }
  }
  for(int c=0;c<Cg_new;c++){
    pi_new(c) = sum(pi0.elem(arma::find(alpha_new_possible_code==alpha_new_possible_index(c))));
    delta_new(c) = sum(delta0.elem(arma::find(alpha_new_possible_code==alpha_new_possible_index(c))));
    for(int j=0;j<J;j++){
      aa = arma::cumprod(vecpow((alpha_new_possible_all.row(c)).t(),(Q_new.row(j)).t()));
      eta_new(j,c)=aa(K-1);
    }
  }
  
  //Ratio part 1: P(Y)
  arma:: mat pr(N,Cg);
  arma:: mat pr_new(N,Cg_new);
  arma::vec pr_ratio(N);
  for(int i=0;i<N;i++){
    for(int c=0;c<Cg;c++){
      //aa = arma::cumprod(vecpow(vecpow((1-s),eta.col(c))%vecpow(g,(1-eta.col(c))),(y.row(i)).t())%
      // vecpow(vecpow(s,eta.col(c))%vecpow((1-g),(1-eta.col(c))),(1-(y.row(i)).t())));
      //pr(i,c) =  aa(J-1);
      pr(i,c) = pYi(eta.col(c),(y.row(i)).t(),s,g);
    }
    
    for(int c=0;c<Cg_new;c++){
      //aa = arma::cumprod(vecpow(vecpow((1-s),eta_new.col(c))%vecpow(g,(1-eta_new.col(c))),(y.row(i)).t())%
      // vecpow(vecpow(s,eta_new.col(c))%vecpow((1-g),(1-eta_new.col(c))),(1-(y.row(i)).t())));
      //pr_new(i,c) =  aa(J-1);
      pr_new(i,c) = pYi(eta_new.col(c),(y.row(i)).t(),s,g);
    }
    pr_ratio(i) = arma::sum((pr_new.row(i)).t()%pi_new)/arma::sum((pr.row(i)).t()%pi);
  }
  aa = arma::cumprod(pr_ratio);
  Ratio1 = aa(N-1);
  //Ratio part 2:P(G)
  //double pr_G = prob_edge;
  double pr_G_new;
  if(action==1){
    arma::uvec index=arma::find(G_new==1);
    int n_edge=index.n_elem;
    pr_G_new = 1.0/n_edge;
    Ratio2 = pr_G_new/pr_G;
    Ratio3 = p1/p2;
  }else{
    R_new = R_new+R_new.t();
    arma::uvec index=arma::find(R_new==0);
    int n_edge=index.n_elem;
    pr_G_new = 1.0/n_edge;
    Ratio2 = pr_G_new/pr_G;
    Ratio3 = p2/p1;
  }
  return Ratio1*Ratio2*Ratio3;
  //return Rcpp::List::create(Rcpp::Named("Ratio1") =Ratio1);
  
}

// [[Rcpp::export]]
double Trans_2to10(arma::vec x,int K){
  double r=0.0;
  for(int k=0;k<K;k++){
    r=r+pow(2,k)*x(K-k-1);
  }
  return r;
}

// [[Rcpp::export]]
arma::vec Trans_2to10_mat(arma::mat x,int K){
  int n=x.n_rows;
  arma::vec code(n);
  for(int i=0;i<n;i++){
    code(i)=Trans_2to10((x.row(i)).t(),K);
  }
  return code;
}
// [[Rcpp::export]]
Rcpp::List Reduced_alpha1(arma::mat alpha_all,arma::mat R,int K)//all possible alpha
{
  int L= pow(2,K);
  arma::mat alpha_current = alpha_all;
  arma::vec delta_new = arma::zeros(L,1);
  arma::uvec prerequisite;
  arma::uvec if_prerequisite;
  arma::mat alpha_prerequisite;
  arma::vec aa=arma::ones(L,1);
  for(int k=0;k<K;k++){
    aa=arma::ones(L,1);
    prerequisite = arma::find(R.col(k)==1);
    alpha_prerequisite=alpha_current.cols(prerequisite);
    if((prerequisite).n_elem > 0){
      if_prerequisite = arma::find(arma::sum(alpha_prerequisite,1) < (prerequisite).n_elem);
      aa.elem(if_prerequisite) = arma::zeros(if_prerequisite.n_elem,1);
      alpha_current.col(k)=aa;
    }
  }
  arma::vec alpha_current_binary;
  arma::vec index;
  alpha_current_binary = Trans_2to10_mat(alpha_current,K);
  index=sort(unique(alpha_current_binary));
  
  return Rcpp::List::create(Rcpp::Named("alpha_current") = alpha_current,
                            Rcpp::Named("alpha_current_binary") = alpha_current_binary,
                            Rcpp::Named("index") = index);
}


// [[Rcpp::export]]
arma::mat AHM_alpha(arma::mat y,arma::vec s,arma:: vec g,arma::mat Q,arma::mat G,
                    arma::mat alpha_all,arma::vec pi0,int N,int J,int K,int L)
{
  arma::mat R;
  arma::vec alpha_current_possible_index;
  arma::vec alpha_current_possible_code;
  //arma::mat alpha_current_possible;
  arma::mat alpha_all_current;
  Rcpp::List rr;
  R=Reachability(G,K);
  rr = Reduced_alpha1(alpha_all,R,K);
  alpha_current_possible_index = Rcpp::as<arma::vec>(rr["index"]);
  alpha_current_possible_code = Rcpp::as<arma::vec>(rr["alpha_current_binary"]);
  //alpha_current_possible = Rcpp::as<arma::mat>(rr["alpha_current"]);
  alpha_all_current = Trans_10to2_mat(K,alpha_current_possible_index);
  
  int Cg=alpha_current_possible_index.n_elem;
  //int K=Q.n_cols;
  arma::mat pi_mat(N,Cg,arma::fill::ones);
  arma::vec pi_vec(N,arma::fill::zeros);
  arma::mat eta(N,J);
  arma::mat pr_mat(N,J);
  arma::vec pi0_new(Cg,arma::fill::zeros);
  arma::mat alpha;
  arma::mat alpha_new(N,K);
  arma::ivec alpha_class(N);
  arma::mat Q_j;
  arma::mat ee;
  //arma::mat ee2;
  arma::vec aa;
  arma::vec bb;
  for(int c=0;c<Cg;c++){
    alpha=arma::ones(N,1)*alpha_all_current.row(c);
    for(int j=0;j<J;j++){
      Q_j=arma::ones(N,1)*Q.row(j);
      ee=cumprod(maipow(alpha,Q_j),1);
      eta.col(j)=ee.tail_cols(1);
      pr_mat.col(j) = ((1-s(j))*eta.col(j)+g(j)*(1-eta.col(j)))%y.col(j)+
        (s(j)*eta.col(j)+(1-g(j))*(1-eta.col(j)))%(1-y.col(j));
      pi_mat.col(c)=pi_mat.col(c)%pr_mat.col(j);
    }
    aa=pi0.elem(arma::find(alpha_current_possible_code==alpha_current_possible_index(c)));
    pi0_new(c)=sum(aa);
    //ee2=cumprod(pr_mat,1);
    pi_mat.col(c)=pi_mat.col(c)*pi0_new(c);
    pi_vec= pi_vec+pi_mat.col(c);
    
  }
  arma::mat iden(1,Cg,arma::fill::ones);
  pi_mat =pi_mat/(pi_vec*iden);
  for(int i=0;i<N;i++){
    alpha_class(i) = rgen::rmultinomial((pi_mat.row(i)).t());
    alpha_new.row(i) = alpha_all_current.row(alpha_class(i));
  }
  return alpha_new;
}

// [[Rcpp::export]]
arma::vec rDirichlet(const arma::vec& deltas){
  unsigned int C = deltas.n_elem;
  arma::vec Xgamma(C);
  
  //generating gamma(deltac,1)
  for(unsigned int c=0;c<C;c++){
    Xgamma(c) = R::rgamma(deltas(c),1.0);
  }
  return Xgamma/sum(Xgamma);
}

// [[Rcpp::export]]
Rcpp::List AHM_alpha_pi(arma::mat y,arma:: vec s,arma:: vec g,arma::mat Q,arma::mat G,
                        arma::mat alpha_all,arma::vec pi0,arma::vec delta0,int N,int J,int K,int L)
{
  arma::mat R;
  arma::vec alpha_current_possible_index;
  arma::vec alpha_current_possible_code;
  //arma::mat alpha_current_possible;
  arma::mat alpha_all_current;
  Rcpp::List rr;
  R=Reachability(G,K);
  rr = Reduced_alpha1(alpha_all,R,K);
  alpha_current_possible_index = Rcpp::as<arma::vec>(rr["index"]);
  alpha_current_possible_code = Rcpp::as<arma::vec>(rr["alpha_current_binary"]);
  //alpha_current_possible = Rcpp::as<arma::mat>(rr["alpha_current"]);
  alpha_all_current = Trans_10to2_mat(K,alpha_current_possible_index);
  
  int Cg=alpha_current_possible_index.n_elem;
  //int K=Q.n_cols;
  arma::mat pi_mat(N,Cg,arma::fill::ones);
  arma::vec pi_vec(N,arma::fill::zeros);
  arma::mat eta(N,J);
  arma::mat pr_mat(N,J);
  arma::vec pi0_new(Cg,arma::fill::zeros);
  arma::vec delta0_new(L,arma::fill::zeros);
  arma::mat alpha;
  arma::mat alpha_new(N,K);
  arma::ivec alpha_class(N);
  arma::ivec pi_class(L);
  arma::mat Q_j;
  arma::mat ee;
  //arma::mat ee2;
  arma::vec aa;
  arma::vec bb;  
  arma::vec count(L,arma::fill::zeros);
  for(int c=0;c<Cg;c++){
    alpha=arma::ones(N,1)*alpha_all_current.row(c);
    for(int j=0;j<J;j++){
      Q_j=arma::ones(N,1)*Q.row(j);
      ee=cumprod(maipow(alpha,Q_j),1);
      eta.col(j)=ee.tail_cols(1);
      pr_mat.col(j) = ((1-s(j))*eta.col(j)+g(j)*(1-eta.col(j)))%y.col(j)+
        (s(j)*eta.col(j)+(1-g(j))*(1-eta.col(j)))%(1-y.col(j));
      pi_mat.col(c)=pi_mat.col(c)%pr_mat.col(j);
    }
    aa=pi0.elem(arma::find(alpha_current_possible_code==alpha_current_possible_index(c)));
    pi0_new(c)=sum(aa);
    aa=delta0.elem(arma::find(alpha_current_possible_code==alpha_current_possible_index(c)));
    delta0_new(alpha_current_possible_index(c))=sum(aa);
    //ee2=cumprod(pr_mat,1);
    pi_mat.col(c)=pi_mat.col(c)*pi0_new(c);
    pi_vec= pi_vec+pi_mat.col(c);
    
  }
  arma::mat iden(1,Cg,arma::fill::ones);
  pi_mat =pi_mat/(pi_vec*iden);
  int kk;
  for(int i=0;i<N;i++){
    alpha_class(i) = rgen::rmultinomial((pi_mat.row(i)).t());
    kk=alpha_class(i);
    alpha_new.row(i) = alpha_all_current.row(alpha_class(i));
    count(alpha_current_possible_index(kk))++;
  }
  arma::vec pi_new=rDirichlet(count+delta0_new);
  return Rcpp::List::create(Rcpp::Named("alpha_new") = alpha_new,
                            Rcpp::Named("pi_new") = pi_new,
                            Rcpp::Named("count") = count,
                            Rcpp::Named("delta0_new") = delta0_new,
                            Rcpp::Named("alpha_class") = alpha_class,
                            Rcpp::Named("alpha_current_possible_index")=alpha_current_possible_index);
}


// [[Rcpp::export]]
arma::mat AHM_Q(arma::mat y,arma::vec s,arma::vec g,arma::mat alpha,arma::mat Q_all,
                      arma::mat Q,int N,int J,int K,arma::vec alpha_current_possible_index,
                      arma::mat eta,arma::vec alpha_code)
{
  int D=pow(2,K)-1;
  int CC=alpha_current_possible_index.n_elem;
  //arma::mat pr_mat_j(N,D);
  //arma::mat eta_j(N,D);
  arma::vec pr_j(D,arma::fill::ones);
  arma::vec yj;
  arma::vec yy;
  arma::vec Q_class(J);
  arma::mat Q_new(J,K);
  arma::vec xi_j(K);
  arma::vec p_xi0(D);
  arma::vec p_xi(D);
  arma::vec pr_ratio(D);
  arma::vec aa;
  arma::vec bb;
  arma::vec dd;
  arma::vec cc;
  arma::mat Q_j;
  arma::mat ee;
  for(int j=0;j<J;j++){
    yj=y.col(j);
    p_xi0.zeros();
    p_xi.zeros();
    for(int k=0;k<K;k++){
      xi_j(k)=R::rbeta(1+Q(j,k),2-Q(j,k));
    }
    for(int c=0;c<D;c++){
      aa = cumprod(xi_j%(Q_all.row(c)).t()+(1-xi_j)%(1-Q_all.row(c)).t());
      p_xi0(c)=aa(K-1);
    }
    p_xi0 = p_xi0/sum(p_xi0);
    for(int c=0;c<D;c++){
      Q_j = arma::ones(N,1)*Q_all.row(c); 
      ee=cumprod(maipow(alpha,Q_j),1);
      pr_j(c)=1;
      for(int cg=0;cg<CC;cg++){
        yy=yj.elem(arma::find(alpha_code==alpha_current_possible_index(cg)));
        pr_j(c)=pr_j(c)*pow( pow(1-s(j),eta(alpha_current_possible_index(cg),c))*
          pow(g(j),1-eta(alpha_current_possible_index(cg),c)),sum(yy))*
          pow( pow(s(j),eta(alpha_current_possible_index(cg),c))*
          pow(1-g(j),1-eta(alpha_current_possible_index(cg),c)),yy.n_elem-sum(yy));
      }
      
    }
    for(int c=0;c<D;c++){
      for(int k=0;k<D;k++){
        pr_ratio(k) = pr_j(k)/pr_j(c)*(p_xi0(k))/(p_xi0(c));
      }
      p_xi(c) = 1/arma::sum(pr_ratio);
    }
    Q_class(j) = rgen::rmultinomial(p_xi);
    Q_new.row(j) = Q_all.row(Q_class(j));
  }
  //return Rcpp::List::create(Rcpp::Named("Q_new") = Q_new);
  return Q_new;
}


// [[Rcpp::export]]
arma::mat update_G(arma::mat y,arma::vec s,arma::vec g,arma::mat alpha,arma::mat alpha_all,arma::mat Q,arma::mat G,
                              arma::mat R,int N,int J,int K,int L,double p1,double p2,arma::vec pi,arma::vec delta0,arma::mat G_new,
                              double prob_edge,int action)
{
  
  
  Rcpp::List rr;
  arma::mat R_new;
  arma::vec alpha_code;
  arma::vec alpha_new_possible_index;
  arma::vec alpha_new_possible_code;
  arma::mat alpha_new_possible;
  arma::vec alpha_possible_index;
  arma::vec alpha_possible_code;
  arma::mat alpha_possible;
  arma::mat alpha_new(N,K);
  
  
  double accept=1.0;
  R_new = Reachability(G_new,K);
  alpha_code = Trans_2to10_mat(alpha, K);
  rr = Reduced_alpha1(alpha_all,R_new,K);
  alpha_new_possible_index = Rcpp::as<arma::vec>(rr["index"]);
  alpha_new_possible_code = Rcpp::as<arma::vec>(rr["alpha_current_binary"]);
  alpha_new_possible = Rcpp::as<arma::mat>(rr["alpha_current"]);
  for(int i=0;i<N;i++){
    alpha_new.row(i) = alpha_new_possible.row(alpha_code(i));
  }
  
  rr = Reduced_alpha1(alpha_all,R,K);
  alpha_possible_index = Rcpp::as<arma::vec>(rr["index"]);
  alpha_possible_code = Rcpp::as<arma::vec>(rr["alpha_current_binary"]);
  alpha_possible = Rcpp::as<arma::mat>(rr["alpha_current"]);
  
  double rand=R::runif(0,1);
  accept = Accept_ratio_piprior(y,s,g,pi,delta0,alpha,alpha_new,Q,Q,G,G_new,alpha_possible_index,alpha_new_possible_index,
                                alpha_possible_code,alpha_new_possible_code,alpha_possible,alpha_new_possible,prob_edge,action,p1,p2);
  arma::vec mm(2);
  mm(0)=1.0;
  mm(1)=accept;
  accept = arma::min(mm);
  if(rand > accept){
    G_new=G;
    action = 0;
  }
  
  
  
  return G_new;
}

// [[Rcpp::export]]
arma::mat AHM_G(arma::mat y,arma::vec s,arma::vec g,arma::mat alpha,arma::mat alpha_all,arma::mat Q,arma::mat G,
                        arma::mat R,int N,int J,int K,int L,double p1,double p2,arma::vec pi0,arma::vec delta0)
{
  int action=0;
  Rcpp::List operate;
  arma::mat G_new(K,K);
  int change;
  double prob_edge;
  //double dd;
  double u=R::runif(0,1);
  G_new=G;
  if(u < p1){
    //action = 1;
    operate = add_path_sample1(R,G,K);
    change = operate["change"];
    if(change == 1){
      G_new = Rcpp::as<arma::mat>(operate["G_new"]);
      prob_edge = operate["prob"];
      G_new=update_G(y,s,g,alpha,alpha_all,Q, G,R, N, J, K, L, p1,p2, pi0, delta0,G_new,prob_edge, action=1);
    }
  }
  if(u > p2){
    //action=2;
    operate = reduce_path_sample1(G,K);
    change = operate["change"];
    if(change == 1){
      G_new = Rcpp::as<arma::mat>(operate["G_new"]);
      prob_edge = operate["prob"];
      G_new=update_G(y,s,g,alpha,alpha_all,Q, G,R, N, J, K, L, p1,p2, pi0, delta0,G_new,prob_edge, action=2);
    }
  }
  
  return G_new;
}  
// [[Rcpp::export]]
arma::vec generate_sequence(int N) {
  arma::vec sequence = arma::linspace<arma::vec>(1, N, N);
  return sequence;
}

// [[Rcpp::export]]
Rcpp::IntegerVector sample_int(int N, int N1) {
  Rcpp::IntegerVector pool = Rcpp::seq(1, N);
  std::random_shuffle(pool.begin(), pool.end());
  return pool[Rcpp::Range(0, N1-1)]-1;
}


// [[Rcpp::export]]
arma::mat random_Q(unsigned int J,unsigned int K){
  
  //Generate identity matrices
  arma::vec one_K = arma::ones<arma::vec>(K);
  arma::mat I_K = arma::diagmat(one_K);
  //arma::mat Two_I_K = arma::join_cols(I_K,I_K);
  
  //generate Q1
  unsigned int JmK = J-K;
  unsigned int J1max = K;
  if(K>JmK){
    J1max = JmK;
  }
  unsigned int J1 = arma::conv_to< unsigned int >::from(arma::randi<arma::vec>(1,arma::distr_param(1,J1max) ) );
  arma::mat U1 = arma::randu<arma::mat>(J1,K);
  arma::mat Q1 = arma::zeros<arma::mat>(J1,K);
  
  //fix elements so columns are nonzero
  arma::vec row_ks = arma::randi<arma::vec>(K,arma::distr_param(0,J1-1) );
  for(unsigned int k=0;k<K;k++){
    Q1(row_ks(k),k) = 1;
  }
  
  Q1.elem(arma::find(Q1 > .5) ).fill(1.0);
  
  arma::mat Q = arma::join_cols(I_K,Q1);
  
  //Generating the remaining elements of Q in Q2 
  unsigned int JmKmJ1 = JmK - J1;
  arma::mat Q2 = arma::zeros<arma::mat>(JmKmJ1,K);
  if(JmKmJ1>0){
    arma::mat U2 = arma::randu<arma::mat>(JmKmJ1,K);
    Q2.elem(arma::find(U2 > .5) ).fill(1.0);
    Q = arma::join_cols(Q,Q2);
  }
  
  //Q
  arma::uvec P = arma::uvec(J);
  for(unsigned int j=0;j<J;j++){
    P(j)=j;
  }
  P = arma::shuffle(P);
  return Q.rows(P);
}

// [[Rcpp::export]]
void AHM_update(arma::mat Y,arma::mat& alpha,arma:: vec& s,arma:: vec& g,arma::mat& Q,arma::mat& G,arma::vec& pi,
                               double p1,double p2,arma::mat alpha_all,arma::mat Q_all,int N,int J,int K,int L,int N1,
                               double a_s0,double a_g0,double b_s0,double b_g0,arma::vec delta0,arma::mat eta)
{
  arma::mat y(N1,J);
  //update alpha, pi
  Rcpp::List r_alpha;
  r_alpha=AHM_alpha_pi(Y,s,g,Q,G,alpha_all,pi,delta0,N,J,K,L);
  alpha=Rcpp::as<arma::mat>(r_alpha["alpha_new"]);
  pi=Rcpp::as<arma::vec>(r_alpha["pi_new"]);
  arma::vec alpha_current_possible_index=Rcpp::as<arma::vec>(r_alpha["alpha_current_possible_index"]);
  //sample mini-batch
  arma::vec x = generate_sequence(N)-1;
  arma::ivec sampled_values = sample_int(N,N1);
  arma::uvec sample_index = arma::conv_to<arma::uvec>::from(sampled_values);
  y=Y.rows(sample_index);
  arma::mat alpha_mini=alpha.rows(sample_index);
  //update s,g
  Rcpp::List rr;
  rr = AHM_sg(y, Q,alpha_mini, s,a_s0,b_s0,a_g0,b_g0);
  s = Rcpp::as<arma::vec>(rr["ss_new"]);
  g = Rcpp::as<arma::vec>(rr["gs_new"]);
  //update Q
  arma::vec alpha_code1=Trans_2to10_mat(alpha,K);
  arma::vec alpha_code=alpha_code1.elem(sample_index);
  arma::mat R = Reachability(G,K);
  Q = AHM_Q(y,s,g,alpha_mini,Q_all,Q,N1,J,K,alpha_current_possible_index,eta,alpha_code);
  //updata G
  G = AHM_G(y,s,g,alpha_mini,alpha_all,Q,G,R,N1,J,K,L,p1,p2,pi,delta0);
  
}


// [[Rcpp::export]]
Rcpp::List AHMQ_MCMC(arma::mat Y,int K,arma::vec& s,arma::vec& g,arma::vec& pi,arma::mat& G,arma::mat& Q,
                                    arma::mat& alpha,int N1=128,int chain_length=20000,int burn_in=10000,
                                    double a_s0=1.0,double a_g0=1.0,double b_s0=1.0,double b_g0=1.0,double p1=0.5,double p2=0.5)
{
  int N=Y.n_rows;int J=Y.n_cols;int L=pow(2,K);int D=L-1;
  int T = chain_length-burn_in;
  arma::cube QQ(J,K,T);
  arma::cube GG(K,K,T);
  arma::cube AA(N,K,T);
  arma::mat SLIP(J,T);arma::mat GUESS(J,T);
  arma::mat PIS(L,T);
  arma::mat alpha_all(L,K);
  arma::mat Q_all(L-1,K);
  arma::vec classvec=arma::linspace<arma::vec>(0, L-1, L);
  alpha_all=Trans_10to2_mat(K,classvec);
  Q_all=alpha_all.rows(1,L-1);
  // prior
  arma::vec delta0=arma::ones(L,1);
  arma::vec pi0=1.0/L*arma::ones(L,1);
  
  arma::mat eta(L,D);
  int tmburn;
  for(int l=0;l<L;l++){
    for(int d=0;d<D;d++){
      eta(l,d)=arma::prod(vecpow((alpha_all.row(l)).t(),(Q_all.row(d)).t()));
    }
  }
  for(int t = 0; t < chain_length; t++){
    //update
    AHM_update(Y,alpha,s,g,Q,G,pi,p1,p2,alpha_all,Q_all,N,J,K,L,N1,a_s0,a_g0,b_s0,b_g0,delta0,eta);    
    if(t>burn_in-1){
      tmburn = t-burn_in;
      SLIP.col(tmburn)  = s;
      GUESS.col(tmburn) = g;
      AA.slice(tmburn)  = alpha;
      GG.slice(tmburn)  = G;
      QQ.slice(tmburn)  = Q;
      PIS.col(tmburn)   = pi;
      
    }
    int length=chain_length/10;
    if ((t+1)%length == 0) {
      // std::cout << t << std::endl;
      float percent = (float)(t+1) / (float)chain_length * 100; 
      Rcout << "\rComplete: "<< percent <<'%';
    }
    
// Rcout << "Iteration = "<< t <<'\r';
    
  }
  
  
  
  return Rcpp::List::create(Rcpp::Named("alpha") = AA,Rcpp::Named("s") = SLIP,Rcpp::Named("g") = GUESS,
                            Rcpp::Named("pi") = PIS,Rcpp::Named("Q") = QQ,Rcpp::Named("G") = GG);
}
