AHMQ_single <- function(Y,K,N1=128,chain_length=20000,burn_in=10000,a_s0=1.0,a_g0=1.0,b_s0=1.0,b_g0=1.0,p1=0.5,p2=0.5){
  #############
  N <- dim(Y)[1];J=dim(Y)[2];L=2^K
  #initial value
  s <- runif(J,0,.4)
  g <-  runif(J,0,.4)
  pi <- runif(L);pi=pi/sum(pi)
  Q <- random_Q(J,K)
  G <- matrix(0,K,K)
  alpha <- matrix(rbinom(N*K,1,0.5),N,K)
  #############################
  r <- AHMQ_MCMC(Y,K,s,g,pi,G,Q,alpha,N1,chain_length,burn_in,a_s0,a_g0,b_s0,b_g0,p1,p2)
  
  r
}
#' @export 
AHMQ <- function(Y,K,N1=128,chain_length=20000,burn_in=10000,chain_num=4,a_s0=1.0,a_g0=1.0,b_s0=1.0,b_g0=1.0,p1=0.5,p2=0.5){
  result <- vector("list",chain_num)
  for(i in 1:chain_num){
    r=AHMQ_single(Y=Y,K=K,N1=N1,chain_length = chain_length,burn_in = burn_in)
    G=r$G;s=r$s;g=r$g;Q=r$Q;pi=r$pi;alpha=r$alpha
    data <-list(Y=Y,K=K,N1=N1,chain_length=chain_length,burn_in=burn_in,chain_num=chain_num,
                a_s0=a_s0,a_g0=a_g0,b_s0=b_s0,b_g0=b_g0,p1=p1,p2=p2)
    chain_sample <- list(G=G,s=s,g=g,Q=Q,pi=pi,alpha=alpha)
    result[[i]] <- list(data=data,chain_sample=chain_sample)
    cat("chain ", i, "\n")
  }
  result
}