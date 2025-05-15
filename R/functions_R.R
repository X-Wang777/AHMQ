eta_onesubject_dina <- function(alpha,Q)
{
  J <- nrow(Q)
  K <- ncol(Q)
  alphamatrix <- matrix(rep(alpha,J),J,K,2)^Q
  eta <- apply(alphamatrix,1,prod)
  return(eta)
}
#' @export 
DINA_data <- function(N,J,K,s,g,Q,alpha)
{
  
  r <- array(NA,c(N,J))
  for(i in 1:N){
    eta <- eta_onesubject_dina(alpha[i,],Q)
    p <- (1-s)^eta*g^(1-eta)
    for(j in 1:J){
      r[i,j] <- stats::rbinom(1,1,p[j])
    }
  }
  r
}

#' @export 
Restricted_Q <- function(Q_all,R,J,K)
{
  Q_restrict  <- Q_all
  for(k in 1:K){
    pre_index <- which(R[,k]==1)
    index <- which(Q_restrict[,k]==1)
    Q_restrict[index,pre_index]=1
  }
  Q_restrict_code <- apply(Q_restrict,1,twoToten)
  index <- sort(unique(Q_restrict_code))
  r <- list(Q_restrict=Q_restrict,Q_restrict_code=Q_restrict_code,Q_restrict_index=index)
  r
}

samevec <- function(Q,q){
  J=dim(Q)[1]
  num <- 0
  for(j in 1:J){
    if(all(Q[j,]==q)){num <- num+1}
  }
  return(num)
}
#' @export 
Est_fun <- function(result,burn_in=NULL,chain_length=NULL,cut_value=NULL){
  Y =result[[1]]$data$Y
  N = dim(Y)[1]; J=dim(Y)[2];K=dim(result[[1]]$chain_sample$G)[1];L=2^K
  permutation=permn(c(1:K))
  permutation_length=length(permutation)
  alpha_all=Trans_10to2_mat(K,c(0:(L-1)))
  Q_all=alpha_all[-1,]
  
  if(is.null(burn_in)){
    chain_length <- length(result[[cl]]$chain_sample$s[1,])
    burn_in=0
  }
  sample_size=chain_length-burn_in
  min=c()
  ll <- c()
  DIC=c()
  for(cl in 1:length(result)){
    s_sample <- result[[cl]]$chain_sample$s[,(burn_in+1):chain_length]
    g_sample <- result[[cl]]$chain_sample$g[,(burn_in+1):chain_length]
    Q_sample <- result[[cl]]$chain_sample$Q[,,(burn_in+1):chain_length]
    G_sample <- result[[cl]]$chain_sample$G[,,(burn_in+1):chain_length]
    alpha_sample <- result[[cl]]$chain_sample$alpha[,,(burn_in+1):chain_length]
    pi_sample <- result[[cl]]$chain_sample$pi[,(burn_in+1):chain_length]
    DIC[cl] <- compute_DIC(Y,N,J, K, sample_size,s_sample,  g_sample,alpha_sample, Q_sample)
  }
  min_ll <- which(DIC == min(DIC))
  DIC_value <- DIC[min_ll]
  s <- apply(result[[min_ll]]$chain_sample$s[,(burn_in+1):chain_length],1,mean)
  g <- apply(result[[min_ll]]$chain_sample$g[,(burn_in+1):chain_length],1,mean)
  s_sample <- result[[min_ll]]$chain_sample$s[,(burn_in+1):chain_length]
  g_sample <- result[[min_ll]]$chain_sample$g[,(burn_in+1):chain_length]
  Q_sample <- result[[min_ll]]$chain_sample$Q[,,(burn_in+1):chain_length]
  G_sample <- result[[min_ll]]$chain_sample$G[,,(burn_in+1):chain_length]
  alpha_sample <- result[[min_ll]]$chain_sample$alpha[,,(burn_in+1):chain_length]
  pi_sample <- result[[min_ll]]$chain_sample$pi[,(burn_in+1):chain_length]
  ####################################
  Q_re_est=Q_sample
  alpha_re_est=alpha_sample
  G_re_est=G_sample
  pi_re_est <- pi_sample
  Q_est=array(0,c(J,K,permutation_length))
  Q_current=Q_sample;Q_pre=Q_sample
  G_current=G_sample;alpha_current=alpha_sample
  d <- c()
  e=1;kkk=0
  while(e>0){
    kkk=kkk+1
    QQ=apply(Q_current,c(1,2),mean)
    Q_pre=Q_current
    for(i in 1:sample_size){
      for(j in 1:permutation_length){
        tem=Q_current[,,i]
        ind <- permutation[[j]]
        for(k in 1:K){
          Q_est[,k,j]=Q_current[,ind[k],i]
        }
        d[j] <- sqrt(sum((Q_est[,,j]-QQ)^2))
      }
      min_ind=which(d==min(d))[1]
      Q_re_est[,,i]=Q_est[,,min_ind]
      Q_current[,,i]=Q_re_est[,,i]
      ind <- permutation[[min_ind]]
      alpha_all_re_est=alpha_all
      if(sum(abs(ind-c(1:K)))>0){
        G_est=G_current[,,i]
        alpha_est=alpha_current[,,i]
        for(k in 1:K){
          G_re_est[,k,i]=G_est[,ind[k]]
          alpha_re_est[,k,i]=alpha_est[,ind[k]]
          alpha_all_re_est[,k]=alpha_all[,ind[k]]
        }
        G_est=G_re_est[,,i]
        for(k in 1:K){
          G_re_est[k,,i]=G_est[ind[k],]
        }
        new_oder=Trans_2to10_mat(alpha_all_re_est,K)
        for(l in 1:L){
          iii <- which(new_oder==l-1)
          pi_re_est[l,i]=pi_sample[iii,i]
        }
      }
    }
    e=sum(abs(Q_pre-Q_current))
    cat("\rkkk=",kkk,";e=",e)
  }
  QQ=apply(Q_re_est,c(1,2),mean) 
  QQ1 <- matrix(0,J,K)
  QQ1[which(QQ>0.5)]=1
  Q <- QQ1
  GG=apply(G_re_est,c(1,2),mean);print(GG)
  if(is.null(cut_value)){
    cut_value <- readline((prompt="Please input a cut_value:"))
    cut_value=as.numeric(cut_value)
  }
  GG1 <- matrix(0,K,K)
  GG1[which(GG>cut_value)]=1
  G <- GG1
  AA <- apply(alpha_re_est,c(1,2),mean)
  AA1 <- matrix(0,N,K)
  AA1[which(AA>0.5)]=1
  alpha=AA1
  pi <- apply(pi_re_est,1,mean);pi=pi/sum(pi)
  G =Transitive(G ,K);G
  #################################
  Est_s <- s
  Est_g <- g
  Est_pi <- pi
  Est_G <- G
  Est_GG <- GG
  R=Reachability(G ,K)
  Est_Q  <- Restricted_Q(Q,R,J,K)$Q_restrict
  Est_alpha <- alpha
  ################
  Est <- list(Est_s=Est_s,Est_g=Est_g,Est_pi=Est_pi,Est_G=Est_G,Est_Q=Est_Q,Est_alpha=Est_alpha,Est_GG=Est_GG,DIC=DIC_value)
  return(Est)
}

