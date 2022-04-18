# generate some data
#'@description generate binomial nuggest data
gen_data = function(ni,alpha,sigma2){
  lp = alpha + rnorm(length(alpha),0,sqrt(sigma2))
  p = 1/(1+exp(-lp))
  rbinom(length(alpha),ni,p)
}


#'@description Normal approximation of log odds
logodds_normal = function(ns,nf){
  n = length(ns)
  alpha_hat = c()
  for(i in 1:n){
    if(ns[i]==0){
      alpha_hat[i] = log((ns[i]+1/2)/(nf[i]+1/2))-1/2
    }else if(nf[i]==0){
      alpha_hat[i] = log((ns[i]+1/2)/(nf[i]+1/2))+1/2
    }else{
      alpha_hat[i] = log((ns[i])/(nf[i]))
    }
  }
  N = ns+nf
  V3 = (ns+nf+1)/(ns+nf)*(1/(ns+1)+1/(nf+1))
  Vs = V3*(1-2/N+V3/2)
  alpha_hat_se = sqrt(Vs-V3^2/2*(V3-4/N))
  return(list(alpha_hat=alpha_hat,alpha_hat_se=alpha_hat_se))
}

library(ashr)

#'@description Fit normal mean nuggest using variational inference
#'@param x,s data and known standard errors
#'@param sigma2_init init value of nugget effect
#'@param fix_sigma2 whether fix sigma2 at its init value. For testing purpose.
nm_nug = function(x,s,
                  max_iter = 1000,
                  tol=1e-2,
                  sigma2_init = NULL,
                  fix_sigma2 = FALSE,
                  printevery = 10){
  n = length(x)

  sigma2_func = function(sigma2,x,s){
    return(sum(x^2/(sigma2+s^2)^2)-sum(1/(sigma2+s^2)))
  }

  if(is.null(sigma2_init)){
    sigma2 = try(uniroot(sigma2_func,interval=c(0,0.5),x=x,s=s)$root,silent=FALSE)
    if(class(sigma2)[1]=='try-error'){
      sigma2 = var(x)
    }
  }else{
    sigma2 = sigma2_init
  }

  sigma2_trace = c(sigma2)

  grid_min = min(s)/10
  grid_max = 2*sqrt(max(x^2-s^2))
  ash_grid = exp(seq(log(grid_min),log(grid_max),by=log(sqrt(2))))
  mu_bar = ash(x,sqrt(sigma2+s^2),method='fdr',mixsd = ash_grid,mixcompdist = 'normal')$result$PosteriorMean

  elbo = c(-Inf)

  for(i in 1:max_iter){
    eps_bar = (x-mu_bar)/(1+s^2/sigma2)
    tau2 = 1/(1/s^2+1/sigma2)
    if(!fix_sigma2){
      sigma2 = mean(tau2+eps_bar^2)
      sigma2_trace[i+1] = sigma2
    }
    res = ashr::ash(x-eps_bar,s,method='shrink',mixsd = ash_grid,mixcompdist = 'normal')
    mu_bar = res$result$PosteriorMean
    mu2_bar = (res$result$PosteriorMean)^2 + (res$result$PosteriorSD)^2
    # calc elbo
    Eloglik = sum(-1/2/s^2*(-2*x*mu_bar+mu2_bar+eps_bar^2+tau2-2*x*eps_bar+2*mu_bar*eps_bar))
    Eloggq_mu = res$loglik + 1/2*sum(1/s^2*((x-eps_bar)^2-2*(x-eps_bar)*mu_bar+mu2_bar))
    Eloggq_eps = -n/2*log(sigma2) - sum(tau2+eps_bar^2)/2/sigma2 + 1/2*sum(log(tau2))
    elbo[i+1] = Eloglik+Eloggq_mu+Eloggq_eps
    if((elbo[i+1]-elbo[i]) <=tol){
      print('Converged')
      break
    }
    if(i%%printevery==0){
      print(sprintf("Done %i iterations, objective = %f", i, elbo[i+1]))
    }
  }
  return(list(mu=res$result$PosteriorMean,mu_var = mu2_bar-mu_bar^2,
              sigma2=sigma2,elbo=elbo,sigma2_trace=sigma2_trace))
}

#'@description Fit x_i = mu_i + e_i + epsilon_i, e_i~N(0,v), epsilon_i~N(0,s_i^2) using EM algorithm
#'@param x,s data and known standard errors
#'@param grid fixed grid of variances
#'@param gridmult default sqrt(2)
#'@param w init of weights
#'@param v init of variance
K_group_model = function(x,s,grid=NULL,gridmult = sqrt(2),maxiter = 1000,w = NULL,v=NULL,printevery = 10,tol=1e-5){
  n = length(x)
  if(is.null(s)){
    s = rep(0,n)
  }
  if(is.null(grid)){
    grid_min = min(s)/10
    grid_max = 2*sqrt(max(x^2-s^2))
    grid = exp(seq(log(grid_min),log(grid_max),by=log(gridmult)))
    grid = c(0,grid^2)
  }
  K = length(grid)
  if(is.null(w)){
    w = rep(1/K,K)
  }
  if(is.null(v)){
    v = var(x) - mean(s^2)
  }
  # a is a matrix
  gradient_v = function(v,x,s,grid,a){
    n = length(s)
    K = length(grid)
    M = 1/(v+rep(1,n)%*%t(grid)+s^2%*%t(rep(1,K)))
    sum(a*(M-x^2%*%t(rep(1,K))/M^2))
  }
  compute_L = function(x,v,grid,s){
    L = matrix(nrow=n,ncol=K)
    for(k in 1:K){
      L[,k] = dnorm(x,0,sqrt(v+grid[k]+s^2))
    }
    L
  }

  iter = 1
  obj = -Inf
  L = compute_L(x,v,grid,s)
  while (iter<=maxiter) {
    if(iter%%printevery==0){
      print(paste('Done iter',iter))
    }
    # E step
    A = L%*%diag(w)
    A = A/rowSums(A)

    # M step
    w = colMeans(A)
    v = try(uniroot(gradient_v,interval=c(0,5),x=x,s=s,grid=grid,a=A,extendInt = 'no'),silent = T)
    if(class(v)=='try-error'){
      v = 0
    }else{
      v = v$root
    }

    # check convergence
    L = compute_L(x,v,grid,s)
    obj[iter+1] = sum(log(L%*%w))
    if((obj[iter+1]-obj[iter])<tol){
      break
    }
    iter  = iter + 1
  }
  return(list(w=w,v=v,A=A,obj=obj))
}

##########################################
########## test K group model ##########
##########################################
pi0 = 0.8
pi1 = 0.2
n = 1000
sigma2_0 = 1
sigma2_1 = 3
s = abs(rnorm(n,0,1))
#s = 0
x = c(rnorm(pi0*n,0,sqrt(sigma2_0+s^2)),rnorm(pi1*n,0,sqrt(sigma2_1+s^2)))
fit = K_group_model(x,s,maxiter = 1000)
##########################################
##########################################
##########################################

two_group_model = function(x,s,tol = 1e-8, maxiter = 10000,v1 = 1, v2 = 2, w = c(0.5,0.5),printevery=10){

  gradient_v = function(v,x,s,a){
    sum(a/(v+s^2) - a*x^2/(v+s^2)^2)
  }

  if(is.null(s)){
    s = rep(0,n)
    gmm = TRUE
  }else{
    gmm = FALSE
  }

  # EM for fitting x_i ~ pi_0 N(0,sigma^2+s_i^2) + pi_1 N(0, eta^2 + s_i^2), s_i known.
  n = length(x)
  iter = 1
  obj = -Inf
  while(iter<=maxiter){
    if(iter%%printevery==0){
      print(paste('Done iter',iter))
    }
    # E step
    a1 = w[1]*dnorm(x,0,sqrt(s^2+v1))
    a2 = w[2]*dnorm(x,0,sqrt(s^2+v2))
    a1 = a1/(a1+a2)
    a2 = 1-a1
    # M step
    w[1] = mean(a1)
    w[2] = 1 - w[1]
    if(gmm){
      v1 = sum(a1*x^2)/sum(a1)
      v2 = sum(a2*x^2)/sum(a2)
    }else{
      v1 = uniroot(gradient_v,c(0,5),x=x,a=a1,s=s,extendInt = 'yes')$root
      v2 = uniroot(gradient_v,c(0,5),x=x,a=a2,s=s,extendInt = 'yes')$root
    }

    obj[iter+1] = sum(log(w[1]*dnorm(x,0,sqrt(v1+s^2))+w[2]*dnorm(x,0,sqrt(v2+s^2))))
    #obj[iter+1] = sum(a1*(log(w[1])+dnorm(x,0,sqrt(v1+s^2),log = T))) + sum(a2*(log(w[2])+dnorm(x,0,sqrt(v2+s^2),log = T)))
    if((obj[iter+1]-obj[iter])<tol){
      break
    }
    iter = iter + 1
  }
  return(list(w=w,v1=v1,v2=v2,a1=a1,a2=a2,obj=obj))
}

##########################################
########## test two group model ##########
##########################################
pi0 = 0.3
pi1 = 1 - pi0
n = 5000
sigma2_0 = 1
sigma2_1 = 3
s = runif(n,0.1,0.5)
#s = 0
x = c(rnorm(pi0*n,0,sqrt(sigma2_0+s^2)),rnorm(pi1*n,0,sqrt(sigma2_1+s^2)))
fit = two_group_model(x,s=NULL,maxiter = 1000)

fit2 = Mclust(x,G = 2,model = 'V')

fit3 = normalmixEM(x,k=2,mean.constr = c(0,0))

fit4 = ash(x,s,pointmass=FALSE,prior='uniform',mixcompdist = 'normal')


x = rnorm(1000)
library(mixtools)
fit = normalmixEM(x,k=2,mean.constr = c(0,0),maxit = 1e5,sd.constr = c('a','a'))
fit = Mclust(x,G=1:2)

library(mr.ash.alpha)
pi0 = 0.1
pi1 = 1 - pi0
n = 500
sigma2_0 = 1
sigma2_1 = 2
x = c(rnorm(pi0*n,0,sqrt(sigma2_0)),rnorm(pi1*n,0,sqrt(sigma2_1)))
fit = mr.ash(diag(n),x,intercept = FALSE,max.iter = 5000,method_q = 'sigma_indep_q')

##################### test Jason's optimize function ###############

optimize.noisy <- function(R2, S2, wts = 1) {
  interval.max <- max((R2 - S2) / wts)
  if (interval.max <= 0)
    return(0)
  opt.res <- optimize(function(v) {
    sum(log(wts * v + S2)) + sum(R2 / (wts * v + S2))
  }, interval = c(0, interval.max), tol = sqrt(.Machine$double.eps))
  return(opt.res$minimum)
}

optimize.noisy.root <- function(R2, S2) {
  interval.max <- max((R2 - S2))
  if (interval.max <= 0)
    return(0)
  opt.res <- uniroot(function(v) {
    sum(1/(v + S2)) - sum(R2 / (v + S2)^2)
  }, interval = c(0, interval.max), tol = sqrt(.Machine$double.eps))
  return(opt.res$root)
}

v = 10
s = runif(100,0.1,0.5)
x = rnorm(100,0,sqrt(s^2+v))
optimize.noisy(x^2,s^2)
optimize.noisy.root(x^2,s^2)



#####################################################################

# set.seed(12345)
# nn = 100
# datax = gen_data(nn,c(rep(0,90),rep(1,10)),0.1)
# est = logodds_normal(datax,nn-datax)
# res = nm_nug(est$alpha_hat,est$alpha_hat_se,max_iter = 500,sigma2_init=0.5,fix_sigma2 = F)




# res$sigma2
#
# plot(est$alpha_hat)
# lines(res$mu)
# lines(ash(est$alpha_hat,est$alpha_hat_se)$result$PosteriorMean,col=4)
# plot(res$elbo,type='l')
#
# plot(res$sigma2_trace,type='l')
#
#
#
#
# set.seed(12345)
# n = 200
# mu = c(rep(-1,n-10),rep(1,10))
# s = runif(n,0.01,0.2)
# sigma2 = 0.05
# x = rnorm(n,mu,sqrt(s^2+sigma2))
#
# ash.fit = ash(x,sqrt(s^2+sigma2))
# plot(x)
# lines(ash.fit$result$PosteriorMean)
#
# res = nm_nug(x,s)
# res$sigma2
#
# plot(x)
# lines(res$mu)
#
# mu = c(1.5,-1.5)
# s = runif(2,0.01,0.2)
# sigma2 = 0.05
# x = rnorm(2,mu,sqrt(s^2+sigma2))
# res = nm_nug(x,s)
# res$sigma2




