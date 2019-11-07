# generate matrix X with 3 topics
set.seed(123)
n = 120
p = 300
K= 3
L = matrix(0, nrow=n, ncol=K)
FF = matrix(0, nrow=K, ncol=p)

L[1:(n/3),1] = 1
L[((n/3)+1):(2*n/3),2] = 1
L[((2*n/3)+1):n,3] = 1


#x_seq = seq(0,1,length.out = p)
#FF[1,] = 10*sin(pi*x_seq)+1
#FF[2,] = c(rep(2,p/3),rep(6,p/6),rep(9,p/3),rep(2,p/6))
#FF[3,] = c(rep(6,100),rep(2,50),rep(2,100),rep(9,50))

FF[1,1:(p/3)] = 1+10
FF[2,((p/3)+1):(2*p/3)] = 1+10
FF[3,((2*p/3)+1):p] = 1+10

lambda = L %*% FF
X = matrix(rpois(n=length(lambda),lambda),nrow=n)
image(X)

fit_scd = NNLM::nnmf(X,K,method='lee',loss='mse')

for(i in 1:K){
  plot(fit_scd$W[,i],main=paste0('scd:estimate of loading',i))
}

for(i in 1:K){
  plot(fit_scd$H[i,],main=paste0('scd:estimate of factor',i))
}


library(mixsqp)
fit_ebmpmf = EBMPMF(X,K,init='lee')

for(i in 1:K){
  plot(fit_ebmpmf$ql$El[,i],main=paste0('ebmpmf:estimate of loading',i))
}

for(i in 1:K){
  plot(fit_ebmpmf$qf$Ef[i,],main=paste0('ebmpmf:estimate of factor',i))
}

fit_ebmpmf = EBMPMF(X,K,init='uniform')
for(i in 1:K){
  plot(fit_ebmpmf$ql$El[,i],main=paste0('ebmpmf:estimate of loading',i))
}

for(i in 1:K){
  plot(fit_ebmpmf$qf$Ef[i,],main=paste0('ebmpmf:estimate of factor',i))
}
#

set.seed(123)
n = 120
p = 300
K= 4
mfac = 2 # controls PVE of dense factor
L = matrix(0, nrow=n, ncol=K)
FF = matrix(0, nrow=K, ncol=p)

L[1:(n/3),1] = 1
L[((n/3)+1):(2*n/3),2] = 1
L[((2*n/3)+1):n,3] = 1
L[,4] = 1+mfac*runif(n)

FF[1,1:(p/3)] = 1+9
FF[2,((p/3)+1):(2*p/3)] = 1+10
FF[3,((2*p/3)+1):p] = 1+11
FF[4,]= 1+mfac*runif(p)

lambda = L %*% FF
X = matrix(rpois(n=length(lambda),lambda),nrow=n)
image(X)


fit_scd = NNLM::nnmf(X,K,method='lee',loss='mse')

for(i in 1:K){
  plot(fit_scd$W[,i],main=paste0('scd:estimate of loading',i))
}

for(i in 1:K){
  plot(fit_scd$H[i,],main=paste0('scd:estimate of factor',i))
}


fit_ebmpmf = EBMPMF(X,K,init='uniform',smooth_f = F,smooth_l = F)
for(i in 1:K){
  plot(fit_ebmpmf$ql$El[,i],main=paste0('ebmpmf:estimate of loading',i))
}

for(i in 1:K){
  plot(fit_ebmpmf$qf$Ef[i,],main=paste0('ebmpmf:estimate of factor',i))
}


##

set.seed(123)
n = 120
p = 300
K= 4
k=4
mfac = 2 # controls PVE of dense factor
L = matrix(0, nrow=n, ncol=k)
FF = matrix(0, nrow=p, ncol=k)
L[1:(n/3),1] = 1
L[((n/3)+1):(2*n/3),2] = 1
L[((2*n/3)+1):n,3] = 1
L[,4] = 1+mfac*runif(n)
FF[1:(p/3),1] = 1+10*runif(p/3)
FF[((p/3)+1):(2*p/3),2] = 1+10*runif(p/3)
FF[((2*p/3)+1):p,3] = 1+10*runif(p/3)
FF[,4]= 1+mfac*runif(p)
FF = t(FF)
lambda = L %*% FF
X = matrix(rpois(n=length(lambda),lambda),nrow=n)
image(X)


fit_ebmpmf = EBMPMF(X,K,init='scd',smooth_f = T,smooth_l = T,
                    scalel = 'estimate',nullweightl = 1e6,shapel = 0.001,point_massl = F,rounding=F)
summary_ebmpmf(fit_ebmpmf,X,lambda,K)

for(i in 1:K){
  plot(fit_ebmpmf$ql$El[,i],main=paste0('ebmpmf:estimate of loading',i))
}

for(i in 1:K){
  plot(fit_ebmpmf$qf$Ef[i,],main=paste0('ebmpmf:estimate of factor',i))
}

