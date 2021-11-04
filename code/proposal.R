source('code/smooth_flash.R')
source('code/wave_ebmf.R')
source('code/wave_ebmf_ndwt.R')

ploter = function(EF,main){
  par(mfrow=c(3,1))
  for(k in 1:ncol(EF)){
    plot(EF[,k],ylab=paste('f',k,sep=''),main=ifelse(k==1,main,""),type='l')
  }
}

ploter2 = function(EF,main){
  p = nrow(EF)
  par(mfrow=c(3,1))
  for(k in 1:ncol(EF)){
    plot(EF[1:(p/2),k],type='l',ylim = range(EF[,k]),ylab=paste('f',k,sep=''),xlab='',main=main)
    lines(EF[(p/2+1):(p),k],type='l',col=4)

  }
}


n = 120
p = 256
K= 3

FF = matrix(0, nrow=p, ncol=K)
f1 = 3
f2 = 1
FF[(p/8*1):(p/8*2),1] = f1
FF[(p/8*3):(p/8*4),2] = f2
FF[(p/8*5):(p/8*7),3] = f1

FF2 = matrix(0, nrow=p, ncol=K)
FF2[(p/8*1):(p/8*2),1] = f2
FF2[(p/8*3+10):(p/8*4-10),2] = f1
FF2[(p/8*5):(p/8*7),3] = f2

FFF = cbind(c(FF[,1],FF2[,3]),c(FF[,2],FF2[,2]),c(FF[,3],FF2[,1]))

par(mfrow=c(3,1))
for(k in 1:K){
  plot(FFF[1:p,k],type='l',ylim = range(FFF),ylab='',xlab='')
  lines(FFF[(p+1):(2*p),k],type='l',col=4)
}

par(mfrow=c(3,1))
plot(FFF[,1],type='l')
abline(v=p,lty=2)
plot(FFF[,2],type='l')
abline(v=p,lty=2)
plot(FFF[,3],type='l')
abline(v=p,lty=2)


l0 = 0
l1 = 3
L = cbind(c(rep(l1,n/3),rep(l0,n/3*2)),
          c(rep(l0,n/3),rep(l1,n/3),rep(l0,n/3)),
          c(rep(l0,n/3*2),rep(l1,n/3)))




s = 3
set.seed(12345)
y = tcrossprod(L,FFF) + matrix(rnorm(n*p*2,0,s),nrow=n,ncol=p*2)



library(tictoc)
tic();fit.flash = flash(y,var_type = 'by_row');toc()
tic();fit.dwt = wave_ebmf(y);toc()
#tic();fit.dwt2 = smooth_flash(y[,1:p],init_fn = 'udv_svd',type='station');toc()

ploter2(fit.flash$ldf$f,main='flash')
ploter2(fit.dwt$ldf$f,main='wave_flash')
#ploter(fit.dwt2$ldf$f,main='smooth_flash')

###########################
###########################
###### Poisson case #######
library(stm)

ploter2 = function(EF,main){
  p = nrow(EF)
  par(mfrow=c(3,1))
  for(k in 1:ncol(EF)){
    plot(EF[1:(p/2),k],type='l',ylim = range(EF[,k]),ylab=paste('f',k,sep=''),xlab='',main=main,col=2)
    lines(EF[(p/2+1):(p),k],type='l',col=3)

  }
}

K = 3
n=45
p=256
l0 = 0.1
l1 = 1
L = cbind(c(rep(l1,n/3),rep(l0,n/3*2)),
          c(rep(l0,n/3),rep(l1,n/3),rep(l0,n/3)),
          c(rep(l0,n/3*2),rep(l1,n/3)))

par(mfrow=c(3,1))
for(k in 1:K){
  plot(L[,k])
}

FF = matrix(0, nrow=p, ncol=K)
f1 = 3
f2 = 1
FF[(p/8*1):(p/8*2),1] = f1
FF[(p/8*3):(p/8*4),2] = f2
FF[(p/8*5):(p/8*7),3] = f1

FF2 = matrix(0, nrow=p, ncol=K)
FF2[(p/8*1):(p/8*2),1] = f2
FF2[(p/8*3+10):(p/8*4-10),2] = f1
FF2[(p/8*5):(p/8*7),3] = f2

FFF = cbind(c(FF[,1],FF2[,3]),c(FF[,2],FF2[,2]),c(FF[,3],FF2[,1]))


par(mfrow=c(3,1),tcl=-0.5,mai=c(0.3,0.6,0.1,0.3))
for(k in 1:K){
  plot(FFF[1:p,k],type='l',ylim = range(FFF),ylab=paste('f',k,sep=''),xlab='',lwd=2,col=2)
  lines(FFF[(p+1):(2*p),k],type='l',col=3,lwd=2)
}

set.seed(12345)
y = matrix(rpois(n*p*2,tcrossprod(L,FFF)),nrow=n)

set.seed(12345)
kmeans.init = kmeans(y, K, nstart = 20)
L0 = rep(1, n) %o% (as.vector(table(kmeans.init$cluster)))
F0 = t(kmeans.init$centers)
fit.stm5 = stm(y,K=3,nugget = F,init = list(L_init=L0,F_init=F0),tol=1e-5)
ploter2(fit.stm5$EF,main=NULL)

set.seed(12345)
fit.nmf = NNLM::nnmf(y,k=3,method = 'lee',loss='mkl',init = list(W=L0,H=t(F0)))
ploter2(t(fit.nmf$H),main=NULL)








