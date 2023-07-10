datax = read.csv('/project2/mstephens/cfbuenabadn/SRSF3.Counts.csv.gz')
rownames(datax)=datax[,1]
datax[1:5,1:5]
datax = datax[,-1]
dim(datax)
datax = as.matrix(datax)

library(smashr)
library(ebpmf)
library(Matrix)

# source('https://github.com/DongyueXie/ebpmf/blob/master/void/cluster_mix.R')
K = 5
fit_cluster_mix = cluster.mix(datax,K=5,tol=1e-3,maxit=50,smooth = T)
saveRDS(fit_cluster_mix,file='output/sgom_srsf3.rds')
#init = kmeans(datax,5,nstart = 5)
datax = Matrix(datax,sparse=T)
set.seed(12345)
for(i in 1:5){
  init = fastTopics::fit_poisson_nmf(datax,k=K,init.method = 'random')
  fit_ebpmf = ebpmf_identity(datax,K=5,tol = 1e-3,maxiter = 100,init = list(L_init = init$L,F_init = init$F))
  saveRDS(list(init=init,fit_ebpmf=fit_ebpmf),file=paste('output/ebpmf_srsf3_run',i,'.rds'))
}

par(mfrow=c(5,1))
for(k in 1:K){
  plot(init$F[,k],type='l')
}

par(mfrow=c(5,1))
for(k in 1:K){
  plot(fit_ebpmf$EF[,k],type='l')
}

par(mfrow=c(5,1))
for(k in 1:K){
  plot(fit_ebpmf$EF_smooth[,k],type='l')
}


for(k in 1:K){
  plot(fit_cluster_mix$phi[k,],type='l')
}

for(k in 1:K){
  plot(init$centers[k,],type='l')
}

library(funflash)
datax2 = datax[,colSums(datax)!=0]
Y_tilde = biwhitening(datax2)
fit_sf = scaledflash(Y_tilde$Y,Y_tilde$u,Y_tilde$v,
                     S2 = NULL,
                     var.type = 'by_column',
                     Kmax=30,
                     tol=0.01,
                     maxiter = 1000,
                     ebnm_fn = 'ebnm_pe',
                     init_fn = 'nnmf_r1',
                     ebnm_param=NULL,
                     verbose=TRUE,
                     nullcheck=TRUE,
                     sigma2 = NULL,
                     seed=12345)

fit_tt = fastTopics::fit_poisson_nmf(datax2,5)



