EBMPMF = function(X,K,init = 'scd',niter=100,seed=123,tol = 1e-4,
                  shapef=NULL,point_massf=TRUE,nullweightf=100,
                  scalel='estimate',shapel=1,point_massl=FALSE,nullweightl=1000,
                  smooth_f=TRUE,smooth_l=FALSE,rounding=FALSE){
  set.seed(seed)
  #initialize q(Z)

  n = dim(X)[1]
  p = dim(X)[2]

  if(init == 'scd'){
    X_init_fit = NNLM::nnmf(X,K,method='scd',loss='mse',show.warning = F,verbose = F)
    L_init = X_init_fit$W
    F_init = X_init_fit$H
    ratio = mean(L_init)/mean(F_init)
    L_init = L_init/sqrt(ratio)
    F_init = F_init*sqrt(ratio)
  }
  if(init == 'lee'){
    X_init_fit = NNLM::nnmf(X,K,method='lee',loss='mse',show.warning = F,verbose = F)
    L_init = X_init_fit$W
    F_init = X_init_fit$H
    ratio = mean(L_init)/mean(F_init)
    L_init = L_init/sqrt(ratio)
    F_init = F_init*sqrt(ratio)
  }
  if(init == 'uniform'){
    L_init = matrix(runif(n*K),nrow=n,ncol=K)
    F_init = matrix(runif(K*p),nrow=K,ncol=p)
    ratio = mean(X)/(mean(L_init)*mean(F_init))
    L_init = L_init*sqrt(ratio)
    F_init = F_init*sqrt(ratio)
  }

  ql_hat = list(El = L_init, Elogl = log(L_init))
  qf_hat = list(Ef = F_init, Elogf = log(F_init))

  gl_hat = list()
  gR_hat = list()

  EZ = array(dim = c(n,p,K))

  EZ = Calc_EZ(X,K,EZ,ql_hat,qf_hat)

  KL = c()

  KL[1] = mKL(X,L_init%*%F_init)


  for(iter in 1:niter){

    if(rounding){
      EZ = round(EZ)
    }

    for(k in 1:K){


      # Update L

      EZk_rowsum = rowSums(EZ[,,k])
      Sum_Efk = sum(qf_hat$Ef[k,])
      if(smooth_l){
        lk_hat = update_smooth(EZk_rowsum, Sum_Efk, shape=shapef,point_mass=point_massf,nullweight=nullweightf)
        ql_hat$El[,k] = lk_hat$E
        ql_hat$Elogl[,k] = lk_hat$Elog
      }else{
        lk_hat = update_nsmooth(EZk_rowsum,Sum_Efk,scalel,point_massl,nullweightl,shapel)
        ql_hat$El[,k] = lk_hat$posterior$mean
        ql_hat$Elogl[,k] = lk_hat$posterior$mean_log
      }

      # Update R

      EZk_colsum = colSums(EZ[,,k])
      Sum_Elk = sum(ql_hat$El[,k])

      if(smooth_f){
        fk_hat = update_smooth(EZk_colsum, Sum_Elk, shape=shapef,point_mass=point_massf,nullweight=nullweightf)
        qf_hat$Ef[k,] = fk_hat$E
        qf_hat$Elogf[k,] = fk_hat$Elog
      }else{
        fk_hat = update_nsmooth(EZk_colsum,Sum_Elk,scalel,point_massl,nullweightl)
        qf_hat$Ef[k,] = fk_hat$posterior$mean
        qf_hat$Elogf[k,] = fk_hat$posterior$mean_log
      }

      #gR_hat[k] = fk_hat$pi_weight


      #gl_hat[k] = lk_hat$fitted_g

      # Update Z




    }

    EZ = Calc_EZ(X,K,EZ,ql_hat,qf_hat)

    KL[iter+1] = mKL(X,ql_hat$El%*%qf_hat$Ef)
    ########
    print(KL[iter+1])
    ########

    if(abs(KL[iter+1]-KL[iter])<=tol){
      break
    }

  }
  if(iter==niter){
    warning('Reached maximum iterations')
  }

  lambda_hat = ql_hat$El%*%qf_hat$Ef
  lambda_init = L_init%*%F_init

  return(list(qf=qf_hat,ql=ql_hat,gR=gR_hat,gl=gl_hat,KL=KL,Lambda=lambda_hat,
              L_init = L_init,F_init=F_init,Lambda_init = lambda_init,EZ=EZ))

}


Calc_EZ = function(X,K,EZ,ql_hat,qf_hat){
  n = nrow(X)
  p = ncol(X)
  for(k in 1:K){
    EZ[,,k] = outer(ql_hat$Elogl[,k], qf_hat$Elogf[k,], "+")
  }
  EZ = softmax3d(EZ)
  EZ = as.vector(EZ)*as.vector(X)
  dim(EZ) = c(n,p,K)
  EZ
}

softmax3d <- function(x){
  score.exp <- exp(x)
  probs <-as.vector(score.exp)/as.vector(rowSums(score.exp,dims=2))
  probs[is.na(probs)] = 0 ##  since softmax((0,0,...,0)) = (NA, NA,...,NA), I manually set them to be 0. But it is not safe!!!
  dim(probs) <- dim(x)
  return(probs)
}


mKL = function(A,B){
  D = A*log(A/B)-A+B
  D[is.nan(D)] = 0
  mean(D)
}

rmse = function(x,y){
  sqrt(mean((x-y)^2))
}

summary_ebmpmf = function(fit_obj,X,lambda,K,makeplot=T,typel='p',typef='l'){
  print(paste0('Mean KL:',round(mKL(lambda,fit_obj$Lambda),5)))
  #print(paste0("Sum of EZ:",sum(fit_obj$EZ)))
  if(makeplot){

    par(mfrow=c(1,2))
    image(lambda,main='True lambda')
    image(fit_obj$Lambda,main='Estimated lambda')
    nr = floor(sqrt(K))
    nc = ceiling(K/nr)

    par(mfrow=c(nr,nc))
    for(i in 1:K){
      plot(fit_obj$ql$El[,i],main=paste0('ebmpmf:estimate of loading',i),ylab = '',type=typel)
    }
    for(i in 1:K){
      plot(fit_obj$qf$Ef[i,],main=paste0('ebmpmf:estimate of factor',i),ylab = '',type=typef)
    }

  }
}

