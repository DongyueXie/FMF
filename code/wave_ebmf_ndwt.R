##########################################
##########################################
##### wave ebmf ###########################
##### perform flash on DWTed matrix #######
##########################################
##########################################


## to do:

# 1. add other wavelet basis
# 2. add NDWT option

library(wavethresh)
devtools::load_all('~/Documents/Rpackages/flashr')
#'@param est_var 'mad': use mad estimator on the wavelet coefficients. 'va': update within variational approx
#'@param nlevel.prior number of levels that are assumed prior. If NULL, all levels
#'@param sum.prior whether put prior on the summation term, after DWT
wave_ebmf = function(Y,Kmax=100,tol=0.01,
                     s = NULL,
                     maxiter = 1000,
                     ebnm_fn = 'ebnm_pn',
                     init_fn = 'udv_si',
                     ebnm_param=NULL,
                     verbose=TRUE,
                     nullcheck=TRUE,
                     #est_var = 'va',
                     greedy=TRUE,
                     family="DaubExPhase",
                     filter.number=1,
                     type='wavelet',
                     #nlevels.prior = NULL,
                     #sum.prior = TRUE,
                     seed=12345){

  set.seed(seed)
  # check if the number of columns is a power of 2;
  # if not, extend it to be a power of 2.
  Y = power_of_2(Y)
  orig.idx = Y$idx
  Y = Y$Y

  n = nrow(Y)
  p = ncol(Y)
  nlevel = log2(p)

  # rows are wavelet coefficients. finest level comes first.
  Yw = get_dwt_wc(Y,filter.number,family,type)
  pw = ncol(Yw)
  # energy is the scaled summation of each row of Y
  energy = rowSums(Y)/sqrt(p)
  # initialize empty result object
  res = init_res(n,pw,Kmax)
  res$energy = energy

  for(k in 1:Kmax){
    res_old = res
    Rk = get_Rk(Yw,res$EL,res$EF)
    if(verbose){
      print(paste("Fitting dimension ", k))
    }
    res = wave_ebmf_rank1(Rk,res,k,
                          ebnm_fn=ebnm_fn,
                          ebnm_param=ebnm_param,
                          maxiter=maxiter,tol=tol,
                          verbose=verbose,
                          #nlevels.prior=nlevels.prior,
                          #sum.prior=sum.prior,
                          init_fn=init_fn,
                          type=type,
                          Yw=Yw,
                          nlevel=nlevel)
                          #est_var=est_var)


    if(nullcheck){
      if(verbose){
        print("Performing nullcheck")
      }
      res = null_check(res,k,verbose,Yw)
    }

    if(is_tiny_fl(res,k)){
      if(k>1){
        res=res_old
      }
      break
    }
    # check to stop

  }

  # transform back to data space
  # also deal with non-power 2: transform it back
  #browser()
  out = construct_object(res,orig.idx,p,filter.number,family,Yw,type)


  out
}

#'@description get the NDWT wavelet coefficients
#'@param Y input data matrix
get_dwt_wc = function(Y,filter.number,family,type){
  Yw = apply(Y,1,function(z){
    w.d = wavethresh::wd(z,filter.number = filter.number,family = family,type=type)
    w.d$D
  })
  return(t(Yw))
}

#'@description get the NDWT wavelet coefficients
#'@param Y input data matrix
inv_dwt_wc = function(EF,f_e,p,filter.number,family,type){
  K = ncol(EF)
  FF = matrix(nrow=p,ncol=K)
  if(type=='wavelet'){

    for(k in 1:K){
      temp0 = wavethresh::wd(rep(0,p),filter.number=filter.number,family=family,type='wavelet')
      temp0$D = EF[,k]
      temp0 = putC(temp0,0,f_e[k])
      FF[,k] = wr(temp0)
    }

  }else if(type == 'station'){

    for(k in 1:K){
      temp0 = wavethresh::wd(rep(0,p),filter.number=filter.number,family=family,type='station')
      temp0$D = EF[,k]
      temp0 = putC(temp0,0,rep(f_e[k],p))
      FF[,k] = AvBasis(convert(temp0))
    }

  }else{
    stop('type should be either wavelet or station')
  }
  return(FF)
}


#'@description evaluate objective function
calc_objective = function(Yw,res){
  obj = Eloglik(Yw,matrix((res$s)^2,nrow=length(res$s),ncol=ncol(Yw)),res$EL,res$EF,res$EL2,res$EF2) + sum(res$KL_L) + sum(res$KL_F)
  obj
}

delete_factor = function(res,k,Yw){
  res$EL[,k] = 0
  res$EF[,k] = 0
  res$EL2[,k] = 0
  res$EF2[,k] = 0
  res$gL[[k]] = list(NULL)
  res$gF[[k]] = list(NULL)
  res$KL_L[k] = 0
  res$KL_F[k] = 0
  res$s = sqrt(rowMeans(get_R2(Yw,res$EL,res$EF,res$EL2,res$EF2)))
  res
}

null_check = function(res,k,verbose,Yw){
  # delete the kth factor
  res0 = delete_factor(res,k,Yw)
  # compare the objective function
  obj0 = calc_objective(Yw,res0)
  obj1 = calc_objective(Yw,res)

  if(obj0>obj1){
    res = res0
  }
  if(verbose){
    if(obj1>obj0){
      print(paste('Deleting factor ',k, ' decreases objective by ', round(obj1-obj0,3)))
    }else{
      print(paste('Deleting factor ',k, ' increases objective by ', round(obj0-obj1,3)))
    }

  }

  res
}

#'@description construct object for returning results.
construct_object = function(res,orig.idx,p,filter.number,family,Yw,type){
  # remove factor and loadings are zero
  rm.idx = which(colSums(res$EL2)==0)
  if(length(rm.idx)>0){
    res$EL = res$EL[,-rm.idx,drop=FALSE]
    res$EF = res$EF[,-rm.idx,drop=FALSE]
    res$EL2 = res$EL2[,-rm.idx,drop=FALSE]
    res$EF2 = res$EF2[,-rm.idx,drop=FALSE]
    res$KL_L = res$KL_L[-rm.idx]
    res$KL_F = res$KL_F[-rm.idx]
  }


  # take the energy into account
  # here E(energy|EL) = EL*f_e; var(energy|EL) = diag(s^2)
  # perform a weighted least square to estimate f_e
  #browser()
  f_e = lm(y~.+0,data.frame(y=res$energy,X = res$EL),weights = 1/(res$s)^2)
  f_e = f_e$coefficients
  # transform F back to data space
  FF = inv_dwt_wc(res$EF,f_e,p,filter.number,family,type)
  if(!is.null(orig.idx)){
    FF = FF[orig.idx,]
  }
  d = sqrt(colSums(res$EL^2) * colSums(FF^2))
  fitted_values = tcrossprod(res$EL,FF)

  FF = scale(FF,scale = sqrt(colSums(FF^2)),center=FALSE)
  L = scale(res$EL,scale = sqrt(colSums(res$EL^2)),center=FALSE)

  ldf = list(d=d,l=L,f=FF)

  nfactors = length(d)
  pve = d^2/(sum(d^2)+sum((res$s)^2)*nrow(FF))

  objective = calc_objective(Yw,res)



  return(list(ldf = ldf,
              nfactors = nfactors,
              pve = pve,
              fitted_values = fitted_values,
              fit = res,
              objective = objective))

}

#'@description get the residual for model fitting, kth factor and loadings
get_Rk = function(Y,EL,EF){
  if(is.null(EL)&is.null(EF)){
    return(Y)
  }else{
    Y - tcrossprod(EL,EF)
  }
}

#'@description initialize the output list
init_res = function(n,p,Kmax){
  return(list(EL = matrix(0,nrow=n,ncol=Kmax),
              EF = matrix(0,nrow=p,ncol=Kmax),
              EL2 = matrix(0,nrow=n,ncol=Kmax),
              EF2 = matrix(0,nrow=p,ncol=Kmax),
              gL = list(),
              gF = list(),
              KL_L = rep(0,Kmax),
              KL_F = rep(0,Kmax),
              s = NULL))
}

#'@description  update the output list
update_res = function(res,k,El,El2,Ef,Ef2,
                      gl,gf,
                      KL.l,KL.f,s=NULL){
  res$EL[,k] = El
  res$EL2[,k] = El2
  res$EF[,k] = Ef
  res$EF2[,k] = Ef2
  res$gL[[k]] = gl
  res$gF[[k]] = gf
  res$KL_L[k] = KL.l
  res$KL_F[k] = KL.f
  if(!is.null(s)){
    res$s = s
  }
  res
}

init_fl = function(Y,init_fn){
  out = do.call(init_fn,list(Y,1))
  out
}

#'@description  fit a rank 1 wave-ebmf model
#'@param Y data matrix
#'@param res current fit of the wave-ebmf model
#'@param k current k to be added
#'@return res
wave_ebmf_rank1 = function(Y,res,k,ebnm_fn,ebnm_param,
                           maxiter=100,tol=0.01,
                           verbose,init_fn,type,Yw,nlevel){

  #s = res$s
  # initialize l and f
  init = init_fl(Y,init_fn)
  El = init$u*sqrt(init$d[1])
  Ef = init$v*sqrt(init$d[1])
  El2 = El^2
  Ef2 = Ef^2
  gl = list()
  gf = list()
  KL.l = 0
  KL.f = 0
  res = update_res(res,k,El,El2,Ef,Ef2,gl,gf,KL.l,KL.l)

  obj = -Inf
  for(i in 1:maxiter){

    # update s

    s = c(sqrt(rowMeans(get_R2(Yw,res$EL,res$EF,res$EL2,res$EF2))))

    # update l

    ## formulate the x and s for ebnm

    #xl = rowSums(rep(1,n)%*%t(Ef) * Y)/sum(Ef2)
    xl = colSums(t(Y)*c(Ef))/sum(Ef2)
    sl = sqrt(s^2/sum(Ef2))
    a = do.call(ebnm_fn, list(xl, sl, ebnm_param))

    El = a$postmean
    El2 = a$postmean2
    gl = a$fitted_g

    KL.l = a$penloglik - NM_posterior_e_loglik(xl, sl,El, El2)

    if(sum(El2)==0){
      res = update_res(res,k,El,El2,Ef,Ef2,gl,gf,KL.l,KL.f,s)
      print('loading zeroed out')
      break
    }

    # update f

    # update wavelet coefficients
    KL.f = 0
    for(t in (0):(nlevel-1)){
      idx = get_wave_index(t,nlevel,type)
      yt = Y[,idx,drop=FALSE]
      #xf = colSums(El%*%t(rep(1,2^t))*yt/s^2)/sum(El2/s^2)
      xf = colSums(yt*c(El)/s^2)/sum(El2/s^2)
      sf = sum(El2/s^2)^(-0.5)
      a = do.call(ebnm_fn, list(xf, sf, ebnm_param))
      Ef[idx] = a$postmean
      Ef2[idx] = a$postmean2
      gf[[t+1]] = a$fitted_g
      KL.f = KL.f + a$penloglik - NM_posterior_e_loglik(xf, sf,a$postmean, a$postmean2)
    }

    if(sum(Ef2)==0){
      res = update_res(res,k,El,El2,Ef,Ef2,gl,gf,KL.l,KL.f,s)
      print('factor zeroed out')
      break
    }

    res = update_res(res,k,El,El2,Ef,Ef2,gl,gf,KL.l,KL.f,s)
    # evaluate objective function
    obj[i+1] = calc_objective(Yw,res)
    if(verbose){
      print(paste('Iteration ', i, ': obj ', round(obj[i+1],3)))
    }

    if((obj[i+1]-obj[i])<0){
      print('An iteration decreased the objective')
    }

    if((obj[i+1]-obj[i])<=tol){
      break
    }


  }

  return(res)


}


# get_R2 = function(Y,res){
#   (Y-tcrossprod(res$EL,res$EF))^2+tcrossprod(res$EL2,res$EF2) - tcrossprod(res$EL^2,res$EF^2)
# }

get_R2 = function(Y,EL,EF,EL2,EF2){
  (Y-tcrossprod(EL,EF))^2+tcrossprod(EL2,EF2) - tcrossprod(EL^2,EF^2)
}

#'@description This function calculates the expected log likelihood wrt q(L,F).
#'@param Y the data matrix
#'@param Sigma a matrix of variance of random errors. [sigma^2 sigma^2 ... sigma^2]
Eloglik = function(Y,Sigma,EL,EF,EL2,EF2){
  -0.5 * sum(log(2 * pi * Sigma) + 1/Sigma * get_R2(Y,EL,EF,EL2,EF2))
}

#'@param t current wavelet level, take value 0:(nt-1)
#'@param nt total number of level
get_wave_index = function(t,nt,type){
  # find level index before the current level t
  levels = 0:(nt-1)
  if(t<nt){

    if(type=='wavelet'){
      if(t==levels[nt]){
        return(1:2^t)
      }else{
        idx = levels[levels>t]
        return((sum(2^idx)+1):(sum(2^idx)+2^t))
      }
    }
    if(type=='station'){
      p = 2^nt
      temp = p*(nt-t-1)
      (temp+1):(temp+p)
    }

  }else{
    stop('wavelet levels are from 0 to (#levels-1)')
  }





}



haar = function(x,scale= sqrt(2)){
  if(length(x)==1){
    return(x)
  }
  else{
    x = matrix(x,nrow=2)
    diff = (x[1,]-x[2,])/scale
    sum = (x[1,]+x[2,])/scale
    return(c(diff, haar(sum)))
  }
}

haar_inv = function(x,scale=sqrt(2)){
  n=length(x)
  if(n==1){
    return(x)
  }
  x = matrix(scale*x,nrow=2,byrow=TRUE)
  smoothed = haar_inv(x[2,])
  return(as.vector(rbind(smoothed+x[1,], smoothed-x[1,]))/2)
}

#'@description   If a sequence is not a power of two, make it a power of two.
power_of_2 = function(Y){
  n = ncol(Y)
  J = log2(n)
  if((J%%1) == 0){
    return(list(Y=Y, idx = NULL))
  }else{
    n.ext = 2^ceiling(J)
    lnum = round((n.ext - n)/2)
    rnum = n.ext - n - lnum
    if (lnum == 0) {
      Y.lmir = NULL
    } else {
      Y.lmir = Y[,lnum:1]
    }
    if (rnum == 0) {
      Y.rmir = NULL
    } else {
      Y.rmir = Y[,n:(n - rnum + 1)]
    }
    Y.ini = cbind(Y.lmir, Y, Y.rmir)
    return(list(Y=Y.ini,idx = (lnum + 1):(lnum + n)))
  }
}

#' @title Reflect and extend a vector.
#'
#' @description Extends the vector to have length a power of 2 (if not already a power of 2) and then
#'   reflects it about its right end.
#'
#' @details The vector x is first reflected about both its left and right ends, by (roughly) the same
#' amount each end, to make its length a power of 2 (if the length of x is already a power of 2 this step is skipped).
#' Then the resulting vector is reflected about its right end to create a vector
#' that is both circular (left and right ends are the same) and a power of 2.
#'
#' @param x An n-vector.
#'
#' @return A list with two list elements: \code{"x"} containing the
#'   extended-and-reflected signal; and \code{"idx"} containing the indices of the
#'   original signal.
reflect <- function (x) {
  n = length(x)
  J = log2(n)
  if ((J%%1) == 0) {

    # if J is an integer, i.e. n is a power of 2.
    x = c(x, x[n:1])
    return(list(x=x, idx = 1:n))
  } else {
    n.ext = 2^ceiling(J)
    lnum = round((n.ext - n)/2)
    rnum = n.ext - n - lnum
    if (lnum == 0) {
      x.lmir = NULL
    } else {
      x.lmir = x[lnum:1]
    }
    if (rnum == 0) {
      x.rmir = NULL
    } else {
      x.rmir = x[n:(n - rnum + 1)]
    }
    x.ini = c(x.lmir, x, x.rmir)
    x.mir = x.ini[n.ext:1]
    x = c(x.ini, x.mir)
    return(list(x = x, idx = (lnum + 1):(lnum + n)))
  }
}
