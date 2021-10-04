
source('code/smash_wave.R')

#
#
# sparseWS2 = function(y,sigma,sigma0=sqrt(0.01),
#                      maxiter=1000,tol = 0.01,fix_pi1 = TRUE,pi1 = 0.1,
#                     filter.number=10,family="DaubLeAsymm"){
#   # init
#   n = length(y)
#   #fit0 = smash.gaus(y,sigma=sigma,filter.number=filter.number,family=family,post.var = TRUE)
#   fit0 = smash_wave(y,sigma=sigma,filter.number=filter.number,family=family)
#   m = fit0$mu.est
#   s2 = fit0$mu.est.var+m^2
#   # init phi
#   phi = exp(abs(m))/(1+exp(abs(m)))
#   phi = phi*2-1
#   #m = y
#   #s2 = y^2
#
#   sigma2 = sigma^2
#   sigma02 = sigma0^2
#   elbo = -Inf
#   for(i in 1:maxiter){
#     print(i)
#     phi.old = phi
#
#     #update beta
#     beta1 = sigma02/(1+sigma02)*y + 1/(1+sigma02)*phi*m
#
#     # update phi
#     phi = pi1*exp((beta1*m-1/2*s2)/(sigma2*sigma02))
#     phi = phi/(phi+(1-pi1))
#     phi[is.nan(phi)] = 1
#     # update pi1
#     if(!fix_pi1){
#       pi1 = sum(phi)/n
#     }
#
#     # update m, s2
#     #fit = smash_wave(beta1,sigma = sqrt(sigma2*sigma02/pmax(phi,1e-2)),filter.number=filter.number,family=family)
#     fit = smash.gaus(beta1,sigma = sqrt(sigma2*sigma02/pmax(phi,1e-2)),filter.number=filter.number,family=family,post.var = TRUE)
#     m = fit$mu.est
#     s2 = fit$mu.est.var+m^2
#     # calc elbo
#     #Eloggq = calc_Eloggq(fit$loglik,y,sigma,phi,m,s2)
#     #elbo[i+1] = calc_elbo(y,sigma,m,s2,phi,pi1,Eloggq)
#     # if(abs(elbo[i+1]-elbo[i])<=tol){
#     #   break
#     # }
#
#     if(sqrt(crossprod(phi.old-phi))<=tol){
#       break
#     }
#
#   }
#
#   return(list(pm=beta1,phi=phi,pi1=pi1,m=m,s2=s2))
# }
#


calc_elbo = function(y,sigma,sigma0,log_lik,phi,w,Ebeta,Ebeta2,Emu,Emu2,K){

  sigma2=sigma^2
  sigma02=sigma0^2
  elbo=0
  Eloggq = 0
  for(k in 1:K){

    elbo = elbo +
      sum(phi[,k]*(-1/2*log(sigma2)-1/2/sigma2*(Ebeta2[,k]-2*Ebeta[,k]*y)+log(w[k])
                   -1/2*log(sigma2*sigma02)-1/2/sigma2/sigma02*(Ebeta2[,k]-2*Ebeta[,k]*Emu[,k]+Emu2[,k])-log(phi[,k])))
    Eloggq = Eloggq + calc_Eloggq(log_lik[k],Ebeta[,k],sigma2*sigma02/phi[,k],Emu[,k],Emu2[,k])
  }
  elbo+Eloggq
}

calc_Eloggq = function(lg,y,v,Emu,Emu2){
  lg+1/2*sum(log(v))+1/2*sum(1/v*(Emu2-2*y*Emu+y^2))
}

#'@param smoother either wavelet or gp, gp stands for gaussian process.
#'@param sparser either ash or delta_0

sparse_smooth = function(x,y,sigma,sigma0=sqrt(0.01),
                         maxiter=1000,tol = 0.01,fix_w = F,w = c(0.5,0.5),
                         smoother = 'gp',
                         sparser = 'ash',
                         filter.number=1,family="DaubExPhase",
                         printevery = 10){
  # init
  n = length(y)
  #fit0 = smash.gaus(y,sigma=sigma,filter.number=filter.number,family=family,post.var = TRUE)


  # init phi

  #phi2 = exp(abs(m))/(1+exp(abs(m)))
  #phi2 = phi*2-1
  #m = y
  #s2 = y^2

  sigma2 = sigma^2
  sigma02 = sigma0^2
  elbo = -Inf
  K = 2
  Ebeta = matrix(nrow=n,ncol=K)
  Ebeta2 = matrix(nrow=n,ncol=K)
  phi = matrix(nrow=n,ncol=K)
  Emu = matrix(nrow=n,ncol=K)
  Emu2 = matrix(nrow=n,ncol=K)
  log_lik = c()

  # if(is.null(w)){
  #   fit.temp = ash(y,sigma)
  #   pi1 = fit.temp$fitted_g$pi[1]
  #   w = c(pi1,1-pi1)
  #   Emu[,1] = fit.temp$result$PosteriorMean
  #   Emu2[,1] = fit.temp$result$PosteriorMean^2+(fit.temp$result$PosteriorSD)^2
  # }else{
  #   fit.temp = ash(y,sigma)
  #   Emu[,1] = fit.temp$result$PosteriorMean
  #   Emu2[,1] = fit.temp$result$PosteriorMean^2+(fit.temp$result$PosteriorSD)^2
  # }

  if(sparser=='ash'){
    fit1 = ash(y,sigma)
    Emu[,1] = fit1$result$PosteriorMean
    Emu2[,1] = fit1$result$PosteriorMean^2+(fit1$result$PosteriorSD)^2
    log_lik[1] = fit1$loglik
  }
  if(sparser=='delta_0'){
    Emu[,1] = 0
    Emu2[,1] = 0
    log_lik[1] = sum(dnorm(y,0,sigma,log=TRUE))
  }

  if(smoother=='wavelet'){
    fit2 = smash_wave(y,sigma=sigma,filter.number=filter.number,family=family)
    Emu[,2] = fit2$mu.est
    Emu2[,2] = fit2$mu.est.var+(fit2$mu.est)^2
    log_lik[2] = fit2$loglik
  }
  if(smoother=='gp'){
    fit2 = eb_gp(x,y,nugget=sigma^2)
    Emu[,2] = fit2$pM
    Emu2[,2] = (fit2$pM)^2+diag(fit2$pV)
    log_lik[2] = fit2$loglik
  }


  for(i in 1:maxiter){

    for(k in 1:K){
      #update beta
      Ebeta[,k] = (sigma02*y+Emu[,k])/(1+sigma02)
      Ebeta2[,k] = Ebeta[,k]^2+sigma2*sigma02/(1+sigma02)
      #update phi
      #phi[,k] = (sigma02*y+Emu[,k])*Ebeta[,k]-(1+sigma02)/2*Ebeta[,k]^2-1/2*Emu2[,k]
    }

    phi = Ebeta^2/(2*sigma2*sigma02/(1+sigma02)) - Emu2/(2*sigma2*sigma02)
    phi = phi-apply(phi,1,max)
    phi = exp(phi)%*%diag(w)
    phi = phi/rowSums(phi)



    # update w
    if(!fix_w){
      w = colMeans(phi)
    }


    # update Emu, Emu2

    if(sparser=='ash'){
      fit1 = ash(Ebeta[,1],sqrt(sigma2*sigma02/phi[,1]))
      Emu[,1] = fit1$result$PosteriorMean
      Emu2[,1] = (fit1$result$PosteriorMean)^2 + (fit1$result$PosteriorSD)^2
      log_lik[1] = fit1$loglik
    }
    if(sparser=='delta_0'){
      log_lik[1] = sum(dnorm(Ebeta[,1],0,sqrt(sigma2*sigma02/phi[,1]),log=TRUE))
    }


    if(smoother=='wavelet'){
      fit = smash_wave(Ebeta[,2],sigma = sqrt(sigma2*sigma02/phi[,2]),filter.number=filter.number,family=family)
      Emu[,2] = fit$mu.est
      Emu2[,2] = fit$mu.est.var+(fit$mu.est)^2
      log_lik[2] = fit$loglik
    }
    if(smoother=='gp'){
      fit = eb_gp(x,Ebeta[,2],sigma2*sigma02/phi[,2])
      Emu[,2] = fit$pM
      Emu2[,2] = (fit$pM)^2+diag(fit2$pV)
      log_lik[2] = fit$loglik
    }


    # fit = smash.gaus(Ebeta[,2],sigma = sqrt(sigma2*sigma02/phi2),filter.number=filter.number,family=family,post.var = TRUE)
    # Emu = fit$mu.est
    # Emu2 = fit$mu.est.var+Emu^2

    # fit = ashr::ash(sigma02/(1+sigma02)*(y-Emu),sqrt(sigma2*sigma02/phi2))
    # Emu = Emu + fit$result$PosteriorMean
    # Emu2 = (fit$result$PosteriorSD)^2+(Emu)^2

    # plot(y,col='grey80')
    # lines(mu,col='grey80')
    # lines(rowSums(phi*Ebeta),col=4,lwd=2)
    # plot(phi[,2])


    # calc elbo

    elbo[i+1] = calc_elbo(y,sigma,sigma0,log_lik,phi,w,Ebeta,Ebeta2,Emu,Emu2,K)

    if(i%%printevery==0){print(sprintf("done iteration %d, elbo %f",i,elbo[i+1]))}
    if(abs(elbo[i+1]-elbo[i])<=tol){
      break
    }



    # if(sqrt(crossprod(w.old-w))<=tol){
    #   break
    # }

  }

  return(list(Ebeta=Ebeta,phi=phi,w=w,Emu=Emu,Emu2=Emu2,pm = rowSums(phi*Ebeta),elbo=elbo))
}


sparseWS0 = function(x,y,sigma,maxiter=1000,tol = 0.01,fix_pi1 = TRUE,pi1 = 0.1,
                     smoother = 'gp',printevery = 1,
                     filter.number=10,family="DaubLeAsymm"){


  calc_elbo = function(y,sigma,m,s2,phi,pi1,Eloggq){
    sigma2 = sigma^2
    n = length(y)
    -1/2*log(2*pi*sigma2)*n-1/(2*sigma2)*sum((y^2-2*y*m*phi+s2*phi))+log(pi1)*sum(phi)
    +log(1-pi1)*sum((1-phi))-sum(phi*log(phi))-sum((1-phi)*log(1-phi))+Eloggq
  }

  calc_Eloggq = function(lg,y,sigma,phi,m,s2){
    n = length(y)
    sigma2 = sigma^2
    lg+1/2*sum(log(2*pi*sigma/phi))+1/(2*sigma2)*sum((phi*(y^2-2*y*m+s2)))
    # lg+n/2*log(2*pi*sigma2)-1/2*sum(log(phi))+1/(2*sigma2)*sum((phi*(y^2-2*y*m+s2)))
  }


  # init
  n = length(y)
  if(smoother=='wavelet'){
    fit = smash_wave(y,sigma=sigma,filter.number=filter.number,family=family)
    m = fit$mu.est
    s2 = fit$mu.est.var+m^2
  }else if(smoother=='gp'){
    fit = eb_gp(x,y,nugget=sigma^2)
    m = fit$pM
    s2 = diag(fit$pV)+m^2
  }


  # init phi

  phi = c(rep(0.1,40),rep(0.9,n-40))


  #m = y
  #s2 = y^2

  sigma2 = sigma^2
  elbo = -Inf
  for(i in 1:maxiter){


    if(smoother=='wavelet'){
      fit = smash_wave(y,sigma = sqrt(sigma2/phi),filter.number=filter.number,family=family)
      m = fit$mu.est
      s2 = fit$mu.est.var+m^2
    }else if(smoother=='gp'){
      fit = eb_gp(x,y,nugget=sigma2/phi)
      m = fit$pM
      s2 = diag(fit$pV)+m^2
    }
    Eloggq = calc_Eloggq(fit$loglik,y,sigma,phi,m,s2)


    # update phi
    phi = (y*m/sigma2-1/2/sigma2*s2)
    phi = cbind(rep(0,n),phi)
    phi = phi - apply(phi,1,max)
    phi = exp(phi)%*%diag(c(1-pi1,pi1))
    phi = phi[,2]/rowSums(phi)

    # print('update phi')
    # print(calc_elbo(y,sigma,m,s2,phi,pi1,Eloggq))
    # update pi1
    if(!fix_pi1){
      pi1 = sum(phi)/n
    }
    # print('update pi')
    # print(calc_elbo(y,sigma,m,s2,phi,pi1,Eloggq))

    # update m, s2


    # print('update mu')
    # print(calc_elbo(y,sigma,m,s2,phi,pi1,Eloggq))

    # calc elbo

    elbo[i+1] = calc_elbo(y,sigma,m,s2,phi,pi1,Eloggq)

    if(i%%printevery==0){print(sprintf("done iteration %d, elbo %f",i,elbo[i+1]))}

    if(abs(elbo[i+1]-elbo[i])<=tol){
      break
    }

    # if(sqrt(crossprod(m.old-m))<=tol){
    #   break
    # }

  }
  pm = m*phi
  #pvar = s2*phi-pm^2
  return(list(pm=pm,elbo=elbo,phi=phi,pi1=pi1,m=m,s2=s2))
}




