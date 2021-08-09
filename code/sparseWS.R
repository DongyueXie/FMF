
smash_wave = function(x,sigma,filter.number=10,family="DaubLeAsymm"){
  n = length(x)
  J = log(n,2)
  if(length(sigma==1)){
    sigma = rep(sigma,n)
  }
  tsum = sum(x)
  x.w = wavethresh::wd(x, filter.number = filter.number,
                       family = family, type = "wavelet")
  ashparam = list()
  ashparam$gridmult = 64
  data.var = sigma^2

  W = GenW(n,filter.number,family)
  W = t(W)
  x.w.v =  colSums(t(W)^2*data.var)[-1]

  return.loglr = TRUE
  post.var = TRUE
  logLR.scale = 0
  loglik.scale = c()
  g_hat = list()
  x.w.v.s = rep(0, 2^J-1)
  for (j in 0:(J - 1)) {
    x.pm = rep(0, 2^j)

    #index = (((J - 1) - j) * n + 1):((J - j) * n)
    index = (n-2^(j+1)+1):(n-2^j)
    x.w.j = wavethresh::accessD(x.w, j)
    x.w.v.j = x.w.v[index]
    ind.nnull = (x.w.v.j != 0)
    zdat.ash = shrink.wc(x.w.j[ind.nnull],
                         sqrt(x.w.v.j[ind.nnull]), ashparam, jash = FALSE,
                         df = NULL, SGD = FALSE)
    g_hat[[j+1]] = zdat.ash$fitted_g
    x.pm[ind.nnull] = get_pm(zdat.ash)
    x.pm[!ind.nnull] = 0
    x.w = wavethresh::putD(x.w, j, x.pm)
    if (return.loglr == TRUE) {
      #spins = 2^(J - j)
      loglik.scale[j + 1] = calc_loglik(get_fitted_g(zdat.ash),
                                        set_data(x.w.j[ind.nnull],
                                                 sqrt(x.w.v.j[ind.nnull]), NULL, 0))
      logLR.temp = loglik.scale[j + 1] -
        sum(dnorm(x.w.j[ind.nnull], 0, sqrt(x.w.v.j[ind.nnull]),
                  log = TRUE))
      logLR.scale[j + 1] = logLR.temp
    }
    if (post.var == TRUE) {
      x.w.v.s[index[ind.nnull]] = get_psd(zdat.ash)^2
      x.w.v.s[index[!ind.nnull]] = 0
    }
  }
  mu.est = wr(x.w)
  loglik = sum(loglik.scale)
  x.w.v.s = c(0,x.w.v.s)
  mu.est.var = colSums(W^2*x.w.v.s)
  return(list(mu.est=mu.est,mu.est.var=mu.est.var,loglik = loglik,g_hat=g_hat))
}

sparseWS = function(y,sigma,maxiter=1000,tol = 0.01,fix_pi1 = TRUE,pi1 = 0.1,
                    filter.number=10,family="DaubLeAsymm"){
  # init
  n = length(y)
  fit0 = smash_wave(y,sigma=sigma,filter.number=filter.number,family=family)
  m = fit0$mu.est
  s2 = fit0$mu.est.var+m^2

  sigma2 = sigma^2
  elbo = -Inf
  for(i in 1:maxiter){
    # update phi
    phi = pi1*exp(y*m/sigma2-1/2/sigma2*s2)
    phi = phi/(phi+(1-pi1))
    # update pi1
    if(!fix_pi1){
      pi1 = sum(phi)/n
    }

    # update m, s2
    fit = smash_wave(y,sigma = sqrt(sigma2/phi),filter.number=filter.number,family=family)

    m = fit$mu.est
    s2 = fit$mu.est.var+m^2
    # calc elbo
    Eloggq = calc_Eloggq(fit$loglik,y,sigma,phi,m,s2)
    elbo[i+1] = calc_elbo(y,sigma,m,s2,phi,pi1,Eloggq)
    if(abs(elbo[i+1]-elbo[i])<=tol){
      break
    }
  }
  pm = m*phi
  pvar = s2*phi-pm^2
  return(list(pm=pm,pvar=pvar,elbo=elbo,phi=phi,pi1=pi1,m=m,s2=s2))
}

calc_elbo = function(y,sigma,m,s2,phi,pi1,Eloggq){
  sigma2 = sigma^2
  n = length(y)
  -1/2*log(sigma2)*n-1/(2*sigma2)*sum((y^2-2*y*m*phi+s2*phi))+log(pi1)*sum(phi)
  +log(1-pi1)*sum((1-phi))-sum((phi*log(phi)),na.rm = TRUE)-sum(((1-phi)*log(1-phi)),na.rm = TRUE)+Eloggq
}

calc_Eloggq = function(lg,y,sigma,phi,m,s2){
  n = length(y)
  sigma2 = sigma
  lg+n/2*log(sigma2)-1/2*sum(log(phi))+1/(2*sigma2)*sum((phi*(y^2-2*y*m+s2)))
}

set.seed(12345)
filter.number=10
family="DaubLeAsymm"
n=512
x = seq(-2*pi,pi,length.out = n)
mu = 1.5*sin(x) + sin(2*x)
idx=which.min(mu[1:200])
mu = mu-mu[idx]
mu[1:idx] = 0
y = mu + rnorm(n,sd=0.5)
plot(y,col='grey80')
lines(mu)
lines(smash.gaus(y,0.5,filter.number=filter.number,family=family),col=2)
fit.ss = sparseWS(y,0.5,maxiter = 100,tol=1e-5,pi1=0.5)
lines(fit.ss$pm,col=3,lwd=2)



