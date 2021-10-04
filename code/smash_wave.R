
library(devtools)
load_all('~/Documents/Rpackages/flashr')
library(wavethresh)


smash_wave = function(x,sigma,filter.number=10,family="DaubLeAsymm",ebnm_fn = 'ebnm_ash',ebnm_param=list()){
  n = length(x)
  J = log(n,2)

  # ebnm_param = list(gridmult = 64)
  # ashparam = list()
  # ashparam$gridmult = 64

  tsum = sum(x)
  x.w = wavethresh::wd(x, filter.number = filter.number,
                       family = family, type = "wavelet")

  data.var = sigma^2
  if(length(data.var==1)){
    data.var = rep(data.var,n)
  }

  W = GenW(n,filter.number,family)
  W = t(W)
  if(length(sigma)==1){
    x.w.v = rep(sigma^2,n-1)
    tsum.var = sigma^2
  }else{
    x.w.v =  colSums(t(W)^2*data.var)
    tsum.var = x.w.v[1]
    x.w.v = x.w.v[-1]
  }



  # return.loglr = TRUE
  # post.var = TRUE
  # logLR.scale = 0
  loglik.scale = c()
  x.w.v.s = rep(0, 2^J-1)
  for (j in 0:(J - 1)) {
    x.pm = rep(0, 2^j)
    #index = (((J - 1) - j) * n + 1):((J - j) * n)
    index = (n-2^(j+1)+1):(n-2^j)
    x.w.j = wavethresh::accessD(x.w, j)
    x.w.v.j = x.w.v[index]
    ind.nnull = (x.w.v.j != 0)
    a = do.call(ebnm_fn,list(x.w.j[ind.nnull], sqrt(x.w.v.j[ind.nnull]), ebnm_param))
    # zdat.ash = shrink.wc(x.w.j[ind.nnull],
    #                      sqrt(x.w.v.j[ind.nnull]), ashparam, jash = FALSE,
    #                      df = NULL, SGD = FALSE)
    x.pm[ind.nnull] = a$postmean
    x.pm[!ind.nnull] = 0
    x.w = wavethresh::putD(x.w, j, x.pm)
    loglik.scale[j + 1] = a$penloglik
    x.w.v.s[index[ind.nnull]] = a$postmean2-(a$postmean)^2
    x.w.v.s[index[!ind.nnull]] = 0
    # if (return.loglr == TRUE) {
    #   spins = 2^(J - j)
    #   loglik.scale[j + 1] = calc_loglik(get_fitted_g(zdat.ash),
    #                                     set_data(x.w.j[ind.nnull],
    #                                              sqrt(x.w.v.j[ind.nnull]), NULL, 0))
    #
    #   logLR.temp = loglik.scale[j + 1] -
    #     sum(dnorm(x.w.j[ind.nnull], 0, sqrt(x.w.v.j[ind.nnull]),
    #               log = TRUE))
    #   logLR.scale[j + 1] = logLR.temp
    # }

  }
  mu.est = wr(x.w)
  loglik = sum(loglik.scale)
  x.w.v.s = c(tsum.var,x.w.v.s)
  mu.est.var = colSums(W^2*x.w.v.s)
  return(list(mu.est=mu.est,mu.est.var=mu.est.var,loglik = loglik))
}

