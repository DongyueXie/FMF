devtools::load_all('~/smashr')
source('code/binomial_nuggest.R')

poisson_multiscale_nugget = function(x,lev=0,glm.approx.param=list()){
  n = length(x)
  J = log2(length(x))
  ls = sum(x)
  y = as.vector(t(ParentTItable(x)$parent))

  zdat = withCallingHandlers(do.call(glm.approx, c(list(x = y, g = NULL),
                                                   glm.approx.param)))

  res = list()

  ind = 1:n
  sigma2_1 = nm_nug(zdat[1, ind],zdat[2, ind])$sigma2
  for (j in 1:(J - lev)) {
    #res = getlist.res(res, j, n, zdat, FALSE, TRUE, ashparam)
    sigma2_init = sigma2_1/(2^(j-1))
    res = getlist.res2(res, j, n, zdat, FALSE, TRUE,sigma2_init=sigma2_init)
  }
  # Do not smooth for coarser levels, with everything the same as
  # above but using the estimate and its variance as the posterior
  # mean and variance ie flat prior.
  if (lev != 0) {
    for (j in (J - lev + 1):J) {
      res = getlist.res2(res, j, n, zdat, FALSE, FALSE)
    }
  }

  est.mean = exp(reverse.pwave(log(ls),
                               log(matrix(res$lp.mean, J,n, byrow = TRUE)),
                               log(matrix(res$lq.mean, J, n, byrow = TRUE))))
  est.ms = exp(reverse.pwave(2 * log(ls),
                                log(matrix(res$lp.var, J, n, byrow = TRUE)),
                                log(matrix(res$lq.var, J, n, byrow = TRUE))))
  est.var = pmax(est.ms - est.mean^2, 0)
  return(list(est.mean = est.mean, est.var = est.var))
}



# n = length(x)
# J = log2(length(x))
# ls = sum(x)
# y = as.vector(t(ParentTItable(x)$parent))
# glm.approx.param = setGlmApproxParam(list())
# zdat = withCallingHandlers(do.call(glm.approx, c(list(x = y, g = NULL),
#                                                  glm.approx.param)))
#
# res = list()
# ashparam = setAshParam.poiss(list())
# # Loop through resolutions, smoothing each resolution separately
# lev = 3
# for (j in 1:(J - lev)) {
#   #res = getlist.res(res, j, n, zdat, FALSE, TRUE, ashparam)
#   res = getlist.res2(res, j, n, zdat, FALSE, TRUE)
# }
# # Do not smooth for coarser levels, with everything the same as
# # above but using the estimate and its variance as the posterior
# # mean and variance ie flat prior.
# if (lev != 0) {
#   for (j in (J - lev + 1):J) {
#     res = getlist.res2(res, j, n, zdat, FALSE, FALSE)
#   }
# }
#
# est.mean = exp(reverse.pwave(log(ls),
#                              log(matrix(res$lp.mean, J,n, byrow = TRUE)),
#                              log(matrix(res$lq.mean, J, n, byrow = TRUE))))
#
# plot(x)
# lines(est.mean)


getlist.res2 = function (res, j, n, zdat, log, shrink, sigma2_init = NULL) {
  ind = ((j - 1) * n + 1):(j * n)
  if (shrink == TRUE) {

    # Apply ash to vector of intercept estimates and SEs.
    # zdat.ash = withCallingHandlers(do.call(ash,
    #                                       c(list(betahat = zdat[1, ind], sebetahat = zdat[2, ind]), ashparam)))
    zdat.ash = nm_nug(zdat[1, ind],zdat[2, ind],sigma2_init=sigma2_init)
    # Find mean and variance of alpha.
    alpha.mv = list(mean = zdat.ash$mu,
                    var = zdat.ash$mu_var)
  } else {
    alpha.mv = list(mean = fill.nas(zdat[1, ind]),
                    var = fill.nas(zdat[2, ind])^2)
  }
  res.j = compute.res(alpha.mv, log)
  res = rbindlist(list(res, res.j))
  return(res)
}



##############
##############
set.seed(12345)
n = 1024
mu = c(rep(30,n/2),rep(50,n/4),rep(30,n/4))
x = rpois(n,exp(log(mu)+rnorm(n,0,0.2)))
#x = rpois(n,log(1+exp(mu+rnorm(n,0,10))))
plot(x,col='grey80')
lines(poisson_multiscale_nugget(x,lev=0)$est.mean,col='grey50')
lines(poisson_multiscale_nugget(x,lev = 1)$est.mean,col=1)
lines(poisson_multiscale_nugget(x,lev = 2)$est.mean,col=2)
lines(poisson_multiscale_nugget(x,lev = 3)$est.mean,col=3)
lines(poisson_multiscale_nugget(x,lev = 4)$est.mean,col=4)
lines(poisson_multiscale_nugget(x,lev = 5)$est.mean,col=5)

plot(x,col='grey80')
lines(smash.poiss(x))

