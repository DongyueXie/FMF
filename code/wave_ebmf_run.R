
set.seed(12345)
n = 40
p = 2^6
mu1 = c(rep(1,p/8),rep(5,p/4),rep(1,p/8*5))
mu2 = c(rep(1,p/8*5),rep(5,p/4),rep(1,p/8))
plot(mu1,type='l')
plot(mu2,type='l')

FF = rbind(mu1,mu2)
l = c(rep(1,n/2),rep(0,n/2))
L = cbind(l,1-l)

Y = L%*%FF+matrix(rnorm(n*p,0,sd=0.5),nrow=n,ncol=p)
plot(Y[1,])
plot(Y[n,])
data_flash = flash_set_data(Y)

fit_flash = flash(data_flash,var_type = 'by_row')
fit_flash$nfactors
plot(fit_flash$ldf$f[,1],type='l')
plot(fit_flash$ldf$f[,2],type='l')
fit_flash$fit$tau

fit_wave = wave_ebmf(Y,Kmax=10,nullcheck=TRUE,verbose = TRUE,est_var = 'va')
plot(fit_wave$ldf$f[,1],type='l')
plot(fit_wave$ldf$f[,2],type='l')






