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
    plot(EF[1:(p/2),k],type='l',ylim = range(EF[,k]),ylab=paste('f',k,sep=''),xlab='',main=main,col=2)
    lines(EF[(p/2+1):(p),k],type='l',col=3)

  }
}


n = 45
p = 256
K= 2

f0=0.5
FF = matrix(f0, nrow=p, ncol=3)
f1 = 5
f2 = 1
FF[(p/8*1):(p/8*2),1] = f1
FF[(p/8*3):(p/8*4),2] = f2
FF[(p/8*5):(p/8*7),3] = f1

FF2 = matrix(f0, nrow=p, ncol=3)
FF2[(p/8*1):(p/8*2),1] = f2
FF2[(p/8*3+10):(p/8*4-10),2] = f1
FF2[(p/8*5):(p/8*7),3] = f2

FFF = cbind(c(FF[,1],FF2[,3]),c(FF[,2],FF2[,2]),c(FF[,3],FF2[,1]))

par(mfrow=c(3,1))
for(k in 1:K){
  plot(FFF[1:p,k],type='l',ylim = range(FFF),ylab='',xlab='')
  lines(FFF[(p+1):(2*p),k],type='l',col=4)
}

set.seed(123)
L = matrix(0,nrow=n,ncol=K)
pi0=1/3
while(sum(rowSums(L)==0)>0){
  for(k in 1:K){
    l1 = c(rep(0,n*pi0),
           rnorm(n*(1-pi0)/2,0,sqrt(0.25)),
           rnorm(n*(1-pi0)/2,0,sqrt(1)))
    L[,k] =  l1[sample(1:n)]
  }
}

s = mean(apply(tcrossprod(L,FFF[,1:K]),1,sd))
y = tcrossprod(L,FFF[,1:K]) + matrix(rnorm(n*p*2,0,s),nrow=n,ncol=p*2)

fit.flash = flash(y,var_type = 'by_row')
fit.dwt = wave_ebmf(y,type='wavelet')
#fit.ndwt = wave_ebmf(y,type = 'station',K=3)
#tic();fit.dwt2 = smooth_flash(y,type='station');toc()

ploter2(fit.flash$ldf$f,main=NULL)
ploter2(fit.dwt$ldf$f,main=NULL)


s = mean(apply(tcrossprod(L,FFF[,1:K]),1,sd))*3
y = tcrossprod(L,FFF[,1:K]) + matrix(rnorm(n*p*2,0,s),nrow=n,ncol=p*2)

fit.flash = flash(y,var_type = 'by_row')
fit.dwt = wave_ebmf(y,type='wavelet')
#fit.ndwt = wave_ebmf(y,type = 'station',K=3)
#tic();fit.dwt2 = smooth_flash(y,type='station');toc()

ploter2(fit.flash$ldf$f,main=NULL)
ploter2(fit.dwt$ldf$f,main=NULL)




# library(reshape2)
# library(ggplot2)
# LL=  melt(L)
# ggplot(LL, aes(x = Var2, y = Var1)) + geom_raster(aes(fill=value)) +
#   scale_fill_gradient(low="grey90", high="red") +
#   labs(x="column", y="row", title="L")



#ploter2(fit.dwt2$ldf$f,main='smooth_flash')


###########################
###########################
###### Poisson case #######
library(ebpmf)

# ploter2 = function(EF,main){
#   p = nrow(EF)
#   par(mfrow=c(ncol(EF),1))
#   for(k in 1:ncol(EF)){
#     plot(EF[1:(p/2),k],type='l',ylim = range(EF[,k]),ylab='',xlab='',main=main,col=2)
#     lines(EF[(p/2+1):(p),k],type='l',col=3)
#
#   }
# }

ploter2 = function(EF,main=NULL){
  p = nrow(EF)
  library(ggplot2)
  library(tidyr)

  # Convert the matrix EF to a data frame and gather the columns
  EF_df <- as.data.frame(EF) %>%
    mutate(row_number = 1:n()) %>%
    pivot_longer(cols = -row_number, names_to = "column", values_to = "value")

  # Create the ggplot2 plot
  pp <- ggplot(EF_df, aes(x = row_number, y = value)) +
    geom_line(color = 2, size = 1) +
    facet_wrap(~column, ncol = 1, scales = "free_y") +
    labs(x = "", y = "") +
    theme_minimal() +
    theme(
      strip.text = element_blank(),
      strip.background = element_blank(),
      #axis.text.y = element_blank()  # Remove y-axis numbers
    )+
    ggtitle(main)

  # Print the plot
  print(pp)
}

K = 2
n=45
p=256
# l0 = 0.1
# l1 = 1
# L = cbind(c(rep(l1,n/3),rep(l0,n/3*2)),
#           c(rep(l0,n/3),rep(l1,n/3),rep(l0,n/3)),
#           c(rep(l0,n/3*2),rep(l1,n/3)))

f0 = 0.5

FF = matrix(f0, nrow=p, ncol=3)
f1 = 10
f2 = 5
FF[(p/8*1):(p/8*2),1] = f1
FF[(p/8*3):(p/8*4),2] = f2
FF[(p/8*5):(p/8*7),3] = f1

FF2 = matrix(f0, nrow=p, ncol=3)
FF2[(p/8*1):(p/8*2),1] = f2
FF2[(p/8*3+10):(p/8*4-10),2] = f1
FF2[(p/8*5):(p/8*7),3] = f2

FFF = cbind(c(FF[,1],FF2[,3]),c(FF[,2],FF2[,2]),c(FF[,3],FF2[,1]))
FFFn = exp(log(FFF)+matrix(rnorm(nrow(FFF)*ncol(FFF),0,0.1),nrow=nrow(FFF)))
#
# par(mfrow=c(K,1),tcl=-0.5,mai=c(0.5,1,0.2,0.3))
# for(k in 1:K){
#   plot(FFF[1:p,k],type='l',ylim = range(FFF),ylab=paste('f',k,sep=''),xlab='',lwd=2,col=2)
#   lines(FFF[(p+1):(2*p),k],type='l',col=3,lwd=2)
#
# }

set.seed(12345)

L = matrix(0,nrow=n,ncol=K)
pi0=1/3
while(sum(rowSums(L)==0)>0){
  for(k in 1:K){
    l1 = c(rexp(n*pi0,10),
           rexp(n*(1-pi0)/2,5),
           rexp(n*(1-pi0)/2,1))
    L[,k] =  l1[sample(1:n)]
  }
}

y = matrix(rpois(n*p*2,tcrossprod(L,FFF[,1:K])),nrow=n)
fit.nmf = fastTopics::fit_poisson_nmf(y,k=K)
temp = poisson_to_multinom(fit.nmf$F,fit.nmf$L)
ploter2(temp$FF,main=NULL)

fit.stm5 = ebpmf_identity(y,K=K)
ploter2(fit.stm5$EF,main=NULL)

L_true_scaled = poisson_to_multinom(FFF[,1:K],L)
ploter2(L_true_scaled$FF,main=NULL)



set.seed(2023)
f0 = 0.5
FF = matrix(f0, nrow=p, ncol=3)
f1 = 5
f2 = 1
FF[(p/8*1):(p/8*2),1] = f1
FF[(p/8*3):(p/8*4),2] = f2
FF[(p/8*5):(p/8*7),3] = f1

FF2 = matrix(f0, nrow=p, ncol=3)
FF2[(p/8*1):(p/8*2),1] = f2
FF2[(p/8*3+10):(p/8*4-10),2] = f1
FF2[(p/8*5):(p/8*7),3] = f2

FFF = cbind(c(FF[,1],FF2[,3]),c(FF[,2],FF2[,2]),c(FF[,3],FF2[,1]))
FFFn = exp(log(FFF)+matrix(rnorm(nrow(FFF)*ncol(FFF),0,0.1),nrow=nrow(FFF)))

L = matrix(0,nrow=n,ncol=K)
pi0=1/3
while(sum(rowSums(L)==0)>0){
  for(k in 1:K){
    l1 = c(rexp(n*pi0,10),
           rexp(n*(1-pi0)/2,5),
           rexp(n*(1-pi0)/2,1))
    L[,k] =  l1[sample(1:n)]
  }
}

y = matrix(rpois(n*p*2,tcrossprod(L,FFFn[,1:K])),nrow=n)
# kmeans.init = kmeans(y, K, nstart = 20)
# L0 = rep(1, n) %o% (as.vector(table(kmeans.init$cluster)))
# F0 = t(kmeans.init$centers)

fit.nmf = fastTopics::fit_poisson_nmf(y,k=K)
temp = poisson_to_multinom(fit.nmf$F,fit.nmf$L)
ploter2(temp$FF,main=NULL)

fit.stm5 = ebpmf_identity(y,K=K)
ploter2(fit.stm5$EF,main=NULL)

L_true_scaled = poisson_to_multinom(FFF[,1:K],L)
ploter2(L_true_scaled$FF,main=NULL)

#
# plot(fit.stm5$EL[,2],L_true_scaled$L[,1])
# plot(fit.stm5$EL[,1],L_true_scaled$L[,2])
#
# plot(temp$L[,1],L_true_scaled$L[,1])


