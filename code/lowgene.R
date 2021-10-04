library(stm)
library(NNLM)


sum_base = function(x,base){
  x = x[1:(floor(length(x)/base)*base)]
  colSums(matrix(x,nrow=base))
}

#compare methods like nmf, stm-bmsm, stm-smashgen, stm-smashgen robust, HALS-wavelet, HALS-runmed

gene_splicing_study = function(X,name,K,nreps = 6,seed=12345,base=1,path,label){

  set.seed(seed)

  if(base>1){
    X = t(apply(X, 1, sum_base,base=base))
  }

  nmf_loss = Inf
  stm_loss = Inf
  stm_nugget_loss = Inf
  stm_nugget_robust_loss = Inf
  hals_wavelet_loss = Inf
  hals_runmed_loss = Inf

  for (reps in 1:nreps) {

    print(reps)

    fit_NMF = nnmf(X,k=K,method = 'scd',loss='mkl',verbose = F,max.iter = 800)
    if(min(fit_NMF$mkl)<nmf_loss){
      fit_NMF$label = label
      nmf_loss = min(fit_NMF$mkl)
      save(fit_NMF,file = paste(path,name,'_NMF_mkl_scd_K',K,'_base',base,'.RData',sep = ''))
    }

    fit_stm = stm(X,K,init = list(L_init=fit_NMF$W+0.01,F_init = fit_NMF$H+0.01),
                  return_all = FALSE,tol=1e-4,maxiter=50,printevery = 1e5)
    if(fit_stm$KL[length(fit_stm$KL)]<stm_loss){
      fit_stm$label = label
      stm_loss=fit_stm$KL[length(fit_stm$KL)]
      save(fit_stm,file=paste(path,name,'_stm_bmsm_K',K,'_base',base,'.RData',sep = ''))
    }

    fit_stm_nugget_robust = stm(X,K,init = list(L_init=fit_NMF$W+0.01,F_init = fit_NMF$H+0.01),
                                return_all = FALSE,tol=1e-2,nugget = TRUE,maxiter=50,printevery = 1e5)
    if(fit_stm_nugget_robust$KL[length(fit_stm_nugget_robust$KL)]<stm_nugget_robust_loss){
      fit_stm_nugget_robust$label = label
      stm_nugget_robust_loss=fit_stm_nugget_robust$KL[length(fit_stm_nugget_robust$KL)]
      save(fit_stm_nugget_robust,file=paste(path,name,'_stm_nugget_robust_K',K,'_base',base,'.RData',sep = ''))
    }

    fit_stm_nugget = stm(X,K,init = list(L_init=fit_NMF$W+0.01,F_init = fit_NMF$H+0.01),
                         return_all = FALSE,tol=1e-2,nugget = TRUE,maxiter=50,printevery = 1e5,
                         nug_control_f = list(robust=F))
    if(fit_stm_nugget$KL[length(fit_stm_nugget$KL)]<stm_nugget_loss){
      fit_stm_nugget$label = label
      stm_nugget_loss=fit_stm_nugget$KL[length(fit_stm_nugget$KL)]
      save(fit_stm_nugget,file=paste(path,name,'_stm_nugget_K',K,'_base',base,'.RData',sep = ''))
    }

    fit_hals_wavelet = NMF_HALS(X,K,smooth_method = 'wavelet',printevery = 1e5)
    if(fit_hals_wavelet$loss[length(fit_hals_wavelet$loss)]<hals_wavelet_loss){
      fit_hals_wavelet$label = label
      hals_wavelet_loss = fit_hals_wavelet$loss[length(fit_hals_wavelet$loss)]
      save(fit_hals_wavelet,file=paste(path,name,'_hals_wavelet_K',K,'_base',base,'.RData',sep = ''))
    }

    fit_hals_runmed = NMF_HALS(X,K,smooth_method = 'runmed',printevery = 1e5)

    if(fit_hals_runmed$loss[length(fit_hals_runmed$loss)]<hals_runmed_loss){
      fit_hals_runmed$label = label
      hals_runmed_loss = fit_hals_runmed$loss[length(fit_hals_runmed$loss)]
      save(fit_hals_runmed,file=paste(path,name,'_hals_runmed_K',K,'_base',base,'.RData',sep = ''))
    }
  }
}


path = '/home/dyxie/NMF/YangLi/readcount'
count.files = list.files(path=path)
data.list = lapply(count.files,function(x){
  datax = read.table(paste(path,'/',x,sep='',collapse = ''),header = T)
  quantile(rowSums(datax[,-c(1:3)])/457,0.99)
  #mean(rowSums(datax[,-c(1:3)])/457)
})
# from lowest expression to highest expression
gene.order = order(unlist(data.list),decreasing = F)

K=3
base=10
path.save = '~/SMF/data/lowgene/'
for(gene in gene.order[11:length(gene.order)]){
  X = read.table(paste(path,'/',count.files[gene],sep='',collapse = ''),header = T)
  X = t(X[,-c(1:3)])
  #remove rows with all zeros.
  label = c(rep('Adipose',226),rep('Skin',231))
  rm.idx = which(rowSums(X)==0)
  if(length(rm.idx)>0){
    X = X[-rm.idx,]
    label = label[-rm.idx]
  }
  print(count.files[gene])
  gene_splicing_study(X,count.files[gene],K,base=base,path=path.save,label=label)
}
