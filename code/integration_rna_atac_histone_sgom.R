library(stm)
library(NNLM)
library(Matrix)

files = list.files('/project2/mstephens/dongyue/luis/assays')
genes = c()
for(i in 1:length(files)){
  genes[i] = strsplit(files[i],split = '_')[[1]][1]
}
genes = unique(genes)
genes

library(readr)
merge_len = 10

for(g in 4:6){
  gene = genes[g]
  print(paste('running ',gene))
  fit = readRDS(paste('output/luis/',gene,'_K10_merge10base.rds',sep=''))
  X = fit$X
  fit_sgom = cluster.mix(X,K=10,tol=1e-2,maxit = 30,nugget=TRUE)
  fit_sgom = list(pi = fit_sgom$pi,phi = fit_sgom$phi)

  saveRDS(list(gene=gene,
               merge_len=merge_len,
               X=Matrix(X,sparse=TRUE),
               assays = c('RNA','H3K4me1','ATAC'),
               fit_sgom = fit_sgom),file=paste('output/luis/',gene,'_K10_merge10base_sgom.rds',sep=''))

}




