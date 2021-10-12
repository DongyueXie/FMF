library(stm)
library(NNLM)
library(Matrix)

files = list.files('/project2/mstephens/dongyue/luis/luis')
genes = c()
for(i in 1:length(files)){
  genes[i] = strsplit(files[i],split = '_')[[1]][1]
}
genes = unique(genes)
genes

library(readr)
merge_len = 10

for(g in 1:3){
  gene = genes[g]
  print(paste('running ',gene))
  RNAseq = read_csv(paste("/project2/mstephens/dongyue/luis/luis/",gene,"_RNAseq.csv.gz",sep=''))
  ATACseq = read_csv(paste("/project2/mstephens/dongyue/luis/luis/",gene,"_ATACseq.csv.gz",sep=''))
  H3K4seq = read_csv(paste("/project2/mstephens/dongyue/luis/luis/",gene,"_H3K4me1.csv.gz",sep=''))

  indis_ATAC = (colnames(ATACseq))[-c(1,2)]
  indis_RNA = (colnames(RNAseq))[-c(1,2)]
  indis_H3K4me1 = (colnames(H3K4seq))[-c(1,2)]
  for(i in 1:length(indis_H3K4me1)){
    name_i = strsplit(indis_H3K4me1[i],split = '_')[[1]]
    indis_H3K4me1[i] = paste(name_i[1],name_i[2],sep = '_')
  }
  indis = intersect(intersect(indis_ATAC,indis_RNA),indis_H3K4me1)

  Y_RNA = t(RNAseq[,match(indis,indis_RNA)+2])
  Y_H3K4 = t(H3K4seq[,match(indis,indis_H3K4me1)+2])
  Y_ATAC = t(ATACseq[,match(indis,indis_ATAC)+2])


  if(var(c(ncol(Y_RNA),ncol(Y_H3K4),ncol(Y_ATAC)))!=0){
    l = min(c(ncol(Y_RNA),ncol(Y_H3K4),ncol(Y_ATAC)))
    Y_RNA = Y_RNA[,1:l]
    Y_H3K4 = Y_H3K4[,1:l]
    Y_ATAC = Y_ATAC[,1:l]
  }

  if(ncol(Y_RNA)%%merge_len!=0){
    l = ncol(Y_RNA)
    l = floor(l/merge_len)*merge_len
    Y_RNA = Y_RNA[,1:l]
    Y_H3K4 = Y_H3K4[,1:l]
    Y_ATAC = Y_ATAC[,1:l]
  }

  Y_RNAr = do.call(cbind, by(t(Y_RNA), (seq(ncol(Y_RNA)) - 1) %/% merge_len, FUN = colSums))
  Y_H3K4r = do.call(cbind, by(t(Y_H3K4), (seq(ncol(Y_H3K4)) - 1) %/% merge_len, FUN = colSums))
  Y_ATACr = do.call(cbind, by(t(Y_ATAC), (seq(ncol(Y_ATAC)) - 1) %/% merge_len, FUN = colSums))

  X = cbind(Y_RNAr,Y_H3K4r,Y_ATACr)

  fit_sgom = cluster.mix(X,K=10,tol=1e-3,maxit = 100,nugget=TRUE)

  saveRDS(list(gene=gene,
               merge_len=merge_len,
               X=X,
               assays = c('RNA,H3K4me1','ATAC'),
               fit_sgom = fit_sgom),file=paste('output/luis/',gene,'_K10_merge10base_sgom.rds',sep=''))

}




