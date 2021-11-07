source('/home/dyxie/smash-gen/code/smash_gen_poisson.R')

files = list.files('/project2/mstephens/dongyue/luis/assays')
genes = c()
for(i in 1:length(files)){
  genes[i] = strsplit(files[i],split = '_')[[1]][1]
}
genes = unique(genes)
g = 7
gene = genes[g]
print(paste('running ',gene))
RNAseq = read_csv(paste("/project2/mstephens/dongyue/luis/assays/",gene,"_RNAseq.csv.gz",sep=''))
ATACseq = read_csv(paste("/project2/mstephens/dongyue/luis/assays/",gene,"_ATACseq.csv.gz",sep=''))
H3K4seq = read_csv(paste("/project2/mstephens/dongyue/luis/assays/",gene,"_H3K4me1.csv.gz",sep=''))
indis_ATAC = (colnames(ATACseq))[-c(1,2)]
indis_RNA = (colnames(RNAseq))[-c(1,2)]
indis_H3K4me1 = (colnames(H3K4seq))[-c(1,2)]
for(i in 1:length(indis_H3K4me1)){
  name_i = strsplit(indis_H3K4me1[i],split = '_')[[1]]
  indis_H3K4me1[i] = paste(name_i[1],name_i[2],sep = '_')
}
indis = intersect(intersect(indis_ATAC,indis_RNA),indis_H3K4me1)

# use only NI individuals
indis = indis[seq(2,54,by=2)]

base.idx = (as.numeric(strsplit(RNAseq$locus[1],split=':')[[1]][2])):(as.numeric(strsplit(RNAseq$locus[length(RNAseq$locus)],split=':')[[1]][2]))
# FAM220A, c(6369042,6388598)
# GTF3C6, c(111279909,111289075)
# HDDC2, c(125596496,125623113)
# HLA-DRB5, c(32485130,32498064)
# IRF2, c(185308883,185395704)
# MRPL18, c(160210844,160219461)
start.idx = 160210844-800
end.idx = 160219461+800
idx = (which(base.idx==start.idx)):(which(base.idx==end.idx))

plot(rowSums(RNAseq[idx,match(indis,indis_RNA)+2]),pch=20,col='grey80',type='h',ylab='counts',xlab='Chr6',xaxt = "n",ylim=c(0,1000))
xlab.idx = round(seq(1,length(idx),length.out = 5))
axis(1, at = xlab.idx,labels = (start.idx:end.idx)[xlab.idx])

plot(rowSums(ATACseq[idx,match(indis,indis_ATAC)+2]),pch=20,col='grey80',type='h',ylab='counts',xlab='Chr6',xaxt = "n")
xlab.idx = round(seq(1,length(idx),length.out = 5))
axis(1, at = xlab.idx,labels = (start.idx:end.idx)[xlab.idx])

plot(rowSums(H3K4seq[idx,match(indis,indis_H3K4me1)+2]),pch=20,col='grey80',type='h',ylab='counts',xlab='Chr6',xaxt = "n")
xlab.idx = round(seq(1,length(idx),length.out = 5))
axis(1, at = xlab.idx,labels = (start.idx:end.idx)[xlab.idx])


fit.rna = smash.gen.poiss(rowSums(RNAseq[idx,match(indis,indis_RNA)+2]),smooth_method = 'smash')
fit.atac = smash.gen.poiss(rowSums(ATACseq[idx,match(indis,indis_ATAC)+2]),smooth_method = 'smash')
fit.h = smash.gen.poiss(rowSums(H3K4seq[idx,match(indis,indis_H3K4me1)+2]),smooth_method = 'smash')

normalize = function(x){x/sum(x)}

par(mfrow=c(3,1))
plot((fit.rna$lambda.est),type='h',col='grey50',ylab='',xlab='Chr6',xaxt = "n",main='RNA-seq')
xlab.idx = round(seq(1,length(idx),length.out = 5))
axis(1, at = xlab.idx,labels = (start.idx:end.idx)[xlab.idx])

plot((fit.atac$lambda.est),type='h',col='grey50',ylab='',xlab='Chr6',xaxt = "n",main='ATAC-seq')
xlab.idx = round(seq(1,length(idx),length.out = 5))
axis(1, at = xlab.idx,labels = (start.idx:end.idx)[xlab.idx])

plot((fit.h$lambda.est),type='h',col='grey50',ylab='',xlab='Chr6',xaxt = "n",main='Histone mark')
xlab.idx = round(seq(1,length(idx),length.out = 5))
axis(1, at = xlab.idx,labels = (start.idx:end.idx)[xlab.idx])


fit.rna2 = smash.pois.gaus(rowSums(RNAseq[idx,match(indis,indis_RNA)+2]))
fit.atac2 = smash.pois.gaus(rowSums(ATACseq[idx,match(indis,indis_ATAC)+2]))
fit.h2 = smash.pois.gaus(rowSums(H3K4seq[idx,match(indis,indis_H3K4me1)+2]))


par(mfrow=c(3,1))
plot((fit.rna2$lambda.est),type='h',col='grey50',ylab='',xlab='Chr6',xaxt = "n",main='RNA-seq')
xlab.idx = round(seq(1,length(idx),length.out = 5))
axis(1, at = xlab.idx,labels = (start.idx:end.idx)[xlab.idx])

plot((fit.atac2$lambda.est),type='h',col='grey50',ylab='',xlab='Chr6',xaxt = "n",main='ATAC-seq')
xlab.idx = round(seq(1,length(idx),length.out = 5))
axis(1, at = xlab.idx,labels = (start.idx:end.idx)[xlab.idx])

plot((fit.h2$lambda.est),type='h',col='grey50',ylab='',xlab='Chr6',xaxt = "n",main='Histone mark')
xlab.idx = round(seq(1,length(idx),length.out = 5))
axis(1, at = xlab.idx,labels = (start.idx:end.idx)[xlab.idx])


# k = 31
# plot(runmed(rowSums(RNAseq[idx,match(indis,indis_RNA)+2]),k=k),type='h')
# plot(runmed(rowSums(ATACseq[idx,match(indis,indis_ATAC)+2]),k=k),type='h')
# plot(runmed(rowSums(H3K4seq[idx,match(indis,indis_H3K4me1)+2]),k=k),type='h')

