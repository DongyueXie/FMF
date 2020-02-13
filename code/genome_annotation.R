load_marks_data <- function(cell_type, chr_name, marks_list){
  dir_data <- paste0(dir_home, '/chrom_cluster/ENCODE_data/hg19/meancounts/', chr_name)
  count_marks.df <- readRDS(file.path(dir_data, paste0(cell_type, '_10marks_meanChIP_hg19_200bp_bins.rds')))
  count_marks.m <- as.matrix(count_marks.df[, marks_list])
  return(count_marks.m)
}


dir_home <- '/project2/mstephens/old_project/'
cell_type_list <- c('K562', 'GM12878', 'HUVEC', 'NHLF')
marks_list <- c('CTCF', 'H3K27me3','H3K36me3', 'H4K20me1', 'H3K4me1', 'H3K4me2', 'H3K4me3','H3K27ac', 'H3K9ac', 'Control')
col_marks <- c('blue',  'black', rep('darkgreen', 2), rep('red', 3), rep('orange', 2), 'darkgray')
count_marks_cells.l <- vector('list', length = length(cell_type_list))
names(count_marks_cells.l) <- cell_type_list
cat('cell types:', cell_type_list, '\n')
for(cell_type in cell_type_list){
  # cat(cell_type, '\n')
  count_marks.m <- load_marks_data(cell_type, 'chr1', marks_list)
  count_marks_cells.l[[cell_type]] <- count_marks.m[c(340361:344860),] # chr1:68,072,000-68,972,000
}
count_marks_cells.m <- do.call(rbind, count_marks_cells.l)
data.m <- round(count_marks_cells.m)
rownames(data.m) <- 1:nrow(data.m)

library(stm)

fit = stm(t(data.m),K=10,maxiter = 5,printevery = 1,nugget = TRUE)

load('/home/dyxie/fit.RData')

data.m2 = data.m[-which(rowSums(data.m)==0),]
dim(data.m2)
fit = stm(t(data.m),K=7,maxiter = 5,printevery = 1,nugget = TRUE)

lf = poisson2multinom(t(fit$qf$Ef),fit$ql$El)

barplot(t(lf$L),col=2:(7+1),axisnames = F, space = 0, border = NA, las = 1, ylim = c(0, 1), cex.axis = 1.5, cex.main = 1.4)

plot(lf$FF[,1],type = 'h')
plot(lf$FF[,2],type = 'h')
plot(lf$FF[,3],type = 'h')
plot(lf$FF[,4],type = 'h')
plot(lf$FF[,5],type = 'h')
plot(lf$FF[,6],type = 'h')
plot(lf$FF[,7],type = 'h')

