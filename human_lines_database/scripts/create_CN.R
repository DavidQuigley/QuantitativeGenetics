# This script integrates CCLE and Sanger copy number calls into a single file for major 
# and minor copy number. It does not directly calculate anything. In this iteration 
# CCLE values are prefered

source('/notebook/code/src/R/quantitative_genetics.R')
cl = load.matrix('/datasets/human_lines_database/cell_line_attributes/cell_line_dictionary_2015_02_18.txt')
cn.ccle.minor = load.matrix('/datasets/human_lines_database/copy_number/CN_minor_CCLE.txt.gz')
cn.ccle.total = load.matrix('/datasets/human_lines_database/copy_number/CN_total_CCLE.txt.gz')
cn.sanger.minor = load.matrix('/datasets/human_lines_database/copy_number/CN_minor_sanger.txt.gz')
cn.sanger.total = load.matrix('/datasets/human_lines_database/copy_number/CN_total_sanger.txt.gz')

#cc = load.matrix('/datasets/human_lines_database/copy_number/CN_total_2015_03_03.txt')

m = match.idx(names(cn.sanger.total), names(cn.ccle.total) )
rhos = rep(0, dim(m)[1])
for(i in 1:dim(m)[1]){
    rhos[i] = as.numeric(cor.test( cn.sanger.total[,m$idx.A[i]], cn.ccle.total[,m$idx.B[i]], method="pearson")$estimate)
}
n_agree = rep(0, dim(m)[1])
for(i in 1:dim(m)[1]){
    n_agree[i] = sum( cn.sanger.total[,m$idx.A[i]] == cn.ccle.total[ , m$idx.B[i] ], na.rm=TRUE)
}
names( cn.sanger.total)[ m[ which(rhos<0.5) ,]$idx.A ]

confusion = matrix(0, nrow=15, ncol=15)
dimnames(confusion)[[1]] = paste("CC_", 0:14, sep="")
dimnames(confusion)[[2]] = paste("SG_", 0:14, sep="")
idx_ccle=384
idx_sanger=385
for(i in 1:dim(cn.ccle.total)[1]){
    rr = cn.ccle.total[i,idx_ccle] + 1
    cc = cn.sanger.total[i,idx_sanger] + 1
    confusion[ rr, cc ] = confusion[ rr, cc ] + 1
}
print(paste( "CCLE",idx_ccle, names(cn.ccle.total)[idx_ccle], 
      ", Sanger", idx_sanger, names(cn.sanger.total)[idx_sanger] ) )
confusion
confusion = matrix(0, nrow=15, ncol=15)
cn.to.test = cn.sanger.total[,idx_sanger]
n_agree_test = rep(NA, dim(cn.ccle.total)[2])
for(i in 1:dim(cn.ccle.total)[2]){
    n_agree_test[i] = sum( cn.to.test == cn.ccle.total[ , i ], na.rm=TRUE)
}
idx_ccle_best_fit = which(n_agree_test==max(n_agree_test))
for(i in 1:dim(cn.ccle.total)[1]){
    rr = cn.ccle.total[i,idx_ccle_best_fit] + 1
    cc = cn.sanger.total[i,idx_sanger] + 1
    confusion[ rr, cc ] = confusion[ rr, cc ] + 1
}
print(paste( "CCLE best fit",idx_ccle_best_fit, names(cn.ccle.total)[idx_ccle_best_fit], 
      ", Sanger", idx_sanger, names(cn.sanger.total)[idx_sanger] ) )
confusion

# Use CCLE copy number call where available; otherwise Sanger
only_sanger = setdiff( names(cn.sanger.total), names(cn.ccle.total) )
only_ccle = setdiff(names(cn.ccle.total),  names(cn.sanger.total) )
m = match.idx(names(cn.sanger.total), only_sanger)
cn.total = data.matrix(cbind( cn.ccle.total, cn.sanger.total[,m$idx.A]))
cn.minor = data.matrix(cbind( cn.ccle.minor, cn.sanger.minor[,m$idx.A]))

cn.total.final = matrix(NA, nrow=dim(cn.total)[1], ncol=dim(cl)[1])
cn.minor.final = matrix(NA, nrow=dim(cn.minor)[1], ncol=dim(cl)[1])

m = match.idx( dimnames(cl)[[1]], dimnames(cn.total)[[2]])
cn.total.final[,m$idx.A] = cn.total[,m$idx.B]
cn.minor.final[,m$idx.A] = cn.minor[,m$idx.B]

dimnames(cn.total.final)[[1]] = dimnames(cn.total)[[1]]
dimnames(cn.total.final)[[2]] = dimnames(cl)[[1]]
dimnames(cn.minor.final)[[1]] = dimnames(cn.minor)[[1]]
dimnames(cn.minor.final)[[2]] = dimnames(cl)[[1]]

cn.total.final = data.frame(cn.total.final)
cn.minor.final = data.frame(cn.minor.final)
names(cn.total.final) = dimnames(cl)[[1]]
names(cn.minor.final) = dimnames(cl)[[1]]

date=format( Sys.Date(), "%Y_%m_%d")

write.matrix( cn.minor.final, paste('/datasets/human_lines_database/copy_number/CN_minor_',date,'.txt', sep=""))
write.matrix( cn.total.final, paste('/datasets/human_lines_database/copy_number/CN_total_',date,'.txt',sep=""))