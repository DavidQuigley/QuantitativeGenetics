# This script integrates CCLE and Sanger copy number calls into a single file for major 
# and minor copy number. It does not directly calculate anything. In this iteration 
# CCLE values are prefered because we do not yet have access to the raw Sanger data.

source('/notebook/code/quantitative_genetics.R')
cl = load.matrix('/datasets/human_lines_database/cell_line_attributes/cell_line_dictionary_2015_02_18.txt')
cn.ccle.minor = load.matrix('/datasets/human_lines_database/copy_number/CN_minor_CCLE.txt.gz')
cn.ccle.total = load.matrix('/datasets/human_lines_database/copy_number/CN_total_CCLE.txt.gz')
cn.sanger.minor = load.matrix('/datasets/human_lines_database/copy_number/CN_minor_sanger.txt.gz')
cn.sanger.total = load.matrix('/datasets/human_lines_database/copy_number/CN_total_sanger.txt.gz')

only_sanger = setdiff( names(cn.sanger.total), names(cn.ccle.total) )
m = match.idx(names(cn.sanger.total), only_sanger)
cn.ccle.total = data.matrix(cbind( cn.ccle.total, cn.sanger.total[,m$idx.A]))
cn.ccle.minor = data.matrix(cbind( cn.ccle.minor, cn.sanger.minor[,m$idx.A]))

cn.ccle.total.final = matrix(NA, nrow=dim(cn.ccle.total)[1], ncol=dim(cl)[1])
cn.ccle.minor.final = matrix(NA, nrow=dim(cn.ccle.minor)[1], ncol=dim(cl)[1])

m = match.idx( dimnames(cl)[[1]], dimnames(cn.ccle.total)[[2]])
cn.ccle.total.final[,m$idx.A] = cn.ccle.total[,m$idx.B]
cn.ccle.minor.final[,m$idx.A] = cn.ccle.minor[,m$idx.B]

dimnames(cn.ccle.total.final)[[1]] = dimnames(cn.ccle.total)[[1]]
dimnames(cn.ccle.total.final)[[2]] = dimnames(cl)[[1]]
dimnames(cn.ccle.minor.final)[[1]] = dimnames(cn.ccle.minor)[[1]]
dimnames(cn.ccle.minor.final)[[2]] = dimnames(cl)[[1]]

cn.ccle.total.final = data.frame(cn.ccle.total.final)
cn.ccle.minor.final = data.frame(cn.ccle.minor.final)
names(cn.ccle.total.final) = dimnames(cl)[[1]]
names(cn.ccle.minor.final) = dimnames(cl)[[1]]

date=format( Sys.Date(), "%Y_%m_%d")

write.matrix( cn.ccle.minor.final, paste('/datasets/human_lines_database/copy_number/CN_minor_',date,'.txt', sep=""))
write.matrix( cn.ccle.total.final, paste('/datasets/human_lines_database/copy_number/CN_total_',date,'.txt',sep=""))