if(!require("pd.mogene.1.1.st.v1")){
    source("http://bioconductor.org/biocLite.R")
    biocLite("pd.mogene.1.1.st.v1")
}
if(!require("oligo")){
    source("http://bioconductor.org/biocLite.R")
    biocLite("oligo")
}
if( !require("pd.mogene.1.1.nosprfvbsnps") ){
    install.packages('results/MoGene_FVB_SPR/pd.mogene.1.1.nosprfvbsnps', repos = NULL, type='source')
}
library("pd.mogene.1.1.st.v1")
library("pd.mogene.1.1.nosprfvbsnps")
library(oligo)

load.matrix = function( fn ){
    read.table(file=fn, sep="\t",row.names=1, check.names=F, header=T, na.strings=c('NA', '-99999'), stringsAsFactors=F)
}
write.matrix = function(o, fn_out, na.string="NA"){
    # Write data frame in standard format used by load.matrix
    rownames = append( 'IDENTIFIER', names(o) )
    write( rownames, fn_out, sep='\t', ncolumns=length(rownames) )
    write.table(o, fn_out, quote=F, sep='\t', row.names=T, na=na.string, col.names=F, append=T)
}
    
get.split.col = function(v, string, col=0, last=F, first=F){
    if( last & first )
        stop("Cannot request both last and first column")
    if( col==0 & !last & !first)
        stop("Must request either a column by index, first, or last")
        
    for(i in 1:length(v)){
        x = strsplit( v[i], string, fixed=T)[[1]]
        if(last){
            v[i] = x[length(x)]
        }
        else if(first){
            v[i] = x[1]
        }
        else{
            v[i] = x[col]
        }
    }
    v
}

match.idx = function(A, B, allow.multiple.B=F){
    # return dataframe of indices into A and B restricted to perfect matches
    # between A and B, where idx.A[i] == idx.B[i] for each i in matched pairs
    if( allow.multiple.B ){
        idx.B = which(B %in% A)
        idx.A = match(B[idx.B], A)
    }
    else{
        in.both = intersect(A,B)
        idx.A = match(in.both, A)
        idx.B = match(in.both, B)
    }
    C= data.frame(idx.A, idx.B)
    if( sum( A[ C$idx.A ] != B[ C$idx.B] )>0 )
        stop("ERROR! At least one in idx.A not the same as matched item in idx.B")
    C
}

cT <- list.celfiles("data/CEL", full.names = TRUE, listGzipped=TRUE)
FS <- read.celfiles(cT, pkgname = "pd.mogene.1.1.st.v1")
expr = data.frame(round( exprs(rma(FS) ), 3))
ga = load.matrix('data/annotation/MoGene-1_1/gene_attributes_MoGene-1_1-st-v1.MM10.txt')
ga = ga[ga$probe.type=="main",]
m = match.idx(rownames(ga), rownames(expr))
ga = ga[m$idx.A,]
expr = expr[m$idx.B,]
sa = load.matrix('data/sample_attributes_mammary.txt')
nn = get.split.col( get.split.col(names(expr), "_", last=TRUE), ".", first=TRUE)
nn = paste("RU109", nn, "mammary", sep="_")
names(expr) = nn
m = match.idx(rownames(sa), names(expr))
expr = expr[,m$idx.B]
keep = rowMeans(data.matrix(expr))>5 & ga$is.refseq=="yes"
# write all high-quality, above.bg probes
write.matrix( expr[ keep,], 'results/expression/mammary/MoGene_expr_above_bg_refseq.txt')
write.matrix( ga[ keep,], 'results/expression/mammary/MoGene_gene_attributes_above_bg_refseq.txt')

rm(expr)
rm(ga)
rm(FS)

#------------------------------------------------------------
# MoGene Mammary corrected annotation
#------------------------------------------------------------
cT <- list.celfiles("data/CEL", full.names = TRUE, listGzipped=TRUE)
FS <- read.celfiles(cT, pkgname = "pd.mogene.1.1.nosprfvbsnps")
expr = data.frame(round( exprs(rma(FS) ), 3))
ga = load.matrix('data/annotation/MoGene-1_1/gene_attributes_MoGene-1_1-st-v1.MM10.txt')
ga = ga[ga$probe.type=="main",]
m = match.idx(rownames(ga), rownames(expr))
ga = ga[m$idx.A,]
expr = expr[m$idx.B,]
sa = load.matrix('data/sample_attributes_mammary.txt')
nn = get.split.col( get.split.col(names(expr), "_", last=TRUE), ".", first=TRUE)
nn = paste("RU109", nn, "mammary", sep="_")
names(expr) = nn
m = match.idx(rownames(sa), names(expr))
expr = expr[,m$idx.B]
keep = rowMeans(data.matrix(expr))>5 & ga$is.refseq=="yes"
# write all high-quality, above.bg probes
write.matrix( expr[ keep,], 'results/expression/mammary/MoGene_nosprfvbsnps_expr_above_bg_refseq.txt')
write.matrix( ga[ keep,], 'results/expression/mammary/MoGene_nosprfvbsnps_gene_attributes_above_bg_refseq.txt')



