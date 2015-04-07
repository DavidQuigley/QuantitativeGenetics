# CCLE and Sanger have good correlation overall
# For cases where we have an array from both, we prefer CCLE because they used the 
# HG-U133 Plus 2 array, which covers about 9300 additional genes

library(matrixStats)
source('/notebook/code/quantitative_genetics.R')
fn_expr_gray = '/datasets/human_lines_gray/expression/expr_heiser_breast_cell_lines.txt'
fn_expr_ccle = '/datasets/human_lines_ccle/CCLE_Expression_2012-09-29.txt'
fn_expr_sanger = '/datasets/human_lines_sanger/expU133A_correct.txt'
fn_ga_ccle = '/notebook/annotations/gene_attributes_HG-U133_Plus_2.na34.txt'
fn_ga_sanger = '/notebook/annotations/gene_attributes_HG-U133A.na34.txt'
fn_ga_gray = '/datasets/human_lines_gray/expression/gene_attributes_heiser_breast_cell_lines.txt'
fn_database = '/datasets/human_lines_database/cell_line_attributes/cell_line_dictionary_2015_02_18.txt'
fn_genes = '/datasets/human_lines_database/annotations/refseq_genes_hg19_2015_02_20.BED'

genes = read.table(fn_genes, header=FALSE, sep='\t', stringsAsFactors=FALSE)
names(genes) = c("chrom", "txStart", "txEnd", "symbol")

expr_ccle_raw = load.matrix(fn_expr_ccle)
expr_sanger_raw = load.matrix(fn_expr_sanger)

# Sanger and CCLE contain technical replicates with identical names.
# pick an arbitrary replicate, as they are very strongly correlated
expr_sanger_raw = expr_sanger_raw[, match.idx( unique(names(expr_sanger_raw)), names(expr_sanger_raw))$idx.A]
expr_ccle_raw = expr_ccle_raw[, match.idx( unique(names(expr_ccle_raw)), names(expr_ccle_raw))$idx.A]

# match sanger and CCLE probes with annotation; 
# some probes have lost annotation since arrays were designed.
ga_sanger = load.matrix(fn_ga_sanger)
ga_ccle = load.matrix(fn_ga_ccle)
m = match.idx(rownames(ga_sanger), rownames(expr_sanger_raw))
ga_sanger = ga_sanger[m$idx.A,]
expr_sanger_raw = expr_sanger_raw[m$idx.B,]
m = match.idx(rownames(ga_ccle), rownames(expr_ccle_raw))
ga_ccle = ga_ccle[m$idx.A,]
expr_ccle_raw = expr_ccle_raw[m$idx.B,]

expr_gray = load.matrix(fn_expr_gray)
ga_gray = load.matrix(fn_ga_gray)

DC = data.matrix(expr_ccle_raw)
DS = data.matrix(expr_sanger_raw)
DG = data.matrix(expr_gray)
    
# Choose single best probe for each symbol.
# prefer to use probes that are present in both CCLE and sanger.
probes_in_both = intersect(rownames(ga_ccle), rownames(ga_sanger))
#probes_in_both = hsh_from_vectors(probes_in_both, 1:length(probes_in_both))

symbols.c = sort( unique( ga_ccle$symbol ))
probes.c = rep("", length(symbols.c))
DC.compact = matrix(0, nrow=length(symbols.c), ncol=dim(expr_ccle_raw)[2])
for(i in 1:length(symbols.c)){
    probes = rownames(ga_ccle)[ which(ga_ccle$symbol==symbols.c[i]) ]
    probes_shared_ccle_sanger = intersect(probes, probes_in_both)
    if( length(probes_shared_ccle_sanger)>=1 )
        probes = probes_shared_ccle_sanger
    idx = which( rownames( ga_ccle ) %in% probes )
    ee = DC[idx,]
    if( length(idx)>1 ){
        mus = rowMeans(ee, na.rm=TRUE)
        idx_max_mu = which(mus==max(mus))
        ee = ee[idx_max_mu,]
        probes.c[i] = names(idx_max_mu)
    }
    else{
        probes.c[i] = rownames(ga_ccle)[ idx ]
    }
    DC.compact[i,] = ee
}   
s2p = hsh_from_vectors( symbols.c, probes.c)
rownames(DC.compact) = symbols.c
dimnames(DC.compact)[[2]] = names(expr_ccle_raw)
DC = DC.compact
ga_ccle = ga_ccle[match.idx(symbols.c, ga_ccle$symbol)$idx.B,]
rownames(ga_ccle) = rownames(DC)

symbols.s = sort( unique( ga_sanger$symbol ))
DS.compact = matrix(0, nrow=length(symbols.s), ncol=dim(expr_sanger_raw)[2])
probes_sanger = rep("", dim(DS.compact)[1])
for(i in 1:length(symbols.s)){
    # if probe already picked for CCLE, use that one
    probe_for_this_symbol_ccle = hsh_get( s2p, symbols.s[i], na.if.not.found=TRUE ) 
    idx_sanger = c()
    if( !is.na( probe_for_this_symbol_ccle ) ){
        idx_sanger = which( rownames(ga_sanger)==probe_for_this_symbol_ccle)
    }
    if( length(idx_sanger)==1){
        DS.compact[i,] = DS[idx_sanger,]
        probes_sanger[i] = probe_for_this_symbol_ccle
    }
    else{
        idx = which(ga_sanger$symbol==symbols.s[i])
        ee = DS[idx,]
        if( length(idx)>1 ){
            mus = rowMeans(ee, na.rm=TRUE)
            ee = ee[mus == max(mus),]
            idx = idx[mus==max(mus)]
        }
        probes_sanger[i] = rownames(ga_sanger)[idx]
        DS.compact[i,] = ee
    }
}   
rownames(DS.compact) = symbols.s
dimnames(DS.compact)[[2]] = names(expr_sanger_raw)
DS = DS.compact
ga_sanger = ga_sanger[match.idx(symbols.s, ga_sanger$symbol)$idx.B,]
rownames(ga_sanger) = rownames(DS)

# Choose single best probe for each symbol: GRAY

symbols.g = sort( unique( ga_gray$symbol ))
DG.compact = matrix(0, nrow=length(symbols.g), ncol=dim(expr_gray)[2])
for(i in 1:length(symbols.g)){
    idx = which(ga_gray$symbol==symbols.g[i])
    ee = DG[idx,]
    if( length(idx)>1 ){
        mus = rowMeans(ee, na.rm=TRUE)
        ee = ee[mus == max(mus),]
    }
    DG.compact[i,] = ee
}   
rownames(DG.compact) = symbols.g
dimnames(DG.compact)[[2]] = names(expr_gray)
DG = DG.compact
ga_gray = ga_gray[match.idx(symbols.g, ga_gray$symbol)$idx.B,]
rownames(ga_gray) = rownames(DG)


cl = load.matrix(fn_database)
m.ccle = match.idx(cl$name_ccle, names(expr_ccle_raw))
m.sanger = match.idx(cl$name_sanger, names(expr_sanger_raw))
m.gray = match.idx(cl$name_gray, names(expr_gray))

ids = sort(unique( c( rownames( cl )[m.ccle$idx.A], rownames( cl )[m.sanger$idx.A], 
                      rownames( cl )[m.gray$idx.A]) ) )
symbols = sort( unique( c( symbols.s, symbols.c, symbols.g ) ) )


m = match.idx(ids, rownames(cl) )
site = cl$site[m$idx.B]
site[is.na(site)] = "unknown"
site[ site %in% names( table(site))[ which( as.numeric(table(site))<10 ) ] ] = "other"

array_ccle = cl$id_ccle_expr[m$idx.B]
has_array_sanger = cl$has_expr_sanger[m$idx.B]
ids_gray = cl$name_gray[m$idx.B]
ids_ccle = cl$name_ccle[m$idx.B]
ids_sanger = cl$name_sanger[m$idx.B]

m_symbol_ccle = match.idx( symbols, symbols.c)
m_symbol_sanger= match.idx( symbols, symbols.s)
m_symbol_gray = match.idx( symbols, symbols.g)

D = matrix(NA, nrow=length(symbols), ncol=length(ids))
dimnames(D)[[1]] = symbols
dimnames(D)[[2]] = ids
array_source = rep("none", dim(D)[2])
for( i in 1:length(ids) ){
    id=ids[i]
    if( !is.na( array_ccle[i]  ) ){
        D[m_symbol_ccle$idx.A,i] = DC[m_symbol_ccle$idx.B, which( dimnames(DC)[[2]] == ids_ccle[i] ) ]
        array_source[i] = "ccle"
    }
    else if( !is.na( ids_gray[i]  ) ){
        D[m_symbol_gray$idx.A,i] = DG[m_symbol_gray$idx.B, which( dimnames(DG)[[2]] == ids_gray[i] ) ]
        array_source[i] = "gray"
    }
    else if( has_array_sanger[i] ){
        D[m_symbol_sanger$idx.A,i] = DS[m_symbol_sanger$idx.B, which( dimnames(DS)[[2]] == ids_sanger[i] ) ]
        array_source[i] = "sanger"
    }
}
nas = rowSums(is.na(D))
D = D[nas<1000,]
nas = rowSums(is.na(D))
table(nas)
#    0    16   251   267 
#11392  1846  4219  5178 

# identify the samples for each class of NA
example_idx_16 = which(nas==16)[1]
example_idx_250 = which(nas==250)[1]
example_idx_267 = which(nas==267)[1]

source('/notebook/code/ComBat.R')

D.adj = D

idx.nona = which(nas==0)
D.nona = D[idx.nona,]
saminfo = data.frame( SampleName=ids, Batch=array_source, site, row.names=ids)
e.adj = ComBat(D.nona, saminfo, write=F, covariates='all', par.prior=T, filter=F,  prior.plots=F)
D.adj[idx.nona,] = e.adj

# items with 16 missing are present in CCLE, Sanger but not genentech
idx.16 = which(nas==16)
idx.samples.16 = which( !is.na(D[example_idx_16,]) )
D.16 = D[idx.16, idx.samples.16]
array_source_16 = array_source[idx.samples.16] # CCLE and sanger
site_16 = site[idx.samples.16] # CCLE and sanger
saminfo.16 = data.frame( SampleName=dimnames(D.16)[[2]], Batch=array_source_16, site_16, 
                         row.names=dimnames(D.16)[[2]])
e.adj.16 = ComBat(D.16, saminfo.16, write=F, covariates='all', par.prior=T, filter=F, prior.plots=F)
D.adj[idx.16,idx.samples.16] = e.adj.16

# items with 251 missing are present in CCLE, genentech but not sanger
idx.250 = which(nas==250)
idx.samples.250 = which( !is.na(D[example_idx_250,]) )
D.250 = D[idx.250, idx.samples.250]
array_source_250 = array_source[idx.samples.250] # CCLE and sanger
site_250 = site[idx.samples.250] # CCLE and sanger
saminfo.250 = data.frame( SampleName=dimnames(D.250)[[2]], Batch=array_source_250, site_250, 
                         row.names=dimnames(D.250)[[2]])
e.adj.250 = ComBat(D.250, saminfo.250, write=F, covariates='all', par.prior=T, filter=F, prior.plots=F)
D.adj[idx.250,idx.samples.250] = e.adj.250

# no need to change values in 267

D.adj.std = standardize(D.adj)
D.adj = round(D.adj, 2)
D.adj.std = round(D.adj.std, 2)

# Adjust D.adj, D.adj.std to include the uniform complement of genes and samples 
D.adj.final = matrix(NA, nrow=dim(genes)[1], ncol=dim(cl)[1])
D.adj.std.final = matrix(NA, nrow=dim(genes)[1], ncol=dim(cl)[1])
idx_symbols_matched = match.idx( genes$symbol, dimnames(D.adj)[[1]] )
idx_samples_matched = match.idx( rownames(cl), dimnames(D.adj)[[2]] )
D.adj.final[idx_symbols_matched$idx.A, idx_samples_matched$idx.A] = 
           D.adj[idx_symbols_matched$idx.B, idx_samples_matched$idx.B]
D.adj.std.final[idx_symbols_matched$idx.A, idx_samples_matched$idx.A] = 
           D.adj.std[idx_symbols_matched$idx.B, idx_samples_matched$idx.B]

D.adj.final = data.frame(D.adj.final)
D.adj.std.final = data.frame(D.adj.std.final)

names(D.adj.final) = dimnames(cl)[[1]]
names(D.adj.std.final) = dimnames(cl)[[1]]

rownames(D.adj.final) = genes$symbol
rownames(D.adj.std.final) = genes$symbol

date=format( Sys.Date(), "%Y_%m_%d")

write.matrix( D.adj.final, paste('/datasets/human_lines_database/expression/expr_CCLE_Sanger_Gray_',date,'_ComBat.txt', sep=''))
write.matrix( D.adj.std.final, paste('/datasets/human_lines_database/expression/expr_CCLE_Sanger_Gray_',date,'_ComBat_standardized.txt', sep=''))

