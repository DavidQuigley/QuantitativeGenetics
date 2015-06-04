#equalizer_root = '/notebook/equalizer/deploy'
equalizer_root = getwd()

args <- commandArgs(trailingOnly = TRUE)
if( length(args) != 2 )
    stop("Two arguments are required: dir_bedtools_bin, fn_equalizer_source")
library(findpython)

fn_python = find_python_cmd()
dir_bedtools_bin = args[2]
fn_equalizer_source = args[3]

if(!file.exists(fn_python))
    stop(paste("Cannot find python at", fn_python))
if(!file.exists(dir_bedtools_bin))
    stop(paste("Cannot find dir_bedtools_bin at", dir_bedtools_bin))
if(!file.exists(fn_equalizer_source))
    stop(paste("Cannot find fn_equalizer_source at",fn_equalizer_source))

#----------------------------------------------------------------
# Functions
#----------------------------------------------------------------

load.matrix = function( fn ){
    read.table(file=fn, sep="\t",row.names=1, check.names=F, header=T, 
    na.strings=c('NA', '-99999'), stringsAsFactors=F)
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

#----------------------------------------------------------------
# Libraries
#----------------------------------------------------------------

if(!require("pd.mogene.1.1.st.v1")){
    source("http://bioconductor.org/biocLite.R")
    biocLite("pd.mogene.1.1.st.v1")
}
if(!require("oligo")){
    source("http://bioconductor.org/biocLite.R")
    biocLite("oligo")
}

if(!require("pdInfoBuilder")){
    source("http://bioconductor.org/biocLite.R")
    biocLite("pdInfoBuilder")
}

library( pdInfoBuilder )
library("pd.mogene.1.1.st.v1")
library(oligo)

install.packages(fn_equalizer_source)
library(equalizer)

#----------------------------------------------------------------
# Call equalizer
#----------------------------------------------------------------
chip = "MoGene-1_1-st-v1"
chip_version = "1.1"
package_name="mogene-1_1_nosprfvbSnps"
fn_vcf_list = c("cdc26.vcf")

ANNOT = "data/annotation/"
affy_bed=paste(ANNOT,"probe_bed/MoGene_probes_MM10.bed", sep="")
annot_version="na34"
annot_revision="1.1"
genome_version_UCSC="mm10"
genome_version_NCBI="38"
species="Mus musculus"
organism="Mouse"
author="David Quigley"
email="David.Quigley@ucsf.edu"
fn_probeset_csv = paste(ANNOT, "MoGene-1_1/MoGene-1_1-st-v1.na34.mm10.probeset.csv", sep="")
fn_transcript_csv = paste(ANNOT, "MoGene-1_1/MoGene-1_1-st-v1.na34.mm10.transcript.csv", sep="")
fn_pgf = paste(ANNOT, "MoGene-1_1/MoGene-1_1-st-v1.r4.pgf", sep="")
fn_mps = paste(ANNOT, "MoGene-1_1/MoGene-1_1-st-v1.r4.mps", sep="")
fn_clf = paste(ANNOT, "MoGene-1_1/MoGene-1_1-st-v1.r4.clf", sep="")
dir_out = "results/MoGene_FVB_SPR"

equalize_gene( fn_python, dir_bedtools_bin,
        package_name, fn_vcf_list, affy_bed,
        annot_version, annot_revision,
        genome_version_UCSC, genome_version_NCBI,
        chip, chip_version, species, organism,
        author, email,
        fn_probeset_csv, fn_transcript_csv, fn_pgf, fn_mps, fn_clf,
        dir_out)

#----------------------------------------------------------------
# Build new package with equalized chip definition
#----------------------------------------------------------------
fn_pgf_nosnps = 'results/MoGene_FVB_SPR/mogene-1_1_nosprfvbSnps.pgf'
fn_mps_nosnps = 'results/MoGene_FVB_SPR/mogene-1_1_nosprfvbSnps.mps'
fn_clf_nosnps = 'data/annotation/MoGene-1_1/MoGene-1_1-st-v1.r4.clf'
fn_probes_nosnps = 'results/MoGene_FVB_SPR/mogene-1_1_nosprfvbSnps.probeset.csv'
fn_trans_nosnps = 'results/MoGene_FVB_SPR/mogene-1_1_nosprfvbSnps.transcript.csv'

seed <- new('AffyGenePDInfoPkgSeed', chipName = package_name, version = chip_version, 
             pgfFile = fn_pgf_nosnps, clfFile = fn_clf_nosnps, coreMps = fn_mps_nosnps, 
             transFile = fn_trans_nosnps, probeFile = fn_probes_nosnps, 
             author = author, email = email, biocViews = 'AnnotationData', 
             genomebuild = genome_version_UCSC, organism = organism, species = species, 
             url ='')

unlink( paste(dir_out, package_name, sep="/pd.", collapse=""), recursive = TRUE)
makePdInfoPackage(seed, destDir = 'results/MoGene_FVB_SPR', unlink=TRUE)
install.packages( paste(dir_out, "pd.mogene.1.1.nosprfvbsnps", sep='/', collapse=""), repos = NULL, type='source')
library("pd.mogene.1.1.nosprfvbsnps")

#----------------------------------------------------------------
# Normalize microarrays
#----------------------------------------------------------------

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

system2( "bin/eqtl", c("-dresults/expression/mammary/MoGene_expr_above_bg_refseq.txt",
         "-fdata/sample_attributes_mammary.txt",
         "-gresults/expression/mammary/MoGene_gene_attributes_above_bg_refseq.txt",
         "-edata/genotypes/calls.txt",
         "-hdata/genotypes/sample_attributes.txt",
         "-sdata/genotypes/gene_attributes.txt",
         "-ysymbol",
         "-cmammary.id",
         "-n1000",
         "-r10513583",
         "-oresults/eqtl/eqtl_mammary_above_bg_refseq_1000.txt") )

system2( "bin/eqtl", c("-dresults/expression/mammary/MoGene_nosprfvbsnps_expr_above_bg_refseq.txt",
         "-fdata/sample_attributes_mammary.txt",
         "-gresults/expression/mammary/MoGene_nosprfvbsnps_gene_attributes_above_bg_refseq.txt",
         "-edata/genotypes/calls.txt",
         "-hdata/genotypes/sample_attributes.txt",
         "-sdata/genotypes/gene_attributes.txt",
         "-ysymbol",
         "-cmammary.id",
         "-n1000",
         "-r10513583",
         "-oresults/eqtl/eqtl_mammary_CORRECTED_above_bg_refseq_1000.txt") )

eqtl_un = read.table( "results/eqtl/eqtl_mammary_above_bg_refseq_1000.txt",
                      sep='\t', stringsAsFactors=F)
names(eqtl_un) = c('raw.p','symbol','probe','snp','perm.p','aa','ab','bb')
eqtl_cor = read.table( "results/eqtl/eqtl_mammary_CORRECTED_above_bg_refseq_1000.txt",
                      sep='\t', stringsAsFactors=F)
names(eqtl_cor) = c('raw.p','symbol','probe','snp','perm.p','aa','ab','bb')                      

print( eqtl_un$raw.p )  # 0.007616683
print( eqtl_cor$raw.p ) # 4.302593e-14 
