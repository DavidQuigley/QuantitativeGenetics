print("equalizer example code")
print("usage: equalizer_example.R <dir_bedtools_bin>")
print(" <dir_bedtools_bin> path to bedtools2/bin")

#-------------------------------------------------------------------------------
# Functions
#-------------------------------------------------------------------------------

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
# accept arguments
# If you are running this code line by line rather than as a script, 
# set the values of dir_bedtools_bin, fn_equalizer_source, and 
# ANNOT by hand.
#----------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if( length(args) != 1 ){
    stop("One argument is required: dir_bedtools_bin")
}
dir_bedtools_bin = args[1]

initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
ANNOT = paste(dirname(script.name) , "data/annotation", sep='/')
fn_equalizer_source = paste(dirname(script.name) , "bin/equalizer_0.3.tar.gz", sep='/')
dir_output = paste(dirname(script.name) , "results", sep='/')
if(!file.exists(ANNOT))
    stop(paste("Cannot find dir_annotation_root at", ANNOT))
if(!file.exists(dir_bedtools_bin))
    stop(paste("Cannot find dir_bedtools_bin at", dir_bedtools_bin))
if(!file.exists(fn_equalizer_source))
    stop(paste("Cannot find fn_equalizer_source at", fn_equalizer_source))
if(!file.exists(dir_output))
    stop(paste("Cannot find dir_output at", dir_output))

dir_CEL = paste(dirname(script.name), "/data/CEL", sep="")
fn_sa = paste(dirname(script.name) , "/data/sample_attributes_mammary.txt", sep='')
dir_geno = paste(dirname(script.name), "/data/genotypes", sep="")

#----------------------------------------------------------------
# Libraries
#----------------------------------------------------------------

if(!require("findpython")){
    install.packages("findpython", repos="http://cran.cnr.Berkeley.edu/")   
}
library(findpython)
fn_python = find_python_cmd()
if(!file.exists(fn_python))
    stop(paste("Cannot find python at", fn_python))

if(!require(pdInfoBuilder)){
    source("http://bioconductor.org/biocLite.R")
    biocLite("pdInfoBuilder")
}
library(pdInfoBuilder)

if(!require(pd.mogene.1.1.st.v1)){
    source("http://bioconductor.org/biocLite.R")
    biocLite("pd.mogene.1.1.st.v1")
}
library(pd.mogene.1.1.st.v1)

if(!require(equalizer)){
    install.packages(fn_equalizer_source)
}
library(equalizer)



#----------------------------------------------------------------
# Call equalizer
#----------------------------------------------------------------
chip = "MoGene-1_1-st-v1"
chip_version = "1.1"
package_name="mogene-1_1_nosprfvbSnps"
fn_vcf_list = c( paste(ANNOT,"/cdc26.vcf", sep="") )
affy_bed=paste(ANNOT, "/probe_bed/MoGene_probes_MM10.bed", sep="")
annot_version="na34"
annot_revision="1.1"
genome_version_UCSC="mm10"
genome_version_NCBI="38"
species="Mus musculus"
organism="Mouse"
author="David Quigley"
email="David.Quigley@ucsf.edu"
fn_probeset_csv = paste(ANNOT, "/MoGene-1_1/MoGene-1_1-st-v1.na34.mm10.probeset.csv", sep="")
fn_transcript_csv = paste(ANNOT, "/MoGene-1_1/MoGene-1_1-st-v1.na34.mm10.transcript.csv", sep="")
fn_pgf = paste(ANNOT, "/MoGene-1_1/MoGene-1_1-st-v1.r4.pgf", sep="")
fn_mps = paste(ANNOT, "/MoGene-1_1/MoGene-1_1-st-v1.r4.mps", sep="")
fn_clf = paste(ANNOT, "/MoGene-1_1/MoGene-1_1-st-v1.r4.clf", sep="")

equalize_gene( fn_python, dir_bedtools_bin,
               package_name, fn_vcf_list, affy_bed,
               annot_version, annot_revision,
               genome_version_UCSC, genome_version_NCBI,
               chip, chip_version, species, organism,
               author, email,
               fn_probeset_csv, fn_transcript_csv, fn_pgf, fn_mps, fn_clf,
               paste(dir_output, "/MoGene_FVB_SPR", sep=""))
#----------------------------------------------------------------
# Build new package with equalized chip definition
#----------------------------------------------------------------
fn_pgf_nosnps = paste(dir_output, '/MoGene_FVB_SPR/mogene-1_1_nosprfvbSnps.pgf',sep="")
fn_mps_nosnps = paste(dir_output, '/MoGene_FVB_SPR/mogene-1_1_nosprfvbSnps.mps',sep="")
fn_probes_nosnps = paste(dir_output, '/MoGene_FVB_SPR/mogene-1_1_nosprfvbSnps.probeset.csv',sep="")
fn_trans_nosnps = paste(dir_output, '/MoGene_FVB_SPR/mogene-1_1_nosprfvbSnps.transcript.csv',sep="")

seed <- new('AffyGenePDInfoPkgSeed', chipName = package_name, version = chip_version, 
             pgfFile = fn_pgf_nosnps, clfFile = fn_clf, coreMps = fn_mps_nosnps, 
             transFile = fn_trans_nosnps, probeFile = fn_probes_nosnps, 
             author = author, email = email, biocViews = 'AnnotationData', 
             genomebuild = genome_version_UCSC, organism = organism, species = species, 
             url ='')

dir_new_package=paste(dir_output, "/MoGene_FVB_SPR/pd.mogene.1.1.nosprfvbsnps", sep="")
print(paste("Installing to", dir_new_package))
unlink( dir_new_package, recursive = TRUE)
makePdInfoPackage(seed, destDir = paste(dir_output, "/MoGene_FVB_SPR", sep="") , unlink=TRUE)
install.packages( dir_new_package, repos = NULL, type='source')
library("pd.mogene.1.1.nosprfvbsnps")

#----------------------------------------------------------------
# Normalize microarrays
#----------------------------------------------------------------

cT <- list.celfiles( dir_CEL, full.names = TRUE, listGzipped=TRUE)
FS <- read.celfiles(cT, pkgname = "pd.mogene.1.1.st.v1")
expr = data.frame(round( exprs(rma(FS) ), 3))
ga = load.matrix( paste(ANNOT, "/gene_attributes_MoGene-1_1-st-v1.MM10.txt", sep="") )
ga = ga[ga$probe.type=="main",]
m = match.idx(rownames(ga), rownames(expr))
ga = ga[m$idx.A,]
expr = expr[m$idx.B,]
sa = load.matrix( fn_sa )
nn = get.split.col( get.split.col(names(expr), "_", last=TRUE), ".", first=TRUE)
nn = paste("RU109", nn, "mammary", sep="_")
names(expr) = nn
m = match.idx(rownames(sa), names(expr))
expr = expr[,m$idx.B]
keep = "10513583"
write.matrix( expr[ keep,], paste( dir_output, '/expression/mammary/MoGene_expr_above_bg_refseq.txt', sep=""))
write.matrix( ga[ keep,], paste( dir_output, '/expression/mammary/MoGene_gene_attributes_above_bg_refseq.txt', sep=""))

rm(expr)
rm(ga)
rm(FS)
rm(m)
gc()

#------------------------------------------------------------
# MoGene Mammary corrected annotation
#------------------------------------------------------------
cT <- list.celfiles( dir_CEL, full.names = TRUE, listGzipped=TRUE)
FS <- read.celfiles(cT, pkgname = "pd.mogene.1.1.nosprfvbsnps")
expr = data.frame(round( exprs(rma(FS) ), 3))
ga = load.matrix( paste(ANNOT, "/gene_attributes_MoGene-1_1-st-v1.MM10.txt", sep="") )
ga = ga[ga$probe.type=="main",]
m = match.idx(rownames(ga), rownames(expr))
ga = ga[m$idx.A,]
expr = expr[m$idx.B,]
sa = load.matrix( fn_sa )
nn = get.split.col( get.split.col(names(expr), "_", last=TRUE), ".", first=TRUE)
nn = paste("RU109", nn, "mammary", sep="_")
names(expr) = nn
m = match.idx(rownames(sa), names(expr))
expr = expr[,m$idx.B]
keep = "10513583"
write.matrix( expr[ keep,], paste( dir_output, '/expression/mammary/MoGene_nosprfvbsnps_expr_above_bg_refseq.txt', sep=""))
write.matrix( ga[ keep,], paste( dir_output, '/expression/mammary/MoGene_nosprfvbsnps_gene_attributes_above_bg_refseq.txt', sep=""))
