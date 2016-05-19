print("equalizer example code")
print("usage: equalizer_minimal_example.R <dir_bedtools_bin>")
print(" <dir_bedtools_bin> path to bedtools2/bin")

#----------------------------------------------------------------
# accept arguments
# If you are running this code line by line rather than as a script, 
# set the values of dir_bedtools_bin and 
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
fn_equalizer_source = paste(dirname(script.name) , "bin/equalizer_0.3.1.tar.gz", sep='/')
dir_output = paste(dirname(script.name) , "results", sep='/')
if(!file.exists(ANNOT))
    stop(paste("Cannot find dir_annotation_root at", ANNOT))
if(!file.exists(dir_bedtools_bin))
    stop(paste("Cannot find dir_bedtools_bin at", dir_bedtools_bin))
if(!file.exists(fn_equalizer_source))
    stop(paste("Cannot find fn_equalizer_source at", fn_equalizer_source))
if(!file.exists(dir_output))
    stop(paste("Cannot find dir_output at", dir_output))

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
               dir_output)

probes = read.table( paste( dir_output, '/probe_count_changes.txt', sep=""), 
                    stringsAsFactors=FALSE, header=TRUE)

cdc26 = probes[probes$transcript.id=="10513583",]
print(paste("Original number of probes for Cdc26:", cdc26$orig.n.probes[1]))
print(paste("Probes remaining after equalizer (should be 3):", 
            sum(cdc26$final.probeset.n.probes)))
