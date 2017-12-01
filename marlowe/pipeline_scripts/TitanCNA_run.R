suppressPackageStartupMessages(
    require(TitanCNA, quietly=TRUE, warn.conflicts=FALSE)
)
suppressPackageStartupMessages(
    require(HMMcopy, quietly=TRUE, warn.conflicts=FALSE)
)
suppressPackageStartupMessages(
    require(doMC, quietly=TRUE, warn.conflicts=FALSE)
)
suppressPackageStartupMessages(
    require(optparse, quietly=TRUE, warn.conflicts=FALSE)
)

option_list <- list(
	make_option(c("--sample_id"), type = "character",
            help = "Sample identifier. [Required]"),
	make_option(c("--fn_het_counts"), type = "character",
            help = "File containing allelic read counts at germline heterozygous sites. [Required]"),
	make_option(c("--fn_wig_normal"), type = "character",
            help = "BigWig of normal read depth. [Required]"),
	make_option(c("--fn_wig_tumor"), type = "character",
            help = "BigWig of tumor read depth. [Required]"),
	make_option(c("--fn_wig_gc"), type = "character",
            help = "BigWig of gc percentage genome-wide. [Required]"),
	make_option(c("--fn_wig_mappable"), type = "character",
            help = "BigWig of mappable regions. [Required]"),
	make_option(c("--dir_out"), type = "character",
            help = "directory where results are written. [Required]"),
	make_option(c("--n_clusters"), type = "character", default=1,
            help = "number of clones, default [Default: %default]"),
    make_option(c("--ploidy"), type = "numeric", default = 2,
            help = "Initial ploidy value; float [Default: %default]"),
    make_option(c("--genomestyle"), type = "character", default = "UCSC",
            help = "genome style [Default: %default]"),            
    make_option(c("--n_iter"), type = "numeric", default = 3,
            help = "Number of EM iterations; numeric [Default: %default]")
)
parseobj = OptionParser(option_list=option_list, usage = "usage: Rscript %prog [options]")
opt = parse_args(parseobj)

sample_id = opt$sample_id
fn_het_counts = opt$fn_het_counts
fn_wig_tum = opt$fn_wig_tumor
fn_wig_nor = opt$fn_wig_normal
fn_wig_gc = opt$fn_wig_gc
fn_wig_map_counter = opt$fn_wig_mappable
dir_out = opt$dir_out
ploidy = opt$ploidy
n_clusters = opt$n_clusters
n_iter = opt$n_iter
genomestyle = opt$genomestyle

fn_results = paste( dir_out, "/", sample_id, "_results.txt", sep="")
fn_segs_out = paste( dir_out, "/", sample_id, "_segments.txt", sep="")
fn_igv_seg_out = paste( dir_out, "/", sample_id, "_segments.seg", sep="")
fn_model_parameters = paste( dir_out, "/", sample_id, "_model_parameters.txt", sep="")

message("Loading allele counts")
data = loadAlleleCounts(fn_het_counts, 
                        genomeStyle = genomestyle)
message("Correcting read depth")
cnData = correctReadDepth(fn_wig_tum, 
                          fn_wig_nor, 
                          fn_wig_gc, 
                          fn_wig_map_counter, 
                          genomeStyle=genomestyle)
logR = getPositionOverlap(data$chr, 
                          data$posn, 
                          cnData) 
data$logR <- log(2^logR) 
chromosomes_permitted = paste( "chr", c(1:22, "X", "Y"), sep='')
data = filterData(data, chromosomes_permitted,                 
                  minDepth = 10, 
                  maxDepth = 2000,
                  positionList = NULL, 
                  centromere = NULL, 
                  centromere.flankLength = 10000)
params = loadDefaultParameters(copyNumber = 8, 
                               numberClonalClusters = n_clusters, 
                               symmetric = TRUE, 
                               hetBaselineSkew = 0, 
                               alleleEmissionModel = "binomial", 
                               data = data) 
#params$ploidyParams$phi_0 
K = length(params$genotypeParams$alphaKHyper)
params$genotypeParams$alphaKHyper = rep(500, K)
params$ploidyParams$phi_0 = 1.5
message("Running expectation maximization")
convergeParams = runEMclonalCN(data, 
                               params=params,
                               maxiter = n_iter, 
                               maxiterUpdate = 50,
                               useOutlierState = FALSE, 
                               txnExpLen = 1e15,
                               txnZstrength = 5e5,
                               normalEstimateMethod = "map",
                               estimateS = TRUE, 
                               estimatePloidy = TRUE)
message("Computing optimal path")
optimalPath = viterbiClonalCN(data, convergeParams)
message("Computing results")
results = outputTitanResults(data, 
                             convergeParams, 
                             optimalPath,
                             filename = NULL, 
                             posteriorProbs = FALSE,
                             subcloneProfiles = TRUE,
                             correctResults = TRUE, 
                             proportionThreshold = 0.05, 
                             proportionThresholdClonal = 0.05,
                             is.haplotypeData = FALSE)
convergeParams = results$convergeParams
results = results$corrResults

# Convert chromosomes to numeric values, required for plotting functions
results_n = results
for(i in 1:22){
    cc = paste("chr",i,sep="")
    results_n$Chr[results_n$Chr==cc] = i   
}
results_n$Chr[results_n$Chr=="chrX"] = "X"
results_n$Chr[results_n$Chr=="chrY"] = "Y"
results = results_n

message(paste("Writing results to",dir_out))
segs = outputTitanSegments(results, 
                           id = sample_id, 
                           convergeParams, 
                           filename = fn_segs_out, 
                           igvfilename = fn_igv_seg_out)
write.table(results, file = fn_results, col.names = TRUE, row.names = FALSE, 
            quote = FALSE, sep = "\t")
outputModelParameters(convergeParams, results, fn_model_parameters) 

