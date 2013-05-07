# Package: carmen
# Version: 0.1-1
# Date: 2011-04-23
# Title: Functions for Quantitative Genetics
# Author: David Quigley <dquigley@cc.ucsf.edu>
# Maintainer: David Quigley <dquigley@cc.ucsf.edu>
# Depends: R (>= 1.8.0), survival, aroma.affymetrix, affy
# Suggests:
# Description: Functions for Quantitative Genetics
# License: Unlimited
# URL: http://www.davidquigley.com
#
############################
# USEFUL REFERENCES        #
############################
CHR.LENGTHS.MOUSE = c(197195432,181748087,159599783,155630120,152537259,149517037,152524553,131738871,124076172,129993255,121843856,121257530,120284312,125194864,103494974,98319150,95272651,90772031,61342430,166650296,15902555)
#
############################
# LOADING AND WRITING DATA #
############################
#
# load.matrix = function( fn )
#   Load file in "carmen" format into data frame
#
# write.matrix = function(o, fn_out, na.string="NA" )
#   writes a data frame in the format expected by load.matrix
#
# read.eqtl = function(fn){
#   Reads table in Quigley eQTL format
#
# read.spear=function(fn){
#   Read output of Quigley spear program
#
# write.stanford = function( o, fn.out, col.annotation, row.annotation){
#   Write data frame in Stanford microarray format
#
# match.identifiers = function( SA.G, SA.P, translation, SA.G.valid.rows=NULL, SA.P.valid.rows=NULL)
#   Identify sample names that are matched between two data sets
#   When there is a column in one data set that matches the identifier of the other data set
#
# match.identifiers.eQTL = function( sa, expr, sa.col.matching.expr.identifiers)
#     use this to match identifiers in sa to identifiers in expr using a column in sa.
#    Genotype data and expression data usually have different identifier for the same samples
#
# match.identifiers.common.col = function( SA.A, SA.B, common.col.A, common.col.B, SA.A.valid.rows=NULL, SA.B.valid.rows=NULL){
#   Identify sample names that are matched between two data sets
#   When two data sets have a column that allows them to be paired up
#
# match.identifiers.sa.expr = function( sa, expr, valid.ids )
#   use this to match identifiers in sa to identifiers in expr.
#   convenience method to allow use of functions requiring this matchup.
#   valid.ids should be a list of ids, not a T/F vector.
#
# match.indices = function(ids, A, B)
#   Identify sample indices (i.e. column indices) that are matched between two data sets.
#
# ids.for.limit = function(sa, limit)
#   parses a comma-delimited list of limits
#   each limit is FOO=BAR or FOO!BAR or FOO*BAR (* is synonym for !)
#   returns vector of sa rownames where all limits are met
#
# vals.in.col = function(names.ids, sa, sa.col.name){
#   names.ids is a vector of identifiers.  sa is a data.frame, sa.col.name is a string.
#   This returns the value of sa at column sa.col.name for each identifier in names.ids.
#   Pass in data frames; return value in sa.col in the order specified by ids.col.
#
# vals.in.row = function(names.ids, E, row.name=NULL)
#   names.ids is a vector of identifiers. E is a data frame
#   if no rowname is passed, the first row is used.
#   returns the value of E for each identifier in names.ids.
#
# expr.idx.from.id = function( expr, ids )
#   returns a vector of integers that returns column(id) for id in names(id)
#
# idx.in.region = function( chr, loc_start, range, ga )
#  returns vector of idx for genes in ga within loc_start +/- range
#
# write.dataset=function(sa, gene.symbols, E, fn)
#   writes out dataset in a format readable by TMEV or suitable for Excel.
#
# do.sam = function(E, ga, A, B, med.FDR){
#   Perform a turn-key SAM analysis returning values with FDR<med.FDR
#
# do.sam.paired = function(E, ga, A, B, med.FDR)
#   Perform a turn-key paired SAM analysis returning values with FDR<med.FDR
# 
# write.SAM.siggenes(siggenes, fn)
#   Given output of SAM analysis, write siggenes to file. Returns data frame written.
#
# SAM.convert.siggenes = function(T){
#   combine nasty SAM output into a single dataframe
#
# write.output = function(fn_out, o)
#  write vector o to file fn_out
#
# vector.to.file = function( V, file="/users/dquigley/temp/test.txt" )
#  dump vector to file one value per line
#
# find.dupes = function(fn){
#  identify the duplicates in column one
#
#####################################
# Protein quantification
#####################################
#
# quantify.protein = function(d)
#  determine concentration of protein from BCA assay
#
#
#####################################
# RNA-seq code
#####################################
# collapse.rs.to.genes=function( counts.rs, GCs, lengths, symbols ){
#   Given data frame of counts for transcripts, GC content, lengths, symbol of each transcript,
#   consolidate to consensus values.
#
# strain.linear.model = function( FF, STRAINS ){
#   Calculate linear model of strain on expression. Not truly RNA-seq specific
#
# calculate.FPKM = function( counts, lengths ){
#   N = counts * 1,000,000 / (sample.total * gene.length.in.kb)
#
#   
#####################################
# AROMA AND AFFYMETRIX HELPER TOOLS #
#####################################
#
# stripNormalize.from.XY = function (XY, useTarget = TRUE, verbose = TRUE)
#   Modified version of quantile normalization wrapper used when we only have XY and not raw iDAT
#
# genotype.CRLMM.from.XY = function(fn.X, fn.Y, fn.confs, fn.calls, fn.sa)
#   Wrapper for CRLMM genotyping used when we only have X and Y matrixes
#
#
################
# IMPUTE TOOLS #
################
#
# impute.prepare = function(calls, snp.ids, 
#                  snp.chr, snp.pos, rs.ids, strand, snp.alleleA, snp.alleleB, 
#                  target.chrom, target.pos.start, target.pos.end, 
#                  dir.output, fn.out)
#   convert 0/1/2 genotype format into impute format.
#
# impute.parse.results = function( target.chrom, sample.ids, dir.output, fn.out )
#   Read the result of a call to impute and convert into 0/1/2 genotype format
#
# impute.generate.command = function(target.chrom, target.pos.start, target.pos.end, bin, dir.haplotypes, dir.output, fn.out)
#   Generate the command to run impute (binary located at bin)
# 
#####################################
# AROMA AND AFFYMETRIX HELPER TOOLS #
#####################################
#
# prepare.affy430_2.dataset = function( dir.CEL, dir.OUT, fn_ga, ga.id )
#   convert CEL files into dataset files (full, above_bg, above_bg_refseq)
#
# prepare.affy133.dataset = function( dir.CEL, dir.OUT, fn.source.ga, ga.id ){
#   convert CEL files into dataset files (full, above_bg, above_bg_refseq)
# affymetrix.probe.dataset = function( Data, probe.id, cdf )
#   returns a list with expr, sa, ga
#
# describe.affymetrix.probes = function( Data, probe_id, cdf )
#   Requires library affy and CDF library
#
# aroma.load.chip.effect = function(path, dataset, combineAlleles)
#   Load and normalize all data in dataset at path.
#
# aroma.plot.region = function(cesN, cesNT, filename, chromosome, region)
#   Given matched chip estimate filesets cesN and cesNT, plot sample
#   filename on chromosome between region[1] and region[2] bases.
#
# aroma.extract.array.log2 = function(cesN, cesT, filename, chromosome=NULL)
#   given normals and tumors, extract log2 of SNP copy for filename.
#   If chromosome is passed, extract only that chromosome.
#
#
##############################
# PLOTTING AND DISPLAY TOOLS #
##############################
#
# plot.discretization = function(M, lbls=NULL, method="SD", bounds=0.5, x_lbl=NULL, y_lbl="expression")
#    plot a single gene, discretized using standard deviation
#
# plot.percent.altered.by.locus= function(M, chr, upper.bound=0.5, lower.bound=-0.5)
#    plot whole genome aCGH or SNP data, percentage amplified or deleted
#
# plot.percent.altered.by.tumor = function(M, upper.bound=0.3, lower.bound=-0.3)
#    plot percentage of genome altered for each tumor (as opposed to across genome)
#
# plot.alterations = function(M, name, chr.list, upper.bound=0.3, lower.bound=-0.3, color=T)
#    Plot the whole genome, label x axis with chromosomes
#
# plot.alterations.on.chr( M, name, chr.size.in.mb, marker.interval.size=1,  upper.bound=0.3, lower.bound=0.3, color=T)
#    Plot a single chromosome, label x axis with Mb
#
# plot.rows = function( M, sa=NULL, legend.labels=NULL, group.by=NULL, sorted=F, plot.type='l', no.plot=F)
#    general multiple probe plotting function
#
# plot.by.identifiers = function( dataset, probe.list, group.by=NULL, sorted=F, plot.type='l')
#    Convenience function to allow us to plot with a simple probe list.
#    dataset is list of {expr, sa, ga}
#
# plot.geno.expr = function(ids, expr.g, expr.e, probe.g, probe.e, main=NULL, show.lm.pval=F)
#    Plot a genotype against a phenotype (e.g. SNP genotype vs. expression of a gene)
#
# plot.geno.expr.bracket = function(ids, expr.g, expr.e, ga.g, probe.g, probe.e, main=NULL)
#   Plot a bracket of three physically consecutive genotypes against a phenotype
#
# plot.probe.paired.data = function(ids, expr.A, expr.B, probe.id, main=NULL, method="barplot")
#   barplot or stripchart plot of probe.id for matched subjects in expr.A, expr.B
#
# plot.chr.CRLMM = function(cS, cCN, s.no, max.y=5, x1=1, x2=0, lblS, lblCN)
#   Separately plots SNP and Copy Number (non-polymorphic) probes for one chromosome
#
# plot.cluster = function(M, labels=names(M), d=NULL, main=NULL, direction="samples", cex=0.5)
#   Plot hierarchical cluster of M 
#
# color.dendrogram.labels = function(n){
#   Helper function for plot.cluster.colors
#
# plot.cluster.colors=function( M, colors, force.flat=F, show.labels=T )
#   Given dataframe M, plot with labels colored according to colors.
#
# plot.clean = function(V, ymin=NA, ymax=NA, colors=c(), cex=0.25 )
#   This is the utility function that just plots what it's given without labels.
#   Intended for use in figures
#
# sorted.heatmap=function(D, ga, target.symbols=NULL, target.probes=NULL, symbol=NULL, scale=F, y.min=NULL, y.max=NULL )
#   Given a matrix and a list of gene names and symbol, plot gene names sorted by symbol
#
# do.pca = function( D, pca=NULL, labels=NULL)
#    Helper function for quick visualization of PCA. Automatically colors based on labels.
#
########################
# STATISTICAL ANALYSIS #
########################
#
### CORRELATION
#
# cor.test.row.vs.genes = function(ids, expr.A, expr.B, ga.genes, probe.id, percent.present=0.9, method="pearson", verbose=T)
#   Calculate all correlations between a row (identified by probe.id) in expr.A and rows in expr.B
#
# cor.test.all.rows.vs.genes(ids, expr.A, expr.B, ga.genes, percent.present=0.9, method="pearson", max.p=1E-05, verbose=T)
#   Calculates all correlations between all rows in expr.A and rows in expr.B.
#   returns data frame with results where p <= max.p
#
# cor.test.row.vs.all = function( ids, expr.row, expr.all, percent.present=0.90, method="spearman", verbose=F){
#    Calculate correlation value of first row of expr.row vs. all rows in expr.all
#
# cor.test.probe.vs.all = function( expr, symbols, probe, percent.present=0.9, method="spearman", verbose=F)
#    Correlation of one probe vs the whole dataset.
#    **This is probably the one you want.**
#
# cor.sa.columns = function( sa, method="spearman" )
#    correlation of all numeric columns in data frame sa
#
# write.spear.to.cytoscape = function( fn.base, DF, attributes=NULL )
#   write spear file to cytoscape
#
# identify.eQTL.networks = function(eqtl, spear, probes.snp, ga.snp.chr, ga.snp.pos, 
#                                  probes.gene, ga.gene.chr, ga.gene.loc, min.trans=1)
#
# write.eqtl.to.cytoscape = function( fn.base, QTL, SPEAR, snps, probes, all.S, S.chr, S.loc, all.G, G.chr, G.loc){
#   write QTL and SPEAR to cytoscape, using snps and probes specified
#
# calculate.min.DC.score = function( n.probes )
#   What is the Z score required for a P value lower than the bonferroni-corrected
#   number of tests for differential correlation
#
# spear = function( fn_spear, fn_expr, fn_ga, fn_sa, min_cor, fn_out, class.a="", class.b="", probe="", y="symbol" )
#   simple wrapper for spear 
#
# overlap.spear.pairs = function(A, B, p2s )
#   find probe pairs in two spear files
#
### ANOVA AND REGRESSION
#
# paired.t.test = function(ids, expr.A, expr.B, method='ttest', percent.present=0.9, verbose=T)
#   two-sided t-test for values in expr.A vs. expr.B matched using ids
#
# unpaired.t.test = function(expr.A, expr.B, method='ttest', percent.present=0.9, verbose=T)
#   two-sided t-test for values in expr.A vs. expr.B
#
# logistic = function(sa, expr, do.plot=F)
#    perform a logistic regression of identifiers in sa with values in expr
#    will use the first column of sa and expr for values, matching sa rownames to expr names
#
# fisher = function(A, B, labels.A = NULL, labels.B = NULL, verbose=T)
#   convert two vectors of {0,1} to four categories and perform fisher exact test
#
### SURVIVAL ANALYSIS
#
# km=function( times, had.events, conditions, main=NULL, legends=NULL, fn_out=NULL, verbose=T, legend.x=3, legend.y=0.25)
#    General Kaplan-meier plotting and p-value calculation
#
# km.multi = function( times, had.events, conditions, legends=NULL, fn_out=NULL)
#    General Kaplan-meier plotting and p-value calculation
#
### QTL TOOLS
#
# calculate.linear.model.se = function(linear.model)
#   Given a solved linear model, return the standard error of the coefficients
#
# calculate.lm.pval=function(linear.model)
#   Given a solved linear model, return the p-values of the coefficients
#
# calculate.EQTL = function(dataset, shared.col, covariate=NULL, gene.idx=NULL, n.perm=0)
#   Calculates eQTLs using a linear model.
#   dataset is a list with members e.gene, s.gene, g.gene, e.snp, s.snp, g.snp, gene.names
#
# eqtl.filter = function(eqtl, snp, ga.chr.snp, ga.loc.snp, gene, ga.chr.gene, ga.loc.gene, type="cis", window=0)
#
# write.rqtl.csv=function( expr, probe.list, sa.geno, sa.col.matching.expr, geno.chrom, calls.geno, fn.rqtl ){
#   write all probes in probe.list as phenotypes for an RQTL CSV file
# 
# generate.rQTL.matrix = function( pheno.ids, pheno, e.snp, ga.snp, sa.snp, fn.out)
#   generate an input file for R/QTL
#
# pheno.from.expr = function( sa.SNP, probe.id, expr, sa.shared.col, sa.valid.ids)
#   Used to match probe expression to genotype ids, push input to generate.rQTL.matrix
#
# genotype.for.expression.ids = function( probe.snp, expr, expr.SNP, sa.SNP, sa.SNP.shared.column )
#   Identify the genotype at probe.snp for expression ids in expr
#
# align.geno.pheno = function( IDLIST, GENO, PHENO, probe.geno, probe.pheno ){
#   Match the value of probe.geno to the value of probe.pheno indexed by IDLIST
#
# means.by.genotype = function(IDLIST, GENO, PHENO, probes_genotype, probes_phenotype)
#   report means by genotype, assumes genotypes are in {0,1}
#
# annotate.eqtl.with.locations = function(E, snp2idx, gene2idx, G.snp.c, G.snp.l, G.expr.c, G.expr.l)
#   Add chromosome and location to output of eQTL analysis.
#
# 
# find.eQTL.overlaps = function(A, B, ga.A, ga.B, OVERLAP.LENGTH=5000000)
#
# validate.eqtl = function( source, ids, ga.source, ga.valid, ga.expr.source, ga.expr.target, geno.target, expr.target ){
#   given a source of eQTL results in Quigley eQTL format, for each eQTL:
#     find the closest SNP in the validation dataset
#     find the set of gene probes in validation with the same entrez ID
#     return the lowest P value for an eQTL in validation
#
### DIFFERENTIAL EXPRESSION (SAM)
#
# SAM = function(expr, ga, valid.p, s1, s2, nperms=100)
#    Calls SAM on expr, using s1 and s2 for labels 1 and 2
#
### CGH DATA
#
# estimate.percents.altered = function( M, chrs, locs, bound=0.3, binsize=1000000 )
#    estimate genomic percentage altered by binning loci and finding nearest marker
#
### DIFFERENTIAL CORRELATION
# 
# generate.simulated.DC.data = function(nA=100, nB=100, n.steps=200)
#   create simulation raw data with DC for every bin between 0 and 2 at intervals of 0.01
#
# DC.within.labels = function( A1, B1, A2, B2, n.perms=1000 ){
#   Given raw values for A and B in two classes, calculate a P value for differential correlation
#   by sampling within each class (as opposed to swapping labels)
#
# DC.between.labels = function( A1, A2, B1, B2, n.perms=1000 ){
#   Given values for probes 1 and 2 in two classes A and B, calculate a P value for differential correlation by sampling between classes
#
#
### PLINK
# calls.to.PLINK = function( calls, snp.ids, sample.ids, sexes, CHR, POS, fn.tped, fn.tfam ){
#   write genotype calls to PLINK TPED and TFAM files
#
##################
# GSCA functions #
##################
#
# (adapted from Choi & Kendziorski Bioinformatics 2009)
#
# singleDC <- function(data, group, probe.lists, nperm, method='spearman', focus.probe=NA)
#   Performs differential correlation analysis.
#
# compute.test.stat = function(d, group, m='spearman', focus.probe=NA){
#   utility function for singleDC
#
# get.idx = function(data, probe.lists){
#   utility function for singleDC
#
#
###################
# PSCBS functions #
###################
#
# psCBS.filter = function(s, min.mark=15, TC.range=c(0,30), min.range=c(0,30), max.range=c(0,30), min.diff=0, ga=NULL){
#   remove results from psCBS run that do not match filter settings. If ga !=NULL, count loci that match settings.
# psCBS.dissect = function(r, ga, chr=NULL, loc=NULL )
#   details of psCBS results by chromosome, or spanning a given locus on a chromosome
# psCBS.plot = function( ga, chr, loc.begin=NULL, loc.end=NULL )
#   plot counts in ga for whole chromosome or subset specified by bases
# psCBS.find.het = function(ga, s.names, s.idx, chrom, loc_begin, loc_end, is.het, A.t, B.t, A.n, B.n ){
#   reports raw values for hets
#
#
################
# HASH WRAPPER #
################
# hsh_new = function()
#   Create new
#
# hsh_in = function(H, key)
#   Does key exist in H
#
# hsh_get = function( H, key )
#   Return value for key in H. No error checking.
#
# hsh_set = function( H, key, value )
#   Set value for key in H to value. Will clobber existing value.
#
# hsh_keys = function( H )
#   Returns sorted vector of keys(H)
#
# hsh_keys_values = function( H )
#   Returns data frame of keys and values(H)
#
# hsh_from_vectors = function( key, value )
#   create hash key -> value
#
#
######################
# hashgraph class
######################
#
# hashgraph = ()
#   constructor
#
# addEdge = function(x, e1, e2, ...)
#   add edge e1, e2 to x
#
# hasEdge(x, e1, e2)
#   return boolean for whether e1,e2 exists in x 
#
# nodes = function(x)
#   return vector of nodes in x
#
# edges=function(x)
#  return dataframe e1, e2 of all edges in x
#
# neighbors = function(x, node)
#   return vector of neighbors of node in x
#
# subgraph = function(x, nlist)
#   restrict x to nodes in nlist
#
# length = function(x)
#   number of nodes in x
#
# intersection = function(x, y)
#   return new graph where edges found in both x and y
#
# difference = function(x, y)
#   return new graph where edges found in x not y
#
# triangles = function(x)
#   return subgraph of x where each node is in a clique of N=3
# 
# minimum.degree = function(G, min.n)
#   return G such that all nodes have degree >= min.n
#
#####################
# UTILITY FUNCTIONS #
#####################
#
# get.split.col = function(v, string, col=0, last=F, first=F)
#
# standardize = function(D)
#   standardize matrix by rows
#
# compress.probes = function( D, idx, min.cor=0.8 )
#   Given a set of probe indexes idx, look for those with correlation >= min.cor
#   Report mean values across all probes with correlation >= min.cor
#   If no genes pairs meet these criteria, default to mean value
#
# count.appearances = function(V, order.by="values")
#   return dataframe how how many times each item in V appears
#
# match.idx = function(A, B)
#   return dataframe of indices into A and B restricted to perfect matches between A and B, 
#   where idx.A[i] == idx.B[i] for each i in matched pairs. 
#
# match.idx.first = function( A, B )
#   for intersection, returns first index of B in first index of A
#
#  mean.by.symbol = function( expr, all.symbols, symbols )
#   find all occurences of each element in symbols in all.symbols; return mean of those values in expr
#
# create.meta.ga = function( ga.1, ga.2 ){
#   Given two gene attributes data frames, creates a new ga with meta-rownames.
#
# create.meta.expr = function( ga, expr.1, expr.2 ){
#   Given gene attributes created by create.meta.ga, combine the 
#   expression data in expr.1 and expr.2 to generate a matched expression dataset
# 
# probe.paired.data = function(ids, expr.A, expr.B, probe.id)
#   reports probe.id values for matched subjects in expr.A, expr.B.  Returns a data frame.
#
# ids.with.attribute = function(ga, col.value, col.name="Gene Name")
#   returns vector of identifiers where data frame ga has col.name equal to col.value
#
# set.difference = function(A, B)
#   returns the set difference set(A) - set(B), not vector of T/F
#
# set.intersection = function(A, B)
#   returns the set intersection set(A, B), not vector of T/F
#
# contains = function(A, B){
#   is there a value in B for each value in A?
#
# nearest.ten = function( x )
#   find the nearest factor of 10 for x between 10 and 100, rounding up
#
# plot_multiple = function(V, lbls, x_lbl, probe_label_list, plot.type='l')
#   utility function to actually draw plots
#
# convert.to.ranks = function(M)
#   Convert columns in M to their rank
#
#
# make_list = function(keys, vals)
#   Convert pairs to a list
#
# abi.chr = function(snp)
#   extract the Chromosome from an ABI-style SNP id (e.g. E16.018.064_10)
#
# abi.mb = function(snp, is.mouse=T)
#   extract the MB from an ABI-style SNP id
#
# rows.at.threshold = function(expr.A, percent.present=0.9)
#   return vector of rows where A meets threshold for percent present
#
# paired.rows.at.threshold = function(ids, expr.A, expr.B, percent.present=0.9){
#   return vector of rows where both A and B meet threshold for percent present
#
# prepend.label = function(list.of.targets, label){
#   prepend label to a vector of characters
#
# describe.affymetrix.probes = function( Data, probe_id, cdf )
#   describes low-level (median, mean, position on chip) details of a probe
#
# numeric.values = function(samples, row.id)
#   extracts value for each sample at row, returns numeric vector
#
# smooth = function(x, chrom, smooth.region = 2, outlier.SD.scale = 4, smooth.SD.scale = 2, trim=0.025)
#
#####################################################################

##################
### BEGIN LOAD ###
##################

load.matrix = function( fn ){
    read.table(file=fn, sep="\t",row.names=1, check.names=F, header=T, na.strings=c('NA', '-99999'), stringsAsFactors=F)
}

read.eqtl = function(fn){
    x = read.table(fn, sep='\t', stringsAsFactors=F)
    if( dim(x)[2] == 9)
        names(x) = c('raw.p','symbol','probe','snp','perm.p','aa','ab','bb','t.stat')
    else
        names(x) = c('raw.p','symbol','probe','snp','perm.p','aa','ab','bb')
    x
}

read.spear=function(fn){
    # Read output of Quigley spear program
    x = read.table(fn, sep='\t', stringsAsFactors=F)
    if(dim(x)[2]==9)
        names(x) = c('symbol.1','probe.1','symbol.2', 'probe.2','rho.a', 'rho.b', 'diff','perm.p', 'z.score')
    else
        names(x) = c('symbol.1','probe.1','symbol.2', 'probe.2','rho.a', 'rho.b', 'diff','perm.p')
    x
}

read.coverage=function(fn){
    # Read output of BAM coverage table
    x=read.table(fn, sep='\t', header=F, stringsAsFactors=F)
    names(x) = c('chrom', 'unknown', 'type', 'start', 'stop', 'a', 'strand', 'b', 'id', 'count', 'n.bases', 'feature.size', 'fraction')
    x
}

write.matrix = function(o, fn_out, na.string="NA"){
    # Write data frame in standard format used by load.matrix
    rownames = append( 'IDENTIFIER', names(o) )
    write( rownames, fn_out, sep='\t', ncolumns=length(rownames) )
    write.table(o, fn_out, quote=F, sep='\t', row.names=T, na=na.string, col.names=F, append=T)
}

write.stanford = function( o, fn.out, col.annotation, row.annotation){
    # Write data frame in Stanford format. This format has a row of sample names, a row of annotation, 
    # then the data starts. There is a column of gene names, a column of gene annotation (e.g. probes) 
    # and then the data starts. The upper-left corner of data is 3,3
    tt = rbind( row.annotation, o)
    o.stan = data.frame( c( '', col.annotation), tt )
    names(o.stan)[1] = "annotation"
    rownames(o.stan)[1] = 'SSP'        
    rn = append( 'CLID', names(o.stan) )
    write( rn, fn.out, sep='\t', ncolumns=length(rn) )
    write.table(o.stan, fn.out, quote=F, sep='\t', row.names=T, col.names=F, append=T)
}
    
match.identifiers = function( SA.G, SA.P, translation, SA.G.valid.rows=NULL, SA.P.valid.rows=NULL){
    # Use this when we have a two datasets without a shared column, but we want a
    # one-way translation from one data sets to another
    # EXAMPLE: Matching genotypes to expression data
    # SA.G is a data frame of sample attributes for the genotype measurements
    # SA.P is a data frame of sample attributes for the phenotype measurements
    # translation is a vector that translates identifiers from SA.G to SA.P
    #    if NA, the row is not valid.
    # if SA.G.valid.rows is not null, T/F vector for SA.G rows to include
    # if SA.P.valid.rows is not null, T/F vector for SA.P rows to include
    # returns a data frame of matched, valid identifiers from SA.G and SA.P

    SA.G.len = dim(SA.G)[1]
    SA.P.len = dim(SA.P)[1]
    if( is.null( translation) ){ stop("translation cannot be NULL") }
    if( dim(SA.G)[1] != length(SA.G.valid.rows)){ stop("SA.G rows not same as length of SA.G.valid.rows.  Check parameter order.")}
    if( dim(SA.P)[1] != length(SA.P.valid.rows)){ stop("SA.P rows not same as length of SA.P.valid.rows.  Check parameter order.")}

    if( length(translation) != SA.G.len ){
        stop("translation and SA.G do not have equal numbers of elements ")
    }

    if( is.null( SA.G.valid.rows) ){
        SA.G.valid.rows = !is.na(translation) # all are valid except those missing translation
    }
    else{
        SA.G.valid.rows = SA.G.valid.rows & !is.na(translation)
    }
    if( is.null( SA.P.valid.rows) ){
        SA.P.valid.rows = vector(mode='logical', SA.P.len)
        SA.P.valid.rows = !SA.P.valid.rows
    }
    SA.P.valid.rows[ is.na(SA.P.valid.rows) ] = F
    SA.G.valid.rows[ is.na(SA.G.valid.rows) ] = F

    ids.G = rownames(SA.G)
    ids.P = rownames(SA.P)
    genotype.id = vector(mode='character')
    phenotype.id = vector(mode='character')
    cnt=0
    for( i in 1:SA.G.len ){
        if( SA.G.valid.rows[i] ){
            find.in.P = translation[i]
            for( j in 1:SA.P.len ){
                if( SA.P.valid.rows[j] ){
                    if( ids.P[j] == find.in.P){
                        cnt = cnt+1
                        genotype.id[cnt] = ids.G[i]
                        phenotype.id[cnt] = ids.P[j]
                        break
                    }
                }
            }
        }
    }
    data.frame(genotype.id, phenotype.id)
}

match.identifiers.eQTL = function( sa, expr, sa.col.matching.expr.identifiers){
    ids.expr = names(expr)
    ids.sa = rownames(sa)
    ids.match = sa.col.matching.expr.identifiers
    if( length(ids.sa) != length(ids.match) ){ stop("sa.col.matching.expr.identifiers must have same length as rownames(sa)" ) }
    idx.in.expr = vector(mode="integer")
    idx.in.snps = vector(mode="integer")
    id.s = vector(mode="character")
    id.e = vector(mode="character")
    for(i in 1:length(ids.match)){
        idx = which( ids.expr == ids.match[i] )
        if(length(idx)==1){
            id.s[ length(id.s)+1 ] = ids.sa[i]
            id.e[ length(id.e)+1 ] = ids.expr[idx]
            idx.in.expr[length(idx.in.expr)+1 ] = idx[1]
            idx.in.snps[length(idx.in.snps)+1 ] = i
        }
    }
    data.frame(id.s, id.e, idx.in.expr, idx.in.snps, stringsAsFactors=F)
}

match.identifiers.common.col = function( SA.A, SA.B, common.col.A, common.col.B, SA.A.valid.rows=NULL, SA.B.valid.rows=NULL){
    # Use this when two data sets have a column that allows them to be paired up
    # EXAMPLE: Matching normal samples to tumor samples from a shared patient ID
    # SA.A is a data frame of sample attributes for the genotype measurements
    # SA.B is a data frame of sample attributes for the phenotype measurements
    # common.col.A and common.col.B are vectors of some value we wish to match
    #   between A and B.  if NA, the row is not valid.
    # if SA.A.valid.rows is not null, T/F vector for SA.A rows to include
    # if SA.B.valid.rows is not null, T/F vector for SA.B rows to include
    # returns a data frame of matched, valid identifiers from SA.A and SA.B

    SA.A.len = dim(SA.A)[1]
    SA.B.len = dim(SA.B)[1]
    if( is.null( common.col.A) ){ stop("common.col.A cannot be NULL") }
    if( is.null( common.col.B) ){ stop("common.col.B cannot be NULL") }

    if( is.null( SA.B.valid.rows ) ){
        SA.B.valid.rows = vector(mode='logical', SA.B.len)
        SA.B.valid.rows = !SA.B.valid.rows
    }
    if( is.null( SA.A.valid.rows ) ){
        SA.A.valid.rows = vector(mode='logical', SA.A.len)
        SA.A.valid.rows = !SA.A.valid.rows
    }
    SA.B.valid.rows[ is.na(SA.B.valid.rows) ] = F
    SA.A.valid.rows[ is.na(SA.A.valid.rows) ] = F

    if( dim(SA.A)[1] != length(SA.A.valid.rows)){ stop("SA.A rows not same as length of SA.A.valid.rows.  Check parameter order.")}
    if( dim(SA.B)[1] != length(SA.B.valid.rows)){ stop("SA.B rows not same as length of SA.B.valid.rows.  Check parameter order.")}

    if( length(common.col.A) != SA.A.len ){
        stop("common.col.A and SA.A do not have equal numbers of elements ")
    }
    if( length(common.col.B) != SA.B.len ){
        stop("common.col.B and SA.B do not have equal numbers of elements ")
    }

    ids.A = rownames(SA.A)
    ids.B = rownames(SA.B)
    A.id = vector(mode='character')
    B.id = vector(mode='character')
    shared = vector(mode='character')
    cnt=0
    for( i in 1:SA.A.len ){
        if( SA.A.valid.rows[i] ){
            find.in.B = common.col.A[i]
            for( j in 1:SA.B.len ){
                if( SA.B.valid.rows[j] ){
                    if( common.col.B[j] == find.in.B){
                        cnt = cnt+1
                        A.id[cnt] = ids.A[i]
                        B.id[cnt] = ids.B[j]
                        shared[cnt] = find.in.B
                        break
                    }
                }
            }
        }
    }
    data.frame(A.id, B.id, shared)
}


match.identifiers.sa.expr = function( sa, expr, valid.ids){
    # use this to match identifiers in sa to identifiers in expr.
    # convenience method to allow use of functions requiring this matchup.
    # valid.ids should be a list of ids, not a T/F vector.
    ids.sa = rownames(sa)
    ids.sa = ids.sa[ contains(ids.sa, valid.ids)]
    ids.expr = names(expr)
    N = length(ids.sa)
    if( length(N)==0 ){ stop("No valid identifiers found to match between sa and valid.ids; check parameters") }
    tf = contains(ids.sa, ids.expr)
    if( sum(tf) != N ){ stop("not every identifier in rownames(sa) has a match in names(expr)") }
    data.frame(A.id = ids.sa, B.id = ids.sa)
}


match.indices = function(ids, A, B){
    # once we have matched identifiers between two data sets, we may wish to know
    # which columns in a data frame (e.g. expression samples) correspond to those
    # matched pairs.  Pass in two data frames A and B and the ids generated through
    # match.identifiers() or a similar function.  Returns a data frame with indices
    # matched in A and B
    #
    N = dim(ids)[1]
    samples.A = names(A)
    samples.B = names(B)
    idx.A = vector(mode="integer", N)
    idx.B = vector(mode="integer", N)
    for(i in 1:N){
        iA = which( samples.A==ids[i,1] )
        iB = which( samples.B==ids[i,2] )
        if(length(iA)==0){ stop( paste("identifier '",ids$A.id[i],"' not found in expr.A at id row",i)) }
        if(length(iB)==0){ stop( paste("identifier '",ids$B.id[i],"' not found in expr.B at id row",i)) }
        idx.A[i] = iA
        idx.B[i] = iB
    }
    data.frame(A = idx.A, B = idx.B)
}


vals.in.col = function(names.ids, sa, sa.col.name){
    # names.ids is a vector of identifiers.  sa is a data.frame, sa.col.name is a string.
    # This returns the value of sa at column sa.col.name
    # for each identifier in names.ids.
	#
    v = vector()
    idx.sa = which(names(sa)==sa.col.name)
    for( i in 1:length(names.ids) ){
        x = sa[ which(rownames(sa)==names.ids[i]), idx.sa ]
        v[i] = x
    }
    v
}


vals.in.row = function(names.ids, E, row.name=NULL){
    # names.ids is a vector of identifiers.
    # E is a data frame; if no rowname is passed, the first row is used.
    # This returns the value of E for each identifier in names.ids.
    #
    v = vector()
    if( is.null(row.name) )
        e = E[1,]
    else
        e = E[ which(rownames(E)==row.name),]
    for( i in 1:length(names.ids) ){
        x = e[ 1, which(names(E)==names.ids[i])]
        v[i] = x
    }
    v
}


ids.for.limit = function(sa, limit.str){
    # parses a comma-delimited list of limits
    # each limit is FOO=BAR (equals) or FOO!BAR or FOO*BAR (! and * are not equals)
    # returns vector of sa rownames where all limits are met
	# zero length limit.str returns all identifiers.
	#
	ids = rownames(sa)
	if(nchar(limit.str)>0){
		limits = strsplit(limit.str, ",", fixed=T)[[1]]
		for(i in 1:length(limits)){
			limit=limits[i]
			if( length(grep("=",limit))>0 ){
				match.type="equals"
				symbol="="
			}
			else if( length(grep("!",limit))>0 ){
				match.type="not equals"
				symbol="!"
			}
			else if( length(grep("\\*",limit))>0 ){
				match.type="not equals"
				symbol="*"
			}
			else{
				stop("badly formed limit: ", limit)
			}
			colname = strsplit(limit, symbol, fixed=T)[[1]][1]
			value = strsplit(limit, symbol, fixed=T)[[1]][2]
			if( length( which(names(sa)==colname) )==0 )
				stop("column ", colname, " not found")
			vals = vals.in.col(ids, sa, colname)
			if( match.type=="equals")
				ids = ids[vals==value]
			else
				ids = ids[vals!=value]
		}
	}
	ids
}


expr.idx.from.id = function( E, ids ){
    # returns a vector of integers that returns column(id) for id in names(id)
    idx = rep(-1, length(ids) )
    for( i in 1:length(ids) ){
        idx[i] = which( names(E) == ids[i] )
    }
    idx
}


idx.in.region = function( chr, loc_start, range, ga.chr, ga.loc ){
    # returns vector of idx for genes in ga within loc_start +/- range
	in.range = ( ga.chr==as.character(chr) & ( ga.loc < (loc_start+range) & ga.loc > (loc_start-range) ) )
	which(in.range==T)
}

align.geno.pheno = function( IDLIST, GENO, PHENO, probe.geno, probe.pheno ){
    # Match the value of probe.geno to the value of probe.pheno indexed by IDLIST
    # IDLIST is a data frame where columns are valid paired identifiers
    #    {genotype, phenotype}
    # GENO is a data frame of genotype data (numeric)
    # PHENO is a data frame of phenotype data (numeric)
    # returns a data frame of matched data from GENO and PHENO
    if( dim(IDLIST)[1] == 0 ){
        print("WARNING: idlist contains no entries")
    }
    row.G = which(rownames(GENO)==probe.geno)
    row.P = which(rownames(PHENO)==probe.pheno)
    if( length(row.P)==0 ){ stop(paste("WARNING: no matches for ",probe.pheno," in phenotype.")) }
    if( length(row.G)==0 ){ stop(paste("WARNING: no matches for ",probe.geno," in genotype.")) }
    samples.G = colnames(GENO)
    samples.P = colnames(PHENO)
    N = dim(IDLIST)[1]
    val.G = vector(mode='numeric', N)
    val.P = vector(mode='numeric', N)
    for(i in 1:N){
        col.G = which(samples.G==IDLIST[i,1])
        col.P = which(samples.P==IDLIST[i,2])
        if( length(col.G)==0 ){
            print( paste( "Cannot find genotype ID matching ", IDLIST[i,1],"in ", samples.G))
        }
        if( length(col.P)==0 ){
            print( paste( "Cannot find phenotype ID matching ", IDLIST[i,2],"in ", samples.P))

        }
        val.G[i] = GENO[row.G, col.G]
        val.P[i] = PHENO[row.P, col.P]
    }
    data.frame( val.G, val.P )
}

means.by.genotype = function(IDLIST, GENO, PHENO, probes_genotype, probes_phenotype){
    # Same input as align.geno.pheno
    if( length(probes_genotype) != length(probes_phenotype) ){
        stop("probes_genotype and probes_phenotype must have same length")
    }
    means.0 = rep(0, length(probes_genotype))
    means.1 = rep(0, length(probes_genotype))	
    for(i in 1:length(probes_genotype) ){
        vals = align.geno.pheno( IDLIST, GENO, PHENO, probes_genotype[i], probes_phenotype[i] )
        means.0[i] = mean(vals$val.P[vals$val.G==0], na.rm=T)
        means.1[i] = mean(vals$val.P[vals$val.G==1], na.rm=T)
    }
    data.frame(probes_phenotype, probes_genotype, means.0, means.1)
}

annotate.eqtl.with.locations=function(E, snp2idx, gene2idx, G.snp.c, G.snp.l, G.expr.c, G.expr.l){
    chr.snp = rep(0, dim(E)[1])
    chr.expr = rep(0, dim(E)[1])
    loc.snp = rep(0, dim(E)[1])
    loc.expr = rep(0, dim(E)[1])
    for(i in 1:dim(E)[1]){
        chr.snp[i] = G.snp.c[ hsh_get( snp2idx, E$snp[i] ) ]
        loc.snp[i] = G.snp.l[ hsh_get( snp2idx, E$snp[i] ) ]
        chr.expr[i] = G.expr.c[ hsh_get( gene2idx, as.character(E$probe[i] ) ) ]
        loc.expr[i] = G.expr.l[ hsh_get( gene2idx, as.character(E$probe[i] ) ) ]                        
    }
    same.chr=chr.snp==chr.expr
    chr.snp.num = chr.snp
    chr.snp.num[chr.snp.num=="X"] = 23
    chr.snp.num = as.numeric(chr.snp.num)
    E$aa = round(E$aa,3)
    E$ab = round(E$ab,3)
    E$bb = round(E$bb,3)
    E$raw.p = signif(E$raw.p,3)
    E$t.stat = round(E$t.stat,3)
    data.frame( E, chr.snp, chr.expr, loc.snp, loc.expr, same.chr, chr.snp.num, stringsAsFactors=F )
}
    
   
annotate.eqtl.old = function(E, fn.ga.snp, fn.ga.expr, cis.window.size=1000000){
    ga.snp = load.matrix(fn.ga.snp)
    ga.expr = load.matrix(fn.ga.expr)
    cc = ga.snp$chromosome; cc[cc=='X'] = 23; cc[cc=='Y'] = 24
    rs2chr = hsh_from_vectors(rownames(ga.snp), as.numeric(cc))
    rs2loc = hsh_from_vectors(rownames(ga.snp), as.numeric(ga.snp$loc_start))
    cc = ga.expr$chromosome; cc[cc=='X'] = 23; cc[cc=='Y'] = 24
    g2chr = hsh_from_vectors(ga.expr$symbol, as.numeric(cc))
    g2loc = hsh_from_vectors(ga.expr$symbol, as.numeric(ga.expr$loc_start))
    gene.chr = hsh_get(g2chr, E$symbol, T )
    gene.loc = hsh_get(g2loc, E$symbol, T )
    snp.loc = hsh_get(rs2loc, E$probe.snp, T)
    snp.chr = hsh_get(rs2chr, E$probe.snp, T)
    cis.trans = rep('unknown', length(snp.chr) )
    cis.trans[ snp.chr == gene.chr & abs( snp.loc - gene.loc ) < cis.window.size ] = 'cis'
    cis.trans[ snp.chr != gene.chr | (snp.chr == gene.chr & abs( snp.loc - gene.loc ) >= cis.window.size) ] = 'trans'
    cbind(E, gene.chr, gene.loc, snp.chr, snp.loc, cis.trans)
}


find.eQTL.overlaps = function(A, B, ga.A, ga.B, OVERLAP.LENGTH=5000000){
    if( length( which( names( A )=="entrez") )==0  | length( which( names( B )=="entrez"))==0 )
        stop("A and B must both have a column named 'entrez'")
    
    # requires that A and B both have "entrez" column for entrez ID
    entrez.A = unique( A$entrez )
    entrez.B = unique( B$entrez )   
    g2c.B = hsh_from_vectors(rownames(ga.B), ga.B$Chr)
    g2c.A = hsh_from_vectors(rownames(ga.A), ga.A$chromosome)
    g2l.B = hsh_from_vectors(rownames(ga.B), ga.B$loc_start)
    g2l.A = hsh_from_vectors(rownames(ga.A), ga.A$loc_start)
    genes.both = intersect(entrez.A, entrez.B )    
    chr.A = hsh_get( g2c.A, A$snp)
    chr.B = hsh_get( g2c.B, B$snp)
    loc.A = hsh_get( g2l.A, A$snp )
    loc.B = hsh_get( g2l.B, B$snp )    
    
    # Create index from A to B where value indicates that the line in 
    # A corresponds to a matching line in B. Match is defined as 
    # 1) same entrez 
    # 2) SNP on same chromosome
    # 3) locus within small interval
    # Note that one row in A could match several rows in B or vice versa
    idx.match.A = rep(-1, 10000)
    idx.match.B = rep(-1, 10000)
    ctr=1
    print("Looking for matches...")
    for(i in 1:dim(A)[1] ){
        idx=which(B$entrez==A$entrez[i])
        if(length(idx)>0){
            for(j in 1:length(idx) ){
                if( !is.na(chr.A[i]) & !is.na(chr.B[idx[j]]) ){
                    if( chr.A[i]==chr.B[ idx[j] ] ){
                        if( abs( loc.A[i] - loc.B[ idx[j] ]) < OVERLAP.LENGTH ){
                            idx.match.A[ctr] = i
                            idx.match.B[ctr] = idx[j]
                            ctr = ctr+1
                        }
                    }
                }
            }
        }
    }
    idx.match.A = idx.match.A[1:(ctr-1)]
    idx.match.B = idx.match.B[1:(ctr-1)]
    
    A.match = A[idx.match.A,]
    B.match = B[idx.match.B,]
    
    # Annotate matches with gene and SNP loci
    ucsc = load.matrix('/notebook/annotations/UCSC_all_known_human_genes_Feb_2009_simplified.txt')
    symbol2c = hsh_from_vectors(rownames(ucsc), ucsc$chr)
    symbol2loc = hsh_from_vectors(rownames(ucsc), ucsc$loc_start)
    gene.chr = hsh_get( symbol2c, A.match$symbol)
    gene.loc = hsh_get( symbol2loc, A.match$symbol)    
    A.match = cbind(A.match, gene.chr, gene.loc )
        
    chr.snp = hsh_get( g2c.A, A.match$snp )
    loc.snp = hsh_get( g2l.A, A.match$snp )
    A.match = cbind( A.match, chr.snp, loc.snp )
    chr.snp = hsh_get( g2c.B, B.match$snp )
    loc.snp = hsh_get( g2l.B, B.match$snp )
    B.match = cbind( B.match, chr.snp, loc.snp )
    B.match = cbind(B.match, gene.chr, gene.loc )
    
    list( A=A.match, B=B.match )
}

find.closest.snp = function( target.chr, target.loc, ga){
    ga.close = ga[ga$chromosome==target.chr & ga$loc_start > (target.loc-10000000) & ga$loc_start < (target.loc+10000000),]
    distance = abs(target.loc - ga.close$loc_start)
    smallest = min(distance)
    rownames(ga.close)[ which(distance==smallest)[1] ]
}
    
find.nearby.snp = function( target.chr, target.loc, ga, max.dist=10000){
   rownames(ga)[ga$chromosome==target.chr & ga$loc_start >= (target.loc-max.dist) & ga$loc_start <= (target.loc+max.dist)]
}
    
validate.eqtl = function( source, ids, ga.source, ga.valid, ga.expr.source, ga.expr.target, geno.target, expr.target, max.dist=10000 ){
    # given a source of eQTL results in Quigley eQTL format, for each eQTL:
    #    find the closest SNP in the validation dataset
    #    find the set of gene probes in validation with the same entrez ID
    #    return the lowest P value for an eQTL in validation
    valid.probe = rep('', dim(source)[1])
    valid.snp = rep('', dim(source)[1])
    valid.raw.p = rep(1, dim(source)[1])
    valid.r = rep(1, dim(source)[1])
    for(i in 1:dim(source)[1]){
        snp.idx = which(rownames(ga.source)==source$snp[i])
        #rs.valid = find.closest.snp( ga.source$chromosome[snp.idx], ga.source$loc_start[snp.idx], ga.valid )
        rs.valid = find.nearby.snp( ga.source$chromosome[snp.idx], ga.source$loc_start[snp.idx], ga.valid, max.dist )
        entrez = ga.expr.source$entrez[ which(rownames(ga.expr.source)==source$probe[i] ) ]
        probes.valid = rownames( ga.expr.target )[ ga.expr.target$entrez==entrez ]
        pval.min=1
        rval.max=0
        probe.min = 'NO_MATCH'
        rs.valid.min = 'NO_MATCH'
        print( paste(i, "of", dim(source)[1], " testing ", length(rs.valid), "snps" ) )
        if(length(probes.valid)>0 & length(rs.valid)>0 ){
            for(j in 1:length(probes.valid)){
                probe.valid = probes.valid[j]
                for(s in 1:length(rs.valid)){
                    rs = rs.valid[s]
                    vals = align.geno.pheno(ids, geno.target, expr.target, rs, probe.valid)
                    L = lm(as.numeric(vals$val.P) ~ as.numeric(vals$val.G) )
                    pval = signif( (summary( L )$coefficients[2,4]), 2)
                    rval = signif( (summary( L )$adj.r.squared), 2)
                    if( pval<pval.min){
                        pval.min=pval
                        rval.max=rval
                        probe.min = probe.valid
                        rs.valid.min = rs.valid[s]
                    }
                }
            }
        }
        valid.probe[i] = probe.min
        valid.snp[i] = rs.valid.min
        valid.raw.p[i] = pval.min
        valid.r[i] = rval.max
    }
    v=data.frame(valid.probe, valid.snp, valid.raw.p, valid.r, stringsAsFactors=F)
    v = v[order(v$valid.raw.p),]
    v
}

write.dataset=function(sa, gene.symbols, E, fn){
    col.1 = c("PROBE_ID", rep('', (dim(sa)[2])+1), rownames(E))
    col.2 = c("Name", "IDENTIFIER", names(sa), gene.symbols)
    S = data.frame( t(sa.fc), stringsAsFactors=F )
    E = rbind(names(E),names(E),S, E)
    E = cbind( col.1, col.2, E )
    write.table( E, fn, row.names=F, col.names=F, sep='\t', quote=F )
}


write.SAM.siggenes = function(siggenes.table, fn){
	# Given output of SAM analysis, write siggenes to file. Returns data frame written.
    has.up=F
    has.dn=F
    if(!is.null(siggenes.table$genes.up)){
        gu = siggenes.table$genes.up[,c(2,3,4,7,8)]
        if(length(dim(gu)) == 0)
          gu = data.matrix(t(gu))
        gu = cbind(gu, rep('up',dim(gu)[1]))
        has.up = T
    }
    if(!is.null(siggenes.table$genes.lo)){
        gd = siggenes.table$genes.lo[,c(2,3,4,7,8)]
        if(length(dim(gd)) == 0)
          gd = data.matrix(t(gd))
        gd = cbind(gd, rep('dn',dim(gd)[1]))
        has.dn = T
    }
    if(has.up & has.dn)
        g.all = rbind(gu,gd)
    else if(has.up)
        g.all = gu
    else if(has.dn)
        g.all = gd
    else
        stop("No significant genes; cannot write file.")
    g.all[,3] = round(as.numeric(g.all[,3]),2)
    g.all[,4] = round(as.numeric(g.all[,4]),2)
    g.all[,5] = round(as.numeric(g.all[,5]),2)
    g.all[,4] = signif(as.numeric(g.all[,4]),4)
    g.all = g.all[order(as.numeric(g.all[,3])),]
    g.all = data.frame(g.all)
    rownames(g.all) = g.all[,2]
    names(g.all) = c('gene.name', 'probe', 'score', 'fold.change', 'q.val', 'dir')
    write.table( g.all[,c(1,3,4,5,6)], fn, quote=F, sep='\t')
    g.all[,c(1,3,4,5,6)]
}


write.output = function(fn_out, o){
	# write values in vector o to fn_out, one per line
	for(i in 1:length(o) ){
		if(i==1)
			write(o[i], file=fn_out)
		else
			write(o[i], file=fn_out, append=T)
	}
}


vector.to.file = function( V, file="/users/dquigley/temp/test.txt" ){
    write.table( data.frame(V), file=file, row.names=F, col.names=F, quote=F)
}

find.dupes = function(fn){
    x = read.table(fn, sep='\t', header=T, stringsAsFactors=F)
    freq = count.appearances(x[,1])
    freq$keys[freq$values>1]
}


################################
# Begin protein code           #
################################

quantify.protein = function(d){
    # Input: one row per sample
    # First 7 rows should be the standard dilutions
    
    means = rowMeans(d)
    x = c(0, 1, 2, 5, 10, 20, 50)
    y = means[1:7]
    sds.dil = rep(0, 7)
    for(i in 1:7){
        sds.dil[i] = sd( as.numeric( d[i,]) )
    }
    sds.dil = round(sds.dil,3)
    means.dil = round(means[1:7], 3)
    print( data.frame(means.dil, sds.dil) )
    
    df = data.frame(x,y)
    L=lm(y~x, data=df)
    plot(x,y, pch=19)
    abline(L)
    prots = means[8:length(means)]
    y.int = L$coefficients[1]
    slope = L$coefficients[2]
    sds = rep(0, length(prots))
    for(i in 1:length(sds)){
        sds[i] = sd( as.numeric( d[i+7,]) )
    }
    sds = round(sds,3)
    um.per.2ul = round((prots - y.int)/slope, 3)
    um.per.ul = round( um.per.2ul / 2, 3)
    means = round(means[8:length(means)],3)
    thirty.ug = round( 30/um.per.ul, 2 )
    twenty.ug = round( 20/um.per.ul, 2 )
    h20.for.30 = round(30 - thirty.ug, 2)
    h20.for.20 = round(20 - twenty.ug, 2)
    points(um.per.2ul, means)
    data.frame(means, sds, um.per.2ul, um.per.ul, thirty.ug, h20.for.30, twenty.ug, h20.for.20 )
}


#######################
# BEGIN RNA-seq CODE  #
#######################

collapse.rs.to.genes=function( counts.rs, GCs, lengths, symbols ){
    # Given data frame of counts for transcripts, GC content, lengths, symbol of each transcript,
    # consolidate to consensus values.
    if(length(GCs)!=length(lengths) | length(GCs)!=dim(counts.rs)[1] | length(GCs) != length(symbols) ){
        stop("Inconsistent data passed; data fram and three vectors must have same length")
    }
    unique.symbols = sort( unique( symbols ) )
    CNT = matrix(0, nrow=length(unique.symbols), ncol=dim(counts.rs)[2] )
    LEN = rep(0, length(unique.symbols))
    GC = rep(0, length(unique.symbols))
    D = data.matrix(counts.rs)
    for(i in 1:length(unique.symbols) ){
        symbol = unique.symbols[i]
        idx = which( symbols==symbol )
        if( length(idx)==1 ){
            CNT[i,] = D[idx,]
            LEN[i] = lengths[idx] 
            GC[i] = GCs[idx]
        }
        else{
            CNT[i,] = colMeans( D[idx,] )
            LEN[i] = max(lengths[idx])
            GC[i] = mean(GCs[idx])
        }
    }
    counts = data.frame( round(CNT) )
    rownames(counts) = unique.symbols
    names(counts) = names(counts.rs)

    list(counts=counts, GC=GC, lengths=LEN)
}
    
calculate.FPKM = function( CNT, LEN ){
    # N = counts * 1,000,000 / (sample.total * gene.length.in.kb)
    total.reads = colSums(CNT)
    fpkm = matrix(0, nrow=dim(CNT)[1], ncol=dim(CNT)[2])
    seq.lengths = matrix(0, nrow=dim(CNT)[1], ncol=dim(CNT)[2])
    totals = matrix(0, nrow=dim(CNT)[1], ncol=dim(CNT)[2])
    for(i in 1:dim(seq.lengths)[2] ){
        seq.lengths[,i] = LEN
        totals[,i] = total.reads[i]
    }
    seq.lengths = seq.lengths/1000
    fpkm = CNT / (totals*seq.lengths) * 1000000
    #for(i in 1:dim(CNT)[1]){
    #    fpkm[i,] = as.numeric( CNT[i,]) / (total.reads * LEN[i] / 1000 )
    #}
    #fpkm = fpkm * 1000000
    fpkm = data.frame(fpkm)
    names(fpkm) = names(CNT)
    rownames(fpkm) = rownames(CNT)
    fpkm
}

strain.linear.model = function( FF, STRAINS ){
    # Calculate linear model of strain on expression
    # Not truly RNA-seq specific
    D = data.matrix(FF)            
    lm.p = rep(-1, dim(D)[1])
    r.p = rep(0, dim(D)[1])
    sums = rowSums(D)
    for(i in 1:dim(D)[1]){
        if( i %% 200==0 )
            print(i)
        if(sums[i]>0){
            L = lm(STRAINS~D[i,] )
            if( dim(summary(L)$coefficients)[1] == 2  ){
                # vector if linear model failed; no second row
                lm.p[i] = summary( L )$coefficients[2,4]
                r.p[i] = summary(L)$r.squared
            }
        }
    }
    idx.FVB = which(STRAINS==0)
    idx.F1 = which(STRAINS==1)
    idx.SPR = which(STRAINS==2)
    FVB = round( rowMeans(D[,idx.FVB], na.rm=T), 3)
    F1 =  round( rowMeans(D[,idx.F1], na.rm=T), 3)
    SPR = round( rowMeans(D[,idx.SPR], na.rm=T), 3) 
    mu10 = round( F1-FVB, 3) 
    mu21 = round( SPR-F1, 3) 
    diff = round( SPR-FVB, 3) 
    mus = data.frame(FVB, F1, SPR, diff, r.p, lm.p)
    mus = mus[ !( mus$FVB==0 & mus$F1==0 & mus$SPR==0),]
    mus$r.p = round(mus$r.p,3)
    mus$lm.p = signif(mus$lm.p,5)
    mus = mus[order(mus$lm.p),]
    mus
}

#######################
# BEGIN CRLMM FROM XY #
#######################

stripNormalize.from.XY = function (XY, useTarget = TRUE, verbose = TRUE){
    print("Modified stripnormalize function")
    if (useTarget) {
        objectsNeeded <- c("stripnum", "reference")
    }
    else{
        objectsNeeded <- "stripnum"
    }
    needToLoad <- !all(sapply(objectsNeeded, crlmm:::isLoaded))
    if (needToLoad) {
        pkgname = crlmm:::getCrlmmAnnotationName(annotation(XY))
        if (!require(pkgname, character.only = TRUE, quietly = !verbose)) {
            suggCall = paste("library(", pkgname, ", lib.loc='/Altern/Lib/Loc')", sep = "")
            msg = paste("If", pkgname, "is installed on an alternative location, please load it manually by using", suggCall)
            message(strwrap(msg))
            stop("Package ", pkgname, " could not be found.")
            rm(suggCall, msg)
       }
       if (verbose)
           message("Loading strip and reference normalization information.")
       crlmm:::loader("preprocStuff.rda", crlmm:::.crlmmPkgEnv, pkgname)
    }
	samplenms = sampleNames(XY)
	snpnms = featureNames(XY)
    stripnum = crlmm:::getVarInEnv("stripnum")
    if (useTarget)
       targetdist = crlmm:::getVarInEnv("reference")
    if (verbose) {
        message("Quantile normalizing ", ncol(XY), " arrays by ", max(stripnum), " strips.")
        if( useTarget ) 
            message("Using target")
        if (getRversion() > "2.7.0"){
            pb = txtProgressBar(min = 0, max = max(stripnum), style = 3)
        }
   }
   for (s in 1:max(stripnum)) {
       if (verbose) {
           if (getRversion() > "2.7.0")
               setTxtProgressBar(pb, s)
           else cat(".")
       }
        sel = stripnum == s
        subX <- as.matrix(assayData(XY)[["X"]][sel, ])
        subY <- as.matrix(assayData(XY)[["Y"]][sel, ])
        if (useTarget){
            tmp = normalize.quantiles.use.target(cbind(subX, subY), targetdist[[s]])
            #print(tmp[1:10,1:10])
            #normalized.X = matrix(as.integer(round(tmp[,1:(ncol(tmp)/2)]+16)), nrow(tmp), ncol(tmp)/2)
            #normalized.Y = matrix(as.integer(round(tmp[,(ncol(tmp)/2+1):ncol(tmp)]+16)), nrow(tmp), ncol(tmp)/2)
            #print(normalized.X[1:10,1:10])
        }
        else{
            tmp = normalize.quantiles(cbind(subX, subY))
        }   
        XY@assayData$X[sel, ] = matrix(as.integer(round(2^(tmp[,1:(ncol(tmp)/2)]+16))), nrow(tmp), ncol(tmp)/2)
        XY@assayData$Y[sel, ] = matrix(as.integer(round(2^(tmp[,(ncol(tmp)/2+1):ncol(tmp)]+16))), nrow(tmp), ncol(tmp)/2)                   
        rm(subX, subY, tmp, sel)
        gc()
    }
    if (verbose)
        cat("\n")
    XY@assayData$X = matrix(as.integer(XY@assayData$X), nrow(XY), ncol(XY))
    XY@assayData$Y = matrix(as.integer(XY@assayData$Y), nrow(XY), ncol(XY))
    sampleNames(XY) = samplenms
    featureNames(XY) = snpnms
    XY
}

genotype.CRLMM.from.XY = function(X, Y, fn.confs, fn.calls, fn.sa){
    library(Biobase)
    library(crlmm)
    library(preprocessCore)
    library(human660quadv1aCrlmm)    
#    X = read.table(fn.X, header=T, row.names=1, sep='\t')
#    Y = read.table(fn.Y, header=T, row.names=1, sep='\t')
#    X = as.matrix(X)
#    Y = as.matrix(Y)
    zeroes = X==0 | is.na(X) | Y==0 | is.na(Y)
    stopifnot(colnames(X)==colnames(Y))
    stopifnot(colnames(X)==colnames(zeroes))
    stopifnot(rownames(X)==rownames(Y))
    stopifnot(rownames(X)==rownames(zeroes))
    samplenames = colnames(X)
    X[is.na(X)] = 0
    Y[is.na(Y)] = 0
    sex=rep(2, length(samplenames) )
    samplesheet = data.frame( samplenames, sex, stringsAsFactors=F)
    # form NChannel set object
    print("Forming NChannelSet object...")
    XY =  new("NChannelSet", X=X, Y=Y, zero=zeroes, annotation="human660quadv1a", storage.mode = "environment")
    rm(X, Y, zeroes)
    gc()
    print("Completed NChannelSet object")    
    # quantile normalize the data using the crlmm:::stripNormalize2() function
    print("Quantile normalizing...")    
    XYsn = stripNormalize.from.XY(XY)
    rm(XY)
    print("Quantile normalizing complete")    
    # finally genotype
    print("Beginning genotyping...")    
    output = crlmmIllumina(XY=XYsn, stripNorm=FALSE, gender=samplesheet$sex,
                          seed=1, mixtureSampleSize=10^5, verbose=TRUE,
                          cdfName="human660quadv1a", sns=samplenames,
                          returnParams=TRUE)
    print("genotyping is complete")
    crlmmcalls = data.frame(calls(output)) # Calls
    crlmmconfs = data.frame(confs(output)) # Call confidence values
    names(crlmmconfs) = samplesheet$samplenames
    names(crlmmcalls) = samplesheet$samplenames    
    write.matrix(round(crlmmconfs, 3), fn.confs)
    sa = load.matrix(fn.sa)
    percent.conf.lt.95 = round( as.numeric(colSums(crlmmconfs<0.95) / dim(crlmmconfs)[1]), 2) 
    sa = cbind(sa, percent.conf.lt.95)
    write.matrix(sa, fn.sa)
    crlmmcalls[crlmmconfs<0.95] = NA
    crlmmcalls=crlmmcalls-1
    write.matrix(crlmmcalls, fn.calls)
    output
}


################
# BEGIN IMPUTE #
################

impute.prepare=function(calls, snp.ids, 
                  snp.chr, snp.pos, rs.ids, strand, snp.alleleA, snp.alleleB, 
                  target.chrom, target.pos.start, target.pos.end, 
                  dir.output, fn.out){
    # convert 0/1/2 genotype format into impute format.
    # SNP 1 : AA AA
    # SNP 2 : GG GT
    # is
    # SNP1 rs1 1000 A C 1 0 0 1 0 0
    # SNP2 rs2 2000 G T 1 0 0 0 1 0    
    # generate strand file for impute
    
    if( dim(calls)[1] != length(snp.ids) )
        stop("number of rows in calls != length of snp.ids")
    if( dim(calls)[1] != length(rs.ids) )
        stop("number of rows in calls != length of rs.ids")
    if( dim(calls)[1] != length(strand) )
        stop("number of rows in calls != length of strand")
    if( dim(calls)[1] != length(snp.alleleA) )
        stop("number of rows in calls != length of snp.alleleA")
    if( dim(calls)[1] != length(snp.alleleB) )
        stop("number of rows in calls != length of snp.alleleB")
    fn.calls = paste(dir.output, '/calls_', fn.out, '_impute.txt', sep='')
    fn.strand = paste(dir.output, '/strand_', fn.out, '_impute.txt', sep='')

    idx=which(snp.chr==target.chrom & (snp.pos>=target.pos.start & snp.pos<=target.pos.end))
    ga.target = data.frame( rs.id=rs.ids[idx], pos=snp.pos[idx], allele.a=snp.alleleA[idx], allele.b=snp.alleleB[idx], strand=strand[idx], stringsAsFactors=F)
    calls.target = calls[idx,]
    rownames(ga.target) = rownames(calls.target)    
    N.rs = dim(calls.target)[1]
    N.samp = dim(calls.target)[2]
    print(paste("Identified", N.rs, "SNPs in region."))
    print(paste("Identified", N.samp, "samples."))    
    calls.i = ga.target[,c("rs.id", "pos", "allele.a", "allele.b")]
    rownames(calls.i) = snp.ids[idx]
    vals.i = matrix(0, N.rs, N.samp*3)
    
    for( col in 1:N.samp){
        v = calls.target[,col]==1
        vals.i[ ,(col*3) - 1]=v
        v = calls.target[,col]==2
        vals.i[ ,(col*3) ]=v
    }   
    
    calls.i = cbind(calls.i, vals.i)
    write.table(calls.i, fn.calls, quote=F,sep=' ', col.names=F)
    write.table(ga.target[,c("pos", "strand")],  fn.strand, quote=F, sep=' ', col.names=F, row.names=F)
    print( paste( "Wrote calls to ", fn.calls) )
    print( paste( "Wrote strand to ", fn.strand) )
}


impute.parse.results = function( target.chrom, sample.ids, dir.output, fn.out ){
    # Read the result of a call to impute and convert into 0/1/2 genotype format 
    fn.imputed = paste(dir.output, '/calls_', fn.out, '_imputed.txt', sep='')
    fn.ga.imputed = paste(dir.output, '/gene_attributes_imputed_', fn.out, '.txt', sep='')
    fn.calls.imputed = paste(dir.output, '/calls_imputed_', fn.out, '.txt', sep='')
    
    print(paste("Reading calls from", fn.imputed))
    imp = read.table(fn.imputed, sep=' ')
    print(paste("Identified", dim(imp)[1], "imputed or directly sampled loci"))
    
    r = imp[,6:dim(imp)[2]]
    N.samp = length(r)/3
    
    r.rows = dim(r)[2]
    AA = r[, seq(1, r.rows-2, by=3)]
    AB = r[, seq(2, r.rows-1, by=3)]
    BB = r[, seq(3, r.rows, by=3)]
    calls.r = matrix(0, nrow=dim(imp)[1], ncol=N.samp)
    calls.r[AB>AA & AB>BB] = 1
    calls.r[BB>AA & BB>AB] = 2

    ga.imputed = imp[,c(3,4,5)]
    Chrom = rep( target.chrom, dim(ga.imputed)[1])
    names(ga.imputed) = c("pos", "allele.a", "allele.b")
    ga.imputed = cbind(Chrom, ga.imputed)
    rownames(ga.imputed) = imp$V2
    calls.r = data.frame(calls.r)
    rownames(calls.r) = rownames(ga.imputed)
    names(calls.r) = sample.ids
    n.zero = rowSums(calls.r==0)
    n.one = rowSums(calls.r==1)
    n.two = rowSums(calls.r==2)
    rows.keep = n.zero!=N.samp & n.one!=N.samp & n.two!=N.samp
    ga.imputed = ga.imputed[rows.keep,]
    calls.r = calls.r[rows.keep,]
    write.matrix(ga.imputed,  fn.ga.imputed)
    write.matrix(calls.r,  fn.calls.imputed)
    print(paste("Wrote attributes to ", fn.ga.imputed))
    print(paste("Wrote calls to ", fn.calls.imputed))
}


impute.generate.command = function(target.chrom, target.pos.start, target.pos.end, bin, dir.haplotypes, dir.output, fn.out){
    # Generate the command to run impute (binary located at bin)
    c.haploc = paste("-m ", dir.haplotypes, "/chr", target.chrom, "_genetic_map.txt", sep='')
    c.haps = paste("-h ", dir.haplotypes, "/chr", target.chrom, "_impute.haps", sep='')
    c.legend = paste("-l ", dir.haplotypes, "/chr", target.chrom, "_impute.legend", sep='')
    c.genos = paste("-g ", dir.output, "/calls_", fn.out, "_impute.txt", sep='')
    c.strand = paste("-strand_g ", dir.output, "/strand_", fn.out, "_impute.txt", sep='' )
    c.int = paste("-int ", target.pos.start, " ", target.pos.end, " -Ne 20000", sep='')
    c.o = paste("-o ", dir.output, "/calls_", fn.out, "_imputed.txt", sep='' )
    cmd = paste(bin, c.haploc, c.haps, c.legend, c.genos, c.strand, c.int, c.o )
    print(cmd)
}


####################
# BEGIN AFFYMETRIX #
####################

prepare.affy.dataset = function( dir.CEL, dir.OUT, fn.source.ga, affinity ){
    library(gcrma)
    affinity <- compute.affinities(affinity)
    setwd(dir.CEL)
    eset <- gcrma(ReadAffy(),affinity.info=affinity)
    sn=sampleNames(eset)
    for( i in 1:length(sn)){
       sn[i] = strsplit(sn[i],'.CEL.gz',fixed=T)[[1]][1]
    }
    sampleNames(eset)=sn
    M = exprs(eset)
    M = data.frame(M)
    names(M) = sn
    rownames(M) = featureNames(eset)
    M = M[,sort(names(M))]
    ga = load.matrix(fn.source.ga)
    expr.all.probes = M[order(rownames(M)),]
    ga = ga[order(rownames(ga)),]
    m = match(rownames(ga), rownames(expr.all.probes))
    expr.all.probes = expr.all.probes[m,]
    write.matrix(ga, paste(dir.OUT,'gene_attributes_HG_U133_Plus_2.na30.txt', sep='/'), F)
    write.matrix(round(M, 3), paste(dir.OUT,'expr_all_probes_raw.txt', sep='/') )
    #------------------------------------------------------------------------------------
    # Create expr_above_bg and expr_above_bg_refseq
    #------------------------------------------------------------------------------------
    means = rowMeans(expr.all.probes, na.rm=T)
    expr.mat = data.matrix(expr.all.probes)
    maxes = rep(0, length(means) )
      for(i in 1:length(maxes)){
      maxes[i] = max(expr.mat[i,], na.rm=T)
    }
    is.above.bg = means>3 | maxes>5
    is.above.bg[is.na(is.above.bg)] = F
    expr.a = round(expr.all.probes[is.above.bg, ], 3)
    ga.a = ga[is.above.bg, ]
    write.matrix(expr.a, paste(dir.OUT,'expr_above_bg.txt', sep='/')) 
    write.matrix(ga.a, paste(dir.OUT,'gene_attributes_above_bg.txt', sep='/'), F)
    is.refseq = ga.a$refseq.probes>0
    expr.ar = expr.a[is.refseq, ]
    ga.ar = ga.a[is.refseq, ]
    write.matrix(expr.a, paste(dir.OUT,'expr_above_bg_refseq.txt', sep='/')) 
    write.matrix(ga.a, paste(dir.OUT,'gene_attributes_above_bg_refseq.txt', sep='/'), F)
}


prepare.affy430_2.dataset = function( dir.CEL, dir.OUT, fn.source.ga, ga.id ){
    library(affy)
    setwd( dir.CEL )
    eset <- justRMA()
    fn.expr.all.probes = paste(dir.OUT,'/','expr_all_probes.txt', sep='')
    fn.ga.all.probes = paste(dir.OUT,'/','gene_attributes_', ga.id, '_all_probes.txt', sep='')    
    fn.expr.above.bg = paste(dir.OUT,'/','expr_above_bg.txt', sep='')
    fn.ga.above.bg = paste(dir.OUT,'/','gene_attributes_', ga.id, '_above_bg.txt', sep='')
    fn.expr.above.bg.refseq = paste(dir.OUT,'/','expr_above_bg_refseq.txt', sep='')
    fn.ga.above.bg.refseq = paste(dir.OUT,'/','gene_attributes_', ga.id, '_above_bg_refseq.txt', sep='')    
    sn=sampleNames(eset)
    for( i in 1:length(sn)){
       sn[i] = strsplit(sn[i],'.CEL',fixed=T)[[1]][1]
    }
    sampleNames(eset)=sn
    exprs(eset) = round(exprs(eset),3)
    write.exprs(eset, file=fn.expr.all.probes, sep='\t', quote=F) 
    expr = read.table(fn.expr.all.probes, check.names=F)
    ga = load.matrix(fn.source.ga)
    expr = expr[order(rownames(expr)),]
    ga = ga[order(rownames(ga)),]
    write.matrix(expr, fn.expr.all.probes) 
    write.matrix(ga, fn.ga.all.probes)
    means = rowMeans(expr)
    expr.mat = data.matrix(expr)
    maxes = rep(0, length(means) )
      for(i in 1:length(maxes)){
      maxes[i] = max(expr.mat[i,])
    }
    is.above.bg = means>4 | maxes>5
    expr = expr[is.above.bg, ]
    ga = ga[is.above.bg, ]
    write.matrix(expr, fn.expr.above.bg) 
    write.matrix(ga, fn.ga.above.bg)
    is.refseq.11 = ga$refseq.probes==11
    expr = expr[is.refseq.11, ]
    ga = ga[is.refseq.11, ]
    write.matrix(expr, fn.expr.above.bg.refseq) 
    write.matrix(ga, fn.ga.above.bg.refseq)
}

prepare.affy133.dataset = function( dir.CEL, dir.OUT, fn.source.ga, ga.id, filenames=NA ){
    library(affy)
    if( dir.CEL != "" ){
        setwd( dir.CEL )
        eset <- justRMA()
    }
    else{
        if( is.vector(filenames) & length(filenames)>0 )
            eset <- justRMA(filenames = filenames)
        else
            stop("Must pass either dir.CEL or vector of filenames")
    }
    fn.expr.all.probes = paste(dir.OUT,'/','expr_all_probes.txt', sep='')
    fn.ga.all.probes = paste(dir.OUT,'/','gene_attributes_', ga.id, '_all_probes.txt', sep='')    
    fn.expr.above.bg = paste(dir.OUT,'/','expr_above_bg.txt', sep='')
    fn.ga.above.bg = paste(dir.OUT,'/','gene_attributes_', ga.id, '_above_bg.txt', sep='')
    fn.expr.above.bg.refseq = paste(dir.OUT,'/','expr_above_bg_refseq.txt', sep='')
    fn.ga.above.bg.refseq = paste(dir.OUT,'/','gene_attributes_', ga.id, '_above_bg_refseq.txt', sep='')    
    sn=sampleNames(eset)
    for( i in 1:length(sn)){
       sn[i] = strsplit(sn[i],'.CEL',fixed=T)[[1]][1]
    }
    sampleNames(eset)=sn
    exprs(eset) = round(exprs(eset),3)
    write.exprs(eset, file=fn.expr.all.probes, sep='\t', quote=F) 
    expr = read.table(fn.expr.all.probes, check.names=F, header=T, row.names=1)
    ga = load.matrix(fn.source.ga)
    m = match.idx(rownames(expr), rownames(ga) )
    expr = expr[m$idx.A,]
    ga = ga[m$idx.B,]
    write.matrix(expr, fn.expr.all.probes) 
    write.matrix(ga, fn.ga.all.probes)
    means = rowMeans(expr)
    expr.mat = data.matrix(expr)
    maxes = rep(0, length(means) )
      for(i in 1:length(maxes)){
      maxes[i] = max(expr.mat[i,])
    }
    is.above.bg = means>4 | maxes>5
    expr = expr[is.above.bg, ]
    ga = ga[is.above.bg, ]
    write.matrix(expr, fn.expr.above.bg) 
    write.matrix(ga, fn.ga.above.bg)
    is.refseq.11 = ga$refseq.probes==11
    expr = expr[is.refseq.11, ]
    ga = ga[is.refseq.11, ]
    write.matrix(expr, fn.expr.above.bg.refseq) 
    write.matrix(ga, fn.ga.above.bg.refseq)
}


affymetrix.probe.dataset = function( Data, probe.id, cdf ){
    # returns a list with expr, sa, ga
    idx.probes = indexProbes(Data, 'pm', genenames=probe.id)[[1]]
    expr = data.frame( intensity(Data)[idx.probes, ] )
    n.probes= dim(expr)[1]
    n = names(expr)
    for( i in 1:length(n)){
       n[i] = strsplit(n[i],'.CEL',fixed=T)[[1]][1]
       n[i] = sub('.', '-', n[i], fixed=T)
    }
    names(expr) = n
    valid = rep(1, dim(expr)[2])
    sa = data.frame(valid)
    rownames(sa) = names(expr)
    xs = rep(0, n.probes);  ys = rep(0, n.probes)
    medians = rep(0, n.probes)
    for(i in n.probes){
        medians[i] = median( as.numeric(expr[i,]))
    	 xy=indices2xy(idx.probes[i], cdf=cdf)
	 xs[i] = xy[1]; ys[i] = xy[2];
    }
    ga = data.frame(medians, round(rowMeans(expr),2), xs, ys )
    rownames(ga) = idx.probes
    list(expr=expr, ga=ga, sa=sa)
}

describe.affymetrix.probes = function( Data, probe_id, cdf ){
    # Requires library affy and CDF library
    idx.probes = indexProbes(Data, 'pm', genenames=probe_id)[[1]]
    raw = intensity(Data)[idx.probes, ]
    n.probes= dim(raw)[1]
    medians = rep(0, n.probes )
    xs = rep(0, n.probes)
    ys = rep(0, n.probes)
    for(i in 1:dim(raw)[1]){
        medians[i] = median( as.numeric(raw[i,]))
    	 xy=indices2xy(idx.probes[i], cdf=cdf)
	 xs[i] = xy[1]; ys[i] = xy[2];
    }
    d = data.frame(medians, round(rowMeans(raw),2), xs, ys )
    rownames(d) = idx.probes
    names(d) = c('median', 'mean', 'x', 'y')
    d
}

### Aroma functions ###

aroma.load.chip.effect = function(path, dataset, combineAlleles){
    # Load and normalize all data in dataset at path.
    # Loads but does not re-process if data have already been processed.
    #memory.limit(3000)
    setwd(path)
    log <- Arguments$getVerbose(-8, timestamp=TRUE)
    cdf <- AffymetrixCdfFile$byChipType("GenomeWideSNP_6", tags="Full")
    csR <- AffymetrixCelSet$byName(dataset, cdf=cdf)
    acc <- AllelicCrosstalkCalibration(csR, model="CRMAv2")
    csC <- process(acc, verbose=log)
    bpn <- BasePositionNormalization(csC, target="zero")
    csN <- process(bpn, verbose=log)
    plm <- AvgCnPlm(csN, mergeStrands=TRUE, combineAlleles=combineAlleles)
    if (length(findUnitsTodo(plm)) > 0) {
      fitCnProbes(plm, verbose=log)
      fit(plm, verbose=log)
    }
    ces <- getChipEffectSet(plm)
    fln <- FragmentLengthNormalization(ces, target="zero")
    process(fln, verbose=log)
}

aroma.plot.region = function(cesN, cesNT, filename, chromosome, region){
    # Given matched chip estimate filesets cesN and cesNT, plot sample
    # filename on chromosome between region[1] and region[2] bases.
    # Confirm that files match
    n.files.N = length(cesN$files)
    n.files.T = length(cesNT$files)
    print(paste("#MESSAGE: plotting file",filename,"chr",chromosome,":",region[1],"-",region[2]))
    if(n.files.N != n.files.T){
        stop("Number of files for normal and tumor are not the same.")
    }
    idx = -1
    for(i in 1:n.files.N ){
        if( getName(cesN$files[[i]])[1] != getName(cesNT$files[[i]])[1] )
            stop(paste( "Mismatch of file names for file index",i ))
        if( getName(cesN$files[[i]])[1] == filename )
            idx = i
    }
    if(idx==-1)
        stop( "Requested filename not found")
    print(paste("#MESSAGE: Extracting chromosome from",n.files.N,"samples"))
    cdf = getCdf(cesN)
    gi <- getGenomeInformation(cdf)
    units <- getUnitsOnChromosome(gi, chromosome=chromosome, region=region)
    pos <- getPositions(gi, units=units)

    thetaTumor = extractTheta(getFile(cesNT, idx), units=units)
    thetaRef = extractTheta(getFile(cesN, idx), units=units)
    C <- log2( thetaTumor/thetaRef )
    par(mar=c(3,4,2,1)+0.1)
    plot(pos/1e6, C, pch=".", cex=3, ylim=c(-2,4))
    stext(side=3, pos=0, getName(cesN$files[[idx]]))
    stext(side=3, pos=1, paste("Chr",chromosome))
}

aroma.extract.array.log2 = function(cesN, cesT, filename=NULL, chromosome=NULL){
    n.files.N = length(cesN$files)
    n.files.T = length(cesT$files)
    if(n.files.N != n.files.T){
        stop("Number of files for normal and tumor are not the same.")
    }
    indices = vector(mode="numeric")
    col.names = vector(mode="character")
    for(i in 1:n.files.N ){
        if( getName(cesN$files[[i]])[1] != getName(cesT$files[[i]])[1] )
            stop(paste( "Mismatch of file names for file index",i ))
        cur.name = getName(cesN$files[[i]])[1]
        if( is.null(filename) ){
            indices[ length(indices)+1 ] = i
            col.names[length(col.names)+1 ] = cur.name
        }
        else if(cur.name == filename){
            indices[ length(indices)+1 ] = i
            col.names[length(col.names)+1 ] = cur.name
        }
    }
    if(length(indices)==0)
        stop( "Requested filename not found")
    cdf = getCdf(cesN)
    gi <- getGenomeInformation(cdf)
    if( is.null(chromosome) )
        units = (1:nbrOfUnits(cdf))
    else
        units <- getUnitsOnChromosome(gi, chromosome)
    unit.types = getUnitTypes(cdf, units)
    values = list()
    for(i in 1:length(indices) ){
        thetaTumor = extractTheta(getFile(cesT, indices[i]), units=units)
        thetaRef = extractTheta(getFile(cesN, indices[i]), units=units)
        values[[i]] = round(log2(thetaTumor/thetaRef)[unit.types==2],3)
    }
    df = data.frame(values)
    rownames(df) = getUnitNames(cdf, units)[unit.types==2]
    names(df) = col.names
    df
}

##############################
# BEGIN PLOTTING AND DISPLAY #
##############################

plot.discretization = function(M, lbls=NULL, method="SD", bounds=0.5, x_lbl=NULL, y_lbl="expression"){
    # M should be a data frame with one row

    if( method != "SD"){
        stop("method must be SD")
    }
    if( !is.data.frame(M) || dim(M)[1] != 1 ){
        stop("M is not a data frame with one row")
    }
    if(is.null(x_lbl)){
        x_lbl = rownames(M)[1]
    }
    V = as.numeric(M)
    N = length(V)
    dis = vector(mode="integer",length=N)
    color_v = vector(mode="character",length=N)
    if( is.null(lbls) )
        lbls = (1:N)
    if( method=="SD" ){
        mu = mean(V, na.rm=T)
        std = sd(V, na.rm=T)
        limit.lower = mu - (bounds*std)
        limit.upper = mu + (bounds*std)

        dis[which(V>limit.upper)] = 1
        dis[which(V<limit.lower)] = -1
        plot.title = paste("Mean", round(mu,digits=2), " ", method,"+/-",bounds  )
    }
    color_v[ dis==-1 ] = "green"
    color_v[ dis==1 ] = "red"
    color_v[ dis==0 ] = "gray"
    plot(V, col=color_v, pch=19,main=plot.title, cex=1.25, axes=F, xlab=x_lbl, ylim=c(1,ceiling(max(V, na.rm=T))), ylab=y_lbl)
    lines(c(0,N), c(limit.lower, limit.lower), col='gray')
    lines(c(0,N), c(limit.upper, limit.upper), col='gray')
    par(ps=8)
    axis(2, at=(1:15) )
    axis(1, at=(1:N), labels=lbls, las=2)
}

plot.percent.altered.by.locus= function(M, chr.list, upper.bound=0.5, lower.bound=-0.5, y.min=0, y.max=0, color=T){
    # M is a data frame with some sort of value where we want to count the percentage
    # that exceed a given bounds (typically SNP amplification or aCGH data)
    # chr is the chromosome for each row in M
    # We assume that rows are sorted along the genome
    # Simple version that plots each probe as a bar, regardless of space intervening
    N.genes = dim(M)[1]
    N.samples = dim(M)[2]
    if( length(chr.list) != N.genes ){
        stop("length of chr not equal to rows in M")
    }
    if(color){
        color.up = "blue"
        color.dn = "red"
    }
    else{
        color.up = "black"
        color.dn="black"
    }
    amps = rowSums(M>=upper.bound,na.rm=T) / N.samples * 100
    dels = rowSums(M<=lower.bound,na.rm=T) / N.samples * 100
    amps[is.na(amps)]=0
    dels[is.na(dels)]=0
    
    if( y.max==0 )
        y.max = nearest.ten( max(amps) )
    if( y.min==0 )
        y.min = nearest.ten( max(dels) )
    #colors = vector(mode="integer", N.genes)
    #cur.color = 4
    cur.chr = chr.list[1]
    chr.breaks = vector()
    chr.lbl.idx = vector()
    for( i in 1:N.genes ){
        if( chr.list[i] != cur.chr ){
            #if(cur.color==4)
            #    cur.color=9
            #else
            #    cur.color=4
            cur.chr = chr.list[i]
            chr.breaks[length(chr.breaks)+1] = i
        }
        #colors[i]=cur.color
    }
    chr.lbl.idx[1] = chr.breaks[1] / 2
    for(i in 2:length(chr.breaks) ){
        prev.idx = chr.breaks[i-1]
        curr.idx = chr.breaks[i]
        chr.lbl.idx[i] = prev.idx + ( (curr.idx-prev.idx) / 2 )
    }
    chr.lbl.idx[i+1] = chr.breaks[i] + (( length(chr.list) - chr.breaks[length(chr.breaks)]) /2 )
    main.lbl = paste("SNP Arrays Percent Amplified, +/- ", upper.bound,",",lower.bound,sep="")
    plot(c(0,N.genes), c(-y.min,y.max), xaxs="i", pch='.',col='white', main=main.lbl, xlab="SNP index", ylab="Percentage", axes=F, xlim=c(0, N.genes))
    for(i in 1:N.genes){
        segments(i,0,i,amps[i], col=color.up)
        segments(i,0,i,-1*dels[i],col=color.dn)
    }
    abline(h=0, col='black')
    max.chr = max(unique(chr.list))
    chr.labels = as.character(sort(unique(chr.list)))
    if( length(chr.labels)==20 )
        chr.labels[20] = 'X'
    axis(side=1, xlim=c(0, N.genes), at=chr.lbl.idx, labels=chr.labels, tick=F, cex.axis=1.25,las=1 )
    axis(side=2, at=seq(-y.min,y.max,by=10), labels=seq(-y.min,y.max,by=10), tick=T, cex.axis=1.5,las=1 )
    chr.breaks = append(0, chr.breaks)
    chr.breaks = append(chr.breaks, length(chr.list)) 
    for(i in 1:length(chr.breaks)){
        segments(chr.breaks, rep(-1*y.min, length(chr.breaks)), y1=y.max, col='black')
    }
    data.frame(amps,dels)
}


plot.percent.altered.by.tumor = function(M, upper.bound=0.3, lower.bound=-0.3){
    # for each sample (column) plot percentage of all in M gt or lt bounds
    N.genes = dim(M)[1]
    N.samples = dim(M)[2]
    amps = colSums(M>=upper.bound, na.rm=T) / N.genes * 100
    dels = colSums(M<=lower.bound, na.rm=T) / N.genes * 100
    y.max = nearest.ten( max(amps,na.rm=T) )
    y.min = nearest.ten( max(dels,na.rm=T) )
    d = data.frame(amps, dels)
    o = order(d$amps + d$dels)
    d = d[o,]
    barplot(amps[o], ylim=c(-1*y.min, y.max), col='blue', main='% Genomic alteration by sample', axes=F, names.arg=rep('', length(amps)))
    barplot(-1*dels[o], add=T, col='red', axes=F, names.arg=rep('', length(amps)))
    axis(side=2, at=seq(-y.min,y.max,by=10), labels=seq(-y.min,y.max,by=10), tick=T, cex.axis=0.75,las=1 )
    for(n in seq(from=-1*y.min,to=y.max,by=5) )
        abline(h=n, col='black')
    d
}

interpolate.alterations = function(M.in, chromosomes, loc_begin, binsize=1000000){
    # given aCGH matrix M and locations of markers chrs, locs, estimate 
    # values for markers at evenly-spaced intervals
    n.bins = 0
    unique.chromosomes = unique(chromosomes)
    for(c.idx in 1:length(unique.chromosomes)){
        chromosome = unique.chromosomes[c.idx]
        max.loc = max(loc_begin[chromosomes==chromosome], na.rm=T)
        n.bins = n.bins + floor(max.loc/binsize) + 1
    }
    M = matrix(0, nrow=n.bins, ncol=dim(M.in)[2])
    chrs.inter = rep(0, n.bins )
    locs.inter = rep(0, n.bins )
    for(s.idx in 1:dim(M)[2]){
        print(s.idx)
        bin.idx=1
        for(c.idx in 1:length(unique.chromosomes)){
            chromosome = unique.chromosomes[c.idx]
            vals.chr = M.in[ chromosomes==chromosome, s.idx ]
            locs.chr = loc_begin[ chromosomes==chromosome ]
            max.loc = max(locs.chr, na.rm=T)
            cur.loc=0
            while(cur.loc<max.loc){
                d.to.cur = abs(cur.loc - locs.chr) # distances from all markers to current locus
                closest.bac.idx = which(d.to.cur==min(d.to.cur) )[1] # could be a tie
                M[bin.idx,s.idx] = vals.chr[closest.bac.idx]
                chrs.inter[bin.idx] = chromosome
                locs.inter[bin.idx] = cur.loc
                bin.idx = bin.idx+1 
                cur.loc = cur.loc + binsize
            }
        }
    }
    list( 'values'=M, 'chromosome'=chrs.inter, 'location'=locs.inter )
}

estimate.percents.altered = function( M.in, chrs, locs, bound=0.3, binsize=1000000 ){
    # binsize is the length in bases of each bin, not the number of bins
    M = matrix(0, nrow=dim(M.in)[1], ncol=dim(M.in)[2])
    M[M.in>bound] = 1
    M[M.in < (-1*bound)] = -1
    chrs = unique(chrs)
    percent.altered = rep(0, dim(M)[2] )
    n.bins = 0
    for(s.idx in 1:dim(M)[2]){
        print(s.idx)
        for(chromosome in 1:length(chrs)){
            vals.chr = M[ chrs==chromosome, s.idx ]
            locs.chr = locs[chrs==chromosome]
            max.loc = max(locs.chr, na.rm=T)
            cur.loc=0
            while(cur.loc<max.loc){
                if( s.idx==1 )
                    n.bins = n.bins + 1
                d.to.cur = abs(cur.loc - locs.chr) # distances from all markers to current locus
                closest.bac.idx = which(d.to.cur==min(d.to.cur) )[1] # could be a tie
                if( !is.na( vals.chr[closest.bac.idx] ) && vals.chr[closest.bac.idx] != 0 ){
                    percent.altered[s.idx] = (percent.altered[s.idx])+1
                }
                cur.loc = cur.loc + binsize
            }
        }
    }
    percent.altered = percent.altered / n.bins
    round(percent.altered,2)
}



generate.simulated.DC.data = function(nA=100, nB=100, step.interval=0.01, max.sd=100){
    # Given the size of simulated data sets nA and NB, create simulation raw
    # data that has differential correlation for every bin between 0 and 2 
    # at intervals of 0.01. If the data generation fails to completely fill the
    # bins, it stops with an error. Otherwise, it returns the observed correlations
    # and the data that will generate those correlations
    steps = seq(from=1, to=max.sd, by=step.interval)
    obs.diff = rep(0, length(steps) )
    max.diff = rep(0, length(steps) )
    data.A1 = matrix( 0, nrow=length(steps), ncol=nA )
    data.A2 = matrix( 0, nrow=length(steps), ncol=nA )
    data.B1 = matrix( 0, nrow=length(steps), ncol=nB )
    data.B2 = matrix( 0, nrow=length(steps), ncol=nB )
    for( i in 1:length(steps) ){
        A1 = (1:nA)
        A2 = A1+rnorm(nA, mean=0, sd=steps[i])
        B1 = (1:nB)
        B2 = B1[nB:1]
        data.A1[i,] = A1
        data.B1[i,] = B1
        data.A2[i,] = A2
        data.B2[i,] = B2
        obs1 = as.numeric(cor.test(A1,A2, method="spearman")$estimate)
        obs2 = as.numeric(cor.test(B1,B2, method="spearman")$estimate)
        obs.diff[i] = abs(obs1-obs2)
    }
    targets = seq(from=0.01,to=2,by=0.01)
    idx = rep(-1, length(targets) )
    for(i in 1:length(targets) ){
        target = targets[i]
        hits = which(obs.diff > (target-0.01) & obs.diff <= target )
        if(length(hits)>0 ){
            idx[i] = hits[1]
        } 
    }
    if( sum( idx==-1 ) ){
        # Try to fill in obs from other direction
        for(i in seq(from=1, to=1000, by=0.1) ){
            if( nA>=nB ){
                A1 = (1:nA) 
                A2 = A1+rnorm(nA, mean=0, sd=1)
                B1 = A1[1:nB]
                B2 = ( A2+rnorm(nA, mean=0, sd=i) )[1:nB]
            }
            else{
                B1 = (1:nB) 
                B2 = B1+rnorm(nB, mean=0, sd=1)
                A1 = B1[1:nA]
                A2 = ( B2+rnorm(nB, mean=0, sd=i) )[1:nA]
            }
            obs1 = as.numeric(cor.test(A1,A2, method="spearman")$estimate)
            obs2 = as.numeric(cor.test(B1,B2, method="spearman")$estimate)
            diff = abs(obs1-obs2)
            for(t in 1:length(targets) ){
                if( idx[t]==-1 & diff > (targets[t]-0.01) & diff <= targets[t] ){
                    data.A1 = rbind( data.A1, A1 )
                    data.A2 = rbind( data.A2, A2 )
                    data.B1 = rbind( data.B1, B1 )
                    data.B2 = rbind( data.B2, B2 )                                                            
                    obs.diff = append( obs.diff, diff )
                    idx[t] = length(obs.diff)
                }
            }
        }
    }
    if( sum( idx==-1 ) ){
        data.A1 = rbind( data.A1, rep(NA, nA) )
        data.A2 = rbind( data.A2, rep(NA, nA) )
        data.B1 = rbind( data.B1, rep(NA, nB) )
        data.B2 = rbind( data.B2, rep(NA, nB) )
        obs.diff = append(obs.diff, NA)
        idx[idx==-1] = dim(data.A1)[1]
        warning("Failed to fill all targets")
    }
    out = list( "obs.diff" = obs.diff, "A1"=data.A1,"A2"=data.A2,"B1"=data.B1,"B2"=data.B2 )
    out$A1 = data.A1[idx,]
    out$B1 = data.B1[idx,]
    out$A2 = data.A2[idx,]
    out$B2 = data.B2[idx,]
    out$obs.diff = out$obs.diff[idx]
    out
    
    #plot(out$obs.diff, seq(from=0.01, to=2, by=0.01), pch='.')
}


DC.within.labels = function( A1, A2, B1, B2, n.perms=1000 ){
    # Given values for probes 1 and 2 in two classes A and B, calculate a P value for differential correlation
    # by sampling within each class (as opposed to swapping labels)
    obs1 = as.numeric(cor.test(A1,A2, method="spearman")$estimate)
    obs2 = as.numeric(cor.test(B1,B2, method="spearman")$estimate)
    obs.diff = abs(obs1-obs2)
    stats = rep(0, n.perms)
    for(perm in 1:n.perms){
        A1s = sample(A1)
        B1s = sample(B1)
        A2s = sample(A2)
        B2s = sample(B2)
        C1 = cor.test(A1s,A2s, method="spearman")
        C2 = cor.test(B1s,B2s, method="spearman")
        stats[perm] = abs(C1$estimate-C2$estimate)
    }
    perm.diff = round(stats[ order(stats, decreasing=T) ], 3)
    perm.gt.obs = which(perm.diff>obs.diff)
    if( length(perm.gt.obs)>0 )
        pval = max(perm.gt.obs) / length(perm.diff) 
    else
        pval = 0
    list( pval=pval, perm.diff=perm.diff, obs.diff=obs.diff )
}


DC.between.labels = function( A1, A2, B1, B2, n.perms=1000 ){
    # Given values for probes 1 and 2 in two classes A and B, calculate a P value for differential correlation
    # by sampling between classes
    obs1 = as.numeric(cor.test(A1,A2, method="spearman")$estimate)
    obs2 = as.numeric(cor.test(B1,B2, method="spearman")$estimate)
    obs.diff = abs(obs1-obs2)
    stats = rep(0, n.perms)
    for(perm in 1:n.perms){
        probe1 = sample(append(A1, B1))
        probe2 = sample(append(A2, B2))
        A1s = probe1[ 1:length(A1) ]
        B1s = probe1[ (length(A1)+1) : length(probe1) ]
        A2s = probe2[ 1:length(B1) ]
        B2s = probe2[ (length(B1)+1) : length(probe2) ]
        C1 = cor.test(A1s,A2s, method="spearman")
        C2 = cor.test(B1s,B2s, method="spearman")
        stats[perm] = abs(C1$estimate-C2$estimate)
    }
    perm.diff = round(stats[ order(stats, decreasing=T) ], 3)
    perm.gt.obs = which(perm.diff>obs.diff)
    if( length(perm.gt.obs)>0 )
        pval = max(perm.gt.obs) / length(perm.diff) 
    else
        pval = 0
    list( pval=pval, perm.diff=perm.diff, obs.diff=obs.diff )
}


calls.to.PLINK = function( calls, snp.ids, sample.ids, sexes, CHR, POS, fn.tped, fn.tfam ){
    # write genotype calls to PLINK TPED and TFAM files
    # PLINK sex codes for are male 1, females 2
    #
    plink.tped = data.frame(CHR, snp.ids, rep(0, length(CHR)), POS, stringsAsFactors=F)
    tped.odd = matrix('0', nrow=dim(calls)[1], ncol=dim(calls)[2] )
    tped.even = matrix('0', nrow=dim(calls)[1], ncol=dim(calls)[2] )
    tped.odd[ calls==0 | calls==1  ] = 'A'
    tped.odd[ calls==2 ] = 'B'
    tped.even[ calls==1 | calls==2 ] = 'B'
    tped.even[ calls==0 ] = 'A'
    tped = matrix('0', nrow=dim(calls)[1], ncol=(2*dim(calls)[2]) )
    even.idx = seq(2,dim(tped)[2], by=2)
    odd.idx = seq(1,dim(tped)[2], by=2)
    tped[,even.idx] = tped.even
    tped[,odd.idx] = tped.odd
    plink.tped = cbind( plink.tped, tped )
    write.table(plink.tped, fn.tped, sep='\t', quote=F, col.names=F, row.names=F)
    
    # WRITE TFAM  
    FAMILY_ID = rep(0, length(sample.ids))
    FATHER_ID = FAMILY_ID
    MOTHER_ID = FAMILY_ID
    DUMMY_PHENO = FAMILY_ID
    match = match.idx( names(calls), sample.ids )
    if( dim(match)[1] != length(sample.ids) )
        stop("mismatch between matches sample.ids and calls")
    plink.tfam = data.frame( FAMILY_ID, sample.ids[ match$idx.B ], FATHER_ID, MOTHER_ID, sexes, DUMMY_PHENO, stringsAsFactors=F )
    write.table(plink.tfam, fn.tfam, sep='\t', quote=F, col.names=F, row.names=F)
}


plot.alterations.on.chr = function( M, name, marker.interval.size=0.1,  upper.bound=0.3, lower.bound=0.3, y.min=NA, y.max=NA, color=T){
    # Plot one chromosome, label X axis with Mb.
    # marker.interval.size is the number of Mb. separating each marker.
    if( !is.vector(M) )
        stop("M should be a vector")
    N.genes = length(M)
    vals = as.numeric(M)
    if(is.na(y.max))
        y.max = ceiling( max(vals,na.rm=T) )
    if(is.na(y.min))
        y.min = floor( min(vals,na.rm=T) )
    if(color){
        cols = rep('lightblue', N.genes)
        cols[vals<0] = 58
        cols[vals > upper.bound] = 'blue'
        cols[vals < lower.bound] = 'red'
    }
    else{
        cols = rep('black', N.genes)
    }
    chr.breaks = vector()
    plot(c(0,N.genes), c(y.min,y.max), xaxs="i", pch='.',col='white', main=name, xlab="Mb.", ylab="Log2", axes=F, xlim=c(0, N.genes))
    chr.in.mb = ceiling( marker.interval.size * length(M) )
    x = 0
    chr.lbl.idx = c( x * (1/marker.interval.size) )
    chr.labels = c( as.character(x) )
    x = x + 10
    while(x <= chr.in.mb ){
        chr.lbl.idx = c( chr.lbl.idx, x * (1/marker.interval.size) )
        chr.labels = c( chr.labels, as.character(x) )
        x = x + 10
    }
    axis(side=1, xlim=c(0, N.genes), at=chr.lbl.idx, labels=chr.labels, tick=T, cex.axis=1.5,las=1 )
    for(i in 1:N.genes){
        segments(i,0,i,vals[i], col=cols[i])
    }
    axis(side=2, at=seq(y.min,y.max,by=0.5), labels=seq(y.min,y.max,by=0.5), tick=T, cex.axis=1.5,las=1 )
}



plot.alterations = function(M, name, chr.list, upper.bound=0.3, lower.bound=-0.3, y.min=NA, y.max=NA, color=T){
    # Plot the whole genome, label x axis with chromosomes
    if( !is.vector(M) )
        stop("M should be a vector")
    N.genes = length(M)
    vals = as.numeric(M)
    if(is.na(y.max))
        y.max = ceiling( max(vals,na.rm=T) )
    if(is.na(y.min))
        y.min = floor( min(vals,na.rm=T) )
    if(color){
        cols = rep('gray', N.genes)
        #cols[vals<0] = 'gray'
        cols[vals > upper.bound] = 'blue'
        cols[vals < lower.bound] = 'red'
    }
    else{
        cols = rep('black', N.genes)
    }
    chr.breaks = vector()
    cur.chr = chr.list[1]
    for( i in 1:N.genes ){
        if( chr.list[i] != cur.chr ){
            chr.breaks[length(chr.breaks)+1] = i
            cur.chr = chr.list[i]
        }
    }
    ylabel='Log2'
    plot(c(0,N.genes), c(y.min,y.max), xaxs="i", pch='.',col='white', main=name, xlab="Chromosome", ylab=ylabel, axes=F, xlim=c(0, N.genes))
    
    max.chr = max(unique(chr.list))
    chr.labels = as.character(rep(1:max.chr))
    chr.lbl.idx = vector()
    chr.lbl.idx[1] = chr.breaks[1] / 2
    for(i in 2:length(chr.breaks) ){
        prev.idx = chr.breaks[i-1]
        curr.idx = chr.breaks[i]
        chr.lbl.idx[i] = prev.idx + ( (curr.idx-prev.idx) / 2 )
    }
    chr.lbl.idx[i+1] = chr.breaks[i] + (( length(chr.list) - chr.breaks[length(chr.breaks)]) /2 )
    chr.labels[length(chr.labels)] = 'X'
    axis(side=1, xlim=c(0, N.genes), at=chr.lbl.idx, labels=chr.labels, tick=F, cex.axis=1.25,las=1 )
    for(i in 1:N.genes){
        segments(i,0,i,vals[i], col=cols[i])
    }
    chr.breaks = append(0, chr.breaks)
    chr.breaks = append(chr.breaks, length(chr.list))
    for(i in 1:length(chr.breaks)){
        segments(chr.breaks, rep(y.min, length(chr.breaks)), y1=y.max, col='gray')
    }
    axis(side=2, at=seq(y.min,y.max,by=0.5), labels=seq(y.min,y.max,by=0.5), tick=T, cex.axis=1,las=1 )
}


plot.rows=function( M, sa=NULL, legend.labels=NULL, group.by=NULL, sorted=F, plot.type='l', no.plot=F, legend.xy = c()){
    # sa is an attribute file data frame that describes the values in M
    # M is a numeric value data frame with one row per gene to plot
    # legend.labels are friendly names for rows; if NULL, rownames(M) is used
    # group.by optionally orders samples by a column pulled from sa, if NULL, sorted by identifier
    #
    samples.M = as.character(names(M))
    if( is.null(sa) ){
        samples.sa = samples.M
        group.by = NULL
    }
    else{
        samples.sa = as.character(rownames(sa))
    }
    sa.in.M = contains( samples.sa, samples.M )
    samples = samples.sa[ sa.in.M ] # samples
    if( is.null(legend.labels) ){
        legend.labels = rownames(M) # default label is rowname
    }
    else{
        if( length(legend.labels) != dim(M)[1] ){
            stop("legend.labels not same length as number of rows in M")
        }
    }
    if( is.null(group.by) ){
        if( !sorted ){
            x.axis.labels = samples
            x.bottom.label = ''
        }
        else{
            sorter = data.frame( vals=as.numeric(M[1,]), ids=samples, stringsAsFactors=F)
            sorter = sorter[ order(sorter$vals),]
            m = match.idx( sorter$ids, samples )
            M = M[,m$idx.B]
            x.axis.labels = as.character(sorter[,2])
            x.bottom.label = ''
        }
    }
    else{
        x.bottom.label = paste( "Grouped by", group.by)
        idx.groupby=which(colnames(sa)==group.by)
        if( length(idx.groupby) == 0){
            stop("group.by value not found in sa")
        }
        values.for.groupby = as.character(sa[,idx.groupby][sa.in.M])
        if(sum(is.na(values.for.groupby))>0){
            values.for.groupby[is.na(values.for.groupby)] = ' NA' # Convert <NA> to string saying "NA"
        }
        if(sorted){
            sorter = data.frame( vals.group=values.for.groupby, vals=as.numeric(M[1,]), ids=samples, stringsAsFactors=F)
            sorter = sorter[ order(sorter$vals.group, sorter$vals),]
            samples.reordered = sorter$ids
        }
        else{
            sorter = data.frame( vals=values.for.groupby, ids=samples, stringsAsFactors=F)
            sorter = sorter[ order(sorter$vals),]        
            samples.reordered = sorter$ids
        }   
        x.axis.labels = as.character(sorter[,1])
        m = match.idx( sorter$ids, samples )
        M = M[,m$idx.B]
    }
    if( no.plot )
        M
    else
        plot_multiple(M, lbls=x.axis.labels, x_lbl=x.bottom.label, legend.labels, plot.type, legend.xy = legend.xy)
}


plot.by.identifiers=function( dataset, probe.list, group.by=NULL, sorted=F, plot.type='l', legend.xy=c() ){
    # Convenience function to allow us to plot with a simple probe list.
    # dataset is list of {expr, sa, ga}
    if( is.null(dataset$expr) | is.null(dataset$ga) | is.null(dataset$sa))
        stop("dataset is incomplete; must contain expr, ga, and sa")
    if( length(rownames(dataset$sa)) != length( names(dataset$expr) ) )
        stop("dimensions of sa and expr not compatible")
    if( sum(rownames(dataset$sa) != names(dataset$expr))>0 )
        stop("rownames of sa != names of expr")
    idx.probes = rep(-1, length(probe.list))
    legend.labels = rep("", length(probe.list))
    idx.symbol = which(names(dataset$ga)=="symbol")    
    for( i in 1:length(probe.list) ){
        idx = which(rownames(dataset$expr)==probe.list[i] )
        if( length(idx) != 1 )
            stop(paste("Probe",probe.list[i],"not found."))
        idx.probes[i] = idx
        legend.labels[i] = paste( dataset$ga[which(rownames(dataset$ga)==probe.list[i]), idx.symbol ], probe.list[i])
    }
    plot.rows( dataset$expr[idx.probes,], dataset$sa, legend.labels, group.by, sorted, plot.type, legend.xy=legend.xy )
    abline(0,0)
}


plot.raw.probe.values2 = function(Data, probe_id){
    idx.probes = indexProbes(Data, 'pm', genenames=probe_id)[[1]]
    raw = intensity(Data)[idx.probes, ]
    pch.pool = c(15,16,17,18,19,0,1,2,3,4,5,6)
    #dfsort = data.frame( as.numeric(raw[1,]), 1:dim(raw)[2] )
    #dfsort = dfsort[order(dfsort[,1]),]
    #idx.order = dfsort[,2]
    #raw = raw[,idx.order]
    n_probes = dim(raw)[1]
    n_samples = dim(raw)[2]
    X = (1:n_samples)
    col.used = c(1)
    pch.used = c(pch.pool[1]);
    ymax = max(c(1,ceiling(max(raw, na.rm=T))))
    par(mar = c(7, 4, 4, 2) + 0.1)
    plot(X, as.numeric(raw[1,]), type='p', axes=F, col=1, pch=pch.pool[1], cex=1.25, xlab=probe_id, ylim=c(1,ymax), ylab='', lwd=2)
    for( i in 2:n_probes){
        new.pch = i %% length(pch.pool)
        lines(X, as.numeric(raw[i,]), col=i, lwd=2)
        points(X, as.numeric(raw[i,]), pch=new.pch, col=i, lwd=2)
        pch.used[length(pch.used)+1] = new.pch
        col.used[length(col.used)+1] = i
    }
    n = names(data.frame(raw))
    for( i in 1:length(n)){
       n[i] = strsplit(n[i],'.CEL',fixed=T)[[1]][1]
    }
    par(ps=8)
    #axis(2, at=(1:ymax),las=1 )
    axis(1, at=(1:length(n)), labels=n, las=2, )
    legend(1, ymax-1, rownames(raw), col=col.used, lty=1, lwd=2, pch=pch.used, cex=1,y.intersp=0.75)
}



plot_multiple=function(V, lbls, x_lbl, probe_label_list, plot.type='l', legend.xy = c()){
    # This is the utility function that just plots what it's given
    # pass plot.type=='l' for lines, 'p' for points
    #pch.pool = c(15,16,17,18,19,0,1,2,3,4,5,6)
    pch.pool = c(19)
    col.pool = c('black', 'chartreuse4', 'darkblue', 'firebrick3', 'darkorange', 'blue', 'gray48', 'burlywood')
    n_probes = dim(V)[1]
    n_samples = dim(V)[2]
    X = (1:n_samples)
    ymax = max(c(1,ceiling(max(V, na.rm=T))))
    ymin = floor( min(V, na.rm=T) )
    if(ymin>1)
        ymin=1
    pch.used = c(pch.pool[1]); 
    col.used = c( col.pool[1] )
    plot(X, as.numeric(V[1,]), col=1, type=plot.type, pch=pch.pool[1], cex=0.5, axes=F, xlab=x_lbl, ylim=c(ymin,ymax), ylab='', lwd=2)
    iter = 2
    if( n_probes>1 ){
        for( i in 2:n_probes){
            if( iter==n_probes+1 ){ iter==1 }
            cur.color = col.pool[iter]
            iter = iter+1
            if(plot.type=='l')
                lines(X, as.numeric(V[i,]), col=cur.color, lwd=2)
            else{
                new.pch = i %% length(pch.pool)
                points(X, as.numeric(V[i,]), pch=19, col=cur.color, lwd=2, cex=0.5)
            }
            col.used[length(col.used)+1] = cur.color
            pch.used[length(pch.used)+1] = 19
        }
    }
    par(ps=8)
    axis(2, at=(ymin:15), las=1)
    axis(1, at=(1:length(lbls)), labels=lbls, las=2)
    legend.x = 1
    legend.y = ymax
    if(length(legend.xy)==2){
        legend.x=legend.xy[1]
        legend.y=legend.xy[2]
    }
    if(plot.type=='l')
        legend(legend.x, y=legend.y, probe_label_list, col=col.used, lty=1, lwd=2, cex=0.75,y.intersp=0.5)
    else
        legend(legend.x, y=legend.y,  probe_label_list, col=col.used, lty=1, lwd=2, pch=pch.used, cex=1,y.intersp=0.75)
}


plot.geno.expr = function(ids, expr.g, expr.e, probe.g, probe.e, main=NULL, show.lm.pval=FALSE){
    # used to plot genotype or SNP copy number vs. real-valued data
    vals = align.geno.pheno(ids, expr.g, expr.e, probe.g, probe.e)
    u = sort(unique( vals$val.G ))
    if(is.null(main))
        main = probe.g
    main = paste(main, '\n', probe.e,'-',probe.g)
    special = F
    significance = NULL
    
    L = lm(as.numeric(vals$val.P) ~ as.numeric(vals$val.G) )
    pval = signif( (summary( L )$coefficients[2,4]), 2)
    rval = signif( (summary( L )$adj.r.squared), 2)
    if( show.lm.pval ){
        significance = paste( "pval", pval, ", adj. r2", rval )
    }
    if( length(u)==2 ){
        if( (u[1]==0 && u[2]==1) || (u[1]==1 && u[2]==2) ){
            special = T
            if( u[1]==0 ){
                gnames = c('Hom', 'Het')
                col1 = as.numeric(vals$val.P[ vals$val.G==0 ])
                col2 = as.numeric(vals$val.P[ vals$val.G==1 ])
            }
            else{
                gnames = c('Het', 'Hom')
                col1 = as.numeric(vals$val.P[ vals$val.G==1 ])
                col2 = as.numeric(vals$val.P[ vals$val.G==2 ])
            }
        }
        means = c( mean(col1,na.rm=T), mean(col2,na.rm=T) )
        SDs = c( sd(col1,na.rm=T), sd(col2,na.rm=T) )
        n.rep = c( sum(!is.na(col1)), sum(!is.na(col2)) )
        pval = rep(pval, 2)
        results = data.frame( means, SDs, n.rep, pval)
    }
    else if( length(u)==3 ){
        if(u[1]==0 && u[2]==1 && u[3]==2){
            special = T
            gnames = c('Hom','Het','Hom')
            col1 = as.numeric(vals$val.P[ vals$val.G==0 ])
            col2 = as.numeric(vals$val.P[ vals$val.G==1 ])
            col3 = as.numeric(vals$val.P[ vals$val.G==2 ])
        }
        means = c( mean(col1,na.rm=T), mean(col2,na.rm=T), mean(col3,na.rm=T) )
        SDs = c( sd(col1,na.rm=T), sd(col2,na.rm=T),sd(col3,na.rm=T) )
        n.rep = c( sum(!is.na(col1)), sum(!is.na(col2)), sum(!is.na(col3)) )
        pval = rep(pval, 3)
        results = data.frame( means, SDs, n.rep, pval)
    }
    if (special){
        if(length(u)==3){
            stripchart(list(col1, col2, col3), group.names=gnames, vertical=T, pch='.', cex=3, method="jitter", col=c('blue', 'blue', 'blue'), main=main, sub=significance)
            boxplot(list( mean(col1, na.rm=T), mean(col2, na.rm=T), mean(col3, na.rm=T)), names=gnames, boxwex=0.25, add=T)
        }
        else{
            stripchart(list(col1, col2), group.names=gnames, vertical=T, pch='.', cex=3, method="jitter", col=c('blue', 'blue'), main=main, sub=significance)
            boxplot(list( mean(col1, na.rm=T), mean(col2, na.rm=T)), boxwex=0.25, add=T, names=gnames)
        }
    }
    else{
        plot(as.numeric(vals$val.G), as.numeric(vals$val.P), ylab=probe.e, xlab=probe.g, pch=19, main=main)
    }
    signif( results, 3 )
}


plot.geno.expr.bracket = function(ids, expr.g, expr.e, ga, probe.g, probe.e, main=NULL){
    #   Plot a bracket of three physically consecutive genotypes against a phenotype
    ga.o = ga[order(ga$Chr, ga$loc_start),]
    ids.o = rownames(ga.o)
    idx = which(ids.o==probe.g)

    par(mfrow = c(1, 3))
    plot.geno.expr(ids, expr.g, expr.e, ids.o[idx-1], probe.e)
    plot.geno.expr(ids, expr.g, expr.e, ids.o[idx], probe.e)
    plot.geno.expr(ids, expr.g, expr.e, ids.o[idx+1], probe.e)
    print( paste(ids.o[idx-1],ga.o$Chr[idx-1], ga.o$loc_start[idx-1] ))
    print( paste(ids.o[idx],ga.o$Chr[idx], ga.o$loc_start[idx] ))
    print( paste(ids.o[idx+1],ga.o$Chr[idx+1], ga.o$loc_start[idx+1] ))
}


create.meta.ga = function( ga.1, ga.2 ){
    # Given two gene attributes data frames, creates a new ga with meta-rownames.
    # Assumes both ga.1 and ga.2 have column called symbol
    # If ga.1 and ga.2 each have one probe for symbol X, only one meta-X is created.
    # Otherwise, one meta-X is created for each combination of ga.1 and ga.2
    # If ga.1 has 2 probes for X and ga.2 has 3 probes, six meta-X probes are created.
    
    common.symbols = sort( intersect(ga.AB$symbol, ga.one$symbol) )
    rn1 = rownames(ga.1)
    rn2 = rownames(ga.2)
    N = 0
    for(i in 1:length(common.symbols) ){
        symbol = common.symbols[i]
        n.1 = length(which(ga.1$symbol==symbol ) )
        n.2 = length(which(ga.2$symbol==symbol ) )
        N = N + (n.1 * n.2)
    }
    meta.ids = paste( rep('meta', N), 1:N, sep='_')
    id.1.orig = rep('', N)
    id.2.orig = rep('', N)
    symbol = rep('', N)
    n=1
    for(i in 1:length(common.symbols) ){
        idx.1 = which(ga.1$symbol==common.symbols[i] )
        idx.2 = which(ga.2$symbol==common.symbols[i] )
        n.1 = length( idx.1 )
        n.2 = length( idx.2 )
        for(j in 1:n.1){
            for(k in 1:n.2){
                id.1.orig[n] = rn1[ idx.1[j] ]
                id.2.orig[n] = rn2[ idx.2[k] ]
                symbol[n] = common.symbols[i]
                n = n+1
            }
        }
    }
    ga = data.frame( symbol, id.1.orig, id.2.orig, stringsAsFactors=F)
    rownames(ga) = meta.ids
    ga
}


create.meta.expr = function( ga, expr.1, expr.2 ){
    # Given gene attributes created by create.meta.ga, combine the 
    # expression data in expr.1 and expr.2 to generate a matched expression dataset
    NN =  c(names(expr.1), names(expr.2) )
    expr.1 = data.matrix(expr.1)
    expr.2 = data.matrix(expr.2)
    expr = matrix(0, nrow=dim(ga)[1], ncol=( dim(expr.1)[2] + dim(expr.2)[2]) )
    for(i in 1:dim(ga)[1]){
        id.1 = ga$id.1.orig[i]
        id.2 = ga$id.2.orig[i]
        e1 = expr.1[which(rownames(expr.1)==id.1),]
        e2 = expr.2[which(rownames(expr.2)==id.2),]
        expr[i,] = c( e1 , e2 )
    }
    expr = data.frame(expr, stringsAsFactors=F)
    names(expr) = NN
    rownames(expr) = rownames(ga)
    expr
} 


plot.probe.paired.data = function(ids, expr.A, expr.B, probe.id, main=NULL, method="barplot"){
    # barplot or stripchart plot of probe.id for matched subjects in expr.A, expr.B
    if(method != 'barplot' && method != 'stripchart'){ stop("method must be barplot or stripchart")}
    N = dim(ids)[1]
    rowA = which( rownames(expr.A) == probe.id )
    rowB = which( rownames(expr.B) == probe.id )
    if( is.null(main) )
        main = probe.id
    matches = match.indices(ids, expr.A, expr.B)
    idx.A = matches$A
    idx.B = matches$B
    vals.A = as.numeric( expr.A[ rowA, idx.A] )
    vals.B = as.numeric( expr.B[ rowB, idx.B] )
    ab= data.frame(vals.A, vals.B)
    sorted.idx = order(ab[,1])
    vals = vector(mode="numeric", N*2)
    cols = vector(mode="character", N*2)
    ctr=1
    for(i in 1:N){
        vals[ctr]=vals.A[sorted.idx[i]]
        vals[ctr+1]=vals.B[sorted.idx[i]]
        cols[ctr]="black"
        cols[ctr+1]="orange3"
        ctr = ctr + 2
    }
    if(method=="barplot"){
        barplot(vals, col=cols, border=0, main=main, sub=paste( 'Black:Normal, Orange:Tumor'))
    }
    else{
        stripchart(list(vals.A, vals.B), group.names=c('Normal', 'Tumor'), vertical=T, pch='.', cex=3, method="jitter", col=c('blue', 'blue'), main=main)
        boxplot(list( mean(vals.A, na.rm=T), mean(vals.B, na.rm=T)) , boxwex=0.25, add=T, names=c('Normal','Tumor'))

    }
    print(paste("Mean A:",mean(vals.A,na.rm=T) ))
    print(paste("Mean B:",mean(vals.B,na.rm=T) ))
}

plot.chr.CRLMM = function(cS, cCN, s.no, max.y=5, x1=1, x2=0, lblS, lblCN){
    # Separately plots SNP and Copy Number (non-polymorphic) probes for one chromosome
    lblS = as.numeric(lblS)
    lblCN = as.numeric(lblCN)
    par(mfrow = c(2,1))
    sn = environ$sns[s.no]    
    if(x2==0){
        x1CN = x1
        x1S = x1
        x2CN = dim(cCN)[1]; 
        x2S = dim(cS)[1];
    }
    else{ 
        idxS = which( lblS>x1 & lblS < x2 )
        idxCN = which( lblCN>x1 & lblCN < x2 )
        x1S = min( idxS );
        x1CN = min( idxCN );
        x2S = max( idxS );
        x2CN = max( idxCN );
    }
    scaleCN = (x2CN-x1CN)/10
    scaleS = round((x2S-x1S)/10)
    plot(cS[x1S:x2S,s.no], pch='.', ylim=c(-2,max.y), main=sn, axes=F)
    lbl.S = lblS[seq(x1S,x2S,scaleS)]
    print(lbl.S)
    print(seq(x1S,x2S,scaleS) )
    axis(1,labels=lbl.S, at=seq(x1S,x2S,scaleS) )
    axis(2)
    abline(2,0,col=3)
    plot(cCN[x1CN:x2CN,s.no]/100,pch='.',ylim=c(-2,max.y),axes=F)
    axis(1,labels=lblCN[seq(x1CN,x2CN,scaleCN)],at=seq(1,x2CN-x1CN+1,scaleCN))
    axis(2)
    abline(2,0,col=3)
}

plot.cluster = function(M, labels=names(M), d=NULL, main=NULL, direction="samples", cex=0.5){
    #plot hieararchical cluster with labels, use distance d is passed
    if( is.null(d) ){
        if( direction=="samples" ){
            d = (dist(as.matrix(t(M))))
            if( is.null( labels) ){
                labels =  paste( "s", 1:(dim(M)[2]), sep='.') 
            }
        }
        else if( direction=="probes" ){
            d = (dist(as.matrix(M)))
            if( is.null( labels ) )
                labels =  paste( "p", 1:(dim(M)[1]), sep='.') 
        }
        else
            stop("direction must be one of: samples, probes")
    }
    plot(hclust(d), hang=-1, cex=cex, labels=as.character(labels), main=main )
}


color.dendrogram.labels = function(n){
    if(is.leaf(n)){
        a = attributes(n)
        color = foo.lookup.1$color[ which(foo.lookup.1$leaf.name==a$label)[1] ]
        if( foo.show.labels.1 ){
            cex.val = 1
            lab.color = color
        }
        else{
            cex.val=0.01
            lab.color="white"   
        }
        attr(n, "nodePar") = c(a$nodePar, list(lab.col = lab.color, lab.cex=cex.val, col=color, pch=15, cex=1 ) ) 
    }
    n
}

    
plot.cluster.colors=function( M, colors, force.flat=F, show.labels=T ){
    # Still a hack, but less awful
    d = (dist(as.matrix(t(M))))
    if( is.null(names(M)) )
        stop("columns must have unique names")
    lookup = data.frame( leaf.name = names(M), color=colors, stringsAsFactors=F )
    if(force.flat)
        hangval = -1
    else
        hangval = 0.1
    D = as.dendrogram(hclust(d), hang=hangval)
    assign("foo.lookup.1", lookup, envir=.GlobalEnv)
    assign("foo.show.labels.1", show.labels, envir=.GlobalEnv)    
    dend_colored = dendrapply(D, color.dendrogram.labels) 
    plot(dend_colored)
}



plot.clean=function( V, ymin=NA, ymax=NA, colors=c(), cex=0.25, y.axis=T ){
    # This is the utility function that just plots what it's given
    n_probes = dim(V)[1]
    n_samples = dim(V)[2]
    X = (1:n_samples)
    if( is.na(ymax) )
        ymax = max(c(1,ceiling(max(V, na.rm=T))))
    if( is.na(ymin) ){
        ymin = floor( min(V, na.rm=T) )
        if(ymin>1)
            ymin=1
    }
    pch.used = c(19)
    if( length(colors) == 0 )
        colors = c('black', 'blue', 'darkgreen', 'red')
    plot(X, as.numeric(V[1,]), col=colors[1], type='p', pch=19, cex=cex,  axes=F, ylim=c(ymin,ymax), ylab='', lwd=2, mar=c(1,0.5,1,0.5))
    if(n_probes>1){
        for( i in 2:n_probes){
            points(X, as.numeric(V[i,]), pch=19, col=colors[i], lwd=2, cex=cex)
        }
    }
    par(ps=8)
    if(y.axis)
        axis(2, at=(ymin:ymax),las=2, cex.axis=1 )
    segments(0,0,x1=n_samples, col='black')
}


sorted.heatmap=function(D, ga, target.symbols=NULL, target.probes=NULL, sort.by=NULL, scale=F, y.min=NULL, y.max=NULL ){
    library(ggplot2)
    # passing a list of probes uses the probe IDs to pick exact targets.
    # If symbol==NULL, sorted by first symbol.
    probes =  c()
    symbols= c()
    if( is.null(target.probes) & is.null(target.symbols) )
        target.probes = rownames(D)
        #stop( "Must pass either target.probes or target.symbols" )
    if( !is.null(target.probes) & !is.null(target.symbols) )
        stop( "Cannot pass both target.probes and target.symbols" )
    if( !is.null(target.probes) ){
        m = match.idx( target.probes, rownames(ga))
        if( dim(m)[1] != length(target.probes) ){
            stop("Not all target probes found in gene attributes.")
        }
        else{
            probes = target.probes
            symbols = ga$symbol[m$idx.B]
        }
    }
    else{
        for(i in 1:length(target.symbols)){
            probes = c(probes, rownames(ga)[ga$symbol==target.symbols[i]][1])
            symbols = c(symbols, ga$symbol[ga$symbol==target.symbols[i]][1])
        }
    }
    DT = data.matrix( D[match.idx(probes, rownames(D))$idx.B,] )
    
    if(!is.null(sort.by)){
        idx.s = which(symbols==sort.by)[1]
        idx.n = which(names(D)==sort.by)[1]
        print(idx.s)
        print(idx.n)
        
        if(!is.na(idx.s)){
            DT = DT[,order(DT[idx.s,])]
        }
        else if(!is.na(idx.n)){
            print(DT[1:4,1:4])
            DT = DT[order(DT[,idx.n]),]
            print(DT[1:4,1:4])            
        }
        else{
            stop("Cannot find sort.by parameter in names or rownames of matrix")
        }
    }
    DT = DT[dim(DT)[1]:1,]

    if( scale ){
        means = rowMeans(DT)
        sds = rep(0, dim(DT)[1])
        for(i in 1:length(sds)){
            sds[i] = sd(DT[i,])
        }
        DT = (DT - means) / sds
    }
    if( !is.null(y.min) ){
        if( sum(DT<y.min, na.rm=T)>0 ){
            print("Truncating at lower bound")
            DT[DT<y.min] = y.min
        }
    }
    if( !is.null(y.max) ){
        if( sum(DT>y.max, na.rm=T)>0 ){
            print("Truncating at upper bound")
            DT[DT>y.max] = y.max
        }
    }    
    if(is.null(y.min) | is.null(y.max)){
        y.min=min(DT, na.rm=T)
        y.max=max(DT, na.rm=T)
    }    
    df = expand.grid(y = 1:dim(DT)[1], x = 1:dim(DT)[2] )
    df = cbind(df, v=as.numeric(DT) )
    y.mid = mean(c(y.max, y.min))
    scg = scale_fill_gradient2(low = "darkblue", high = "red2", midpoint=y.mid, limits=c(y.min,y.max))
    ggplot(df, aes(x, y, fill = v)) + geom_tile() + scg + theme_bw() 
}

do.pca=function( D, pca=NULL, labels=NULL, xlim=NULL, ylim=NULL, show.legend=T, colors=NULL, legend.xy=NULL){
    if(is.null(colors))
        color.wheel = c("black", "blue", "gold", "darkgreen", "slategray1", "gray", "magenta", "darkblue",
        "violetred", "bisque", "chartreuse3", "orange", "darksalmon", "green1", "red","pink")
    else{
        color.wheel = colors
    }
    if(is.null(pca)){
        pca = prcomp(as.matrix(t(D)))
    }
    colors = rep("black", dim(D)[2])
    if( !is.null(labels) ){
        labels[is.na(labels)]="NA"
        unique.labels = sort(unique(labels))
        legend.cex=1.5
        for( i in 1:length(labels)){
            colors[i] = color.wheel[ which(unique.labels==labels[i]) ] 
        }
        
        if(length(unique.labels)>8){
            legend.cex=0.75
        }
    }
    y.min = min( pca$x[,2] )
    y.max = max( pca$x[,2] )    
    x.min = min( pca$x[,1] )
    x.max = max( pca$x[,1] )
    
    if( dim(D)[2]>100 ){
        plot.cex = 0.75
    }
    else{
        plot.cex=1
    }
    if(is.null(xlim))
        xlim = c(x.min, x.max)
    if(is.null(ylim))
        ylim = c(y.min, y.max)
    
    plot(pca$x[,1], pca$x[,2], col=colors, pch=19, cex=plot.cex, xlim=xlim, ylim=ylim )  
    
    if( show.legend & !is.null(labels) ){
        if( is.null(legend.xy) )
            legend.xy = c(xlim[1], ylim[2])
        legend(legend.xy[1], legend.xy[2], unique.labels, col=color.wheel[1:length(unique.labels)], 
               pch=19,cex=legend.cex,box.col="white" )
    }
    pca
}

####################
# BEGIN STATISTICS #
####################

cor.test.row.vs.genes = function(ids, expr.A, expr.B, ga.genes, probe.id, percent.present=0.9, method="pearson", verbose=T){
    # Convenience method for correlation testing of one row vs. all rows
    # ids is data frame where columns are valid paired identifiers {expr.row, expr.matrix}
    # expr.A is data frame
    # expr.B is data frame
    # ga.genes is a data frame with gene attributes for expr.B
    # probe is a rowname from expr.A
    # percent.present is minimum percent of measurements present required in expr.matrix
    # method is spearman or pearson
    # if verbose==T, report progress
    #
    # RETURNS: data frame of index in expr.matrix, r-value, p-value
    row.idx = which(rownames(expr.A)==probe.id)
    if( length(row.idx)==0 ){ stop(paste("probe",probe.id,"not found in expr.A")) }
    if( dim(ga.genes)[1] != dim(expr.B)[1] ){stop("ga.genes has different number of rows from expr.B")}

    expr.row = expr.A[ row.idx,]
    d=cor.test.row.vs.all( ids, expr.row, expr.B, percent.present=0.90, method=method, verbose)
    pvals = d$pvals
    rhos = d$rhos
    d = cbind(ga.genes[d$genes,], pvals, rhos)
    d
}

cor.test.all.rows.vs.genes = function(ids, expr.A, expr.B, ga.genes, percent.present=0.9, method="pearson", max.p=1E-05, verbose=T){
    # Convenience method for correlation testing of one row vs. all rows
    # ids is data frame where columns are valid paired identifiers {expr.row, expr.matrix}
    # expr.A is data frame
    # expr.B is data frame
    # ga.genes is a data frame with gene attributes for expr.B
    # percent.present is minimum percent of measurements present required in expr.matrix
    # method is spearman or pearson
    # if verbose==T, report progress
    #
    # RETURNS: data frame of index in expr.matrix, r-value, p-value
    probe.ids.A = rownames(expr.A)
    d = data.frame()
    for(i in 1:length(probe.ids.A) ){
        if(verbose){
            print( paste("Row ", probe.ids.A[i],"(",i,"of ",length(probe.ids.A),")" ))
        }
        d.i = cor.test.row.vs.genes(ids, expr.A, expr.B, ga.genes, probe.ids.A[i], percent.present=percent.present, method=method, verbose)
        d.i = d.i[d.i$pvals<=max.p & !is.na(d.i$pvals), ]
        if( dim(d.i)[1] > 0 ){
            d.i = cbind(rep(probe.ids.A[i], dim(d.i)[1]), rownames(d.i), d.i)
            names(d.i)[1] = "probe.id.A"
            names(d.i)[2] = "probe.id.B"
            rownames( d.i ) = prepend.label( rownames(d.i), probe.ids.A[i])
            if( i==1 )
                d = d.i
            else
                d = rbind(d, d.i)
        }
    }
    d
}

cor.test.probe.vs.all=function( expr, symbols, probe, percent.present=0.90, method="spearman", verbose=F){
    # Calculate correlation 
    # percent.present is minimum percent of measurements present required in expr.matrix
    # method is spearman or pearson
    # if verbose==T, report progress
    #
    # RETURNS: data frame of index in expr.matrix, r-value, p-value
    if( percent.present > 1 || percent.present < 0 ){ stop("percent.present not between 0 and 1") }
    
    D = data.matrix(expr)
    vals.row = D[ which(rownames(expr)==probe),]
    N.genes = dim(D)[1]
    N = dim(D)[2]
    threshold = floor( N * percent.present )
    n.NA = rowSums(is.na(D))
    if( N - sum( is.na(vals.row) ) < threshold ){
        warning("Skipping row due to insufficient number of present measurements")
        d=data.frame( 1:N.genes, rep(NA,N.genes), rep(NA,N.genes))
        names(d) = c('genes', 'rhos', 'pvals')
        d
    }
    else{
        options(warn=-1)
        genes = vector(mode="integer", N.genes)
        pvals = vector(mode="numeric", N.genes)
        rhos = vector(mode="numeric", N.genes)
        for(i in 1:N.genes){
            if( N - n.NA[i] >= threshold ){
                result = cor.test(vals.row, D[i,], method=method, rm.na=T)
                pval = result$p.value
                genes[i] = i
                pvals[i] = pval
                rhos[i] = as.numeric(result$estimate)
            }
            else{
                genes[i] = i
                pvals[i] = NA
                rhos[i] = NA
            }
            if( verbose && i %% 5000==0 ){
                print(paste("Gene", i, "of", N.genes) )
            }
        }
        options(warn=0)
        d = data.frame(symbols, rho=round(rhos,3), pval=signif(pvals,4))
        rownames(d) = rownames(expr)
        d = d[ order(d$pval, d$rho), ]
        d
    }
}

cor.test.row.vs.all = function( ids, expr.row, expr.all, percent.present=0.90, method="spearman", verbose=F){
    # Calculate spearman value of first row of expr.row vs. all rows in expr.all
    #
    # ids is data frame where columns are valid paired identifiers {expr.row, expr.matrix}
    # expr.row is data frame; if more than one row, only first is used.
    # expr.all is data frame
    # percent.present is minimum percent of measurements present required in expr.matrix
    # method is spearman or pearson
    # if verbose==T, report progress
    #
    # RETURNS: data frame of index in expr.matrix, r-value, p-value
    if( !is.data.frame(expr.row) ){ stop( "expr.row is not a data frame") }
    if( !is.data.frame(expr.all) ){ stop( "expr.matrix is not a data frame") }
    if( dim(expr.row)[1] != 1 ){ warning( "expr.row has more than one row; using first row") }
    if( percent.present > 1 || percent.present < 0 ){ stop("percent.present not between 0 and 1") }
    
    matches = match.indices(ids, expr.row, expr.all)
    idx.row = matches$A
    idx.all = matches$B
    vals.row = as.numeric( expr.row[ 1, idx.row] )
    N.genes = dim(expr.all)[1]
    N = dim(ids)[1]
    threshold = floor( N * percent.present )
    if( N - sum( is.na(vals.row) ) < threshold ){
        warning("Skipping row due to insufficient number of present measurements")
        d=data.frame( 1:N.genes, rep(NA,N.genes), rep(NA,N.genes))
        names(d) = c('genes', 'rhos', 'pvals')
        d
    }
    else{
        options(warn=-1)
        genes = vector(mode="integer", N.genes)
        pvals = vector(mode="numeric", N.genes)
        rhos = vector(mode="numeric", N.genes)
        for(i in 1:N.genes){
            vals.gene = as.numeric(expr.all[i,idx.all])
            if( N - sum( is.na(vals.gene) ) >= threshold ){
                result = cor.test(vals.row, vals.gene, method=method, rm.na=T)
                pval = result$p.value
                genes[i] = i
                pvals[i] = pval
                rhos[i] = as.numeric(result$estimate)
            }
            else{
                genes[i] = i
                pvals[i] = NA
                rhos[i] = NA
            }
            if( verbose && i %% 5000==0 ){
                print(paste("Gene", i, "of", N.genes) )
            }
        }
        options(warn=0)
        d = data.frame(genes, rhos, pvals)
        d = d[ order(d$pvals, d$rhos), ]
        d
    }
}


cor.sa.columns=function( sa, method="spearman" ){
    N = names(sa)
    idx = vector()
    for(i in 1:length(N) ){
        if( is.numeric( sa[,i] ) )
            idx[length(idx)+1] = i
    }
    pvals = vector()
    stats = vector()
    is = vector()
    js=vector()
    n.i = vector(); n.j = vector()
    for(i in 1:(length(idx)-1) ){
        vals.i = sa[,idx[i]]
        for( j in (i+1) : length(idx) ){
            if( i!=j ){
                vals.j = sa[,idx[j]]
                if( sum( !(is.na(vals.i)) & !(is.na(vals.j)) )>2 ){
                    r = cor.test(vals.i, vals.j, method=method, na.rm=T )
                    if( !is.na(r$p.value) ){
                        pvals[length(pvals) + 1] = r$p.value
                        stats[length(stats) + 1] = r$estimate
                        n.i[ length(n.i)+1 ] = N[idx[i]]
                        n.j[ length(n.j)+1 ] = N[idx[j]] 
                        is[ length(is)+1 ] = idx[i]
                        js[ length(js)+1 ] = idx[j]
                    }
                }
            }
        }
    }
    df = data.frame(pvals, stats, n.i, n.j, is, js)
    df[order(df$pvals),]
}


identify.eQTL.networks = function(eqtl, spear, probes.snp, ga.snp.chr, ga.snp.pos, 
                                  probes.gene, ga.gene.chr, ga.gene.loc, min.trans=1){
    # find SNPs that appear more than once
    freq = count.appearances(eqtl$snp)
    m = match.idx(freq$keys, probes.snp)
    freq = cbind(freq, chr=ga.snp.chr[m$idx.B], pos=ga.snp.pos[m$idx.B], stringsAsFactors=F)
    freq = freq[freq$values>1,]
    is.first=T
    net.ctr=1
    for(i in 1:dim(freq)[1]){
        
        # make ee, eQTL network for snp i. Call "CIS" anything within 50 Mb.
        snp = freq$keys[i]
        snp.chr = ga.snp.chr[which(probes.snp==snp)]
        snp.loc = ga.snp.pos[which(probes.snp==snp)]
        ee = eqtl[eqtl$snp==snp,]
        m = match.idx( ee$probe, probes.gene )
        ee = cbind(ee, chrom=ga.gene.chr[m$idx.B], start=ga.gene.loc[m$idx.B],stringsAsFactors=F)
        cis = ee$chrom==snp.chr & abs( ee$start-snp.loc ) < 50000000
        ee = cbind(ee, cis, stringsAsFactors=F)
        idx.cis = which(cis)
        n.matches = dim(ee)[1]
        n.cis = length(idx.cis)
        print(paste(i,"of", dim(freq)[1],"SNP", snp, "matches",n.matches,"eQTL",n.cis,"are cis"))
        # Only continue if there is at least one cis-eQTL
        if(n.cis>0 & n.matches-n.cis){
            for(j in 1:length(idx.cis)){
                # For each cis-eQTL, find probes correlated with it.
                # If any of those are in this eQTL network, write out the network
                probe = ee$probe[ idx.cis[j] ]
                ss = spear[ (spear$probe.1==probe | spear$probe.2==probe) & spear$symbol.1 != spear$symbol.2,]
                neighbors = setdiff( c( ss$probe.1 , ss$probe.2), probe )
                with.eqtl = intersect(neighbors, ee$probe )
                
                if(length(with.eqtl)>0){
                    m = match.idx( c( probe, with.eqtl ), ee$probe )
                    found = ee[m$idx.B,]
                    if( sum(!found$cis)>min.trans ){
                        network = rep(net.ctr, dim(found)[1])
                        net.ctr = net.ctr+1
                        found = cbind(found, network)
                        if(is.first){
                            results=found
                            is.first=F
                        }
                        else{
                            results = rbind(results, found)
                            print(paste("Adding",dim(found)[1],"candidates"))
                        }
                    }
                }
            }
        }
    }
    results
}

write.eqtl.to.cytoscape = function( fn.base, EQTL, SPEAR, snps, probes, all.S, S.chr, S.loc, all.G, G.chr, G.loc){

    CYTOSCAPE = '/Applications/Cytoscape_v2.8.1/cytoscape.sh'
    VIZ = '/notebook/code/release/Correlation.props'
    fn.sif = paste( fn.base, '.sif', sep='')
    fn.type = paste( fn.base, '.noa', sep='')
    fn.symbol = paste( fn.base, '_symbol.noa', sep='')
    fn.chr = paste( fn.base, '_chr.noa', sep='')
    fn.loc = paste( fn.base, '_loc.noa', sep='')
    fn.pval = paste( fn.base, '_pval.eda', sep='')
    fn.rho = paste( fn.base, '_rho.eda', sep='')
    fn.cis = paste( fn.base, '_cis.eda', sep='' )
    fn.edgetype = paste( fn.base, '_type.eda', sep='' )
    
    noa.files = c(fn.type, fn.chr, fn.loc, fn.symbol,fn.cis)
    eda.files = c(fn.pval, fn.rho, fn.edgetype)
    
    print(paste("Passed", length(snps), "SNPs and ", length(probes), "probes"))
    EQTL = EQTL[EQTL$snp %in% snps,]    
    SPEAR = SPEAR[SPEAR$probe.1 %in% probes & SPEAR$probe.2 %in% probes,]
    probes.from.eqtl = EQTL$probe
    probes.from.spear = sort(unique(c(SPEAR$probe.1, SPEAR$probe.2)))
    
    gene.probes = intersect( probes.from.spear, probes.from.eqtl )
    SPEAR = SPEAR[SPEAR$probe.1 %in% gene.probes & SPEAR$probe.2 %in% gene.probes,]
    EQTL = EQTL[EQTL$probe %in% gene.probes,]
    snp.probes = sort(unique(EQTL$snp))
    
    print("Restricting to probes with an eQTL and a correlation" )
    print(paste("Using ", length(snps), "SNPs and", length(gene.probes), "probes" ) )
        
    # Write SIF and rhos
    unlink(fn.sif)
    unlink(fn.rho)
    unlink(fn.edgetype)
    write( "rho (java.lang.Double)" , file=fn.rho )
    write( "FOO = 1.0", file=fn.rho, append=T ) # get around cytoscape bug   
    write( "edgetype (java.lang.String)" , file=fn.edgetype )
    for(i in 1:dim(SPEAR)[1]){
        write(paste( SPEAR$probe.1[i], 'gg', SPEAR$probe.2[i], sep=' '), file=fn.sif, append=T)
        write(paste( SPEAR$probe.1[i], '(gg)', SPEAR$probe.2[i], "=", round(SPEAR$rho.a[i],3), sep=' '), file=fn.rho, append=T )
        write(paste( SPEAR$probe.1[i], '(gg)', SPEAR$probe.2[i], "= correlation", sep=' '), file=fn.edgetype, append=T )
    }
    
    # Write chromosomes, locations
    gene.locs = G.loc[ match.idx( gene.probes, all.G )$idx.B ]
    gene.chrs = G.chr[ match.idx( gene.probes, all.G )$idx.B ]
    snp.locs = S.loc[ match.idx( snp.probes, all.S )$idx.B ]
    snp.chrs = S.chr[ match.idx( snp.probes, all.S )$idx.B ]

    # Write chromosomes and locations
    write( "chr (java.lang.String)" , file=fn.chr )
    write( "loc (java.lang.String)" , file=fn.loc )    
    write( "FOO = X", file=fn.chr, append=T ) # get around cytoscape bug   
    write( "FOO = X", file=fn.loc, append=T ) # get around cytoscape bug   
    for(i in 1:length(gene.probes)){
        write(paste( gene.probes[i], '= ', gene.locs[i], sep=' '), file=fn.loc, append=T)
        write(paste( gene.probes[i], '= ', gene.chrs[i], sep=' '), file=fn.chr, append=T)
    }
    for(i in 1:length(snp.probes)){
        write(paste( snp.probes[i], '= ', snp.locs[i], sep=' '), file=fn.loc, append=T)
        write(paste( snp.probes[i], '= ', snp.chrs[i], sep=' '), file=fn.chr, append=T)
    }
            
    # Write types
    write( "type (java.lang.String)" , file=fn.type )
    write( "FOO = gene", file=fn.type, append=T ) # get around cytoscape bug   
    for(i in 1:length(gene.probes))
        write(paste( gene.probes[i], '= gene', sep=' '), file=fn.type, append=T)
    for(i in 1:length(snp.probes))
        write(paste( snp.probes[i], '= locus', sep=' '), file=fn.type, append=T)
    
    # Write pvalues, symbols
    write( "symbol (java.lang.Double)" , file=fn.symbol )
    write( "pval (java.lang.Double)" , file=fn.pval )
    write( "cis (java.lang.String)" , file=fn.cis )
    write( "FOO = gene", file=fn.symbol, append=T ) # get around cytoscape bug
    write( "FOO (gs) BAR = 1.0", file=fn.pval, append=T ) # get around cytoscape bug
    write( "FOO = CIS", file=fn.cis, append=T ) # get around cytoscape bug
    for(i in 1:dim(EQTL)[1]){
        write(paste( EQTL$probe[i], 'gs', EQTL$snp[i], sep=' '), file=fn.sif, append=T)
        write(paste( EQTL$probe[i], '=', EQTL$symbol[i], sep=' '), file=fn.symbol, append=T )
        write(paste( EQTL$snp[i], '=', EQTL$snp[i], sep=' '), file=fn.symbol, append=T )        
        write(paste( EQTL$probe[i], '(gs)', EQTL$snp[i], "=", EQTL$perm.p[i], sep=' '), file=fn.pval, append=T )
        idx.probe = which(gene.probes==EQTL$probe[i])
        idx.snp = which(snp.probes==EQTL$snp[i])
        write(paste( EQTL$probe[i], '(gs)', EQTL$snp[i], "= eqtl", sep=' '), file=fn.edgetype, append=T )        
        if( is.na( gene.chrs[idx.probe]) | is.na(gene.locs[idx.probe]) | is.na(snp.locs[idx.snp]) )
            write(paste( EQTL$probe[i], '= UNKNOWN', sep=' '), file=fn.cis, append=T )
        else{
             if( gene.chrs[idx.probe]==snp.chrs[idx.snp] & abs(gene.locs[idx.probe]-snp.locs[idx.snp])<70000000)
                write(paste( EQTL$probe[i], '= CIS', sep=' '), file=fn.cis, append=T )
            else
                write(paste( EQTL$probe[i], '= TRANS', sep=' '), file=fn.cis, append=T )
        }
    }
    
    # Write symbols
    for(i in 1:dim(EQTL)[1]){
        write(paste( EQTL$probe[i], '(gs)', EQTL$snp[i], "=", EQTL$perm.p[i], sep=' '), file=fn.pval, append=T )
    }
        
    sh = paste(CYTOSCAPE, "-N", fn.sif, "-V", VIZ, sep=" ")
    for(i in 1:length(noa.files))
        sh = paste(sh, "-n", noa.files[i], sep=" ")
    for(i in 1:length(eda.files))
        sh = paste(sh, "-e", eda.files[i], sep=" ")
    write( sh, file=paste( fn.base, '.sh', sep='') )

    print(paste("Wrote to", paste( fn.base, '.sh', sep='') ) )
}

write.spear.to.cytoscape=function( fn.base, DF, DF.node.attr=NULL ){
    CYTOSCAPE = '/Applications/Cytoscape_v2.8.1/cytoscape.sh'
    VIZ = '/notebook/code/release/Correlation.props'
    fn.sif = paste( fn.base, '.sif', sep='')
    edge.attr = c()
    fn.edge.attr = c()
    numeric.attr = c()
    numeric.node.attr = c()
    sh = paste(CYTOSCAPE, "-N", fn.sif, "-V", VIZ, sep=" ")

    if( length(names(DF))>4 ){
        edge.attr = names(DF)[5:length(names(DF))]
        fn.edge.attr = paste(fn.base, '_', edge.attr, '.eda', sep='')
        for(j in 1:length(fn.edge.attr) ){
            numeric.attr[j] = is.numeric(DF[,4+j])
            sh = paste(sh, "-e", fn.edge.attr[j], sep=" ")
        }
    }
    if( !is.null(DF.node.attr) ){
        node.attr = names(DF.node.attr)
        fn.node.attr = paste(fn.base, '_', node.attr, '.noa', sep='')
        for(j in 1:length(fn.node.attr) ){
            numeric.node.attr[j] = is.numeric(DF.node.attr[,j])
            sh = paste(sh, "-n", fn.node.attr[j], sep=" ")
            if( numeric.node.attr[j] ){
                output = paste( node.attr[j], "(java.lang.Double)"  )
                output = c(output, "FOO = 1.0" )# get around cytoscape bug
            }
            else{
                output = paste( node.attr[j], "(java.lang.String)"  )
                output = c(output, "FOO = yay" )# get around cytoscape bug
            }
            print(output)
            vector.to.file(output, fn.node.attr[j])
            write.table(data.frame(rownames(DF.node.attr), DF.node.attr[, j]), 
                        file=fn.node.attr[j], sep=" = ", 
                        quote=F, append=T, col.names=F, row.names=F)
        }
    }
    write( sh, file=paste( fn.base, '.sh', sep='') )
    unlink(fn.sif)
    N = dim(DF)[1]
    write.table(data.frame(  DF$probe.1, DF$probe.2), 
                        file=fn.sif, sep=" gg ", 
                        quote=F, col.names=F, row.names=F)
    if( length(fn.edge.attr)>0){
        for(j in 1:length(fn.edge.attr) ){
            output = paste( edge.attr[j], "(java.lang.Double)"  )
            output = c(output, "FOO = 1.0" )# get around cytoscape bug
            vector.to.file( output, fn.edge.attr[j] )
            
            col.gg = rep('(gg)', dim(DF)[1])
            col.eq = rep('=', dim(DF)[1])
            
            if( numeric.attr[j] ){
                write.table(data.frame(  DF$probe.1, col.gg, DF$probe.2, col.eq, signif(DF[,j+4], 5) ), 
                        file=fn.edge.attr[j], sep=" ", append=T, 
                        quote=F, col.names=F, row.names=F)
            }
            else{
                write.table(data.frame(  DF$probe.1, col.gg, DF$probe.2, col.eq, DF[,j+4] ), 
                        file=fn.edge.attr[j], sep=" ", append=T, 
                        quote=F, col.names=F, row.names=F)            
            }
        }
    }
    print(paste("Wrote to", paste( fn.base, '.sh', sep='') ) )
}



calculate.min.DC.score = function( n.probes ){
    # What is the Z score required for a P value lower than the bonferroni-corrected
    # number of tests for differential correlation
    bonf = 0.05 / (n.probes*n.probes/2)
    pval=1
    sc = -5.0
    while( pval > bonf ){
        pval=pnorm(sc)
        sc = sc - 0.05
    }
    sc = sc + 0.05
    sc*-1
}


spear = function( fn.spear, fn_expr, fn_ga, fn_sa, min_cor, fn_out, class.a="", class.b="", probe="", y="symbol", score="", min_var="", neighbors=F ){
    # simple wrapper for spear 
    c1 = fn.spear
    c2 = paste( "-d", fn_expr, sep='' )
    c3 = paste( "-f", fn_sa, sep='' )
    c4 = paste( "-g", fn_ga, sep='' )
    c5 = paste( "-y", y, sep='' )
    c6 = "-vT"
    c7 = paste( "-o", fn_out, sep='' )
    c8 = paste( "-s", min_cor, sep='' )
    
    cmd = paste(c1, c2, c3, c4, c5, c6, c7, c8)
    if( class.a != "" ){
        aa = paste( "-a", class.a, sep='' )
        cmd = paste( cmd, aa )
    }
    if( class.b != "" ){
        bb = paste( "-b", class.b, sep='' )
        cmd = paste( cmd, bb )
    }
    if( probe != "" ){
        pp = paste( "-p", probe, sep='' )
        cmd = paste( cmd, pp )
    }
    if( score != "" ){
        zz = paste( "-x", score, sep="" )
        cmd = paste(cmd, zz)
    }
    if( min_var != "" ){
        zz = paste( "-m", min_var, sep="" )
        cmd = paste(cmd, zz)
    }
    if( neighbors ){
        cmd = paste(cmd, "-rT" )
    }
    print( cmd )
    system(cmd)
}


overlap.spear.pairs = function(A, B, p2s ){
    # Assumes that the second condition for A and the first condition for B are identical
    print("Assumes that the second condition for A and the first condition for B are identical")
    print(paste("Original size of A:", dim(A)[1]))
    print(paste("Original size of B:", dim(B)[1]))    
    gA = hashgraph()
    gB = hashgraph()
    addEdge(gA, A$probe.1, A$probe.2)
    addEdge(gB, B$probe.1, B$probe.2)
    AB = intersection(gA, gB)
    ee=edges( AB )
    symbol.1 = hsh_get( p2s, ee$e1 )
    symbol.2 = hsh_get( p2s, ee$e2 )
    probe.1 = ee$e1
    probe.2 = ee$e2
    ee = data.frame(symbol.1, probe.1, symbol.2, probe.2, stringsAsFactors=F )
    keep.A = rep(F, dim(A)[1])
    keep.B = rep(F, dim(B)[1])
    for(i in 1:length(keep.A)){
        if( hasEdge(AB, A$probe.1[i], A$probe.2[i]) ){
            keep.A[i] = T
        }
    }
    for(i in 1:length(keep.B)){
        if( hasEdge(AB, B$probe.1[i], B$probe.2[i]) ){
            keep.B[i] = T
        }
    }
    print(paste( "keeping", sum(keep.A), "probes") )
    A.keep = A[keep.A,]
    B.keep = B[keep.B,]
    A.keep = A.keep[order(A.keep$probe.1, A.keep$probe.2),]
    B.keep = B.keep[order(B.keep$probe.1, B.keep$probe.2),]
    if( sum(A.keep$probe.1!=B.keep$probe.1)!=0 )
        stop("Error synchronizing A and B while merging")
    R = cbind( A.keep[,c(1,2,3,4,5,6)], B.keep[,6], A.keep[,9], B.keep[,9] )
    names(R)[5:9] = c('rho.A', 'rho.B', 'rho.C', 'tstat.AB', 'tstat.BC' )
    R$rho.A = round( R$rho.A, 2)
    R$rho.B = round( R$rho.B, 2)
    R$rho.C = round( R$rho.C, 2)
    R$tstat.AB = round( R$tstat.AB, 2)
    R$tstat.BC = round( R$tstat.BC, 2)
    R = R[order(R$rho.B, decreasing=T),]
    R
}

paired.t.test = function(ids, expr.A, expr.B, method='ttest', percent.present=0.9, verbose=T){
    # Calculate pearson or wilcox statistic for paired data sets
    if( percent.present > 1 || percent.present < 0 ){ stop("percent.present not between 0 and 1") }
    N = dim(ids)[1]
    N.probes = dim(expr.A)[1]
    threshold = floor( N * percent.present )
    probes.A = rownames(expr.A)
    probes.B = rownames(expr.B)
    if( N.probes != length(probes.B) ){ stop("Number of probes in expr.A not equal to number of probes in expr.B") }
    if(method != 'ttest' && method != 'wilcox'){ stop("method must be ttest or wilcox") }
    pvals = vector(mode="numeric", N.probes)
    stats = vector(mode="numeric", N.probes)
    mu.A= vector(mode="numeric", N.probes)
    mu.B = vector(mode="numeric", N.probes)
    probes = vector(mode="character", N.probes)
    matches = match.indices(ids, expr.A, expr.B)
    idx.A = matches$A
    idx.B = matches$B

    oA = order(rownames(expr.A))
    oB = order(rownames(expr.B))

    for( i in (1:N.probes) ){
        if( probes.A[oA[i]] != probes.B[oB[i]] ){
            stop( paste("probe identifier in sorted row",i,"is not the same between expr.A and expr.B"))
        }
        vals.A = as.numeric( expr.A[ oA[i], idx.A] )
        vals.B = as.numeric( expr.B[ oB[i], idx.B] )
        probes[i] = probes.A[oA[i]]
        n.present.A = N - sum( is.na(vals.A) )
        n.present.B = N - sum( is.na(vals.B) )
        p = NA
        s = NA
        ma = NA
        mb = NA
        if( n.present.A >= threshold && n.present.B >= threshold ){
            ma = round(mean(vals.A, na.rm=T),2)
            mb = round(mean(vals.B, na.rm=T),2)
            if(method=='ttest'){
                t.res = try( t.test(vals.A, vals.B, paired=T), silent=T )
                if( !is.null(names(t.res))) {
                    p = signif(t.res$p.value,3)
                    s = t.res$statistic
                }
            }
            else{
                w.res = try( wilcox.test(vals.A, vals.B, paired=T), silent=T )
                if( !is.null(names(w.res))) {
                    p = signif(w.res$p.value,3)
                    s = as.numeric(w.res$statistic)
                }
            }
        }
        pvals[i]=p
        stats[i]=s
        mu.A[i] = ma
        mu.B[i] = mb
        if( verbose && i %% 500==0 ){
            print(paste("Probe", i, "of", N.probes) )
        }
    }
    data.frame( IDENTIFIER = probes, pval=pvals, t.stat=stats, mean.A=mu.A, mean.B = mu.B)
}


unpaired.t.test = function(expr.A, expr.B, method='ttest', percent.present=0.9, verbose=T){
    # Calculate pearson or wilcox statistic for paired data sets
    if( percent.present > 1 || percent.present < 0 ){ stop("percent.present not between 0 and 1") }
    N.probes = dim(expr.A)[1]
    N.A = dim(expr.A)[2]
    N.B = dim(expr.B)[2]
    threshold.A = floor( N.A * percent.present )
    threshold.B = floor( N.B * percent.present )
    probes.A = rownames(expr.A)
    probes.B = rownames(expr.B)
    if(method != 'ttest' && method != 'wilcox'){ stop("method must be ttest or wilcox") }
    pvals = vector(mode="numeric", N.probes)
    stats = vector(mode="numeric", N.probes)
    mu.A= vector(mode="numeric", N.probes)
    mu.B = vector(mode="numeric", N.probes)
    probes = vector(mode="character", N.probes)

    oA = order(rownames(expr.A))
    oB = order(rownames(expr.B))

    for( i in (1:N.probes) ){
        if( probes.A[oA[i]] != probes.B[oB[i]] ){
            stop( paste("probe identifier in sorted row",i,"is not the same between expr.A and expr.B"))
        }
        vals.A = as.numeric( expr.A[ oA[i], ] )
        vals.B = as.numeric( expr.B[ oB[i], ] )
        probes[i] = probes.A[oA[i]]
        n.present.A = N.A - sum( is.na(vals.A) )
        n.present.B = N.B - sum( is.na(vals.B) )
        p = NA
        s = NA
        ma = NA
        mb = NA
        if( n.present.A >= threshold.A && n.present.B >= threshold.B ){
            ma = round(mean(vals.A, na.rm=T),2)
            mb = round(mean(vals.B, na.rm=T),2)
            if(method=='ttest'){
                t.res = try( t.test(vals.A, vals.B), silent=T )
                if( !is.null(names(t.res))) {
                    p = signif(t.res$p.value,3)
                    s = signif(t.res$statistic,5)
                }
            }
            else{
                w.res = try( wilcox.test(vals.A, vals.B), silent=T )
                if( !is.null(names(w.res))) {
                    p = signif(w.res$p.value,3)
                    s = signif(as.numeric(w.res$statistic),5)
                }
            }
        }
        pvals[i]=p
        stats[i]=s
        mu.A[i] = ma
        mu.B[i] = mb
        if( verbose && i %% 500==0 ){
            print(paste("Probe", i, "of", N.probes) )
        }
    }
    data.frame( IDENTIFIER = probes, pval=pvals, t.stat=stats, mean.A=mu.A, mean.B = mu.B)
}


logistic = function(sa, expr, do.plot=F){
    # perform a logistic regression of identifiers in sa with values in expr
    # will use the first column of sa and first row of expr for values, matching sa rownames to expr names
    N = dim(sa)[1]
    ids.sa = rownames(sa)
    ids.expr = names(expr)
    tf = contains(ids.sa, ids.expr)
    if( sum(tf) != N ){ stop("not every identifier in rownames(sa) has a match in names(expr)") }
    idx.expr = vector(mode="integer", N)
    for( i in 1:N ){ idx.expr[i] = which( names(expr)==ids.sa[i] ) }
    real.values = as.numeric( expr[1,idx.expr])
    logistic.values = as.numeric(sa[,1])
    u = sort(unique(logistic.values))
    if( length(u) != 2 || u[1] != 0 || u[2] != 1 ){ stop("logistic values (first row of sa) not exclusivly 0 or 1") }
    d=data.frame(x=real.values, y=logistic.values)
    if( do.plot ){
        plot(real.values, logistic.values)
    }
    glm("y ~ x", binomial(), d)
}

fisher = function(A, B, labels.A=NULL, labels.B=NULL, verbose=T){
    # convert two vectors of {0,1} to four categories and perform fisher exact test
    if( length(A) != length(B) ){ stop("Length of A not equal to Length of B") }
    v = matrix(c(0,0,0,0),nrow=2,ncol=2)
    miss = 0
    for(i in 1:length(A) ){
        if( is.na( A[i] ) | is.na(B[i]) ){
            miss = miss + 1
        }
        else{
            if( A[i]==0 & B[i]==0 ){ v[1,1] = v[1,1]+1 }
            else if( A[i]==0 & B[i]==1 ){ v[1,2] = v[1,2]+1  }
            else if( A[i]==1 & B[i]==0 ){ v[2,1] = v[2,1]+1  }
            else if( A[i]==1 & B[i]==1 ){ v[2,2] = v[2,2]+1  }
        }
    }
    if(verbose){
        print(paste("Excluded for missing data:",miss,"of",length(A)))
    }
    d=data.frame(v)
    if( is.null(labels.A) ){
        rownames(d) = c('zero', 'one')
    }
    else{
        if( length(labels.A) != 2 )
            stop("Must provide exactly two labels for A or pass NULL")
        rownames(d) = labels.A
    }
    if( is.null(labels.B) ){
        names(d) = c('zero', 'one')
    }
    else{
        if( length(labels.B) != 2 )
            stop("Must provide exactly two labels for B or pass NULL")
        names(d) = labels.B
    }    
    if(verbose)
        print(d)
    fisher.test(v)
}

#############
# QTL TOOLS #
#############
calculate.linear.model.se = function(linear.model){
    # Given a solved linear model, return the standard error of the coefficients
    # Called by calculate.lm.pval
    # Adapted from /R-2.11.0/src/library/stats/R/lm.R
    Qr = linear.model$qr
    p = linear.model$rank
    p1 = 1L:p
    R = chol2inv(Qr$qr[p1, p1, drop = FALSE])
    r = linear.model$residuals
    n.vars = dim(r)[2]
    rss = as.numeric(colSums(r^2), na.rm=T)
    n = NROW(Qr$qr)
    rdf = n - p
    resvar = rss/rdf
    se = matrix(nrow=p, ncol=n.vars)
    for(i in 1:n.vars){
        se[,i] = diag(R) * resvar[i]
    }
    sqrt(se)
}

calculate.lm.pval=function(linear.model){
    # Given a solved linear model, return the p-values of the coefficients
    #print(linear.model$coefficients)
    n.obs = dim(linear.model$residuals)[1]
    n.models = dim(linear.model$coefficients)[2]
    cc = matrix(coef(linear.model), nrow=linear.model$rank, ncol=n.models)
    se = calculate.linear.model.se(linear.model)
    t.stat = cc/se
    rdf = n.obs - linear.model$rank
    dm=data.frame( 2*pt(abs(t.stat), rdf, lower.tail = FALSE)  )
    names(dm) = as.vector(dimnames(linear.model$coefficients)[[2]])
    rn = dimnames(linear.model$coefficients)[[1]]
    rownames(dm) = rn[1:dim(dm)[1]]
    dm
}

calculate.lm.pval.from.parts=function(coefficients, residuals, rank, qr){
    # Given a solved linear model, return the p-values of the coefficients
    n.obs = dim(residuals)[2]
    if( is.null( n.obs ) )
        n.obs = length(residuals)
    coefficients = matrix( coefficients, nrow=length(coefficients), ncol=1)
    n.models = dim(coefficients)[2]
    print( paste(n.obs, "observations", n.models, "models") )
    cc = matrix(coefficients, nrow=rank, ncol=n.models)
    Qr = qr
    p = rank
    p1 = 1L:p
    R = chol2inv(Qr$qr[p1, p1, drop = FALSE])
    r = residuals
    rss = as.numeric(colSums(r^2))
    n = NROW(Qr$qr)
    rdf = n - p
    resvar = rss/rdf
    se = matrix(nrow=p, ncol=n.models)
    for(i in 1:n.models){
        se[,i] = diag(R) * resvar[i]
    }
    se = sqrt(se)
    t.stat = cc/se
    rdf = n.obs - rank
    print(rdf)
    dm=data.frame( 2*pt(abs(t.stat), rdf, lower.tail = FALSE)  )
    names(dm) = as.vector(dimnames(coefficients)[[2]])
    rownames(dm) = dimnames(coefficients)[[1]]
    dm
}

calculate.EQTL = function(dataset, shared.col, covariate=NULL, gene.idx=NULL, n.perm=0, is.verbose=T){
    # For each probe in the dataset, identify the SNP whose genotype is most 
    # strongly associated with value of that probe.
    # shared.col must be a vector with length equal to the number of rows in 
    # dataset$s.snp, indicating the identifier in dataset$s.gene that corresponds
    # to that row in dataset$s.snp.
    # Optionally pass a covariate (should be a factor).
    # Optionally pass a vector of indexes to restrict the probesets.
    ids = match.identifiers.eQTL( dataset$s.snp, dataset$e.gene, shared.col )
    G = vector("list", dim(dataset$g.snp)[1] )
    M = data.matrix(dataset$e.snp)
    for(i in 1:dim(dataset$g.snp)[1] ){
        G[[i]] = as.factor( M[i,ids$idx.in.snps] )
    }
    if(is.null(gene.idx)){
        gene.idx = 1:dim(dataset$e.gene)[1]
    }
    n.probes = length(gene.idx)
    n.snps = dim(dataset$g.snp)[1]
    P = matrix(nrow=n.probes, ncol=length(ids$idx.in.expr))
    PHENOS = data.matrix(dataset$e.gene)
    for(i in 1:length(ids$idx.in.expr)){
        P[,i] = PHENOS[gene.idx,ids$idx.in.expr[i]]
    }
    gene.names = dataset$gene.names
    if(!is.null(gene.names)){
        gene.names = gene.names[gene.idx]
    }
    min.pvals.obs = rep(1, n.probes)
    min.idxs = rep(-1, n.probes)
    if(is.verbose){ print(paste("Calculating eQTL for",n.snps,"snps,", n.probes, "phenotypes.")) }
    for(i in 1:n.snps){
        if( i %% 10 == 0 & is.verbose ){ print(paste(date(), ": completed", i , "of", n.snps) ) }
        if( is.null(covariate) ){
            pvals.obs = as.numeric(calculate.lm.pval(lm( t(P) ~ G[[i]] ))[2,])
        }
        else
            pvals.obs = as.numeric(calculate.lm.pval(lm( t(P) ~ covariate + G[[i]] ))[3,])
        is.smaller = pvals.obs < min.pvals.obs
        min.pvals.obs[is.smaller] = pvals.obs[is.smaller]
        min.idxs[is.smaller] = i
    }
    if(n.perm>0){
        if(is.verbose){ print(paste("Calculating", n.perm, "permutations")) }
        permutation.pvals = rep(n.perm, n.probes)
        for(p in 1:n.perm){
            min.pvals.perm = rep(1, n.probes)
            for(i in 1:n.snps){
                geno.perm = sample(G[[i]])
                if( i %% 10 == 0  & is.verbose){ print(paste(date(), ": PERM",p,"of",n.perm,"completed", i , "of", n.snps) ) }
                if( is.null(covariate) )
                    pvals.perm = as.numeric(calculate.lm.pval( lm( t(P) ~ geno.perm ) )[2,])
                else
                    pvals.perm = as.numeric(calculate.lm.pval(lm( t(P) ~ covariate + geno.perm ))[3,])
                is.smaller = pvals.perm < min.pvals.perm
                min.pvals.perm[is.smaller] = pvals.perm[is.smaller]
            }
            is.smaller = min.pvals.obs < min.pvals.perm
            permutation.pvals[is.smaller] = permutation.pvals[is.smaller] - 1
        }
        permutation.pvals = permutation.pvals / n.perm
    }    
    probes.gene = rownames(dataset$g.gene)[gene.idx]
    probes.snp = rownames(dataset$g.snp)[min.idxs]
    df=data.frame(probes.gene, probes.snp, min.pvals.obs)
    if(!is.null(gene.names)){
        df = cbind(gene.names, df)
    }
    if( n.perm>0 ){
        df = cbind(df, pvals.perm=permutation.pvals)
    }
    df
}



eqtl.filter = function(eqtl, snp, ga.chr.snp, ga.loc.snp, gene, ga.chr.gene, ga.loc.gene, type="cis", window=0){
    valid.snp = !is.na(ga.chr.snp) & !is.na(ga.loc.snp) 
    valid.gene = !is.na(ga.chr.gene) & !is.na(ga.loc.gene)
    ga.chr.snp = ga.chr.snp[valid.snp]
    ga.loc.snp = ga.loc.snp[valid.snp]
    ga.chr.gene = ga.chr.gene[valid.gene]
    ga.loc.gene = ga.loc.gene[valid.gene]
    hsh.gene = hsh_from_vectors(as.character(gene), 1:length(gene) ) 
    hsh.snp = hsh_from_vectors(as.character(snp), 1:length(snp) )
    idx.gene = hsh_get(hsh.gene, as.character(eqtl$probe) )
    idx.snp = hsh_get(hsh.snp, as.character(eqtl$snp) )
    e.chr.snp = ga.chr.snp[idx.snp]
    e.loc.snp = ga.loc.snp[idx.snp]
    e.chr.gene = ga.chr.gene[idx.gene]
    e.loc.gene = ga.loc.gene[idx.gene]    
    eqtl = cbind(eqtl, e.chr.snp, e.loc.snp, e.chr.gene, e.loc.gene, stringsAsFactors=F)
    if( type=="cis" ){
        if(window==0){
            eqtl[e.chr.snp == e.chr.gene,]
        }
        else{
            eqtl[e.chr.snp == e.chr.gene & abs(e.loc.snp-e.loc.gene)<=window,]
        }
    }
    else{
        if( window==0){
            eqtl[e.chr.snp != e.chr.gene,]
        }
        else{
            eqtl[e.chr.snp != e.chr.gene | 
                 ( e.chr.snp == e.chr.gene & abs(e.loc.snp-e.loc.gene)>window) ,]
        }
    }   
}


write.rqtl.csv=function( expr, probe.list, sa.geno, sa.col.matching.expr, geno.chrom, calls.geno, fn.rqtl ){
    # write all probes in probe.list as phenotypes for an RQTL CSV file
    sa.match = sa.geno[,which(names(sa.geno)==sa.col.matching.expr)]
    m = match.idx(sa.match, names(expr))
    expr = expr[,m$idx.B]
    sa.geno = sa.geno[m$idx.A,]

    expr.top = t(expr[ match.idx(probe.list, rownames(expr))$idx.B, ])
    rqtl.pheno = rbind( rep("", dim(expr.top)[2] ), expr.top )
    if(sum(rownames(ga.geno)!=rownames(calls.geno))>0)
        stop("WARNING: gene attributes not equal to calls")
    calls.geno = t(calls.geno[,m$idx.A])
    rqtl.geno = rbind( geno.chrom, calls.geno )
    rownames(rqtl.geno) = rownames(rqtl.pheno)
    rqtl.csv = cbind(rqtl.pheno, rqtl.geno)
    rownames = append( 'IDENTIFIER', dimnames(rqtl.csv)[[2]] )
    write( rownames, fn.rqtl, sep=',', ncolumns=length(rownames) )
    write.table(rqtl.csv, fn.rqtl, quote=F, sep=',', row.names=T, na="NA", col.names=F, append=T)
}


generate.rQTL.matrix=function( pheno.ids, pheno, e.snp, ga.snp, sa.snp, fn.out){
    # Generates a matrix that can be written to be read by R/QTL
    # writes that matrix to fn.out and returns the matrix
    #
    # pheno.ids and pheno should have the same order (taken from the same data frame)
    #
    common.ids = set.intersection(rownames(sa.snp), pheno.ids)
    pheno.idx = vector()
    geno.idx = vector()
    sexes = rep("NA", length(common.ids) )
    idx.sex = which(names(sa.snp)=="sex")
    for(i in 1:length(common.ids) ){
        pheno.idx[i] = which(pheno.ids==common.ids[i])
        geno.idx[i] = which(names(e.snp)==common.ids[i])
        if( length(idx.sex)==1 ){
        	sexes[i] = sa.snp$sex[ which(rownames(sa.snp)==common.ids[i]) ]
        }
    }
    if( length(geno.idx)==0 | length(pheno.idx)==0 )
        stop("No genotypes selected; check that e.snp sample names match pheno.ids")
    if( is.vector(pheno) ){
        pheno.names = c('pheno')
        pheno = pheno[pheno.idx]
    }
    else{
        pheno.names = names(pheno)
        pheno = pheno[pheno.idx, ]
    }
    n.col=dim(ga.snp)[1] + 2 + length(pheno.names)
    has.cM = F
    has.sex = F
    offset = 2
    if( "cM" %in% names(ga.snp) ){
        has.cM = T
        offset = 3
    }
    M = matrix( nrow=length(common.ids) + offset, ncol=n.col)
    M[1,1] = 'IDENTIFIER'
    M[1,2] = 'sex'
    last.pheno.col=3
    for(i in 1:length(pheno.names) ){
        M[1,last.pheno.col] = pheno.names[i]
        last.pheno.col = last.pheno.col+1
    }
    M[2:3,1:last.pheno.col] = ''
    M[1, (last.pheno.col):n.col] = rownames(ga.snp) # ROW 1: SNP names
    M[2, (last.pheno.col):n.col] = ga.snp$Chr       # ROW 2: Chromosomes
    if( has.cM ){
        M[3, (last.pheno.col):n.col] = ga.snp$cM    # OPTIONAL ROW 3: cM
    }
    id.row = offset+1                               # First row of genotypes
    geno.rows = id.row:dim(M)[1]                    # row idx from id.row to end
    M[geno.rows,1] = common.ids                     # sample IDs in first column
    M[geno.rows,2] = sexes                          # sexes in second column
    if( length(pheno.names)==1 ){                   
        M[geno.rows,3] = pheno 					    # If only one pheno ID, third column
    }
    else{
        print(dim(M))
        print( 3:(2+dim(pheno)[2]) )
        print( dim(data.matrix(pheno)))
        M[ geno.rows, 3:(2+dim(pheno)[2]) ] = data.matrix(pheno) # else, third through X columns
    }
    genotypes = t(e.snp[,geno.idx])
    
    M[ geno.rows, last.pheno.col: dim(M)[2] ] = genotypes
    write.table(M, fn.out ,quote=F,sep=',',row.names=F,col.names=F)
    M
}


pheno.from.expr = function( sa.SNP, probe.id, expr, sa.shared.col, valid.pheno.ids){
    # Used to match probe expression to genotype ids.
    # Data frame useful for generate.rQTL.matrix()
    #
    ids.eqtl=match.identifiers.eQTL( sa.SNP, expr, sa.shared.col)
#    ids.two=set.intersection(as.character(ids.eqtl$id.e), rownames(sa.SNP)[sa.valid.ids])
	idx.valid = which(ids.eqtl$id.e %in% valid.pheno.ids)
    pheno.ids = as.character(ids.eqtl$id.s[ idx.valid ])
    pheno.idx = ids.eqtl$idx.in.expr[ idx.valid]
    pheno.vals = as.numeric(expr[ which(rownames(expr)==probe.id), pheno.idx ])
    data.frame(pheno.ids, pheno.vals)
}

genotype.for.expression.ids = function( probe.snp, expr, expr.SNP, sa.SNP, sa.SNP.shared.column ){
    # Identify genotype at probe.snp for expression ids in expr
    #
    ids=match.identifiers.eQTL( sa.SNP, expr, sa.SNP.shared.column )
    idx.snp = which(rownames(expr.SNP)==probe.snp)
    genos = rep(0, dim(ids)[1] )
    for(i in 1:dim(ids)[1]){
        genos[ i ] = expr.SNP[ idx.snp, which(names(expr.SNP)==ids$id.s[i]) ]
    }
    out=data.frame( ids$id.e, ids$id.s, genos, stringsAsFactors=F )
    names(out) = c('id.expression', 'id.genotype', 'genotype')
    out
}

do.sam = function(E, ga, A, B, med.FDR=10){
    # Perform a turn-key SAM analysis returning values with FDR<med.FDR
    if( length(A) != length(B) | length(A) != dim(E)[2] )
        stop( "incorrectly formatted input: bad dimensions")
    labels = rep(0, dim(E)[2])
    labels[A] = 1
    labels[B] = 2
    E = E[,A|B]
    labels = labels[A|B]
    data = list(x=data.matrix(E), y=labels, geneid=rownames(ga), genenames = ga$symbol, logged2=TRUE)
    samr.obj = samr(data, resp.type="Two class unpaired", nperms=100)
    delta.table = samr.compute.delta.table(samr.obj)
    print(delta.table)
    deltas = delta.table[,1]
    med.fdrs = delta.table[,5]
    med.fdrs[is.nan(med.fdrs)]=1
    if( med.FDR == 100 ){
        delta=0
    }
    else if( min(med.fdrs) > (med.FDR/100) ){
        delta = min(deltas[med.fdrs==min(med.fdrs)])
    }
    else{
        delta = min(deltas[med.fdrs<=(med.FDR+1)/100])
    }
    if( is.infinite(delta) ){
        NA
    }
    else{
        print(paste("Delta: ", delta) )
        sgt = samr.compute.siggenes.table( samr.obj, delta, data, delta.table)
        sig = SAM.convert.siggenes( sgt )
        sig = sig[sig$q.value.percent<=med.FDR,]
        sig = sig[order(sig$fold.change, decreasing=T),]
        fc = sig$fold.change
        fc[fc<1] = -1/fc[fc<1]
        fc[fc==1]=0
        sig = cbind(sig, fc=round(fc,2), stringsAsFactors=F )
        sig
    }
}  

do.sam.paired = function(E, ga, labels, med.FDR=10){
    # Perform a turn-key SAM analysis returning values with FDR<med.FDR
    if( length(labels) != dim(E)[2] )
        stop( "incorrectly formatted input: bad label dimensions")
    data = list(x=data.matrix(E), y=labels, geneid=rownames(ga), genenames = ga$symbol, logged2=TRUE)
    samr.obj = samr(data, resp.type="Two class paired", nperms=100)
    delta.table = samr.compute.delta.table(samr.obj)
    print(delta.table)
    deltas = delta.table[,1]
    med.fdrs = delta.table[,5]
    med.fdrs[is.nan(med.fdrs)]=1
    if( min(med.fdrs) > (med.FDR/100) ){
        delta = min(deltas[med.fdrs==min(med.fdrs)])
    }
    else{
        delta = min(deltas[med.fdrs<=(med.FDR+1)/100])
    }
    if( is.infinite(delta) ){
        NA
    }
    else{
        print(paste("Delta: ", delta) )
        sgt = samr.compute.siggenes.table( samr.obj, delta, data, delta.table)
        sig = SAM.convert.siggenes( sgt )
        sig = sig[sig$q.value.percent<=med.FDR,]
        sig
    }
}

SAM.data = function(expr, ga, valid.p, s1, s2, nperms=100){
    # Wrapper to create data object for unpaired call to SAM
    library(samr)
    gene.names = as.character(ga[valid.p,1])
    gene.ids = as.character(rownames(ga)[valid.p])
    valid.s = append(s1, s2)
    if( length( set.intersection(s1, s2) )>0 )
        stop("s1 and s2 are not disjoint")
    labels=append( rep(1,length(s1)), rep(2,length(s2)))
    M = as.matrix( cbind( expr[valid.p,s1], expr[valid.p, s2] ) )
    list(x=M,y=labels, geneid=gene.ids, genenames = gene.names, logged2=TRUE)
}

SAM.convert.siggenes = function(T){
   # combine nasty SAM output into a single dataframe
   #     Row     Gene ID    Gene Name      Score(d)            Numerator(r)        Denominator(s+s0)   Fold Change          q-value(%)
   no.lo = is.null(T$genes.lo)
   no.up = is.null(T$genes.up)   
   if( !no.lo ){
      rownames(T$genes.lo) = paste('down', 1:dim(T$genes.lo)[1], sep='' )
   }
   if( !no.up ){
      rownames(T$genes.up) = paste('up', 1:dim(T$genes.up)[1], sep='' )
   }
   if( !no.lo & !no.up ){
      genes = data.frame( rbind( T$genes.lo[,c(2,3)], T$genes.up[,c(2,3)]), stringsAsFactors=F )
      names(genes) = c("symbol", "probe.id")
      values = data.frame( rbind( T$genes.lo[,c(4:8)], T$genes.up[,c(4:8)] ) )
   }
   else if( no.lo & !no.up){
      genes = data.frame( T$genes.up[,c(2,3)], stringsAsFactors=F )
      names(genes) = c("symbol", "probe.id")
      values = data.frame( T$genes.up[,c(4:8)] )
   }   
   else if( no.up & !no.lo){
      genes = data.frame( T$genes.lo[,c(2,3)], stringsAsFactors=F )   
      names(genes) = c("symbol", "probe.id")
      values = data.frame( T$genes.lo[,c(4:8)] )      
   }
   if( !no.up | !no.lo){
      values[,1] = round( as.numeric( as.character(values[,1]) ), 3)
      values[,2] = round( as.numeric( as.character(values[,2]) ), 3)
      values[,3] = round( as.numeric( as.character(values[,3]) ), 3)
      values[,4] = round( as.numeric( as.character(values[,4]) ), 3)
      values[,5] = round( as.numeric( as.character(values[,5]) ), 3)
      names(values) = c('score', 'numerator.r', 'denominator.s.plus.s0', 'fold.change', 'q.value.percent')
      T = cbind(genes, values)
      T
   }
   else{
      data.frame(score=c(), numerator.r=c(), denominator.s.plus.s0=c(), fold.change=c(), q.value.percent=c())
   }
}


km=function( times, had.events, conditions, main=NULL, legends=NULL, fn_out=NULL, verbose=T, legend.x=3, legend.y=0.25){
    # Wrapper for Kaplan-meier analysis (using library "survival")
    # If fn_out is passed, the plot is saved to the specified file as a png
    # returns the coxph object

    library(survival)
    if( is.null(times) ){ stop( "times is NULL") }
    if( is.null(had.events) ){ stop( "had.events is NULL") }
    if( is.null(conditions) ){ stop( "conditions is NULL") }
    if( length(times) != length(had.events) || length(times) != length(conditions) ){
        stop("All of times, had.events, and conditions do not have equal length")
    }

    surv.all = survfit(Surv(times, had.events)~conditions)
    cox = coxph(formula = Surv(times,had.events==1)~conditions)
    p.val = as.numeric(summary(cox)$logtest[3])
    if(p.val>0.01)
        p.val = round(p.val, 3)
    else
        p.val = signif(p.val, 3)
    if( !is.null(fn_out) ){
        png(fn_out)
    }
    color.list =c("black", "blue", "gray", "red", "darkgreen", "orange")
    if( is.null(main) )
        main.title=paste('KM curve')
    else
        main.title=main
    if( length(unique(conditions))== 2 ){
        main.title = paste(main.title,'P =', p.val, sep=' ')
    }

    plot(surv.all, lwd=3, col=color.list,
        main=main.title, las=1,
        xlab="time")
    if( !is.null(legends) ){
        legend(legend.x,legend.y, legends, col=color.list,lty=1,lwd=3)
    }
    else{
        legend(legend.x, legend.y, sort(unique(conditions,na.rm=T)), col=color.list,lty=1,lwd=3)
    }
    if( !is.null(fn_out) ){
        i=dev.off()
    }
    if( verbose ){
        print( summary(cox) )
        print( paste( sum(had.events), " subjects had the event", sep='' ) )
        print( paste( length(times)-sum(conditions), " subjects had condition=0", sep='' ) )
        print( paste( sum(conditions), " subjects had condition=1", sep='' ) )
    }
    cox
}

km.multi = function( times, had.events, conditions ){
    #    General Kaplan-meier plotting and p-value calculation, allowing multiple conditions to be passed
    library(survival)
    if( is.null(times) ){ stop( "times is NULL") }
    if( is.null(had.events) ){ stop( "had.events is NULL") }
    if( is.null(conditions) ){ stop( "conditions is NULL") }
    if( length(times) != length(had.events) || length(times) != dim(conditions)[2] ){
        stop("All of times, had.events, and conditions do not have equal length")
    }
    if( !is.data.frame(conditions) ){ stop("conditions should be a data frame")}
    if( max(conditions, na.rm=T)>1 || min(conditions, na.rm=T) < 0 ){ stop("conditions should only consist of 0 or 1") }
    n.conditions = dim(conditions)[1]
    n.samples = dim(conditions)[2]
    pvals = vector(mode="numeric", n.conditions)
    count_0 = vector(mode="integer", n.conditions)
    count_1 = vector(mode="integer", n.conditions)
    for(i in 1:n.conditions){
        cond = as.numeric(conditions[i,])
        count_1[i] = sum( cond, na.rm=T )
        count_0[i] = sum( cond==0, na.rm=T)
        pval = 1
        if( count_1[i] > 5 & count_0[i] > 5 ){
            surv.all = survfit(Surv(times, had.events)~ as.numeric(cond))
            cox = coxph(formula = Surv(times,had.events==1)~as.numeric(cond))
            pval = as.numeric(summary(cox)$logtest[3])
        }
        pvals[i] = pval
        if( i %% 500==0 ){
            print(paste("Processing KM ",i,"of",n.conditions))
        }
    }
    d = data.frame( pvals, count_0, count_1)
    rownames(d) = rownames(conditions)
    d
}


##############
# BEGIN GSCA #
##############
# (adapted from Choi & Kendziorski Bioinformatics 2009)
# 
get.idx = function(data, probe.lists){
	# Utility function for singleDC, takes a list of named probe vectors and
	# converts those into indices in data.
    p2idx = hsh_from_vectors(rownames(data), seq(from=1,to=length(rownames) ) )
	for(j in 1:length(probe.lists) ){
        probes = probe.lists[[j]]
        gene.vector = vector()
		valid.probes = vector()
        for(i in 1:length(probes)){
			#idx = which(rownames(data)==probes[i])
			if( hsh_in(p2idx, probe[i] ) ){
				idx = hsh_get(p2idx, probe[i] )
				gene.vector[ length(gene.vector) + 1 ] = idx 
				valid.probes[length(valid.probes) + 1 ] = probes[i]
			}
			else{
				warning("Unable to find index for probe ", probes[i])
			}
        }
		names(gene.vector) = valid.probes
        probe.lists[[j]] = gene.vector
    }
    probe.lists
}

singleDC <- function(data, group, probe.lists, nperm, method='spearman', focus.probe=NA){
	# Performs differential correlation analysis.
	# data: matrix of quantitative values 
	# group: vector of 0 or 1 indicating sample group
	# probe.lists: list of named vectors, where vector contains probe IDs or index.
	#   lists will be searched for enrichment
	#   finding index in R is slow, so you can pass index as integers 
    # focus.probe: if passed, look for DC of all genes in lists with this probe only
    print(paste("Calculating with", length(probe.lists), "gene sets."))
    DI=rep(0,length(probe.lists))
    pvalue=rep(1,length(probe.lists))
    permv=matrix(0,length(probe.lists),nperm)
    if( is.numeric(gene.sets[[1]][1]) )
    	GSdefList = probe.lists
	else    
		GSdefList = get.idx(data, probe.lists)
    for(j in 1:length(GSdefList)){
        print(paste("Gene list", j, names(probe.lists)[j] ) )
        probe.idx = GSdefList[[j]]
        if( length(probe.idx)>1 )
            DI[j]=compute.test.stat(data[probe.idx,], group, method, focus.probe)
        else
            DI[j] = NA
    }
    print('Calculated Observed DI')
    if(nperm>0){
        for(i in 1:nperm){
            data.perm=data[,sample(length(data[1,]))]
            for(j in 1:length(GSdefList)){
                probe.idx = GSdefList[[j]]
                if(length(probe.idx)>1)
                    permv[j,i]=compute.test.stat(data.perm[probe.idx,], group, method, focus.probe)
                else
                    permv[j,i]=NA
            }
            print(paste('Perm',i,'completed'))
        }
        pvalue=rowSums(permv>=DI, na.rm=T) / (rowSums( !is.na(permv) ))
    }
    ret=NULL
    names(DI)=names(GSdefList)
    names(pvalue)=names(GSdefList)
    rownames(permv) = names(GSdefList)
    permnames=NULL
    #for(m in 1:nperm){
    #    permnames[m]=paste('P',m,sep='')
    #}
    names(permv) = permnames
    ret$DI=DI
    ret$pvalue=pvalue
    ret$permv=permv
    return(ret)

	# Debugging code
	#xx=c(1,1,4,2,2,3,3,3,2,4,3,1,5,1,5,6,2,2,7,9,3,8,0,1)
	#dim(xx) = c(3,8)
	#xx=data.frame(xx)
	#names(xx) = c('g1','g2','g3','g4','g5','g6','g7','g8')
	#rownames(xx) = c('p1', 'p2', 'p3')
	#> xx
	#   g1 g2 g3 g4 g5 g6 g7 g8
	##p1  1  2  3  4  5  6  7  8
	##p2  1  2  3  3  1  2  9  0
	##p3  4  3  2  1  5  2  3  1
	#singleDC(as.matrix(xx), c(0,0,0,0,1,1,1,1), list(test=c('p1', 'p2', 'p3')) )
}


compute.test.stat = function(d, group, m='spearman', focus.probe=NA){
    # utility function for singleDC
	n.pairs = choose(nrow(d),2)
    corr.mats=NULL
    g1 = d[,group==0]
    g2 = d[,group==1]
    s1 = rowSums(is.na(g1))
    s2 = rowSums(is.na(g2))
    if(is.na(focus.probe)){
        corr.mats[[1]] = cor(t(g1), method=m, use="na.or.complete")
        corr.mats[[2]] = cor(t(g2), method=m, use="na.or.complete")
        if( is.na( corr.mats[[1]][1] ) ){ return( NA ) }
        if( is.na( corr.mats[[2]][1] ) ){ return( NA ) }
    }
    else{
        corr.mats[[1]] = matrix(NA, dim(g1)[1], 1); 
        corr.mats[[2]] = matrix(NA, dim(g2)[1], 1);
        idx.fp = which(rownames(d)==focus.probe)
        focus.c1 = as.numeric( g1[ idx.fp, ] )
        focus.c2 = as.numeric( g2[ idx.fp, ] )
        for(i in 1:dim(g1)[1]){
            if(idx.fp==i){
                corr.mats[[1]][i]=NA
                corr.mats[[2]][i]=NA
            }
            else{
                corr.mats[[1]][i] = cor(focus.c1, as.numeric(g1[i,]), use="na.or.complete", method=m)
                corr.mats[[2]][i] = cor(focus.c2, as.numeric(g2[i,]), use="na.or.complete", method=m)
                if( is.na(corr.mats[[1]][i] ) ){ return( NA ) }
                if( is.na(corr.mats[[2]][i] ) ){ return( NA ) }
            }
        }   
    }
    vals=NULL
    k=1
    n.groups=2
    for(i in 1:(n.groups-1)){
        for(j in 2:n.groups){
            sq.diff.mat=(corr.mats[[i]]-corr.mats[[j]])^2
            if(is.na(focus.probe))
                vals[[k]]=sqrt((1/n.pairs)*sum(sq.diff.mat[upper.tri(sq.diff.mat)],na.rm=T))
            else
                vals[[k]]=sqrt((1/n.pairs)*sum(sq.diff.mat,na.rm=T))
          k=k+1
       }
    }
    test.stat=mean(vals)
    return(test.stat)
}

###################
# PSCBS functions #
###################
#
psCBS.plot = function( S, A, B, S.samples, CN.samples, sample, chrom, loci, chromosome=NULL, loc_start=NULL, loc_end=NULL, Y.MAX=5){
    # S is psCBS output for the sample of interest. 
    # A and B are raw copy number calls. 
    # chrom and loci identify the physical location of rows (SNPs) in rawCN.A and rawCN.B
    if( length(A) != length(B) ){ stop("A and B have different length")  }
    if( length(chrom) != length(loci) ){ stop("chrom has diffferent length from number of loc_start" ) }
    if( !is.null(loc_start) & is.null(chromosome) ){ stop("Must pass chromosome if limiting by locus") }
    if( !is.null(loc_end) & is.null(chromosome) ){ stop("Must pass chromosome if limiting by locus") }
    if( length(chrom) != dim(A)[1] ){ stop("chrom has diffferent length from number of rows in A" ) }       
    if( length(S.samples) != dim(S)[1] ){ stop("S.samples has different length from number of rows in s") }
    if( length(CN.samples) != dim(A)[2] ){ stop("CN.samples has different length from number of columns in A")}
    if( length( which(CN.samples==sample))==0 ){ stop("sample was not found in CN.samples") }
    idx.seg = which(as.character(S.samples)==sample)
    idx.pat = which(as.character(CN.samples)==sample)
    if( length(idx.seg)==0 ){ stop("No segments found for sample")}
    A = as.numeric(A[,idx.pat])
    B = as.numeric(B[,idx.pat])
    S = S[idx.seg,]
    A[A>30]=NA
    B[B>30]=NA
    A[A<0] = NA
    B[B<0] = NA
    ds = data.frame(chrom, loci, A, B)
    ds = ds[order(ds$chrom, ds$loci), ]
    if( !is.null(chromosome) ){
        ds = ds[ds$chrom==chromosome,]
        S = S[S$chrom==chromosome,]
        print(paste("Limiting to chromosome",chromosome))
    }
    if( !is.null(loc_start) & !is.null(loc_end) ){
        ds = ds[ds$loci>=loc_start & ds$loci<=loc_end,]
        S = S[ (S$loc.start<=loc_start & S$loc.end >= loc_start) | (S$loc.start<=loc_end & S$loc.end >= loc_end) | (S$loc.start>loc_start & S$loc.end < loc_end), ]
        print(paste("Limiting to range",loc_start,"-",loc_end))
        print(S)
        if(dim(S)[1]==0){ stop("No segments in range requested") }
    }
    A = ds$A
    B = ds$B
    totalCN = A+B
    totalCN[totalCN>30]=NA
    totalCN[totalCN<=0] = 0    
    A.is.major = which(A>B)
    B.is.major = which(B>=A)
    minor=A
    minor[A.is.major] = B[A.is.major]
    plot( totalCN, col='black', ylim=c(0,Y.MAX), pch='.', xlab="", ylab="Copy Number", main=sample, axes=F)
    par(cex=0.25)
    points( minor, col='blue', ylim=c(0,Y.MAX), pch=20)
    par(cex=1)
    axis(2, at=(0:Y.MAX) )
    x.prev = 0
    if( is.null(chromosome) ){
        for(i in 2:23){
            x.idx = min(which(ds$chrom==i))
            lines( c(x.idx, x.idx), c(0, Y.MAX),col='darkgray')
            text( x.prev + (x.idx - x.prev)/2, Y.MAX, as.character(i-1) )
            x.prev = x.idx
            }
        text( x.prev + (length(ds$chrom) - x.prev)/2, Y.MAX, 'X' )
    }
    else{
        lbl = round( ds$loci / 1000000, 2 )
        lbl.b = 1
        lbl.e = length(lbl)
        scale = round((lbl.e - lbl.b)/10, 2)
        axis(1,labels=lbl[seq(lbl.b, lbl.e, scale)], at=seq(lbl.b, lbl.e, scale) )
    }
    for(i in 1:dim(S)[1]){
        x.begin = which( ds$chrom==S$chrom[i] & ds$loci==S$loc.start[i] )
        x.end = which( ds$chrom==S$chrom[i] & ds$loci==S$loc.end[i] )
        if(length(x.begin)==0){ 
            x.begin = 0 
        }
        if(length(x.end)==0){ 
            x.end = length(ds$loci)
        }
        y.minor = S$min.mean[i]
        y.total = S$total.mean[i]
        if( !is.na(y.minor) ){
            lines( c(x.begin, x.end), c(y.minor, y.minor), col='blue' )
            lines( c(x.begin, x.end), c(y.total, y.total), col='red' )
        }
    }
    box()
}


psCBS.filter = function(s, min.mark=15, TC.range=c(0,30), min.range=c(0,30), max.range=c(0,30), min.diff=0, ga=NULL){
    # Expects output from Adam Olshen's psCBS as s
    # If ga is passed, must have "Chr" and "loc_start" columns. 
    # Count of the number of times a probe is within a range specified by the output
    # is placed in an additional ga column called "counts". 
    # Returns a list with "results" and "ga", where "ga" is sorted in increasing genomic order.
    #
    if( length(TC.range) != 2 ){ stop("TC.range must have two elements" ) }
    if( length(min.range) != 2 ){ stop("min.range must have two elements" ) }
    if( length(max.range) != 2 ){ stop("max.range must have two elements" ) }
    mu.tc.lo = TC.range[1];   mu.tc.hi = TC.range[2]
    mu.min.lo = min.range[1]; mu.min.hi = min.range[2]
    mu.max.lo = max.range[1]; mu.max.hi = max.range[2]
    s = s[!is.na(s$max.mean),]
    s = s[ s$num.mark>=min.mark & s$total.mean<mu.tc.hi & s$total.mean>mu.tc.lo, ]
    s = s[ s$min.mean<mu.min.hi & s$min.mean>mu.min.lo, ] 
    s = s[ s$max.mean<mu.max.hi & s$max.mean>mu.max.lo, ]     
    s = s[(s$max.mean-s$min.mean) >= min.diff, ]
    s$min.mean = round( s$min.mean,2)
    s$max.mean = round( s$max.mean,2)
    s$total.mean = round( s$total.mean,2)
    if(is.null(ga) ){
        print(paste("Found", dim(s)[1],"matches"))
        list('results'=s)
    }
    else{
        print(paste('Counting locations for', dim(s)[1], 'results') )
        counts = rep(0, dim(ga)[1])
        for(i in 1:dim(s)[1]){
            loc.match = ga$loc_start>= s$loc.start[i] & ga$loc_start <= s$loc.end[i]
            idx = which( ga$Chr==s$chrom[i] & loc.match )
            counts[idx] = counts[idx] + 1
            if( i %% 100==0 ){ print(paste("Completed", i ) ) }
        } 
        ga.ord = cbind(ga, counts)
        ga.ord = ga.ord[order(ga.ord$Chr, ga.ord$loc_start),]
        list( 'results'=s, 'ga'=ga.ord )
    }
}

psCBS.dissect = function(r, ga, chr=NULL, loc=NULL, do.plot=F ){
    # Currently assumes results (r) has a patient.no column
    # If only ga is passed, calculates maximum number of events on each chromosome.
    # If chr is also passed, returns events on that chromosome and the bounds of all events on that chromosome
    # If loc is also passed, restricts events to those that overlap that location on chromosome chr
    max.on.chr = rep(0, 23)
    Chr = seq(from=1,to=23)
    for(i in 1:23){
        max.on.chr[i] = max( ga$counts[ ga$Chr==i] )
    }
    max.on.chr = data.frame(max.on.chr)
    rownames(max.on.chr)=Chr
    print( t(max.on.chr) )
    if( is.null(chr) ){
        t(max.on.chr)
    }
    else{
        if( is.null(loc) ){
            region = r[r$chrom==chr,]
            region = region[ order(region$loc.start), ]
            ga.max = ga[ga$Chr==chr & ga$counts==max.on.chr[chr,1],]
            print(paste("From",min(ga.max$loc_start,na.rm=T),"to",max(ga.max$loc_start,na.rm=T)))
        }
        else{
            region = r[r$chrom==chr & r$loc.start <= loc & r$loc.end >= loc,]
            if(dim(region)[1]==0){ stop( "No events found that overlap this bound" )}
            print("Bounds:")
            print(paste("chr",chr,":",max(region$loc.start),"-",min(region$loc.end),sep='') )
        }
        region = data.frame(region$patient.no, region$chrom, region$loc.start, region$loc.end, region$num.mark, region$num.hetero, region$min.mean, region$max.mean, region$total.mean)
        names(region) = c('patient.no', 'Chr', 'loc.start', 'loc.end', 'num.mark', 'num.hetero', 'min.mean', 'max.mean', 'total.mean')
        if( do.plot ){
            par(cex=0.75)
            samples = sort(unique(region$patient.no))
            n.samples = length(samples)
            loc_start = min(region$loc.start)
            loc_end = max(region$loc.end)
            locations = sort(ga$loc_start[ ga$loc_start>=loc_start & ga$loc_start<= loc_end ])
            lbl = round( locations / 1000000, 2 )
            lbl.b = 1
            lbl.e = length(lbl)
            scale = round((lbl.e - lbl.b)/10, 2)
            plot(c(0,lbl.e), c(0, n.samples), col='white', axes=F, ylab="Sample", xlab='Megabases', main=paste('Chromosome',chr) )
            axis(1,labels=lbl[seq(lbl.b, lbl.e, scale)], at=seq(lbl.b, lbl.e, scale) )
            for(i in 1:dim(region)[1]){
                y.idx = which(samples==region$patient.no[i])
                l.start = which(locations==region$loc.start[i])[1]
                l.end = which(locations==region$loc.end[i])[1]
                lines( c(l.start,l.end ), c(y.idx, y.idx) )
            }
            for(i in 1:n.samples){
                text( rep(-0.5,n.samples), seq(from=1, to=n.samples), labels=samples )       
            }
            par(cex=1)
        }
        region
    }
}


psCBS.show.het.alleles = function(sa, ga, is.het, cn.t.A, cn.t.B, cn.n.A, cn.n.B, rs.id){
    # reports the raw values for heterozygotes 
    idx = which(rownames(ga)==rs.id )
    if(length(idx)==0){ stop("rs.id not found") }
    tumor.A = cn.t.A[idx, is.het[idx,]]
    tumor.B = cn.t.B[idx, is.het[idx,]]
    normal.A = cn.n.A[idx, is.het[idx,]]
    normal.B = cn.n.B[idx, is.het[idx,]]    
    vals = data.frame(tumor.A, tumor.B, normal.A, normal.B  )
    patients = sa$patient.no[sa$valid==1 & sa$shared==1]
    rownames(vals) =  patients[ is.het[idx,] ]
    vals
}

psCBS.copies.at.rs = function( ga, patient.list, rs.id, segments ){
    # find the called values for rs.id A and rs.id B for each patient in patient.list
    A = rep(0, length(patient.list) )
    B = rep(0, length(patient.list) )
    for(i in 1:length(A) ){
        
    }
}


psCBS.find.het = function(ga, s.names, s.idx, chrom, loc_begin, loc_end, is.het, A.t, B.t, A.n, B.n ){
    if( dim(A.t)[1] != dim(B.t)[1]){ stop("Number of rows in A.t does not match B.t") }
    if( dim(B.t)[2] != dim(B.t)[2]){ stop("Number of columns in A.t does not match B.t") }
    if( dim(A.n)[1] != dim(B.n)[1]){ stop("Number of rows in A.n does not match B.n") }
    if( dim(B.n)[2] != dim(B.n)[2]){ stop("Number of columns in A.n does not match B.n") }
    if( dim(A.n)[2] != dim(A.t)[2]){ stop("Number of columns in A.n does not match A.t") }
    if( dim(B.n)[2] != dim(B.t)[2]){ stop("Number of columns in B.n does not match B.t") }
    if( dim(A.t)[1] != dim(is.het)[1]){ stop("Number of rows in A.t does not match is.het") }
    if( dim(A.t)[2] != dim(is.het)[2]){ stop("Number of columns in A.t does not match is.het") }
    if( loc_end < loc_begin ){ stop("loc_begin must be smaller than loc_end") }
    if( max(s.idx) > dim(A.t)[2] ){ stop("Sample.idx larger than number of columns in A.t") }
    if( length(s.idx) != length(s.names) ){ stop("Sample.idx not equal in length to s.names") }    
    ga.ord = cbind(seq(from=1, to=dim(ga)[1]), ga)
    names(ga.ord)[1] = 'orig'
    ga.ord = ga.ord[order(ga.ord$Chr, ga.ord$loc_start),]
    A.t = A.t[ ga.ord$orig, ]
    B.t = B.t[ ga.ord$orig, ]
    A.n = A.n[ ga.ord$orig, ]
    B.n = B.n[ ga.ord$orig, ]
    is.het = is.het[ga.ord$orig, ]
    probe.in.range = ga.ord$Chr==chrom & ga.ord$loc_start >= loc_begin & ga.ord$loc_start <= loc_end
    if( sum(probe.in.range)==0 ){ 
        stop("No probes selected by chrom, loc_begin, loc_end")
    }
    else{
        print(paste("Identified",sum(probe.in.range),"total rows"))
    }
    ga.ord = ga.ord[probe.in.range, ]
    A.t = A.t[ probe.in.range, ]
    B.t = B.t[ probe.in.range, ]
    A.n = A.n[ probe.in.range, ]
    B.n = B.n[ probe.in.range, ]
    is.het = is.het[ probe.in.range,]
    n.samples = length(s.idx)
    r = list()
    limit.to.valid=T
    valid.rows = which( rowSums( is.het[ ,s.idx] )==length(s.idx) )
    if( length(valid.rows)==0 ){ 
        print("No SNPs are het in all samples") 
        limit.to.valid=F
        valid.rows = which( rowSums( is.het[ ,s.idx] ) > 1 )
    }
    for(i in 1:n.samples){
        AB = A.t[valid.rows, s.idx[i]] + B.t[valid.rows, s.idx[i]]
        S = cbind( ga.ord[valid.rows,], total=AB, A.t = A.t[valid.rows, s.idx[i]], B.t = B.t[valid.rows, s.idx[i]], A.n = A.n[valid.rows, s.idx[i]], B.n = B.n[valid.rows, s.idx[i] ] )
        r[[ s.names[i] ]] = S
    }
    if( length(valid.rows)==0 )
        print("No SNP in selected range is het in more than one sample.")
    else{
        bigger.A = data.frame(A.t[valid.rows, s.idx] > B.t[valid.rows, s.idx])
        names(bigger.A) = s.names
        rownames(bigger.A) = rownames(ga.ord)[valid.rows]
        bigger.A[ !is.het[valid.rows,s.idx] ] = NA
        r[['bigger.A']] = bigger.A
    }
    r
}


##################
### BEGIN HASH ###
##################

hsh_new = function(){
    new.env(hash=TRUE, parent=emptyenv()) 
}

hsh_in = function(H, key){
    exists(key, H)
}

hsh_get = function( H, key, na.if.not.found=F ){
    if( length(key)==1 ){
        if( na.if.not.found ){
            if( exists(key, H) )
                get(key, H)
            else
                NA
        }
        else{
            get(key, H)
        }
    }
    else{
        results = rep(0, length(key) )
        if( !na.if.not.found ){
            for(i in 1:length(key) ){
                if( exists(key[i], H) ){
                    results[i] = get(key[i], H )
                }
                else{
                    results[i] = NA
                }
            }
        }
        else{
            for(i in 1:length(key) ){
                results[i] = get(key[i], H )
            }
        }
        results
    }
}

hsh_set = function( H, key, value ){
    assign(key, value, envir=H)
}

hsh_keys = function( H ){
  return(sort(ls(H)))
}

hsh_keys_values = function( H ){
  keys = ls(H)
  values = hsh_get(H, keys)
  data.frame( keys, values, stringsAsFactors=F)
}

hsh_from_vectors = function( v1, v2 ){
    # Create a hash from vectors v1, v2 with keys from v1 and values from v2
    if( length(v1) != length(v2) ){
        stop("Length of v1 != length of v2")
    }
    H = hsh_new()
    for( i in 1:length(v1) ){
        hsh_set(H, v1[i], v2[i] )
    }
    H
}



#######################
### BEGIN HASHGRAPH ###
#######################

hashgraph <- function(...){
    G=hsh_new()
    attr(G, "class") = "hashgraph"
    class(G) = "hashgraph"
    return( G )
}

addEdge <- function(x, e1, e2, value){
  if(is.null(attr(x, "class"))){ stop("Must be called on a class") }
  else{ UseMethod("addEdge") }
}

addEdge.hashgraph = function(x, e1, e2, value=NULL, ...){
    if( !is.character(e1) | !is.character(e2) )
        stop("Edge names must be legal variable names in R.")
    if( length(e1) != length(e2) )
        stop("Edge vectors must have identical length")
    if( is.null(value) ){
        value = rep(1, length(e1))    
    }
    else{
        if( length(value) != length(e1) )
            stop("Value vector must have same length as edge vectors")
        if( !is.numeric(value) )
            stop("Value vector must be numeric")
    }
    for(i in 1:length(e1) ){
        ee1 = e1[i]
        ee2 = e2[i]
        if( !hsh_in(x, ee1) ){
            new.hash = hsh_new()
            hsh_set( new.hash, ee2, value[i])
            hsh_set(x, ee1, new.hash )
        }
        else{
            existing.hsh = hsh_get( x, ee1 )
            hsh_set( existing.hsh, ee2, value[i] )
        }
        if( !hsh_in( x, ee2 ) ){
            new.hash = hsh_new()
            hsh_set( new.hash, ee1, value[i])
            hsh_set(x, ee2, new.hash )
        }
        else{
            existing.hsh = hsh_get( x, ee2 )    
            hsh_set( existing.hsh, ee1, value[i] )
        }
        attr(x, "n.edges") = attr(x, "n.edges") + 1
    }
}


hasEdge <- function(x, e1, e2)  {
  if(is.null(attr(x, "class"))){ stop("Must be called on a class") }
  else{ UseMethod("hasEdge") }
}

hasEdge.hashgraph = function(x, e1, e2, ...){
    if( !is.character(e1) | !is.character(e2) )
        stop("Edge names must be legal variable names in R.")
    if( exists( e1, x ) ){
        existing.hsh = hsh_get(x, e1 )
        if( hsh_in( existing.hsh, e2 ) ){
            TRUE
        }
        else{
          FALSE
        }
    }
    else{
        FALSE
    }
}

getValue <- function(x, e1, e2){
  if(is.null(attr(x, "class"))){ stop("Must be called on a class") }
  else{ UseMethod("getValue") }
}

getValue.hashgraph = function(x, e1, e2, ...){
    if( hasEdge( x, e1, e2 )){
        hsh_get( hsh_get(x, e1), e2 )
    }
    else{
        stop( paste( "Edge {",e1,e2,"} not found") )
    }
}

nodes <- function(x)  {
  if(is.null(attr(x, "class"))){ stop("Must be called on a class") }
  else{ UseMethod("nodes") }
}

nodes.hashgraph = function(x){
  hsh_keys( x )
}

edges <- function(x)  {
  if(is.null(attr(x, "class"))){ stop("Must be called on a class") }
  else{ UseMethod("edges") }
}

edges.hashgraph = function(x){
    # Returns a dataframe of edge1, edge2, value.
    # Spins through edges one time to count number of pairs; this is required for large 
    # graphs because dynamically growing vectors is terribly slow
    keys = hsh_keys( x )
    if(length(keys)==0 ){
        return( data.frame() )
    }
    ctr=0
    for(i in 1:length(keys) ){
        k.hsh = hsh_get(x, keys[i] )
        i.keys = hsh_keys(k.hsh)
        for(j in 1:length(i.keys ) ){
            ee1 = keys[i]
            ee2 = i.keys[j]
            if( ee1<ee2 ){
                ctr = ctr + 1
            }
        }
    }
    e1 = rep(0, ctr)
    e2 = rep(0, ctr)
    values = rep(1, ctr)
    ctr = 1
    for(i in 1:length(keys) ){
        k.hsh = hsh_get(x, keys[i] )
        i.keys = hsh_keys(k.hsh)
        for(j in 1:length(i.keys ) ){
            ee1 = keys[i]
            ee2 = i.keys[j]
            vv = hsh_get( k.hsh, ee2 )
            if( ee1<ee2 ){
                e1[ctr] = ee1
                e2[ctr] = ee2
                values[ctr] = vv
                ctr = ctr + 1
            }
        }
    }
    data.frame(e1, e2, values, stringsAsFactors=F)  
}

neighbors <- function(x, node)  {
  if(is.null(attr(x, "class"))){ stop("Must be called on a class") }
  else{ UseMethod("neighbors") }
}

neighbors.hashgraph = function(x, nodes){
    neighbors = c()
    for(i in 1:length(nodes)){
        node = nodes[i]
        if( hsh_in( x, node ) ){
            k.hsh = hsh_get(x, node )
            neighbors = c(neighbors, hsh_keys( k.hsh ) )
        }
    }
    unique(neighbors)
}

degree <- function(x, nlist, do.sort)  {
  if(is.null(attr(x, "class"))){ stop("Must be called on a class") }
  else{ UseMethod("degree") }
}

degree.hashgraph = function(x, nlist=NULL, do.sort=TRUE){
    if( is.null(nlist) ){
        nlist = nodes(x)
    }
    deg = rep(0, length(nlist) )
    for(i in 1:length(nlist) ){
        deg[i] = length( neighbors( x, nlist[i] ) )
    }
    dl = data.frame( nlist, deg, stringsAsFactors=FALSE )
    if( do.sort )
        dl = dl[order(dl$deg, decreasing=T),]
    names(dl) = c('node', 'degree')
    dl
}

length.hashgraph = function(x){
    length(nodes(x))
}

subgraph <- function(x, nodes)  {
  if(is.null(attr(x, "class"))){ stop("Must be called on a class") }
  else{ UseMethod("subgraph") }
}

subgraph.hashgraph = function(x, nodes){
    G.new = hashgraph()
    if( length(nodes) > 0 ){
        node.subset = hsh_from_vectors(nodes, rep(1, length(nodes) ) )
        for(i in 1:length( nodes ) ){
            if( !hsh_in(x, nodes[i] ) )
                stop( paste( "Node", nodes[i], "not found in graph") )
            neighbors = hsh_keys( hsh_get(x, nodes[i] ) )
            for( j in 1:length(neighbors) ){
                if( hsh_in( node.subset, neighbors[j] ) ){
                    addEdge( G.new, nodes[i], neighbors[j] )
                }
            }
        }
    }
    G.new
}

triangles <- function(x)  {
  if(is.null(attr(x, "class"))){ stop("Must be called on a class") }
  else{ UseMethod("triangles") }
}

triangles.hashgraph = function(x){
    all = nodes(x)
    keep = rep(F, length(all) )
    for(i in 1:length(all)){
        if(!keep[i]){
            N = neighbors(x, all[i])
            sub.edges = edges( subgraph(x, N) )
            if( dim(sub.edges)[1]>0 ){
                keep[i] = T
                for(j in 1:dim(sub.edges)[1]){
                    keep[ which(all==sub.edges[j,1] ) ] = T
                    keep[ which(all==sub.edges[j,2] ) ] = T                    
                }
            }
        }
    }
    subgraph(x, all[keep])
}

intersection <- function(x, y)  {
  if(is.null(attr(x, "class"))){ stop("Must be called on a class") }
  else{ UseMethod("intersection") }
}

intersection = function(x, y){
    # return new graph where edges found in both x and y
    e.x = edges(x)
    G = hashgraph()
    for(i in 1:dim(e.x)[1] ){
        if( hasEdge( y, e.x[i,1], e.x[i,2] ) ){
            addEdge( G, e.x[i,1], e.x[i,2] )
        }
    }
    G
}

difference <- function(x, y)  {
  if(is.null(attr(x, "class"))){ stop("Must be called on a class") }
  else{ UseMethod("difference") }
}

difference = function(x, y){
    # return new graph where edges found in x not y
    e.x = edges(x)
    G = hashgraph()
    for(i in 1:dim(e.x)[1] ){
        if( !hasEdge( y, e.x[i,1], e.x[i,2] ) ){
            addEdge( G, e.x[i,1], e.x[i,2], getValue( x, e.x[i,1], e.x[i,2] ) )
        }
    }
    G
}

restrict3cliques = function(x){
    if(is.null(attr(x, "class"))){ stop("Must be called on a class") }
    else{ UseMethod("restrict3cliques") }
}

restrict3cliques.hashgraph = function(G){
    in_clique = hsh_new()
    G.nodes = nodes(G)
    for( i in 1:length(G.nodes)){
        g=G.nodes[i]
        if(!hsh_in(in_clique, g)){
            neighbors = neighbors( G, g)
            N = length(neighbors)
            if( N>=2 ){
                for( i in 1:N ){
                    for( j in (i+1):N ){
                        if( hasEdge(G, neighbors[i], neighbors[j]  ) ){
                            hsh_set( in_clique, g, 1 )
                            hsh_set( in_clique, neighbors[i], 1 )
                            hsh_set( in_clique, neighbors[i], 1 )
                            break
                        }
                    }
                }
            }
        }
    }
    subgraph(G, hsh_keys( in_clique ) )
}

intersect.hashgraph = function(x, y){
    in.both = intersect( nodes(x), nodes(y) )
    subgraph( G, in.both )
}

minimum.degree = function(G, min.n){
    if(is.null(attr(x, "class"))){ stop("Must be called on a class") }
    else{ UseMethod("minimum.degree") }
}

minimum.degree.hashgraph = function(G, min.n){
    is.done=F
    count = 1
    while(! is.done ){
        print(paste("Iteration",count,"N =",length(G)))
        count = count+1
        deg = degree(G, do.sort=F)
        if(min(deg$degree) >= min.n){
            is.done=T
        }
        else{
            G = subgraph( G, deg$node[deg$degree>=min.n] )
        }
        if( length(G)==0 ){
            is.done=T
        }
    }
    G
}

#######################
### BEGIN UTILITIES ###
#######################

set.difference = function(A, B){
    # returns the set difference set(A) - set(B), not vector of T/F
    a.not.b = vector()
    for( i in 1:length(A) ){
        idx = which(B==A[i])
        if( length(idx)==0 )
            a.not.b[ length(a.not.b)+1 ] = A[i]
    }
    a.not.b
}

set.intersection = function(A, B){
   # returns the set intersection set(A, B), not vector of T/F
    a.and.b = vector()
    for( i in 1: length(A) ){
        if( length( which(B==A[i]) >0 ) )
            a.and.b[ length(a.and.b)+1 ] = A[i]
    }
    a.and.b
}

set.union = function(A, B){
    a.or.b = vector()
    ab = append(A,B)
    h = hsh_from_vectors( ab, rep(1, length(ab) ) )
    hsh_keys(h)
}


count.appearances=function(V, order.by="values"){
    h = hsh_new()
    for(i in 1:length(V)){
        if( hsh_in(h, V[i] ) ){
            cnt = hsh_get(h, V[i])
            hsh_set(h, V[i], cnt+1)
        }
        else{
            hsh_set(h, V[i], 1)
        }
    }
    kv = hsh_keys_values(h)
    if(order.by=="values")
        kv[order(kv$values, decreasing=T),]
    else
        kv[order(kv$keys),]
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

standardize = function(D){
    D = data.matrix(D)
    library(matrixStats)
    sds = rep(0, dim(D)[1])
    for(i in 1:length(sds)){
        sds[i] = sd(D[i,], na.rm=T)
    }
    (D-rowMeans(D, na.rm=T))/sds
}

compress.probes = function( D, idx, min.cor=0.8 ){
    # Given a set of probe indexes idx, look for those with correlation >= min.cor
    # Report mean values across all probes with correlation >= min.cor
    # If no genes pairs meet these criteria, default to mean value
    ee = D[ idx, ]
    if( length(idx)==1 ){
        ee
    }
    else{
        ee.strong = cor(t(ee), use="complete") >= min.cor
        keep = rep(F, length(idx) )
        for(row in 1:(dim(ee)[1]-1) ){
            for(col in (row+1):dim(ee)[1] ){
                if(ee.strong[row,col]){
                    keep[row]=T
                    keep[col]=T
                }
            }
        }
        if(sum(keep)==0){
           # No two probes with sufficient correlation; default to mean of all
           keep[!keep]=T
        }
        colMeans(ee[keep,], na.rm=T)
    }
}

match.idx = function(A, B){
    # return dataframe of indices into A and B restricted to perfect matches
    # between A and B, where idx.A[i] == idx.B[i] for each i in matched pairs
    in.both = intersect(A,B)
    idx.A = match(in.both, A)
    idx.B = match(in.both, B)
    C= data.frame(idx.A, idx.B)
    if( sum( A[ C$idx.A ] != B[ C$idx.B] )>0 )
        stop("ERROR! At least one in idx.A not the same as matched item in idx.B")
    C
}

match.idx.first = function( A, B ){
    # for intersection, returns first index of B in first index of A
    in.both = intersect(A, B)
    n.intersect = length(in.both)
    idx.A = rep(0, n.intersect )
    idx.B = rep(0, n.intersect )
    for(i in 1:n.intersect ){
        idx.A[i] = min( which(A==in.both[i]) )
        idx.B[i] = min( which(B==in.both[i]) )    
    }
    data.frame(idx.A, idx.B)
}

mean.by.symbol = function( expr, all.symbols, symbols ){
    # find all occurences of each element in symbols in all.symbols; return mean of those values in expr
    # Error if element in symbols is not found.
    means = matrix(0, nrow=length(symbols), ncol=dim(expr)[2] )
    D = data.matrix(expr)
    for(i in 1:length(symbols)){
        idx = which(all.symbols==symbols[i])
        if(length(idx)==0 ){
           stop(paste("No match found for requested symbol", symbols[i], "in all.symbols" ) )
        }
        else if(length(idx)==1){
            means[i,] = D[idx,] 
        }
        else{
            means[i,] = colMeans( D[idx,], na.rm=T )
        }
    }
    df=data.frame( means )
    names(df)=names(expr)
    rownames(df) = symbols
    df
}
    
contains = function(A, B){
    # is there a value in B for each value in A?
    # returns logical vector
    N = length(A)
    if( length(A)==0 ){ stop("Length of A is 0 in utility function 'contains'") }
    if( length(B)==0 ){ stop("Length of B is 0 in utility function 'contains'") }
    tf = vector(mod="logical", N)
    for( i in 1:length(A) ){
        if( A[i] %in% B )
            tf[i] = T
    }
    tf
}

make_list = function(keys, vals){
    N_K = length(keys)
    N_V = length(vals)
    if( N_K != N_V ){
        stop("length of keys and vals is not the same")
    }
    h = list()
    for(i in 1:N_K){
        h[[ as.character(keys[i]) ]] = as.character(vals[i])
    }
    h
}

abi.chr = function(snp, is.mouse=T){
    if( is.vector( snp ) ){
        chrs = rep(0, length(snp) )
        for(i in 1:length(snp) ){
            snpi = as.character(snp[i])
            s = strsplit(snpi,'.',fixed=T)[[1]][1]
            s = strsplit(s,'E',)[[1]][2]
            chrs[i] = as.numeric(s)
        }
        chrs
    }
    else{
        snp = as.character(snp)
        s = strsplit(snp,'.',fixed=T)[[1]][1]
        s = strsplit(s,'E',)[[1]][2]
        as.numeric(s)
    }
}

abi.mb = function(snp, is.mouse=T){
    snp = as.character(snp)
    s = strsplit(snp,'.',fixed=T)[[1]][2]
    as.numeric(s)
}


nearest.ten = function( x ){
    # find the nearest factor of 10 for x between 10 and 100, rounding up
    rv = 100
    for(i in seq(from=10,to=100,by=10) ){
        if( x<i ){
            rv=i
            break
        }
    }
    rv
}

prepend.label = function(list.of.targets, label){
    # prepend label to a vector of characters
    for(i in 1:length(list.of.targets)){
        list.of.targets[i] = paste(label, list.of.targets[i], sep='_')
    }
    list.of.targets
}

convert.to.ranks = function(M){
    # M is a data frame with real-valued data.  We will convert to ranks by columns
    b = M
    for(i in 1:dim(M)[2]){
        b[,i] = rank(as.numeric(b[,i]),na.last="keep")
    }
    b
}

probes.for.gene = function(ga, probes){
    ids.with.attribute(ga,probes)
}

ids.with.attribute = function(ga, col.value, col.name="Gene Name"){
    col.idx = which( names(ga)==col.name)
    if(length(col.idx)==0)
        stop(paste("Column", col.name, "not present in data frame."))
    rownames(ga)[ which(ga[,col.idx]==col.value) ]
}

probe.paired.data = function(ids, expr.A, expr.B, probe.id){
    # reports probe.id values for matched subjects in expr.A, expr.B.  Returns a data frame.
    N = dim(ids)[1]
    rowA = which( rownames(expr.A) == probe.id )
    rowB = which( rownames(expr.B) == probe.id )
    matches = match.indices(ids, expr.A, expr.B)
    idx.A = matches$A
    idx.B = matches$B
    vals.A = as.numeric( expr.A[ rowA, idx.A] )
    vals.B = as.numeric( expr.B[ rowB, idx.B] )
    ab= data.frame(ids[,1], vals.A, ids[,2], vals.B)
    sorted.idx = order(ab[,2])
    ab[sorted.idx,]
}

rows.at.threshold = function(A, percent.present=0.9){
    # given data frame A and percent.present between 0 and 1, return vector of rows
    # where A meets threshold for perecnt.present
    N.probes = dim(A)[1]
    N = dim(A)[2]
    valid = vector(mode="logical", N.probes)
    threshold = floor( N * percent.present )
    for(i in 1:N.probes){
        vals.A = as.numeric( A[ i,] )
        n.present.A = N - sum( is.na(vals.A) )
        valid[i] = n.present.A >= threshold
        if( i %% 500==0 ){
            print(paste("Probe", i, "of", N.probes) )
        }
    }
    valid
}


paired.rows.at.threshold = function(ids, expr.A, expr.B, percent.present=0.9){
    # given data frames expr.A and expr.B and threshold between 0 and 1
    # return vector of rows where both A and B meet threshold for percent present
    if( dim(expr.A)[1] != dim(expr.B)[1] ){ stop( "expr.A does not have the same number of rows as expr.B") }
    N.probes = dim(expr.A)[1]
    N = dim(ids)[1]
    valid = vector(mode="logical", N.probes)
    matches = match.indices(ids, expr.A, expr.B)
    idx.A = matches$A
    idx.B = matches$B

    threshold = floor( N * percent.present )
    a = expr.A[,idx.A]
    b = expr.B[,idx.B]
    for(i in 1:N.probes){
        vals.A = as.numeric( a[ i,] )
        vals.B = as.numeric( b[ i,] )
        n.present.A = N - sum( is.na(vals.A) )
        n.present.B = N - sum( is.na(vals.B) )
        valid[i] = (n.present.A >= threshold && n.present.B >= threshold)
        if( i %% 500==0 ){
            print(paste("Probe", i, "of", N.probes) )
        }
    }
    valid
}

numeric.values = function(expr, samples, row.id){
    # extracts value for each sample at row, returns numeric vector
	v = vector()
	idx = which(rownames(expr)==row.id)
	for( i in 1:length(samples) ){
		v[length(v)+1] = expr[idx, which( names(expr)==samples[i] )]
	}
	as.numeric(v)
}

smooth = function(x, chrom, smooth.region = 2, outlier.SD.scale = 4, smooth.SD.scale = 2, trim=0.025){
    # Adapted directly from Adam Olshen's DNAcopy software
    nsample = dim(x)[2]
    uchrom = unique(chrom)
    for (isamp in 1:nsample) {
        genomdat <- x[,isamp]
        ina <- which(!is.na(genomdat) & !(abs(genomdat) == Inf))
        X = sort(genomdat[ina])
        n = floor(trim*length(X)) + 1
        trimmed.SD <- sqrt( var( X[n : (length(X) - n)] ) )
        outlier.SD <- outlier.SD.scale * trimmed.SD
        smooth.SD <- smooth.SD.scale * trimmed.SD
        k <- smooth.region
        for (i in uchrom) {
            ina <- which(!is.na(genomdat) & !(abs(genomdat) == Inf) & chrom == i)
            n <- length(genomdat[ina])
            smoothed.data <- sapply(1:n, function(i, x, n, nbhd, oSD, sSD) {
                xi <- x[i]
                nbhd <- i + nbhd
                xnbhd <- x[nbhd[nbhd > 0 & nbhd <= n]]
                if (xi > max(xnbhd) + oSD) 
                  xi <- median(c(xi, xnbhd)) + sSD
                if (xi < min(xnbhd) - oSD) 
                  xi <- median(c(xi, xnbhd)) - sSD
                xi
            }, genomdat[ina], n, c(-k:-1, 1:k), outlier.SD, smooth.SD)
            x[, isamp][ina] <- smoothed.data
        }
    }
    x
}

# REAL USE CASES

do.km = function(){
    sa.km = load.matrix('c:/code/microTools/src/R/km_input.txt')
    km(sa.km$time, sa.km$had.event, sa.km$condition)
}

amplification.vs.expression.lung.tumors = function(){
    sa.g.AMP = load.matrix('C:/code/lung/experiments/snps/log2_results/sample_attributes.txt')
    sa.e.LT = load.matrix('c:/code/lung/data/expression/matched_tumors/sample_attributes.txt')
    expr.g.AMP = load.matrix('C:/code/lung/experiments/snps/log2_results/expr_log2.txt')
    expr.e.LT = load.matrix('c:/code/lung/data/expression/matched_tumors/expr.txt')
    ids.AMP.LT = match.identifers( sa.g.AMP, sa.e.LT, sa.g.AMP$expression_id, sa.g.AMP$valid==1, sa.e.LT$valid==1)
    plot.geno.expr(ids.AMP.LT, expr.g.AMP, expr.e.LT, "rs4947952", "ILMN_1798975")
}

genotype.vs.expression.lung.normal = function(){
    sa.g.MAF_5 = load.matrix('C:/code/lung/data/genotypes/MAF_5/sample_attributes.txt')
    sa.e.LN = load.matrix('c:/code/lung/data/expression/normal/sample_attributes.txt')
    expr.g.MAF_5 = load.matrix('C:/code/lung/data/genotypes/MAF_5/expr.txt')
    expr.e.LN = load.matrix('c:/code/lung/data/expression/normal/expr.txt')
    ids.MAF_5.LN = match.identifers( sa.g.MAF_5, sa.e.LN, sa.g.MAF_5$shared_snp_expr, sa.g.MAF_5$is.valid==1 & sa.g.MAF_5$is.tumor==0, sa.e.LN$valid==1)
    plot.geno.expr(ids.MAF_5.LN, expr.g.MAF_5, expr.e.LN, "rs4646491", "ILMN_1659215")
}

plot.expression = function(){
    plot.rows( expr.e.LN[1,])
    plot.rows( expr.e.LN[1:3,], sa.e.LN, sort.within.group=T)
    plot.rows( expr.e.LN[1:3,], sa.e.LN, group.by='mutant.p53', sort.within.group=T)
}

plot.genomic.alteration.frequency = function(){
    expr.g.AMP = load.matrix('C:/code/lung/experiments/snps/log2_results/expr_log2.txt')
    ga.g.AMP = load.matrix('C:/code/lung/data/snps/gene_attributes.txt')
    plot.percent.altered.by.locus(expr.g.AMP, ga.g.AMP$Chr, upper.bound=0.3, lower.bound=-0.3)
    d=plot.percent.altered.by.tumor(expr.g.AMP, upper.bound=0.3, lower.bound=-0.3)
}

plot.discretization = function(){
    plot.discretization(expr.e.LN[1,], bounds=1)
}


correlation.tests = function(){
    expr.g.AMP = load.matrix('C:/code/lung/experiments/snps/log2_results/expr_log2.txt')
    expr.e.LT = load.matrix('c:/code/lung/data/expression/matched_tumors/expr.txt')
    sa.g.AMP = load.matrix('C:/code/lung/experiments/snps/log2_results/sample_attributes.txt')
    sa.e.LT = load.matrix('c:/code/lung/data/expression/matched_tumors/sample_attributes.txt')
    ga.g.LT = load.matrix('c:/code/lung/data/expression/matched_tumors/gene_attributes_stripped.txt')
    ids.g.e=match.identifers( sa.g.AMP, sa.e.LT, sa.g.AMP$expression_id, sa.g.AMP$valid==1, sa.e.LT$valid==1)
    cor.test.snp.genes(ids.g.e, expr.g.AMP, expr.e.LT, ga.g.LT,  "rs4646491")


    expr.miRNA = load.matrix('c:/code/mouse/miRNA/data/expr.txt')
    expr.tail = load.matrix('c:/code/mouse/tails/data/expr_above_bg.txt')
    ga.tail = load.matrix('c:/code/mouse/tails/data/gene_attributes_above_bg.txt')
    sa.miRNA = load.matrix('c:/code/mouse/miRNA/data/sample_attributes.txt')
    sa.tail = load.matrix('c:/code/mouse/tails/data/sample_attributes.txt')
    ids.mir.tail = match.identifers( sa.miRNA, sa.tail, sa.miRNA$tail_id, sa.miRNA$tail_id!='NA', sa.tail$valid==1)
    d=cor.test.row.vs.genes(ids.mir.tail, expr.miRNA, expr.tail, ga.tail, "hmr-miR-200c", percent.present=0.9, method="pearson", verbose=T)


    plot.geno.expr(ids.AMP.LT, expr.g.AMP, expr.e.LT, "rs13234448", "ILMN_1793118")


    d=cor.test.row.vs.genes(ids.mir.tail, expr.miRNA, expr.tail, ga.tail[,c(1,6,7)], probe.id="hmr-miR-18a", method="spearman")
}

unit.test = function(){
    # TEST CODE
    sa.P = data.frame( c(1,1,0,1,1))
    rownames(sa.P) = c('p1', 'p2', 'p3', 'p4', 'p5')
    names(sa.P) = c('valid')
    sa.G = data.frame( c('p1', 'p2', 'p3', 'p4', 'p5'), c(0,1,1,1,1))
    rownames(sa.G) = c('g4', 'g2', 'g3', 'g1', 'g5')
    names(sa.G) = c('translate', 'valid')
    
    P = data.frame( c(1), c(2), c(3), c(4), c(5) )
    names(P) = c('p1', 'p2', 'p3', 'p4', 'p5')
    rownames(P) = c('pheno1')
    
    G = data.frame( c(0), c(1), c(2), c(1), c(NA) )
    names(G) = c('g1', 'g2', 'g3', 'g4', 'g5')
    rownames(G) = c('geno1')
    ids = match.identifers( sa.G, sa.P, sa.G$translate, sa.G$valid==1, sa.P$valid==1)
    g.p.vals = align.geno.pheno( ids, G, P, 'geno1', 'pheno1')
    h = make_list(rownames(sa.G), sa.G$translate)
    
    A = data.frame(c(1,5), c(2,4), c(3,3), c(4,2) , c(5,1) )
    rownames(A) = c('g1','g2')
    names(A) = c('s1','s2','s3','s4','s5')
    
    B = data.frame(c(2,5), c(6,4), c(8,3), c(7,2) , c(4,1) )
    rownames(B) = c('x1','x2')
    names(B) = c('s1','s2','s3','s4','s5')
    ids = data.frame( c('s1','s2','s3','s4','s5'), c('s1','s2','s3','s4','s5'))
    names(ids) = c('genotype.id', 'phenotype.id')
    ga.genes = data.frame( c('gene1', 'gene2'), c('1', '2'))
    names(ga.genes) = c('Gene Name', 'Chr')
    cor.test.row.vs.all( ids, B[1,], A, percent.present=0.90, method="spearman", verbose=F)
    cor.test.row.vs.genes(ids, B, A, ga.genes, 'x1', percent.present=0.9, method="pearson", verbose=T)
    
    eN = data.frame( c(1),c(2),c(3),c(4),c(5) )
    names(eN) = c('s1', 's2', 's3', 's4', 's5')
    rownames(eN) = c('ID1')
    eT = data.frame( c(3),c(4),c(5),c(6),c(9) )
    names(eT) = c('s1', 's2', 's3', 's4', 's5')
    rownames(eT) = c('ID1')
    ids.test = data.frame(c('s1','s2','s3','s4','s5'), c('s1','s2','s3','s4','s5'), c(1,2,3,4,5))
    names(ids.test) = c('A.id', 'B.id', 'shared')
    paired.t.test(ids.test, eN, eT, percent.present=0.9, verbose=T)
    
    a = data.frame( c(1,2,3,4,5), c(5,4,3,2,1), c(2,4,5,1,3) )
    convert.to.ranks(a)
    
    a = rep(0,100)
    a[51:100]=1
    A = data.frame(a)
    rownames(A) = c(1:100)
    B = t(A)
    B[1:50] = ((1:50)/10)+1
    B[51:100] = ((51:100)/10)-1
    names(B) = c(1:100)
    logistic(A,B,do.plot=T)

}
