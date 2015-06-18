#' equalizer: A package for removing SNP effects from Affymetrix microarrays
#' @docType package
#' @name equalizer
#' 

check_file = function( fn ){
    if( !file.exists( fn) ){
        stop( paste( fn, "parameter not found, passed", fn ) )   
    }
    TRUE
}


#' Run equalizer on IVT-format arrays. 
#' 
#'   \code{equalize_IVT} uses one or more VCF files defining the locations of 
#'   Single Nucleotide Polymorphisms (SNPs) to rewrite Affymetrix array description files 
#'   for a source package so they exclude any probes which overlap SNPs. 
#'   The \code{equalize_IVT} function is a wrapper around a Python script.
#'   
#'   To run this function, \strong{you must have Python and bedtools installed locally}.
#'   
#'   The files where (sourced from Netaffx) is listed are available 
#'   from Affymetrix's Netaffx website. You may also be able to find copies at 
#'   davidquigley.com . 
#'   
#' @param fn_python Path to python executable
#' @param dir_bedtools_bin Path to bin folder of bedtools
#' @param package_name Name of the new package to create
#' @param fn_vcf_list Vector of paths to VCF files
#' @param affy_bed path to Affymetrix probe locations in BED format (sourced from Netaffx)
#' @param annot_version Affymetrix annotation version for the source chip
#' @param annot_revision Affymetrix chip library revision number for the source chip
#' @param genome_version_UCSC UCSC Genome version for the source chip
#' @param genome_version_NCBI NCI Genome version for the source chip
#' @param chip Affymetrix chip identifier for the source chip
#' @param chip_version Affymetrix chip version string for the source chip
#' @param species Scientific species name for the source chip
#' @param organism Common species name for the source chip
#' @param author Author name for new package
#' @param email Author email for new package
#' @param fn_CDF Path to Affymetrix CDF file (sourced from Netaffx)
#' @param fn_probe_tab Path to Affymetrix probe_tab file (sourced from Netaffx)
#' @param fn_CEL Path to an Affymetrix CEL file of the type produced by the source chip
#' @param fn_probeset_csv Path to affymetrix probeset.csv (sourced from Netaffx)
#' @param dir_out Folder where new package files will be created
#' @seealso \code{\link{equalize_gene}} for ST-format arrays (e.g. gene ST 1.1) 
equalize_IVT = function(
        fn_python, dir_bedtools_bin,
        package_name, fn_vcf_list, affy_bed,
        annot_version, annot_revision,
        genome_version_UCSC, genome_version_NCBI,
        chip, chip_version, species, organism,
        author, email,
        fn_CDF, fn_probe_tab, fn_CEL, fn_probeset_csv,
        dir_out){
    # This function is a wrapper for equalizer.py.
    # It is used for IVT-format arrays (e.g. HGU133A, M430)
    # Confirm all required files exist, build arguments, run script 
    check_file( fn_python )
    check_file( dir_bedtools_bin )
    if( length(fn_vcf_list)==0 )
        stop("must pass at least one VCF file in fn_vcf_list")
    for( i in 1:length(fn_vcf_list) )
        check_file( fn_vcf_list[i] )
    check_file( affy_bed )
    check_file( fn_CDF )
    check_file( fn_probe_tab )
    check_file( fn_CEL )
    check_file( fn_probeset_csv )   
    check_file( dir_out )
    
    args = c(fn_equalizer_script, "-f IVT", collapse=" ")
    args = c(args, paste("-p", package_name,        collapse=" ") )
    args = c(args, paste("-v", paste(fn_vcf_list, collapse=",") , collapse=" ") )
    args = c(args, paste("-a", affy_bed,            collapse=" ") )
    args = c(args, paste("-b", dir_bedtools_bin,    collapse=" ") )
    args = c(args, paste("-s", annot_version,       collapse=" ") )
    args = c(args, paste("-e", annot_revision,      collapse=" ") )
    args = c(args, paste("-g", genome_version_UCSC, collapse=" ") )
    args = c(args, paste("-d", genome_version_NCBI, collapse=" ") )
    args = c(args, paste("-c", chip,                collapse=" ") )
    args = c(args, paste("-i", chip_version,        collapse=" ") )
    args = c(args, paste("-r", species,             collapse=" ") )
    args = c(args, paste("-w", organism,            collapse=" ") )
    args = c(args, paste("-u", author,              collapse=" ") )
    args = c(args, paste("-l", email,               collapse=" ") )
    args = c(args, paste("-n", fn_CDF,              collapse=" ") )
    args = c(args, paste("-x", fn_probe_tab,        collapse=" ") )
    args = c(args, paste("-c", chip,                collapse=" ") )
    args = c(args, paste("-y", fn_probeset_csv,     collapse=" ") )
    args = c(args, paste("-z", fn_CEL,              collapse=" ") )
    args = c(args, paste("-o", dir_out,             collapse=" ") )
    
    system2( fn_python, args )
}

#' Run equalizer on ST-format arrays. 
#' 
#'   \code{equalize_gene} uses one or more VCF files defining the locations of 
#'   Single Nucleotide Polymorphisms (SNPs) to rewrite Affymetrix array description files 
#'   for a source package so they exclude any probes which overlap SNPs. 
#'   The \code{equalize_gene} function is a wrapper around a Python script.
#'   
#'   To run this function, \strong{you must have Python and bedtools installed locally}.
#'   
#'   The files Affymetrix where (sourced from Netaffx) is listed are available 
#'   from Affymetrix's Netaffx website. You may also be able to find copies at 
#'   davidquigley.com . 
#'   
#' @param fn_python Path to python executable
#' @param dir_bedtools_bin Path to bin folder of bedtools
#' @param package_name Name of the new package to create
#' @param fn_vcf_list Vector of paths to VCF files
#' @param affy_bed path to Affymetrix probe locations in BED format (sourced from Netaffx)
#' @param annot_version Affymetrix annotation version for the source chip
#' @param annot_revision Affymetrix chip library revision number for the source chip
#' @param genome_version_UCSC UCSC Genome version for the source chip
#' @param genome_version_NCBI NCI Genome version for the source chip
#' @param chip Affymetrix chip identifier for the source chip
#' @param chip_version Affymetrix chip version string for the source chip
#' @param species Scientific species name for the source chip
#' @param organism Common species name for the source chip
#' @param author Author name for new package
#' @param email Author email for new package
#' @param fn_probeset_csv Path to affymetrix probeset.csv (sourced from Netaffx)
#' @param fn_transcript_csv Path to affymetrix probeset.csv (sourced from Netaffx)
#' @param fn_pgf Path to affymetrix PGF file (sourced from Netaffx)
#' @param fn_mps Path to affymetrix MPS file (sourced from Netaffx)
#' @param fn_clf Path to affymetrix CLF file (sourced from Netaffx)
#' @param dir_out Folder where new package files will be created
#' 
#' 
#' @return none
equalize_gene = function(
        fn_python, dir_bedtools_bin,
        package_name, fn_vcf_list, affy_bed,
        annot_version, annot_revision,
        genome_version_UCSC, genome_version_NCBI,
        chip, chip_version,
        species, organism,
        author, email,
        fn_probeset_csv, fn_transcript_csv, fn_pgf, fn_mps, fn_clf,
        dir_out){
    # This function is a wrapper for equalizer.py.
    # It is used for gene and exon-format arrays (e.g. gene ST 1.0)
    # Confirm all required files exist, build arguments, run script 
    check_file( fn_python )
    check_file( dir_bedtools_bin )
    if( length(fn_vcf_list)==0 )
        stop("must pass at least one VCF file in fn_vcf_list")
    for( i in 1:length(fn_vcf_list) )
        check_file( fn_vcf_list[i] )
    check_file( affy_bed )
    check_file( fn_probeset_csv )
    check_file( fn_transcript_csv )
    check_file( fn_pgf )
    check_file( fn_mps )
    check_file( fn_clf )
    check_file( dir_out )
    fn_equalizer_script = paste(installed.packages()["equalizer","LibPath"], 
                                "/equalizer/exec/equalizer.py", collapse="", 
                                sep="")
    args = c(fn_equalizer_script, "-f gene", collapse=" ")
    args = c(args, paste("-p", package_name,        collapse=" ") )
    args = c(args, paste("-v", paste(fn_vcf_list, collapse=",") , collapse=" ") )
    args = c(args, paste("-a", affy_bed,            collapse=" ") )
    args = c(args, paste("-b", dir_bedtools_bin,    collapse=" ") )
    args = c(args, paste("-s", annot_version,       collapse=" ") )
    args = c(args, paste("-e", annot_revision,      collapse=" ") )
    args = c(args, paste("-g", genome_version_UCSC, collapse=" ") )
    args = c(args, paste("-d", genome_version_NCBI, collapse=" ") )
    args = c(args, paste("-c", chip,                collapse=" ") )
    args = c(args, paste("-i", chip_version,        collapse=" ") )
    args = c(args, paste("-r", species,             collapse=" ") )
    args = c(args, paste("-w", organism,            collapse=" ") )
    args = c(args, paste("-u", author,              collapse=" ") )
    args = c(args, paste("-l", email,               collapse=" ") )
    args = c(args, paste("-c", chip,                collapse=" ") )
    args = c(args, paste("-y", fn_probeset_csv,     collapse=" ") )
    args = c(args, paste("-t", fn_transcript_csv,   collapse=" ") )
    args = c(args, paste("-q", fn_pgf,              collapse=" ") )
    args = c(args, paste("-m", fn_mps,              collapse=" ") )
    args = c(args, paste("-k", fn_clf,              collapse=" ") )
    args = c(args, paste("-o", dir_out,             collapse=" ") )
    
    system2( fn_python, args )
}
