source('/notebook/code/src/R/quantitative_genetics.R')
args <- commandArgs(trailingOnly = TRUE)
if( length(args) != 2 )
    stop("Two arguments are required: dir_segment_beds, fn_out")
dir_seg = args[1]
fn_out = args[2]

dir_seg = '/datasets/human_lines_CCLE/CN/segment_beds'
files=list.files(dir_seg, pattern='.*BED')
cell_line = c()
chrom = c()
loc_start = c()
loc_end = c()
minor = c()
major = c()
for(i in 1:length(files)){
    print(paste( i, files[i] ))
    bed = read.table(paste(dir_seg, files[i], sep='/'), stringsAsFactors=FALSE, header=FALSE, sep='\t')
    cell_line = c(cell_line, rep( strsplit( files[i], "___", fixed=T)[[1]][1], dim(bed)[1]))
    chrom = c(chrom, bed$V1) 
    loc_start = c(loc_start, bed$V2)
    loc_end = c(loc_end, bed$V3)
    major = c( major, get.split.col( bed$V4, "major:", last=TRUE) )
    minor = c( minor, get.split.col( get.split.col( bed$V4, "minor:", last=TRUE), ",", first=TRUE) )
}
beds = data.frame( cell_line, chrom, loc_start, loc_end, minor, major, stringsAsFactors=FALSE)
write.table(beds, fn_out, row.names=FALSE, quote=FALSE)