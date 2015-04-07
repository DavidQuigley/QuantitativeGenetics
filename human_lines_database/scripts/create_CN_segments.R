# new code to generate fundamental measure for HCLDB
library(foreach)
library(doParallel)
source('/notebook/code/src/R/quantitative_genetics.R')

fn.segments = '/datasets/human_lines_database/copy_number/CCLE_segments.txt'
segs = read.table( fn.segments, header=TRUE, sep=' ', stringsAsFactors=FALSE)
line_names = sort(unique(segs$cell_line))
n_lines = length(line_names)

chr_lengths = rep(0,24)
for(i in 1:24){
    chr_lengths[i] = max(segs$loc_end[segs$cell_line==line_names[1] & segs$chrom==i])
}

bands = read.table('/datasets/human_lines_database/annotations/cytobands_HG19', header=FALSE, stringsAsFactors=FALSE)
names(bands) = c('chrom', 'loc_start', 'loc_end', 'cytoband', 'type')
bands = bands[bands$type=="acen",]
bands$chrom = get.split.col(bands$chrom, "chr", last=TRUE)
bands$chrom[bands$chrom=="X"] = 23
bands$chrom[bands$chrom=="Y"] = 24
bands$chrom=as.numeric(bands$chrom)
bands = bands[order(bands$chrom, bands$loc_start),]
centromere_start = rep(0, 24)
centromere_end = rep(0, 24)
for(chrom in 1:24){
    centromere_start[chrom] = min(bands$loc_start[bands$chrom==chrom] )
    centromere_end[chrom] = max(bands$loc_end[bands$chrom==chrom] )
}
centromeres = data.frame(start=centromere_start, end=centromere_end)    

N_CORES=8

n_states_by_chr = matrix(NA, nrow=n_lines, ncol=24)
weighted_density_states_by_chr = matrix(NA, nrow=n_lines, ncol=24)
whole_arm_alt_by_chr = matrix(NA, nrow=n_lines, ncol=48)
focal_alt_by_chr = matrix(NA, nrow=n_lines, ncol=24)

cl = makeForkCluster(nnodes=N_CORES)
registerDoParallel(cl, cores=N_CORES)

foreach( idx=1:n_lines ) %dopar% {
#for(idx in 1:n_lines){
    print(idx)
    S = segs[segs$cell_line==line_names[idx],]
    loc_center = S$loc_start + ( (S$loc_end - S$loc_start) / 2 )
    seg_length = S$loc_end - S$loc_start
    for( i in 1:24 ){
        is_chrom = S$chrom==i
        loc_center_chr = loc_center[is_chrom]
        n_states_by_chr[idx,i] = length( loc_center_chr )
        weighted_density_states_by_chr[idx, i ] = sum( 1/dist( loc_center_chr ) )
        focal_alt_by_chr[idx,i] = sum( S$major[is_chrom] != 2 &
                                       seg_length[is_chrom] < 20000000 ) 
        whole_arm_alt_by_chr[idx, (i*2)-1] = sum( S$major!= 2 & S$chrom==i & 
                                            S$loc_start < 1000000 & 
                                            S$loc_end > (centromeres$start[i]-1000000)
                                           ) > 0
        whole_arm_alt_by_chr[idx, i*2] = sum( S$major != 2 & S$chrom==i & 
                                              S$loc_start < (centromeres$start[i]+100000) & 
                                              S$loc_end > (chr_lengths[i]-1000000)) > 0    
    }
}
stopCluster(cl)

n_states_by_chr = matrix(NA, nrow=n_lines, ncol=24)
weighted_density_states_by_chr = matrix(NA, nrow=n_lines, ncol=24)
whole_arm_alt_by_chr = matrix(NA, nrow=n_lines, ncol=48)
focal_alt_by_chr = matrix(NA, nrow=n_lines, ncol=24)

n_states_by_chr = data.frame( n_states_by_chr, row.names=line_names )
weighted_density_states_by_chr = data.frame( signif(weighted_density_states_by_chr,3), row.names=line_names )
whole_arm_alt_by_chr = data.frame( whole_arm_alt_by_chr, row.names=line_names )
focal_alt_by_chr = data.frame( focal_alt_by_chr, row.names=line_names )
names( n_states_by_chr ) = paste("chr", 1:24, sep="")
names( weighted_density_states_by_chr ) = paste("chr", 1:24, sep="")
x=c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,13,14,14,15,15,
    16,16,17,17,18,18,19,19,20,20,21,21,22,22,23,23,24,24)
x = paste( "chr", x, rep( c("p", "q"), 24), sep="")
names( whole_arm_alt_by_chr ) = x
names( focal_alt_by_chr ) = paste("chr", 1:24, sep="")

write.matrix( n_states_by_chr, '/datasets/human_lines_database/CN_n_states_by_chr.txt')
write.matrix( weighted_density_states_by_chr, '/datasets/human_lines_database/CN_weighted_density_states_by_chr.txt')
write.matrix( whole_arm_alt_by_chr, '/datasets/human_lines_database/CN_whole_arm_alt_by_chr.txt')
write.matrix( focal_alt_by_chr, '/datasets/human_lines_database/CN_focal_alt_by_chr.txt')
