# All but 285 loci were lifted over from HG18. However, some of those which were missed 
# are important as they are at the end of a chromosmoe.
# This code finds the closest assignment for deleted in original items. 
del = read.table('/datasets/human_lines_CCLE/CN/PICNIC/annotation/loci_deleted_in_hg19.txt',
      stringsAsFactors=FALSE, header=FALSE, sep='\t')
names(del) = c('chrom', 'pos', 'pos_plus_1', 'IDENTIFIER')

lift = read.table('/datasets/human_lines_CCLE/CN/PICNIC/annotation/hg19_SNP6_liftover_results.bed',
      stringsAsFactors=FALSE, header=FALSE, sep='\t')
names(lift) = c('chrom', 'pos', 'pos_plus_1', 'IDENTIFIER')

chrom_del = get.split.col(del$chrom, 'chr', last=TRUE)
new_ids = rep('', dim(del)[1])
new_pos = rep(0, dim(del)[1])
for(i in 1:dim(del)[1]){
    lift_same_chr = lift[ which( lift$chrom == chrom_del[i] ),]
    distance = abs( del$pos[i] - lift_same_chr$pos)
    idx = which( distance==min(distance) )[1]
    new_ids[i] = lift_same_chr$IDENTIFIER[ idx ]
    new_pos[i] = lift_same_chr$pos[ idx ]
}
cbind(del, new_ids, new_pos)

best_fits = data.frame(chrom_del, new_pos, new_pos+1, del$IDENTIFIER, stringsAsFactors=FALSE)
write.table(best_fits, 
  '/datasets/human_lines_CCLE/CN/PICNIC/annotation/hg19_SNP6_liftover_best_fits_for_deleted_loci.bed',
  row.names=FALSE, col.names=FALSE, sep='\t', quote=FALSE)