mappings = read.table('/datasets/human_lines_CCLE/CN/PICNIC/annotation/old_Snp6FeatureMappings.csv',  
                       sep=',', stringsAsFactors=FALSE)
mappings$V2 = paste("chr", mappings$V2, sep="")
mappings$V4 = mappings$V3+1
mappings$V5 = mappings$V1
mappings = mappings[,2:5]
write.table(mappings, '/datasets/human_lines_CCLE/CN/PICNIC/annotation/features_for_liftover',
            sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
