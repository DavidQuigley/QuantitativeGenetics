fn_gene = '/Volumes/Exp15_sgRNA_script/classifications_DNArepair.txt'
fn_sgrna = '/Volumes/Exp15_sgRNA_script/sgRNAs_dnarepair.txt'
fn_neg = '/Volumes/Exp15_sgRNA_script/negative_controls.txt'
fn_adapt = '/Volumes/Exp15_sgRNA_script/adapter_sequences_final.txt'

genes = read.table(fn_gene, header=TRUE, stringsAsFactors=FALSE)
sgrna = read.table(fn_sgrna, header=TRUE, stringsAsFactors=FALSE, sep='\t')
negs = read.table(fn_neg, header=TRUE, stringsAsFactors=FALSE, sep='\t')
adapt = read.table(fn_adapt, header=TRUE, stringsAsFactors=FALSE, sep='\t')

# restrict sgrna to elements where gene_ID is in gene$gene_ID
sgrna = sgrna[sgrna$gene_ID %in% genes$gene_ID,]

# match genes$Sub_library value to elements in sgrna, accounting for the 
# fact that genes appear more than once in sgrna
sgrna_sub = rep("", dim(sgrna)[1])
for(i in 1:length(sgrna_sub)){
    sgrna_sub[i] = genes$Sub_library[ genes$gene_ID==sgrna$gene_ID[i] ]
}
sgrna = cbind( sgrna, Sub_library=sgrna_sub, stringsAsFactors=FALSE)

# calculate 2% of the number of elements in each sublibrary
#
n_guides_required = ceiling( table(sgrna$Sub_library)*0.02 )

# pick the corresponding number of negative guides
#
is_picked = rep( FALSE, dim(negs)[1] )
pick_length = as.numeric(n_guides_required)
neg_guides = list()
for( i in 1:length(pick_length)){    
    idx_to_pick = which(!is_picked)[1:pick_length[i]]
    negs_selected = negs[idx_to_pick,]
    negs_selected = cbind( negs_selected, 
                           Sub_library=names(n_guides_required)[i], 
                           stringsAsFactors=FALSE )
   neg_guides = c(neg_guides, list( negs_selected ) ) 
   is_picked[idx_to_pick] = TRUE
}
names(neg_guides) = names(n_guides_required)

# create a new column called full.sequence that consists of the Forward, proto,
# and reverse primer for each individual guide library
neg_guides$Check$full.sequence = 
  paste( adapt$Forward_Full[adapt$Matching.library=="DNA_Check"], 
         neg_guides$Check$protospacer.sequence, 
         adapt$Reverse_full[adapt$Matching.library=="DNA_Check"], sep='' )

neg_guides$DSB$full.sequence = 
  paste( adapt$Forward_Full[adapt$Matching.library=="DNA_DSB"], 
         neg_guides$DSB$protospacer.sequence, 
         adapt$Reverse_full[adapt$Matching.library=="DNA_DSB"], sep='' )

neg_guides$Other$full.sequence = 
  paste( adapt$Forward_Full[adapt$Matching.library=="DNA_Other"], 
         neg_guides$Other$protospacer.sequence, 
         adapt$Reverse_full[adapt$Matching.library=="DNA_Other"], sep='' )

neg_guides$SSB$full.sequence = 
  paste( adapt$Forward_Full[adapt$Matching.library=="DNA_SSB"], 
         neg_guides$SSB$protospacer.sequence, 
         adapt$Reverse_full[adapt$Matching.library=="DNA_SSB"], sep='' )


# create a new column called full.sequence that consists of the Forward, proto,
# and reverse primer for each individual guide library

sgrna$full.sequence[ sgrna$Sub_library=="Check"] = 
    paste( adapt$Forward_Full[adapt$Matching.library=="DNA_Check"], 
           sgrna$protospacer.sequence[ sgrna$Sub_library=="Check" ], 
           adapt$Reverse_full[adapt$Matching.library=="DNA_Check"], sep='' )

sgrna$full.sequence[ sgrna$Sub_library=="DSB"] = 
    paste( adapt$Forward_Full[adapt$Matching.library=="DNA_DSB"], 
           sgrna$protospacer.sequence[ sgrna$Sub_library=="DSB" ], 
           adapt$Reverse_full[adapt$Matching.library=="DNA_DSB"], sep='' )

sgrna$full.sequence[ sgrna$Sub_library=="Other"] = 
    paste( adapt$Forward_Full[adapt$Matching.library=="DNA_Other"], 
           sgrna$protospacer.sequence[ sgrna$Sub_library=="Other" ], 
           adapt$Reverse_full[adapt$Matching.library=="DNA_Other"], sep='' )

sgrna$full.sequence[ sgrna$Sub_library=="SSB"] = 
    paste( adapt$Forward_Full[adapt$Matching.library=="DNA_SSB"], 
           sgrna$protospacer.sequence[ sgrna$Sub_library=="SSB" ], 
           adapt$Reverse_full[adapt$Matching.library=="DNA_SSB"], sep='' )

M = neg_guides$Check[,c("sgID", "Sub_library", "full.sequence") ]
M = rbind(M, neg_guides$DSB[,c("sgID", "Sub_library", "full.sequence") ] )
M = rbind(M, neg_guides$SSB[,c("sgID", "Sub_library", "full.sequence") ] )
M = rbind(M, neg_guides$Other[,c("sgID", "Sub_library", "full.sequence") ]  )
M = rbind(M, sgrna[, c("sgID", "Sub_library", "full.sequence") ] )
