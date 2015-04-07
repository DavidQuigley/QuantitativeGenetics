# The purpose of this script is to populate the MUT matrix with values:
#  0 for no known cancer-associated mutation
#  -1 for a truncating or nonsense mutation
#  1 for a missense mutation associated with cancer
# We will pull first from the sanger exomes
# for samples where there is no sanger coverage, pull from the CCLE's targeted sequencing
# for samples lacking other coverage, pull from genentech's RNA-seq derived values

source('/notebook/code/quantitative_genetics.R')

fn_genes = '/datasets/human_lines_database/annotations/refseq_genes_hg19_2015_02_20.BED'

fn_ccle = "/datasets/human_lines_CCLE/CCLE_mutations_hg19_NoCommonSNPs_NoNeutralVariants_CDS_2012.05.07_snp141_annot.maf"

# The actual COSMIC database
fn_cosmic_raw = "/datasets/human_lines_database/COSMIC/CosmicMutantExportCensus_2015_01_21.tsv"
fn_cosmic = "/datasets/human_lines_database/COSMIC/CosmicMutantExportCensus_2015_01_21_compressed.txt"
fn_lines = '/datasets/human_lines_database/cell_line_attributes/cell_line_dictionary_2015_02_18.txt'
# The Sanger cell lines database
fn_sanger_raw = "/datasets/human_lines_sanger/CosmicCLP_MutantExport.tsv"
fn_sanger = "/datasets/human_lines_sanger/CosmicCLP_MutantExport_compressed.txt"
fn_ensembl = '/datasets/human_lines_database/annotations/ensembl_transcripts_hg19.txt'
fn_genentech = '/datasets/human_lines_genentech/140331_cellLineMutations_compressed_by_gene.txt'

cmd_cut_cosmic = paste( "cut -f 1,5,13,14,15,16,17,18,20,22", fn_cosmic_raw, ">", fn_cosmic)
cmd_cut_sanger = paste("cut -f 1,2,5,6,13,15,16,17,20,21", fn_sanger_raw, ">", fn_sanger)

#cut -f 1,50 /datasets/human_lines_genentech/140331_cellLineMutations.txt | uniq | sort -f -k 2 > /datasets/human_lines_genentech/sam_id

#system(cmd_cut_cosmic)
#system(cmd_cut_sanger)

#----------------------------------------------------------------
# LOAD GENE REFERENCE
#----------------------------------------------------------------

genes = read.table(fn_genes, header=FALSE, sep='\t', stringsAsFactors=FALSE)
names(genes) = c("chrom", "txStart", "txEnd", "symbol")
symbols = genes$symbol
hsh_symbol = hsh_from_vectors(symbols, 1:length(symbols))

cl = load.matrix(fn_lines)

MUT = matrix(NA, nrow=length(symbols), ncol=dim(cl)[1])

dimnames(MUT)[[1]] = symbols
dimnames(MUT)[[2]] = rownames(cl)
#----------------------------------------------------------------
# LOAD SANGER EXOMES
#----------------------------------------------------------------
sanger = read.table(fn_sanger, sep='\t', header=TRUE, stringsAsFactors=FALSE)
names(sanger) = c("symbol", "ensembl", "name_sanger", "id_sample", "mut_id", "mut_aa", "mut_desc", 
                  "mut_zyg", "SNP", "mut_somatic_status")

# Sanger double-lists some transcripts for the same gene. Consolidate if the symbol, 
# cell line, and mutation aa are identical
seen = hsh_new()
keep = rep(FALSE, dim(sanger)[1])
token = paste( get.split.col( sanger$symbol, "_", first=TRUE), 
       sanger$name_sanger, sanger$mut_aa, sep="___")
for(i in 1:length(keep)){
    if( !hsh_in( seen, token[i] ) ){
        keep[i] = TRUE
        hsh_set( seen, token[i], 1)
    }
}
sanger = sanger[keep,]
sanger$symbol = get.split.col( sanger$symbol, "_", first=TRUE)

# Remove SNPs and silent mutations
sanger = sanger[sanger$SNP=="n",]
sanger = sanger[sanger$mut_desc!="Substitution - coding silent",]
sanger = sanger[sanger$mut_desc != "Unknown",]
sanger = sanger[sanger$mut_desc != "Insertion - In frame",]
sanger = sanger[sanger$mut_desc != "No detectable mRNA/protein",]
sanger = sanger[!(sanger$mut_desc=="Substitution - Missense" & sanger$mut_somatic_status=="PASSENGER/OTHER"),]
sanger = sanger[!(sanger$mut_desc=="Substitution - Missense" & sanger$mut_somatic_status==""),]

# remove genes that are not present in primary gene list
# attempt to rescue genes with incorrect names by cross-referencing their ensembl id
ensembl = load.matrix(fn_ensembl)
symbol.no.match = setdiff(sanger$symbol, symbols ) 
for(i in 1:length(symbol.no.match)){
    idx = which( sanger$symbol==symbol.no.match[i] )
    ensembl_id = sanger$ensembl[idx[1]]
    idx.in.ensembl = which(rownames(ensembl)==ensembl_id)
    if(length(idx.in.ensembl)==1){
        sanger$symbol[idx] = ensembl$symbol[ idx.in.ensembl ]
    }
}
symbol.no.match = setdiff(sanger$symbol, symbols )

# remove mutation calls from the ~700 genes that still don't match
hsh_symbols = hsh_from_vectors(symbols, 1:length(symbols))
keep = rep(FALSE, dim(sanger)[1])
for(i in 1:dim(sanger)[1]){
    keep[i] = hsh_in( hsh_symbols, sanger$symbol[i])
}
sanger = sanger[keep,]
# At this point sanger mutation set consists only of known genes


hsh_sanger = hsh_from_vectors(cl$name_sanger, 1:dim(cl)[1])
# standardize sample names
sanger$name_sanger[ sanger$name_sanger=="FADU"] = "FaDu"
sanger$name_sanger[ sanger$name_sanger=="NTERA-2_cl_D1"] = "NTERA-S-cl-D1"

# for the 1025 lines with mutation calls, mark their default as 0 rather than 
# NA. This allows us to distinguish "no mutation found" from "no information present"
names_mut = unique(sanger$name_sanger)
m = match.idx( cl$name_sanger, names_mut)
for( i in 1:dim(m)[1]){
    MUT[, m$idx.A[i]] = 0
}
for(i in 1:dim(sanger)[1]){
    row = hsh_get(hsh_symbol, sanger$symbol[i])
    col = hsh_get(hsh_sanger, sanger$name_sanger[i])
    if( sanger$mut_desc[i]=="Substitution - Missense" )
        MUT[row,col] = 1
    else if( sanger$mut_desc[i]=="Deletion - Frameshift" | 
             sanger$mut_desc[i]=="Substitution - Nonsense"  )
        MUT[row, col] = -1
}

#----------------------------------------------------------------
# LOAD CCLE targeted sequencing
# The CCLE database contains calls that are actually SNPs, even in the file annotated as 
# "no common SNPs." I check the CCLE calls against snp141Common, obtained by
# curl http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/snp141Common.txt.gz
# unzip snp141Common.txt.gz
# cut -f 2,3,4,5 snp141Common.txt > snp141Common_locations.txt'
# gzip snp141Common_locations.txt
# mv snp141Common_locations.txt.gz /datasets/human_lines_database/annotations
# python /datasets/human_lines_database/scripts/CCLE_flag_SNPs.py \
# --dbSNP /datasets/human_lines_database/annotations/snp141Common_locations.txt.gz \
# --ccle /datasets/human_lines_ccle/CCLE_mutations_hg19_NoCommonSNPs_NoNeutralVariants_CDS_2012.05.07.maf \
# --ccle_out /datasets/human_lines_ccle/CCLE_mutations_hg19_NoCommonSNPs_NoNeutralVariants_CDS_2012.05.07_snp141_annot.maf
#
#----------------------------------------------------------------
ccle = read.table( fn_ccle, sep='\t', header=TRUE, stringsAsFactors=FALSE)

vc_clean = table(ccle$Variant_Classification[ccle$snp141 == "none"])
vc_snp = table(ccle$Variant_Classification[ccle$snp141 != "none"])

vcs = sort(unique(ccle$Variant_Classification))
vcs_clean = rep(0, length(vcs))
vcs_snp = rep(0, length(vcs))
for(i in 1:length( vcs )){
    vcs_clean[i] = sum( ccle$Variant_Classification==vcs[i] & ccle$snp141=="none")
    vcs_snp[i] = sum( ccle$Variant_Classification==vcs[i] & ccle$snp141!="none")
}
data.frame( clean=vcs_clean, snp=vcs_snp, row.names=vcs, stringsAsFactors=FALSE)
ccle_raw = ccle
ccle = ccle[ccle$snp141 == "none",]
ccle = ccle[ ccle$Variant_Classification %in% 
            c("Frame_Shift_Del", "Missense_Mutation","Nonsense_Mutation","Splice_Site_Del","Splice_Site_SNP"),]
freq = table( paste(ccle$Chromosome,":",ccle$Start_position, sep=''))
freq = hsh_from_vectors( names(freq), as.numeric(freq) )
mutfreq = hsh_get( freq, paste(ccle$Chromosome,":",ccle$Start_position, sep='') )
ccle = cbind(ccle, mutfreq, stringsAsFactors=FALSE)
ccle_keep = ccle$Variant_Classification %in% c("Frame_Shift_Del", "Nonsense_Mutation","Splice_Site_Del","Splice_Site_SNP") | 
              (ccle$Variant_Classification=="Missense_Mutation" & ccle$mutfreq>1 )
ccle = ccle[ccle_keep,]

samples_ccle=unique(ccle$Tumor_Sample_Barcode)
sample2idx = hsh_from_vectors( samples_ccle, match.idx( samples_ccle, cl$name_ccle )$idx.B )
gene2idx = hsh_from_vectors( dimnames(MUT)[[1]], 1:dim(MUT)[1])

columns = hsh_get(sample2idx, ccle$Tumor_Sample_Barcode)
rows = hsh_get( gene2idx, ccle$Hugo_Symbol)

columns_sanger = match.idx( cl$name_sanger, names_mut)$idx.A
rows_in_ccle = sort(unique(rows))

# There are 1472 genes with coverage in CCLE. Mark them as 0 rather than NA
MUT[ rows_in_ccle, setdiff(columns, columns_sanger) ] = 0

# now mark the mutations seen among those cell lines
for(i in 1:dim(ccle)[1]){
    if( !is.na( rows[i] & !columns[i] %in% columns_sanger  )){
        if( ccle$Variant_Classification[i]=="Missense_Mutation" ){
            MUT[rows[i], columns[i] ] = 1
        }
        else{
            MUT[rows[i], columns[i] ] = -1
        }
    }
}

m = data.frame(MUT)
names(m) = rownames(cl)
date=format( Sys.Date(), "%Y_%m_%d")
write.matrix(m, paste('/datasets/human_lines_database/mutation/mutations_mis_del_non_', date, '.txt', sep=''))
