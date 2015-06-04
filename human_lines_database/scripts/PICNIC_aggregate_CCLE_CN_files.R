source('/notebook/code/src/R/quantitative_genetics.R')

args <- commandArgs(trailingOnly = TRUE)
if( length(args) != 4 )
    stop("Four arguments are required: fn_genes, dir_copies, fn_out, allele_type")
fn_genes = args[1]
dir_copies = args[2]
fn_out = args[3]
allele = args[4]

if( allele != "total" & allele != "minor")
    stop("allele argument must be one of {total,minor}")

print(paste("Reading genes from", fn_genes))
genes = read.table(fn_genes, sep='\t', header=FALSE, stringsAsFactors=FALSE)
names(genes) = c('chrom', 'txStart', 'txEnd', 'symbol')
symbols = genes$symbol

print(paste("Read", length(symbols), "symbols"))
print(paste("Reading files from", dir_copies))
files=list.files(dir_copies, pattern='.*expanded.txt')
print(paste("Read", length(files), "files"))
symbol2idx = hsh_from_vectors( symbols, 1:length(symbols))
CN = matrix(NA, nrow=length(symbols), ncol=length(files))
cell_lines = rep("", length(files))

for( column in 1:length(files) ){
    copy_data = read.table(paste(dir_copies,files[column],sep="/"), sep='\t', header=TRUE, 
                           stringsAsFactors=FALSE)
    x=strsplit( files[column], "/", fixed=TRUE)[[1]]
    cell_line_name = strsplit( x[length(x)], "___", fixed=TRUE)[[1]][1]
    cell_lines[column] = cell_line_name
    print(paste(cell_lines[column],column))
        
    if( dim(copy_data)[1] > 0 ){
        # a few genes will appear more than once if they are on the bounds of a copy number 
        # change. We want to use the value that is not 2 if one exists.
    
        m = match.idx(symbols, copy_data$IDENTIFIER, allow.multiple.B=TRUE)
        copy_data = copy_data[m$idx.B,]
        if( allele=="total" ){
            allele_count = copy_data$major_allele
        }
        else{
            allele_count = copy_data$minor_allele
        }
        f=which(table(copy_data$IDENTIFIER)>1)
        has_dupes = hsh_new()
        if( length(f)>0)
            has_dupes = hsh_from_vectors( names(f), rep(1, length(f)))
        for(copy_row in 1:length(copy_data$IDENTIFIER)){
            symbol = copy_data$IDENTIFIER[copy_row] 
            if( hsh_in( has_dupes, symbol ) ){
                idx_dupes = which( copy_data$IDENTIFIER == symbol)
                largest_change = 0
                idx_largest = 1
                for(i in 1:length(idx_dupes)){
                    if( abs( 2 - allele_count[ idx_dupes[i] ]) > largest_change ){
                        idx_largest = idx_dupes[i]
                        largest_change = abs(2 - allele_count[ idx_dupes[i] ])
                    }
                }
                copies = allele_count[idx_largest]
            }
            else{
                copies = allele_count[copy_row] 
            }
            CN[ hsh_get(symbol2idx, symbol), column ] = copies
        }
    }
}

dimnames(CN)[[1]] = symbols
dimnames(CN)[[2]] = cell_lines
copies_to_write = data.frame(CN, stringsAsFactors=FALSE)
names(copies_to_write) = cell_lines
rownames(copies_to_write) = symbols
write.matrix(copies_to_write, fn_out)