# Have to use entrez table to find certain validate genes that are not in the standard 
# refseq table download. Example: TEN1, entrez ID 100134934, shows up as 
# chr17	73975297	73996667	n/a	n/a	100134934 in this download.

# Entrez2symbol comes from
# curl ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2refseq.gz
# gunzip -c gene2refseq.gz | grep -e 'REVIEWED' -e 'VALIDATED' | cut -f 2,16 > /datasets/human_lines_database/annotations/entrez2symbol

# ucsc_all_refseq_genes_2015_02_16 comes from a table download 

entrez2symbol = {}
f = open('/datasets/human_lines_database/annotations/entrez2symbol')
for line in f:
    if line[0] != "#":
        a = line.rstrip('\r\n').split('\t')
        entrez2symbol[a[0]] = a[1]

f.close()

f = open('/datasets/human_lines_database/annotations/ucsc_all_refseq_genes_2015_02_16.txt')
fo = open('/datasets/human_lines_database/annotations/refseq_genes_hg19.BED', 'w')
s2chrom, s2start, s2stop = {}, {}, {}
h = f.readline()
valid_chr = set(["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"] )
for line in f:
    a = line.rstrip('\r\n').split('\t')
    chrom, txStart, txEnd, symbols, refseq, entrez = a
    if symbols=="n/a" or refseq=="n/a":
        if entrez in entrez2symbol:
            symbols = entrez2symbol[entrez]
            refseq = "NM_foo"
        
        else:
            continue
    
    chrom = chrom.replace("chr", "")
    if not chrom in valid_chr:
        continue
    
    if not "NM_" in refseq:
        if not "NR_" in refseq:
            continue
    
    symbols = symbols.rstrip(",").split(",")
    for symbol in symbols:
        if symbol in s2chrom:
            if int(txStart) < int(s2start[symbol]):
                s2start[ symbol ] = txStart
            
            if int(txEnd) < int(s2stop[symbol]):
                s2stop[ symbol ] = txEnd
            
        else:
            s2chrom[ symbol ] = chrom
            s2start[ symbol ] = txStart
            s2stop[ symbol ] = txEnd

for symbol in sorted(s2chrom.keys()):
    fo.write( s2chrom[symbol] + '\t' + s2start[symbol] + '\t' + s2stop[symbol] + '\t' + symbol + '\n' )

fo.close()
f.close()
