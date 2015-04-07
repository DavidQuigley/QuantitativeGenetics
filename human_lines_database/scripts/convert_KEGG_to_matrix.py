fn_msig = "/notebook/annotations/msigdb_gene_sets_2013_05.txt"
fn_matrix = '/datasets/human_lines_database/entities/KEGG_matrix_from_msigdb_2013_05.txt'
f = open( fn_msig )
f.readline() # header
pathways = {}
all_symbols = {}

for line in f:
    a = line.rstrip('\r\n').split('\t')
    if a[0]=="KEGG":
        symbols = a[3].split(' ')
        pth = {}
        for symbol in symbols:
            all_symbols[symbol] = 1    
        
        pathways[ a[1] ] = symbols

f.close()
fo = open(fn_matrix, 'w')
all_symbols = sorted(all_symbols)
fo.write("IDENTIFIER\t" + '\t'.join(all_symbols) + '\n')
for pathway in sorted(pathways):
    o = [pathway]
    for symbol in all_symbols:
        if symbol in pathways[pathway]:
            o.append("1")
        
        else:
            o.append("0")
    
    fo.write('\t'.join(o) + '\n')

fo.close()
