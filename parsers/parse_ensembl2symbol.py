rs2s = {}
f=open('/notebook/annotations/UCSC_all_known_mouse_genes_December_2011_MM10.txt')
h = f.readline()
for line in f:
    a = line.rstrip('\r\n').split('\t')
    rs2s[a[0]] = a[1]

f.close()

e2rs = {}
f = open('/notebook/annotations/gene2ensembl_2012_12_28.txt')
h=f.readline()
h=f.readline()
for line in f:
    a = line.rstrip('\r\n').split('\t')
    e2rs[ a[2] ] = a[3].split('.')[0]

f.close()

f = open('/notebook/annotations/ensembl2symbol.txt', 'w')
f.write( "ensembl\trs_id\tsymbol\n")
for ens in e2rs.keys():
    try:
        rs=e2rs[ens]
        symbol = rs2s[rs]
        f.write( ens + '\t' + rs + '\t' + symbol + '\n' )
    except:
        pass

f.close()

