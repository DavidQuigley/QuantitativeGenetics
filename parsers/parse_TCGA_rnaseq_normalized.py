from optparse import OptionParser

parser = OptionParser()
parser.add_option("-f", "--file_in", dest="fn_in", help="tab-delimited path to rnaseq file from TCGA firehose", default="")
parser.add_option("-s", "--file_sa_in", dest="fn_sa_in", help="tab-delimited path to sample file from TCGA firehose", default="")
parser.add_option("-g", "--file_genes_out", dest="fn_genes_out", help="Path to write output gene file", default="")
parser.add_option("-o", "--file_expr_out", dest="fn_expr_out", help="Path to write output expression file", default="")
parser.add_option("-a", "--file_sa_out", dest="fn_sa_out", help="Path to write output sample file", default="")
parser.add_option("-r", "--require_symbol", dest="require_symbol", help="require gene symbol, {T,F}", default="")

(options, args) = parser.parse_args()
fn_in = options.fn_in
fn_sa_in = options.fn_sa_in
fn_genes_out = options.fn_genes_out
fn_expr_out = options.fn_expr_out
fn_sa_out = options.fn_sa_out
require_symbol = options.require_symbol

if require_symbol=="T":
    require_symbol=True
elif require_symbol=="F":
    require_symbol=False
else:
    raise RuntimeError("require_symbol must be in {T,F}")

print( "Reading " + fn_in )
print( "Writing to " + fn_genes_out )
print( "Writing to " + fn_expr_out )
print( "Writing to " + fn_sa_out )

f = open(fn_in)
fo_expr = open(fn_expr_out, 'w')
fo_ga = open(fn_genes_out, 'w')
header = f.readline().rstrip('\r\n').split('\t')
skip = f.readline()
for idx, h in enumerate(header):
    if idx==0:
        header[idx] = "IDENTIFIER"
    else:
        header[idx] = '-'.join( h.split("-")[0:3] )

fo_expr.write( "\t".join(header) + '\n' )
fo_ga.write( "IDENTIFIER\tsymbol\tgene_number\tTCGA_ID\n")
ctr=0
for line in f:
    a = line.rstrip('\r\n').split('\t')
    symbol, gene_num = a[0].split("|")
    if require_symbol and symbol=="?":
        continue
    else:
        a[0] = symbol
        fo_expr.write( '\t'.join(a) + '\n' )
        fo_ga.write( symbol + '\t' + symbol + '\t' + gene_num + '\t' + symbol + '|' + gene_num + '\n' )
        ctr = ctr+1

fo_expr.close()
fo_ga.close()
print("Wrote " + str(ctr) + " symbols")

f = open(fn_sa_in)
lines = []
for line in f:
    lines.append( line.rstrip('\r\n').split("\t") )

fo = open(fn_sa_out, 'w')

for x in zip(*lines):
    out = []
    is_first=True
    for y in x:
        if y=="Hybridization REF":
            out.append( "IDENTIFIER" )
        elif is_first:
            out.append( y.upper() )
            is_first=False
        else:
            out.append(y)
    
    fo.write( '\t'.join(out) + '\n' )

fo.close()

