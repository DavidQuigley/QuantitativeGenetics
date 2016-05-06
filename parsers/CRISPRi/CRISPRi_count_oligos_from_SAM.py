from optparse import OptionParser
parser = OptionParser()
parser.add_option("-i", "--file_in", dest="fn_in", help="path to SAM in", default="")
parser.add_option("-o", "--file_out", dest="fn_out", help="path to COUNTS base out", default="")
(options, args) = parser.parse_args()
fn_in = options.fn_in
fn_out = options.fn_out

symbol_counts, oligo_counts = {}, {}
f = open(fn_in)
for line in f:
    a = line.rstrip('\r\n').split('\t')
    oligo = a[2]
    symbol= oligo.split("_")[0]
    try:
        oligo_counts[oligo]=oligo_counts[oligo]+1
    
    except KeyError:
        oligo_counts[oligo]=1
    
    try:
        symbol_counts[symbol]=symbol_counts[symbol]+1
    
    except KeyError:
        symbol_counts[symbol]=1

f.close()

fos = open(fn_out + "_symbol.counts", "w")
for symbol in sorted(symbol_counts.keys()):
    fos.write( symbol + '\t' + str(symbol_counts[symbol]) + '\n' )

fos.close()
fos = open(fn_out + "_oligo.counts", "w")
for oligo in sorted(oligo_counts.keys()):
    fos.write( oligo + '\t' + str(oligo_counts[oligo]) + '\n' )

fos.close()
