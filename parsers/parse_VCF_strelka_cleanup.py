# Add a FORMAT specifier to the header 
#    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
# For each result, add
#    GT: to the front of column 8 (zero-indexed), 
#    0/0: to column 9
#    1/1: to column 10

from optparse import OptionParser
parser = OptionParser()
parser.add_option("-i", "--file_in", dest="fn_in", help="path to VCF in", default="")
parser.add_option("-o", "--file_out", dest="fn_out", help="path to VCF out", default="")
(options, args) = parser.parse_args()
fn_in = options.fn_in
fn_out = options.fn_out
if fn_in==fn_out:
    raise RuntimeError("cannot write to file you're reading from")

f = open(fn_in)
fo = open(fn_out, "w")
written_GT = False
for line in f:
    if line[0] == "#":
        if not written_GT and "##FORMAT" in line:
            written_GT= True
            fo.write( '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        
        fo.write(line)
    
    else:
        a = line.rstrip('\r\n').split('\t')
        a[8] = "GT:" + a[8]
        a[9] = "0/0:" + a[9]
        a[10] = "1/1:" + a[10]
        line = "\t".join( a )
        fo.write(line + '\n')

fo.close()
f.close()
