from optparse import OptionParser
parser = OptionParser()
parser.add_option("-i", "--file_in", dest="fn_in", help="path to humandb/refseq/hg19_refGene.txt", default="")
parser.add_option("-o", "--file_out", dest="fn_out", help="output exon file to write", default="")
(options, args) = parser.parse_args()
fn_in = options.fn_in
fn_out= options.fn_out
f = open(fn_in)
fo = open(fn_out, 'w')
for line in f:
    a = line.rstrip('\r\n').split('\t')
    chromosome, refseq = a[2], a[1]
    chromosome = chromosome.replace("chr","")
    exidx = 1
    for exStart, exEnd in zip( a[9].split(","), a[10].split(",") ):
        if len(exStart)>0:
            fo.write( chromosome + '\t' + exStart + '\t' + exEnd + '\t' + refseq + ":" + str(exidx) +  "\n")
        
        exidx += 1

fo.close()
f.close()

