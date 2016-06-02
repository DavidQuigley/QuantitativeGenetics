import vcf
from optparse import OptionParser
parser = OptionParser()
parser.add_option("-i", "--file_in", dest="fn_in", help="path to vcf", default="")
parser.add_option("-o", "--file_out", dest="fn_out", help="output grid file to write", default="")
(options, args) = parser.parse_args()
fn_in = options.fn_in
fn_out = options.fn_out

vcf_rdr = vcf.Reader(open(fn_in, 'r'))
record = next(vcf_rdr)
samples=[]
for s in record.samples:
    samples.append(s.sample)

samples = sorted( samples )

missense = set(['nonsynonymous_SNV'])
nonsense = set(['stopgain','frameshift_insertion', 'frameshift_deletion'])

vcf_rdr = vcf.Reader(open(fn_in, 'r'))
symbol2sample2val = {}
for record in vcf.Reader(open(fn_in, 'r')):
    symbol = record.INFO["Gene.refGene"][0]
    if not symbol in symbol2sample2val:
        sample2val = {}
        for sample in samples:
            sample2val[sample] = 0
        
        symbol2sample2val[symbol] = sample2val
    
    for sample_obj in record.samples:
        gt = sample_obj['GT']
        s_func=record.INFO['ExonicFunc.refGene'][0]
        if not (gt == './.' or gt=="0|0" or gt=="0/0"):
            if s_func in missense:
                mut_code=1
            elif s_func in nonsense:
                mut_code=-1
            else:
                mut_code=0
            symbol2sample2val[symbol][sample_obj.sample] = mut_code

symbols = sorted( symbol2sample2val.keys() )

fo = open(fn_out, 'w')
fo.write( "symbol\t" + '\t'.join( samples ) + '\n' )
for symbol in symbols:
    vals = []
    for sample in samples:
        vals.append( symbol2sample2val[symbol][sample] )
    fo.write( symbol + '\t' + '\t'.join([str(x) for x in vals]) + '\n' )

fo.close()
