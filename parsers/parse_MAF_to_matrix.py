# given a mutation file in the MAF format, parse to a matrix of genes

import argparse
parser = argparse.ArgumentParser(description='convert TCGA copy number data to BED files')
parser.add_argument('--in', dest="fn_maf", help='path to MAF file')
parser.add_argument('--samples', dest="fn_samples", help='path to sample attribute file')
parser.add_argument('--genes', dest="fn_genes", help='path to BED with accepted genes')
parser.add_argument('--out', dest="fn_out", help='file to write')

args = parser.parse_args()
fn_genes=args.fn_genes
fn_in = args.fn_maf

symbols = []
f=open(fn_genes)
for line in f:
    a = line.rstrip('\r\n').split('\t')
    symbols.append( a[3] )

f.close()
symbols = sorted(symbols)

s2r2v = {}
barcode2s = {}
f = open(fn_sa)
h=f.readline()
for line in f:
    a = line.rstrip('\r\n').split('\t')
    barcode2s[a[1]] = a[0]
    s2r2v[ a[0] ] = {}

f.close()


sample_names = sorted( barcode2s.values() )

f=open(fn_in)
h = f.readline().rstrip('\r\n').split('\t')
idx_barcode = h.index("Tumor_Sample_Barcode")
idx_chrom = h.index("Chrom")
idx_start = h.index("Start_Position")
idx_end = h.index("End_Position")
idx_symbol = h.index("Hugo_Symbol")
idx_class = h.index( "Variant_Classification")
idx_type = h.index("Variant_Type")

mut_del = set([ "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site"])

for line in f:
    a = line.rstrip('\r\n').split('\t')
    patient_no = barcode2s[ a[idx_barcode] ]
    varclass = a[idx_class]
    if varclass in mut_del:
        s2r2v[ patient_no ][ a[idx_symbol] ] = "-1"    
    elif varclass == "Missense_Mutation":
        s2r2v[ patient_no ][ a[idx_symbol] ] = "1"

f.close()

fo = open( fn_out, 'w' )
fo.write( "IDENTIFIER\t" + '\t'.join( sample_names ) + '\n' )
for symbol in symbols:
    res = [ s2r2v[sample].get( symbol, "0" ) for sample in sample_names ]
    print(len(res))
    fo.write( symbol + '\t' + '\t'.join( res ) + '\n' )

fo.close()

