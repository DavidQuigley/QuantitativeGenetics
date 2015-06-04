import os
import argparse
parser = argparse.ArgumentParser(description='convert TCGA copy number data to BED files')
parser.add_argument('--dir_source', dest="dir_source", help='path to read files')
parser.add_argument('--out', dest="fn_out", help='output file')

args = parser.parse_args()
dir_segments = args.dir_source
fn_out = args.fn_out

sample_names = []
row_names = {}
s2r2v = {}

for fn_input in os.listdir(dir_segments):
    sample_names.append(fn_input)
    s2r2v[ fn_input ] = {}
    f = open( dir_segments + '/' + fn_input )
    for line in f:
        a = line.rstrip('\r\n').split('\t')
        row_names[a[0]] = 1
        s2r2v[ fn_input ][ a[0] ] = a[1]
    
    f.close()

fo = open( fn_out, 'w' )
sample_names = sorted(sample_names)
row_names = sorted(list(row_names))
fo.write( "IDENTIFIER\t" + '\t'.join( sample_names ) + '\n' )
for row in row_names:
    fo.write( row + '\t' + '\t'.join( [ s2r2v[sample].get( row, "NA" ) for sample in sample_names ] ) + '\n' )

fo.close()

