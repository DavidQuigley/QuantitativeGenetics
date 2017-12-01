import sys
import argparse
parser = argparse.ArgumentParser(description='convert TCGA copy number data to BED files')
parser.add_argument('--in', dest="fn_in", help='path to wig file')
parser.add_argument('--out', dest="fn_out", help='file to write')

args = parser.parse_args()

fn_in = args.fn_in
fn_out = args.fn_out

def wig_to_standard_chrom( fn_in, fn_out ):
    good_chroms = set(['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY'])
    write_ok = False
    f = open( fn_in )
    fo = open( fn_out, 'w')
    for line in f:
        if "chrom" in line:
            write_ok = line.split('=')[1].split(' ')[0] in good_chroms
        
        if( write_ok ):
            fo.write(line)
    
    f.close()
    fo.close()

wig_to_standard_chrom(fn_in, fn_out)
