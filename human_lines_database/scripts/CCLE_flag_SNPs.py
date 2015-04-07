import gzip 
import argparse

#fn_dbsnp = '/datasets/human_lines_database/annotations/snp141Common_locations.txt.gz'
#fn_ccle = '/datasets/human_lines_ccle/CCLE_mutations_hg19_NoCommonSNPs_NoNeutralVariants_CDS_2012.05.07.maf'
#fn_ccle_out = '/datasets/human_lines_ccle/CCLE_mutations_hg19_NoCommonSNPs_NoNeutralVariants_CDS_2012.05.07_snp141_annot.maf'
parser = argparse.ArgumentParser(description='Identify SNPs in CCLE mutation file')

parser.add_argument('ccle', help='CCLE mutation file')
parser.add_argument('dbSNP', help='dbSNP mutations, in BED format')
parser.add_argument('ccle_out', help='CCLE mutation file annotated with SNPs')
args = parser.parse_args()

fn_ccle = args.ccle
fn_pos_new = args.dbSNP

chr_loc = {}

f=gzip.open( fn_dbsnp )
for line in f:
    a = line.rstrip('\r\n').split('\t')
    chr_loc[ a[0].replace('chr', '') + ':' + a[1] ] = a[3]
    chr_loc[ a[0].replace('chr', '') + ':' + a[2] ] = a[3]

f.close()

f = open( fn_ccle )
fo = open( fn_ccle_out, 'w' )
h = f.readline().rstrip('\r\n')
fo.write( h + '\tsnp141' + '\n')

for line in f:
    a = line.rstrip('\r\n').split('\t')
    mstart = a[2] + ":" + a[3]
    mend = a[2] + ":" + a[4]
    if mstart in chr_loc:
        a.append( chr_loc[ mstart ] )
    
    elif mend in chr_loc:
        a.append( chr_loc[ mend ] )
    
    else:
        a.append( "none" )
    
    fo.write( '\t'.join( a ) + '\n' )

fo.close()
