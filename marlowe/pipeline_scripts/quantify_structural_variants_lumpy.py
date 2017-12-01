import sys
import argparse
import os
parser = argparse.ArgumentParser(description='quantify structural variants, filter VCF, create CIRCOS input')
parser.add_argument('--in', dest="fn_in", help='vcf file from Lumpy', default="")
parser.add_argument('--out', dest="fn_out", help='base of file to write', default="")
parser.add_argument('--min_del_distance', dest="min_del_distance", help='minimum deletion distance, default 1000', default="1000")
parser.add_argument('--min_break_distance', dest="min_break_distance", help='minimum break distance, default 30', default="30")
parser.add_argument('--min_tumor_reads', dest="min_reads", help='minimum tumor read support, default 10', default="10")
parser.add_argument('--max_normal_reads', dest="max_reads", help='maximum normal read support, default 0', default="0")

args = parser.parse_args()

fn_in = args.fn_in
fn_out = args.fn_out
if not os.path.exists(fn_in):
    print "ERROR: parameter passed to --in (" + fn_in + ") cannot be found as file"
    sys.exit(1)

fn_out_vcf = fn_out.rstrip('.vcf') + '_filtered.vcf'
fn_out_circos = fn_out.rstrip('.vcf') + '_filtered_circos.txt'
min_tumor_reads = int( args.min_reads )
max_normal_reads = int( args.max_reads )
min_break_distance = int( args.min_break_distance )
min_del_distance = int( args.min_del_distance )

f = open( fn_in )
fo_vcf = open(fn_out_vcf, 'w')
fo_circos = open(fn_out_circos, 'w')

print( "MESSAGE: read " + fn_in )
print( "MESSAGE: writing VCF to " + fn_out_vcf )
print( "MESSAGE: writing Circos to " + fn_out_circos )
print( "MESSAGE: Minimum deletion distance: " + str(min_del_distance))
print( "MESSAGE: Minimum breakpoint distance: " + str(min_break_distance))
print( "MESSAGE: Minimum tumor supporting reads: " + str(min_tumor_reads))
print( "MESSAGE: Maximum normal supporting reads: " + str(max_normal_reads))

n_del = 0
n_dup = 0
n_inv = 0
n_trans = 0

n_lines = 0
for line in f:
    if "#" in line:
        fo_vcf.write(line)
    else:
        n_lines = n_lines + 1
        a = line.rstrip('\r\n').split('\t')
        chr_1=a[0]        
        start_1 = a[1]
        info, info_t, info_n = a[7], a[9], a[10]
        GT_t,SU_t,PE_t,SR_t	= info_t.split(":")
        GT_n,SU_n,PE_n,SR_n	= info_n.split(":")
        info = info.split( ";" )
        if (int(SU_t)>min_tumor_reads or int(PE_t)>min_tumor_reads) and (int(SU_n)==max_normal_reads and int(PE_n)==max_normal_reads ):          
            if "BND" in line:
                muddle = a[4]
                try:
                    chr_2 = "chr" + muddle.split(":")[0].split("chr")[1]
                    start_2 = muddle.split(":")[1].split("[")[0].split("]")[0]
                    if chr_1 != chr_2 or abs(float(start_1) - float(start_2)) > min_break_distance:
                        fo_vcf.write( line )                    
                        fo_circos.write( '\t'.join( [chr_1, start_1, str(int(start_1)+1), chr_2, start_2, str(int(start_2)+1)] ) + '\n')        
                        n_trans = n_trans + 1
                except IndexError:
                    continue
            else: 
                for nugget in info:
                    try:
                        key, val = nugget.split("=")
                    except ValueError:
                        key = nugget
                    if key=="SVLEN":
                        observed_del_length = abs(int(val))
                if "<DUP>" in line:
                    fo_vcf.write( line )
                    n_dup = n_dup + 1
                
                elif "<DEL>" in line and observed_del_length >= min_del_distance:
                    fo_vcf.write( line )
                    n_del = n_del + 1
                
                elif "<INV>" in line:
                    fo_vcf.write( line )
                    n_inv = n_inv + 1

fo_vcf.close()
fo_circos.close()
f.close()

print( "MESSAGE: completed reading " + str(n_lines) + " VCF entries" )
