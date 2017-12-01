import os
import sys
import argparse
parser = argparse.ArgumentParser(description='extract summary data from a Manta run')
parser.add_argument('--in', dest="fn_in", help='vcf file from Lumpy', default="")
parser.add_argument('--out', dest="fn_out", help='file to write', default="")
args = parser.parse_args()
fn_in = args.fn_in
fn_out = args.fn_out

if not os.path.exists(fn_in):
    print "ERROR: parameter passed to --in (" + fn_in + ") cannot be found as file"
    sys.exit(1)

print( "MESSAGE: read " + fn_in )
print( "MESSAGE: writing Circos to " + fn_out )

f = open(fn_in)
fo=open(fn_out, "w")
fo.write( "Chromosome\tchromStart\tchromEnd\tChromosome.1\tchromStart.1\tchromEnd.1\n")
n_lines=0
for line in f:
    if line[0] != "#":
        n_lines = n_lines + 1
    
    if "PASS" in line and "IMPRECISE" not in line and "MantaBND" in line:
        a = line.rstrip('\r\n').split('\t')
        Chromosome=a[0]        
        chromStart = a[1]
        chromEnd = str(int(a[1])+1)
        muddle = a[4]
        if not "chr" in muddle or not "chr" in Chromosome: 
            continue
        n_pr, n_sr = a[9].split(':')
        t_pr, t_sr = a[10].split(':')
        n_pr = n_pr.split(",")[1]
        n_sr = n_sr.split(",")[1]
        t_pr = t_pr.split(",")[1]
        t_sr = t_sr.split(",")[1]
        Chromosome_1 = "chr" + muddle.split(":")[0].split("chr")[1]
        chromStart_1 = muddle.split(":")[1].split("[")[0].split("]")[0]
        chromEnd_1 = str(int(chromStart_1)+1)
        if int(n_pr)==0 and int(n_sr)==0 and int(t_sr)>10 :
            fo.write( '\t'.join( [Chromosome, chromStart, chromEnd, Chromosome_1, chromStart_1, chromEnd_1] ) + '\n')        

fo.close()
f.close()

print( "MESSAGE: completed reading " + str(n_lines) + " VCF entries" )