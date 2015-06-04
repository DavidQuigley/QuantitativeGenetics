# rewrite the header of a file to use a new value that comes from a column in a sample 
# attribute file

import sys
import argparse
parser = argparse.ArgumentParser(description='convert TCGA copy number data to BED files')
parser.add_argument('--in', dest="fn_in", help='path to MAF file')
parser.add_argument('--samples', dest="fn_samples", help='path to sample attribute file')
parser.add_argument('--colname_existing', dest="colname_existing", help='column in sample attribute file')
parser.add_argument('--colname_new', dest="colname_new", help='column in sample attribute file')
parser.add_argument('--out', dest="fn_out", help='file to write')

args = parser.parse_args()

fn_in = args.fn_in
fn_out = args.fn_out
fn_samples = args.fn_samples
colname_existing = args.colname_existing
colname_new = args.colname_new

old2new = {}
f=open(fn_samples)
h = f.readline().rstrip('\r\n').split('\t')
if colname_existing not in h:
    print "ERROR: sample attribute file lacks a column " + colname_existing
    sys.exit(-1)
if colname_new not in h:
    print "ERROR: sample attribute file lacks a column " + colname_new
    sys.exit(-1)

idx_old = h.index(colname_existing)
idx_new = h.index(colname_new)
old2new = {}
for line in f:
    a = line.rstrip('\r\n').split('\t')
    if len(a)==1 and a[0]=="":
        break
    old2new[ a[idx_old] ] = a[idx_new]

f.close()

f = open(fn_in)
fo = open(fn_out, "w")
h = f.readline().rstrip('\r\n').split('\t')
for i, val in enumerate(h):
    if val in old2new:
        h[i] = old2new[val]

fo.write( '\t'.join(h) + '\n' )
for line in f:
    fo.write(line)

fo.close
f.close()
