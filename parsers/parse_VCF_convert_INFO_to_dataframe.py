from optparse import OptionParser
import os
parser = OptionParser()
parser.add_option("-i", "--file_in", dest="fn_in", help="comma-delimited list of input VCF files", default="")
parser.add_option("-s", "--sample_ids", dest="sample_ids", help="comma-delimited list of sample_IDs", default="")
parser.add_option("-o", "--file_out", dest="fn_out", help="path to file to write", default="")

(options, args) = parser.parse_args()
fn_in = options.fn_in
sample_ids = options.sample_ids
fn_out = options.fn_out

file_list = fn_in.split(',')
for fn in file_list:
    if not os.path.exists(fn):
        raise RuntimeError("file not found: " + fn )

sample_ids = sample_ids.split(",")
if len(sample_ids) != len(file_list):
    raise RuntimeError("Length of filename list not equal to length of sample ID list")

# identify columns in INFO
info_cols = []
f = open(file_list[0])
for line in f:
    if line[0] != "#":
        a = line.rstrip('\r\n').split('\t')
        format_cols = a[-2].split(":")
        break
    if "##INFO=" in line:
        ##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
        info = line.rstrip('\r\n')[:-1].split("<")[1]
        info_cols.append( info.split(",")[0].split("=")[1] )

f.close()

info_idx = {}
INFO = [] # to be used for writing out to file
for info in info_cols:
    INFO.append("")
    info_idx[info] = len(info_idx)+1

fo = open(fn_out, 'w' )
fo.write( "SAMPLE\tCHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\t" + '\t'.join(info_cols) + '\t' + '\t'.join(format_cols) + '\n' )
for sample_id, fn in zip( sample_ids, file_list):
    f = open(fn)
    for line in f:
        if line[0] == "#":
            continue
        a = line.rstrip('\r\n').split('\t')
        fo.write( sample_id + '\t' + '\t'.join( a[0:7] ) )
        for i in range(0,len(INFO)):
            INFO[i]="NA"
        for info_item in a[7].split(";"):
            try:
                k,v = info_item.split("=")
            except ValueError:
                k = info_item # flag, no value present
                v = "TRUE"
            INFO[ info_idx[k]-1 ] = v
        fo.write("\t" + "\t".join(INFO) + '\t' )
        fo.write( '\t'.join( a[-1].split(":") ) ) # format line
        fo.write("\n")
        
    f.close()

fo.close()