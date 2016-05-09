import argparse
parser = argparse.ArgumentParser(description='sort and split fastq files')
parser.add_argument('--in', dest="fn_in", help='input file')
parser.add_argument('--out', dest="fn_out", help='output file stub')

args = parser.parse_args()
fn_in = args.fn_in
fn_out = args.fn_out

f = open(fn_in)
idx = 1
reads_1 = {}
reads_2 = {}
while True:
    l1 = f.readline()
    if l1=="":
        break
    l2 = f.readline()
    l3 = f.readline()
    l4 = f.readline()
    a = l1.rstrip('\r\n').split("/")
    if a[1]=="1":
        reads_1[ a[0] ] = l1 + l2 + l3 + l4
        
    else:
        reads_2[ a[0] ] = l1 + l2 + l3 + l4

print("read " + str(len(reads_1)) + " reads with key=1")
print("read " + str(len(reads_2)) + " reads with key=2")

fo1 = open( fn_out + '_1.fq', 'w')
fo2 = open(fn_out + '_2.fq', 'w')
nomatch=0
for key in sorted(reads_1.keys()):
    try:
        fo1.write(  reads_1[key] )
        fo2.write( reads_2[key] )
    
    except KeyError:
        nomatch += 1

fo1.close()
fo2.close()