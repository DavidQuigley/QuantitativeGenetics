from optparse import OptionParser
from 

# identify_homology_arms.py \
#    --file_chopchop=/notebook/human_lines_hr/oligos/design/chopchop_brca1.txt \
#    --fasta=/notebook/human_lines_hr/oligos/design/brca1_chr17_41256016_41276014.txt \
#    --fasta_bounds=chr17:41,256,016-41,276,014 \
#    --primer3_bin=/usr/bin/primer3_core \
#    --product_size_bounds=1850,2250 \


parser = OptionParser()
parser.add_option("-i", "--file_fasta", dest="fn_fasta", help="path to fasta", default="")
parser.add_option("-c", "--file_chopchop", dest="fn_cc", help="path to results from chopchop", default="")
parser.add_option("-b", "--fasta_bounds", dest="fasta_bounds", help="UCSC-formatted location defining fasta bounds", default="")
parser.add_option("-p", "--bin_primer3", dest="fn_primer3", help="location of primer3 binary", default="")
parser.add_option("-m", "--file_primer3_input", dest="fn_primer3_input", help="location of primer3 input file to generate", default="primer3_input_file")
parser.add_option("-r", "--product_size_bounds", dest="product_size_bounds", help="dash-separated low,high bound on product size", default="")

(options, args) = parser.parse_args()
fn_fasta = options.fn_fasta
fn_cc = options.fn_cc
fasta_bounds = options.fasta_bounds
fn_primer3 = options.fn_primer3
product_size_bounds = options.product_size_bounds

# extract guide locations
#
cc_results = []
f = open(fn_cc)
for line in f:
    if line[0] != "#":
        a = line.rstrip('\r\n').split('\t')
        a[1], a[2] = sequence, location
        chrom, pos = location.split(":")
        cc_results.append( [sequence, location, pos])

f.close()

# generate primer3-formatted file
# 
f = open(fn_fasta)
header = f.readline()
fasta_str = ""
for line in f:
    fasta_str += line.rstrip('\r\n')

f.close()

fo = open(fn_primer3_input, 'w')
fo.write("SEQUENCE_TEMPLATE=" + fn_fasta + '\n')
fo.write("PRIMER_PRODUCT_SIZE_RANGE=" + product_size_bounds + '\n')
fo.write("=\n")
fo.close()


# run primer3 using bounds, FASTA
# 
fn_primer3_output = fn_primer3_input + ".output"
fo = open(fn_primer3_output, 'w')
subprocess.call([fn_primer3, fn_primer3_input], stdout=fo)
fo.close()

# parse output
# 
f = open(fn_primer3_output)
idx = 0
idx_token = "_0"
process_line = False

# each result is indexed starting at 0 by primer3
# store the output of the primer matches in a hash and append to a list which can 
# then be retreived using that same index
p3_results = []
result = {}
p3_header = {}
for line in f:
    if line[0] == "=":
        if( len(result)>0 ):
            p3_results.append( result )
    else:
        token, payload = line.split("=")
        if idx_token in token:
            process_line = True
        else:
            if "_" + str(idx+1) in token:
                p3_results.append( result )
                idx += 1
                idx_token = "_" + str(idx) 
                process_line = True
        else:
            process_line = False
        if process_line:
            # PRIMER_PAIR_0_PENALTY
            # PRIMER_LEFT_0
            command, subcommand = token.split(idx_token)
            if len(subcommand)>0:
                command = command + "_" + subcommand
            result[ command ] = payload
            p3_header[command] = 1

f.close()

p3_header = sorted(list(p3_header.keys()))
print( '\t'.join( p3_header ) )
for result in p3_results:
    print( '\t'.join( [result[h] for h in p3_header] ) )
