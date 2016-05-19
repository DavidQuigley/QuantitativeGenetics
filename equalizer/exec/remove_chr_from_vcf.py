from optparse import OptionParser

parser = OptionParser()
parser.add_option("-f", "--file_in", dest="fn_in", help="tab-delimited path to file from TCGA firehose", default="")
parser.add_option("-o", "--file_out", dest="fn_out", help="Path to write output file", default="")
(options, args) = parser.parse_args()
fn_in = options.fn_in
fn_out = options.fn_out

f = open(fn_in)
fo = open(fn_out, "w")

# remove the chr from the first line if it's there
for line in f:
    if len(line)>3:
        if line[0:3]=='chr':
            fo.write( line[3:] )

fo.close()
