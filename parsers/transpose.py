from optparse import OptionParser

parser = OptionParser()
parser.add_option("-f", "--file_in", dest="fn_in", help="tab-delimited path to file from TCGA firehose", default="")
parser.add_option("-o", "--file_out", dest="fn_out", help="Path to write output file", default="")

(options, args) = parser.parse_args()
fn_in = options.fn_in
fn_out = options.fn_out

f = open(fn_in)
lines = []
for line in f:
    lines.append( line.rstrip('\r\n').split("\t") )

fo = open(fn_out, 'w')
for x in zip(*lines):
  fo.write( '\t'.join(x) + '\n')

fo.close()
