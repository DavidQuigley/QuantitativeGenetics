# The sequence listed in Martin's CRISPRi library can be used to create a bowtie library.
# Skip the first nucleotide from the protospacer (which is always G anyway) and add in the 
# first six nucleotides from the constant region.

def RC(x):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
    return ''.join( [complement[base] for base in list( x[::-1].upper() ) ] )

from optparse import OptionParser
parser = OptionParser()
parser.add_option("-i", "--file_in", dest="fn_in", help="path to VNC in", default="")
parser.add_option("-o", "--file_out", dest="fn_out", help="path to VNC out", default="")
parser.add_option("-v", "--halves", dest="halves", help="sublibrary halves (comma-delimited)", default="")
parser.add_option("-s", "--sublibraries", dest="sublibraries", help="sublibraries (comma-delimited)", default="")

(options, args) = parser.parse_args()
fn_in = options.fn_in
fn_out = options.fn_out
requested_halves = options.halves.split(',')
requested_sublibraries = options.sublibraries.split(',')

print "selecting Halves: " + ','.join(requested_halves)
print "selecting Sublibraries: " + ','.join(requested_sublibraries)

f = open(fn_in)
fo = open(fn_out, 'w')

a = f.readline().rstrip('\r\n').split('\t')
idx_proto = a.index("protospacer_sequence")
idx_half = a.index("sublibrary_half")
idx_sub = a.index("sublibrary")
for line in f:
    a = line.rstrip('\r\n').split('\t')
    has_half = False
    has_sublib = False
    for half in a[idx_half].split(','):
        if half in requested_halves:
            has_half=True
    
    for sublib in a[idx_sub].split(','):
        if sublib in requested_sublibraries:
            has_sublib=True
    
    if has_half and has_sublib:
        fo.write(">" + a[0] + "\n" + RC( a[idx_proto][1:20] + "GTTTAA") + "\n")

fo.close()
f.close()
