from optparse import OptionParser

parser = OptionParser()
parser.add_option("-i", "--file_in", dest="fn_in", help="path to VNC in", default="")
parser.add_option("-o", "--file_out", dest="fn_out", help="path to VNC out", default="")
(options, args) = parser.parse_args()
fn_in = options.fn_in
fn_out = options.fn_out
if fn_in==fn_out:
    raise RuntimeError("cannot write to file you're reading from")

n_kept=0
n_total=0

f = open(fn_in)
fo = open(fn_out, "w")

for line in f:
    if line[0] == "#":
        fo.write(line)
    else:
        n_total=n_total+1
        a = line.rstrip('\r\n').split('\t')
        if a[6]=="PASS":
            if "bnd" not in a[2] and "rs" not in a[2]:
                fo.write(line)
                n_kept=n_kept+1

f.close()
fo.close()

print( "Kept " + str(n_kept) + " of " + str(n_total) + " variants, " + str( round( float(n_kept)/float(n_total)*100.0, 1 ) ) + "%" )