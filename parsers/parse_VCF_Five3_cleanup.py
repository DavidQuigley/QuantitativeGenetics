from optparse import OptionParser
parser = OptionParser()
parser.add_option("-i", "--file_in", dest="fn_in", help="path to VNC in", default="")
parser.add_option("-c", "--file_contigs", dest="fn_contigs", help="path to contig file", default="")
parser.add_option("-o", "--file_out", dest="fn_out", help="path to VNC out", default="")
(options, args) = parser.parse_args()
fn_in = options.fn_in
fn_out = options.fn_out
fn_contigs = options.fn_contigs
if fn_in==fn_out:
    raise RuntimeError("cannot write to file you're reading from")

n_kept=0
n_total=0
f = open(fn_in)
fo = open(fn_out, "w")
contigs_written=False
seen_PC=False
seen_VT=False
seen_NP=False
for line in f:
    if line[0] == "#":
        if not contigs_written and "INFO" in line:
            fc = open(fn_contigs)
            for lc in fc:
                fo.write(lc)
            fc.close()
            contigs_written=True
        if "ID=PC" in line:
            seen_PC = True
        if "ID=VT" in line:
            seen_VC = True
        if not seen_PC and "FILTER" in line:
            fo.write( "##INFO=<ID=PC,Number=1,Type=Float,Description=\"PhastCons conservation score\">\n")
            seen_PC = True
        if not seen_VT and "FILTER" in line:
            fo.write( "##INFO=<ID=VT,Number=1,Type=String,Description=\"Variant type from Five3 call\">\n")
            seen_VT = True
        if not seen_NP and "FILTER" in line:
            fo.write( "##INFO=<ID=NP,Number=1,Type=String,Description=\"No clue\">\n")
            seen_NP = True
        
        fo.write(line)
    
    else:
        n_total=n_total+1
        a = line.rstrip('\r\n').split('\t')
        if a[6]=="PASS":
            if "bnd" not in a[2] and "rs" not in a[2]:
                INFO_tumor = a[10].split(":")
                INFO_tumor[0] = "0/1"
                a[10] = ":".join(INFO_tumor)
                INFO_normal = a[9].split(":")
                INFO_normal[0] = "0/1"
                a[9] = ":".join( INFO_normal )
                line = "\t".join( a )
                fo.write(line + '\n')
                n_kept=n_kept+1

f.close()
fo.close()

print( "Kept " + str(n_kept) + " of " + str(n_total) + " variants, " + str( round( float(n_kept)/float(n_total)*100.0, 1 ) ) + "%" )
