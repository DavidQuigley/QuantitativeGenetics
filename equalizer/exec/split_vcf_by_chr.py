from optparse import OptionParser

parser = OptionParser()
parser.add_option("-f", "--file_in", dest="fn_in", help="Path to input VCF file", default="")
parser.add_option("-o", "--base_out", dest="fn_out", help="Path to base name of output VCF files", default="")
(options, args) = parser.parse_args()
fn_in = options.fn_in
fn_out = options.fn_out

header = ""
f = open(fn_in)
cur_chr = ""
for line in f:
    if line[0] == "#":
        header += line
    
    else:
        col1 = line.split("\t")[0]
        if col1 == cur_chr:
            fo.write( line )
        else:            
            if not cur_chr == "":
                fo.close()
            cur_chr = col1
            print( "Began writing chromosome " + cur_chr )
            fo = open( fn_out + "_" + cur_chr + ".vcf", "w" )
            fo.write(header)
            fo.write(line)

fo.close()
