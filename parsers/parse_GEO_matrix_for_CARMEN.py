# Extracts expression data and sample attributes from GEO Matrix files
# Files can be gzipped or unzipped.

import sys
import gzip
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-f", "--file_in", dest="fn_in", help="Comma-delimited Path to input matrix files", default="")
parser.add_option("-o", "--path_out", dest="path_out", help="Path to write output files", default="")

(options, args) = parser.parse_args()
fn_in = options.fn_in
path_out = options.path_out
if path_out !="":
    path_out = pathout + '/'

if fn_in=="":
    print "REQUIRED: -f parameter specifying the name(s) of the matrix files"
    sys.exit(0)

file_idx = 1
fn_expr = []
s2prop2v = {}
probe2v = {}
identifiers = []

for fn in fn_in.split(","):    
    in_expression_data = False
    print "opening " + fn
    if( fn.count(".gz") > 0 ):
        f = gzip.open(fn, 'rb')
    
    else:
        f = open(fn, 'r')
    
    for line in f:
        a = line.replace('"', '').rstrip('\n').split("\t")
        if in_expression_data:
            if line[0] != '!':
                if a[0] in probe2v:
                    probe2v[a[0]] += '\t' + '\t'.join( a[1:])
                else:
                    probe2v[a[0]] =  '\t'.join( a[1:])
        else:
            if line[0] == "!":  
                if a[0] == "!Sample_geo_accession":
                    accessions = a[1:]
                    for accession in accessions:
                        s2prop2v[accession] = {}
                    
                elif a[0] == "!Sample_characteristics_ch1" or a[0] == "!Sample_characteristics_ch2":
                    # GEO is not very consistent in how it encodes this field.
                    # the value of a column here may be blank so we need to check all 
                    # of the values to find the first case where the name of this characteristic
                    # is listed; store that as property
                    # The value may also be "property:" lacking a space.
                    for x in a[1:]:
                        if x.count(": ")>0:
                            property = x.split(": ")[0]        
                            break
                    property = property.replace(" ", "_")
                    print "Read property: " + property
                    values = []
                    for x in a[1:]:
                        value = "NA"
                        if x != "":
                            x = x.split(":")
                            if len(x)==2 and x[1] != "":
                                value= x[1].lstrip()
                        values.append(value)
                    for accession, value in zip(accessions, values):
                        s2prop2v[accession][property] = value
                    
            elif line[0] == '"':
                # header for expression data
                identifiers_this_file = line.replace('"', '').rstrip("\n").split('\t')[1:]
                for identifier in identifiers_this_file:
                    identifiers.append(identifier)
                in_expression_data=True
                file_idx += 1

f.close()

f_expr = open(path_out + "expr_from_matrix.txt", "w")
f_expr.write("IDENTIFIER\t" + '\t'.join(identifiers) + '\n' )
for probe in sorted(probe2v.keys()):
    f_expr.write( probe + "\t" + probe2v[probe] + '\n' )

f_expr.close()

properties = sorted(s2prop2v[accessions[0]].keys())
f_sa = open(path_out + "sample_attributes_from_matrix.txt", "w")
f_sa.write( 'IDENTIFIER' + '\t' + '\t'.join( properties ) + '\n')
for accession in identifiers:
    f_sa.write(accession)
    for property in properties:
        f_sa.write("\t"+s2prop2v[accession][property] )
    f_sa.write("\n")

f_sa.close()
