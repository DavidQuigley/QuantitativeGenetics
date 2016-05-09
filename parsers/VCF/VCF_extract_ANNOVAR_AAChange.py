from optparse import OptionParser

parser = OptionParser()
parser.add_option("-i", "--file_in", dest="fn_in", help="path to output of Annovar protein_changes", default="")
parser.add_option("-f", "--filter", dest="filter", help="filter requirement", default="")
parser.add_option("-r", "--refseq", dest="refseq", help="refseq to extract from AAChange.refGene", default="")
parser.add_option("-s", "--sample_id", dest="sample_id", help="sample_id to write", default="")
(options, args) = parser.parse_args()
fn_in = options.fn_in
refseq = options.refseq
required_filter = options.filter

f = open(fn_in)

for line in f:
    if line[0] != "#":
        a = line.rstrip('\r\n').split('\t')
        filter = a[6]
        info = a[7]
        if required_filter != "":
            if filter==required_filter:
                info_tokens = info.split(";")
                for info_token in info_tokens:
                    if info_token.split("=")[0] == "AAChange.refGene":
                        refseq_tokens = info_token.split("=")[1]
                        for refseq_token in refseq_tokens.split(","):
                            #TP53:NM_000546:exon4:c.C215G
                            if refseq in refseq_token:
                                print( sample_id'\t' + '\t'.join( a[0:7] ) + '\t' + '\t'.join( refseq_token.split(":") ) )

                
    