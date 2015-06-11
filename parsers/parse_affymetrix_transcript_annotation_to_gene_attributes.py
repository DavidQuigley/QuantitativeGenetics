from optparse import OptionParser
parser = OptionParser()
parser.add_option("-i", "--input", dest="fn_in", help="file to convert", default="")
parser.add_option("-o", "--output", dest="fn_out", help="file to write", default="")
(options, args) = parser.parse_args()

f = open(options.fn_in)
fo = open(options.fn_out, 'w')
fo.write("IDENTIFIER\tchrom\tstart\tstop\tn_probes\tsymbol\tprobe_type\tis_refseq\trefseq_id\n")
for line in f:
    if line[0] == "#":
        continue
    else:   
        a = line.rstrip('\r\n')[1:-1].split('","')
        if a[0][0] == 't':
            idx_chr = a.index("seqname")
            idx_start = a.index("start")
            idx_stop = a.index("stop")
            idx_probes = a.index("total_probes")
            idx_gene = a.index("gene_assignment")
            idx_type = a.index("category")
        
        else:
            identifier, chrom, start, stop, probes, symbol, type = a[0], a[idx_chr], a[idx_start], a[idx_stop], a[idx_probes], a[idx_gene], a[idx_type]
            chrom = chrom.replace('chr', '')
            chrom = chrom.replace('---', 'NA')
            is_refseq = "false"
            refseq_id = "NA"
            if symbol=="---":
                symbol = identifier
            else:
                assignments = symbol.split(" /// ")
                for assignment in assignments:
                    parts = assignment.split(" // ")
                    if "NM_" in parts[0]:
                        is_refseq = "true"
                        symbol = parts[1]
                        refseq_id = parts[0]
                        break
                    else:
                        symbol = parts[1]
            out = [identifier,chrom, start, stop, probes, symbol, type, is_refseq, refseq_id]
            fo.write( '\t'.join(out) + '\n')

fo.close()
