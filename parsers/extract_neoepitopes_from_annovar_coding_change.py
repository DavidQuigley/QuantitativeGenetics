from optparse import OptionParser

def extract_nmers( seq_wt, seq_mut, protein_index, nmer_size ):
    # identify all nmers across seq_wt and seq_mut of size nmer_size that include the 
    # nucleotide at position protein_index. require both wt and mut has protein sequence,
    # which could potentially not be the case for deletions where protein index is near 
    # the end of a string
    # returns array of arrays [first nucleotide, nmer_wt, nmer_mut] 
    nmers=[]
    idx = protein_index-nmer_size
    for i in range(1,nmer_size+1):
        print idx + nmer_size
        if idx >= 0 and idx+nmer_size < len(seq_wt) and idx+nmer_size < len(seq_mut):
            nmers.append( [idx, seq_wt[ idx:idx+nmer_size ], seq_mut[ idx:idx+nmer_size ] ] )
        
        idx += 1
    
    return nmers

    
parser = OptionParser()
parser.add_option("-i", "--file_in", dest="fn_in", help="path to output of Annovar protein_changes", default="")
parser.add_option("-g", "--file_genemodel", dest="fn_gene", help="path to gene model file", default="")
parser.add_option("-n", "--nmers", dest="nmer_size", help="sequence length to generate, default 9", default="9")
parser.add_option("-w", "--file_out_wt", dest="fn_out_wt", help="path to wildtype fasta out", default="")
parser.add_option("-m", "--file_out_mut", dest="fn_out_mut", help="path to mutant fasta out", default="")
(options, args) = parser.parse_args()
fn_in = options.fn_in
fn_out_wt = options.fn_out_wt
fn_out_mut = options.fn_out_mut
fn_gene = options.fn_gene
nmer_size = int(options.nmer_size)
if fn_in==fn_out:
    raise RuntimeError("cannot write to file you're reading from")
if nmer_size<8 or nmer_size>14:
    raise RuntimeError("nmer_size must be an integer between 8 and 14, inclusive")

rs2desc = {}
f = open(fn_gene)
for line in f:
    a = line.rstrip('\r\n').split('\t')
    rs2desc[ a[1] ] = a[12] + ":" + a[2] + ":" + a[4] + ":" + a[5] + a[3]

f.close()

f = open(fn_in)
fo_wt = open(fn_out_wt, "w")
fo_mut = open(fn_out_mut, "w")
seq_wt = ""
seq_mut = ""
cur_seq_type = ""

for line in f:
    line = line.rstrip('\r\n')
    if line[-1]=="*": # Stop codon
        seq += line[0:-1] # append last line, excluding stop
        if cur_seq_type == "mutant":
            # completed reading both WT and mutant
            if not is_synonymous:
                nmers = extract_nmers( seq_wt, seq_mut, protein_index, nmer_size )
                if len(nmers_)>0:
                    for idx, epitope_wt, epitope_mut in nmers:
                        fo_wt.write( ">" + refseq + "_wt " + rs2desc[refseq] + ":" + str(protein_index) + "\n")
                        fo_wt.write( epitope_wt "\n")
                        fo_mut.write( ">" + refseq + "_mut " + rs2desc[refseq] + ":" + str(protein_index) + "\n")
                        fo_mut.write( epitope_mut "\n")                        
        else:
            # Finished reading wildtype sequence, next will be mutant
            cur_seq_type="mutant"
    
    elif line[0]==">":
        if "WILDTYPE" in line:
            cur_seq_type="wildtype"
            refseq = line.split(" ")[1]
        else:
            #>line10 NM_004446 c.C2525G protein-altering  (position 842 changed from A to G)
            header = line[ line.index("NM_") : ]
            if header.split(" ")[2] == "synonymous":
                is_synonymous=True
                protein_index=-1
            else:
                is_synonymous=False
                protein_index = header[ header.index("position "):].split(" ")[1]

fo_wt.close()
fo_mut.close()
