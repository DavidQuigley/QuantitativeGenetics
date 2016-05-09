from optparse import OptionParser

def extract_nmers( seq_wt, seq_mut, protein_index, nmer_size ):
    # identify all nmers across seq_wt and seq_mut of size nmer_size that include the 
    # nucleotide at position protein_index. require both wt and mut has protein sequence,
    # which could potentially not be the case for deletions where protein index is near 
    # the end of a string
    # returns array of arrays [first nucleotide, nmer_wt, nmer_mut] 
    nmers=[]
    idx = int(protein_index)-nmer_size
    for i in range(1,nmer_size+1):
        if idx >= 0 and idx+nmer_size < len(seq_wt) and idx+nmer_size < len(seq_mut):
            nmers.append( [idx, seq_wt[ idx:idx+nmer_size ], seq_mut[ idx:idx+nmer_size ] ] )
        idx += 1
    
    return nmers

    
parser = OptionParser()
parser.add_option("-i", "--file_in", dest="fn_in", help="path to output of Annovar protein_changes", default="")
parser.add_option("-v", "--file_variant_func", dest="fn_vf", help="path to output of Annovar variant function prediction", default="")
parser.add_option("-g", "--file_genemodel", dest="fn_gene", help="path to gene model file", default="")
parser.add_option("-n", "--nmers", dest="nmer_size", help="sequence length to generate, default 9", default="9")
parser.add_option("-w", "--file_out_wt", dest="fn_out_wt", help="path to wildtype fasta out", default="")
parser.add_option("-m", "--file_out_mut", dest="fn_out_mut", help="path to mutant fasta out", default="")
(options, args) = parser.parse_args()
fn_in = options.fn_in
fn_vf = options.fn_vf
fn_out_wt = options.fn_out_wt
fn_out_mut = options.fn_out_mut
fn_gene = options.fn_gene
nmer_size = int(options.nmer_size)
if fn_in==fn_out_wt or fn_in==fn_out_mut:
    raise RuntimeError("cannot write to file you're reading from")
if nmer_size<8 or nmer_size>14:
    raise RuntimeError("nmer_size must be an integer between 8 and 14, inclusive")

# To allow tracing a neoepitope back to the originating VCF, use the "line###" in 
# ${OUT_BASE}_protein_changes.txt to look up columns 4 and 5 of 
# ${OUT_BASE}.refGene.exonic_variant_function and add this information into the FASTA that 
# will be generated

# Read in variant function file to make line_number -> chr_pos.
f = open(fn_vf)
line2chrpos = {}
for line in f:
    a = line.rstrip('\r\n').split('\t')
    line2chrpos[ a[0] ] = a[3] + "_" + a[4]

f.close()

rs2desc = {}
f = open(fn_gene)
for line in f:
    a = line.rstrip('\r\n').split('\t')
    rs2desc[ a[1] ] = a[12] + ":" + a[2] + ":" + a[4] + ":" + a[5]

f.close()

f = open(fn_in)
fo_wt = open(fn_out_wt, "w")
fo_mut = open(fn_out_mut, "w")
seq_wt = ""
seq_mut = ""
cur_seq_type = ""
running_idx=1
line_number=""
for line in f:
    line = line.rstrip('\r\n')
    if len(line)==0:
        continue
    if line[-1]=="*": # Stop codon
        if cur_seq_type=="wildtype":
            # Finished reading wildtype sequence, next will be mutant header
            seq_wt += line[0:-1] # append last line, excluding stop
            cur_seq_type="mutant"
            
        elif cur_seq_type == "mutant":
            seq_mut += line[0:-1] # append last line, excluding stop
            # completed reading both WT and mutant
            if not is_synonymous:
                print header
                nmers = extract_nmers( seq_wt, seq_mut, protein_index, nmer_size )
                if len(nmers)>0:
                    for idx, epitope_wt, epitope_mut in nmers:
                        fo_wt.write( ">idx_" + str(running_idx) + ":" + refseq + "_wt " + rs2desc[refseq] + ":p" + str(protein_index) + ":" + line2chrpos[line_number] + "\n")
                        fo_wt.write( epitope_wt + "\n")
                        fo_mut.write( ">idx_" +str(running_idx) + ":" + refseq + "_mut " + rs2desc[refseq] + ":p" + str(protein_index) + ":" + line2chrpos[line_number] + "\n")
                        fo_mut.write( epitope_mut  + "\n")
                        running_idx += 1
            seq_wt = ""
            seq_mut = ""  
    
    elif line[0]==">":
        if "WILDTYPE" in line:
            cur_seq_type="wildtype"
            refseq = line.split(" ")[1]
            line_number = line.split(" ")[0][1:]
        else:
            #>line10 NM_004446 c.C2525G protein-altering  (position 842 changed from A to G)
            header = line[ line.index("NM_") : ]
            line_number = line.split(" ")[0][1:]
            if header.split(" ")[2] == "synonymous":
                is_synonymous=True
                protein_index=-1
            else:
                is_synonymous=False
                protein_index = header[ header.index("position "):].split(" ")[1]
    
    else:
        if cur_seq_type=="wildtype":
            seq_wt += line
        else:
            seq_mut += line
fo_wt.close()
fo_mut.close()
