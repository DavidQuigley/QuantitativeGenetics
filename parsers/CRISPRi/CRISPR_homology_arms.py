# INSTRUCTIONS
# Before using:
# Install the UCSC stack (for in-silico PCR binary isPcr) and Primer3 (for primer pairs)
# Install a 2bit formatted version of that reference (used by isPcr)
# Install samtools (to extract FASTA-formatted pieces of the genome)
# Install and index a genome reference (used by samtools, e.g. hg19)
#
# To use:
# Identify a list of CRISPR guide positions using CHOPCHOP or another tool
# 
# The simplest usage will find PCR primer pairs that can amplify a region of size bounded 
# by --product_size_bounds no more than --distance_to_cut bases from a guide site. 
# Indicate where the cut site is relative to the nearby arm with --cut_position, either 
# left or right. To exclude from considerating any sequence that bears a restriction 
# enzyme binding sequence, pass one or more sequence strings to --restriction_sequences.. 
# If you pass a list of distant_primer_bounds pairs, the code will attempt to place a 
# second homology  arm where the primer nearer to the cut site has its inner endpoint 
# within the bounds.
#
# .genome......\---   distant primers  ___/...............X.\---  near primers  ___/....
#  |                                        |              --
#  |  distant primer bounds                 |          cutsite    
#
# Example call:
# python /home/ubuntu/scripts/CRISPR_homology_arms.py \
#  --gene_symbol=BRCA1_plasmid1_rightarm \
#  --fn_samtools=/opt/samtools-1.2/samtools \
#  --file_guides=/home/ubuntu/temp/chopchop_BRCA1_41256016_41276014.txt \
#  --fn_reference=/datasets/human_sequence_hg19/ucsc.hg19.fasta \
#  --bin_primer3=/usr/bin/primer3_core \
#  --bin_ispcr=/home/ubuntu/bin/x86_64/isPcr \
#  --file_ref_twobit=/datasets/human_sequence_hg19/hg19.2bit \
#  --cut_position=left \
#  --distance_to_cut=100 \
#  --dir_out=/home/ubuntu/temp \
#  --product_size_bounds=1850-2250 \
#  --restriction_sequences=GAAGAC,CTTCTG

import os.path
import subprocess
from optparse import OptionParser

def parse_primer3_output(fn_primer3_output):
    # return a vector of primer attributes (values in a hash)
    f = open(fn_primer3_output)
    idx = 0
    idx_token = "_0"
    process_line = False # skip header lines

    # each result is indexed starting at 0 by primer3
    # store the output of the primer matches in a hash and append to a list which can 
    # then be retreived using that same index
    p3_results = []
    result = {}
    for line in f:
        if line[0] == "=":
            if( len(result)>0 ):
                p3_results.append( result )
        else:
            token, payload = line.rstrip('\r\n').split("=")
            if idx_token in token:
                process_line = True
            elif "_" + str(idx+1) in token:
                p3_results.append( result )
                idx += 1
                idx_token = "_" + str(idx) 
                process_line = True
                result = {}
            else:
                process_line = False
            if process_line:
                result["PRIMER3_IDX"] = str(idx)
                # PRIMER_PAIR_0_PENALTY
                # PRIMER_LEFT_0
                command, subcommand = token.split(idx_token)
                result[ command + subcommand ] = payload
    
    f.close()
    return p3_results

parser = OptionParser()
parser.add_option("-g", "--gene_symbol", dest="symbol", help="symbol identifier for filenames", default="test")
parser.add_option("-s", "--fn_samtools", dest="fn_samtools", help="path to samtools", default="/opt/samtools-1.2/samtools")
parser.add_option("-c", "--file_guides", dest="fn_cc", help="path to guide file, with sequence, location in columns 2,3", default="")
parser.add_option("-e", "--restriction_sequences", dest="res_seq", help="comma-delimited list of restriction sequences to avoid", default="/datasets/human_sequence_hg19/ucsc.hg19.fasta")
parser.add_option("-r", "--fn_reference", dest="fn_ref", help="path to indexed reference genome", default="/datasets/human_sequence_hg19/ucsc.hg19.fasta")
parser.add_option("-p", "--bin_primer3", dest="fn_primer3", help="location of primer3 binary", default="")
parser.add_option("-i", "--bin_ispcr", dest="fn_ispcr", help="location of UCSC in-silico PCR binary", default="")
parser.add_option("-t", "--file_ref_twobit", dest="fn_ref_twobit", help="location of reference in twobit format", default="")
parser.add_option("-u", "--cut_position", dest="cut_position", help="cut position with respect to homology region, {left,right}", default="/")
parser.add_option("-d", "--distance_to_cut", dest="distance_to_cut", help="maximum distance to cut, default 100", default="100")
parser.add_option("-o", "--dir_out", dest="dir_out", help="path to output files", default="")
parser.add_option("-b", "--product_size_bounds", dest="product_size_bounds", help="dash-separated low,high bound on product size", default="")

(options, args) = parser.parse_args()
fn_samtools = options.fn_samtools
fn_ref = options.fn_ref
fn_cc = options.fn_cc
fn_primer3 = options.fn_primer3
fn_ispcr = options.fn_ispcr
fn_ref_twobit = options.fn_ref_twobit
dir_out = options.dir_out
symbol=options.symbol
max_distance_to_cut = options.distance_to_cut
product_size_bounds = options.product_size_bounds
cut_position=options.cut_position
res_seq = options.res_seq

if cut_position != "left" and cut_position !="right":
    raise RuntimeError("Cut position must be either left or right")
if not os.path.exists(dir_out):
    raise RuntimeError("Output folder does not exist:" + dir_out)
if not os.path.exists(fn_primer3):
    raise RuntimeError("Primer3 binary does not exist:" + fn_primer3)
if not os.path.exists(fn_ref):
    raise RuntimeError("Reference genome file (indexed FASTA) does not exist:" + fn_ref)
if not os.path.exists(fn_samtools):
    raise RuntimeError("Samtools binary file does not exist:" + fn_samtools)
if not os.path.exists(fn_ispcr):
    raise RuntimeError("UCSC in silico PCR binary does not exist:" + fn_ispcr)
if not os.path.exists(fn_ref_twobit):
    raise RuntimeError("Reference genome file (2bit format) does not exist:" + fn_ref_twobit)

# split off FASTA sequences to avoid amplifying any restriction sequence that we're going
# to use downstream, e.g. BbsI (GAAGAC)
restriction_sequences = res_seq.split(',')

min_productsize, max_productsize = product_size_bounds.split("-")
min_productsize=int(min_productsize)
max_productsize=int(max_productsize)

fasta_file = dir_out + '/' + symbol + "_search_region.fasta"
fn_isPcr_input = dir_out + '/' + symbol + "_isPcr_input"
fn_isPcr_output = dir_out + '/' + symbol + "_isPcr_output.bed"
fn_primer3_input = dir_out + '/' + symbol + "primer3_input"
fn_primer3_output = dir_out + '/' + symbol + "primer3_output"
fn_out = dir_out + '/' + symbol + '_homology_results.txt'

print "# MESSAGE: reading guide positions from " + fn_cc 
print "# MESSAGE: calling primer3 at " + fn_primer3
print "# MESSAGE: Avoiding restriction sites: " + res_seq

# extract guide locations and write FASTA files for corresponding search regions
#
f_isPcr = open(fn_isPcr_input, "w")
cc_results = []
f = open(fn_cc)
idx = 0
for line in f:
    if line[0] != "#":
        a = line.rstrip('\r\n').split('\t')
        sequence, location = a[1], a[2]
        chrom, pos = location.split(":")
        cc_result = {}
        cc_result["GUIDE_IDX"] = str(idx)
        cc_result["GUIDE_SEQUENCE"] = sequence
        cc_result["GUIDE_LOCATION"] = location
        if cut_position=="left":
            fasta_region = chrom + ":" + str(int(pos)+1) + "-" + str(int(pos)+1+max_productsize)
            cc_result["PRIMER_LOCATION_RESTRICTION"] = "1," + max_distance_to_cut + ",,"
        else:
            fasta_region = chrom + ":" + str(int(pos)-1-max_productsize) + "-" + str(int(pos)-1)
            cc_result["PRIMER_LOCATION_RESTRICTION"] = ",," + str(max_productsize-int(max_distance_to_cut)) + "," + max_distance_to_cut
        fo = open(fasta_file, 'w')
        subprocess.call([fn_samtools, "faidx", fn_ref, fasta_region], stdout=fo)
        fo.close()
        fo = open( fasta_file )
        h = fo.readline()
        cc_result["SEARCH_SEQUENCE"] = ""
        for line in fo:
            cc_result["SEARCH_SEQUENCE"]=cc_result["SEARCH_SEQUENCE"]+line.rstrip('\r\n')
        
        # This isn't perfect as it could in theory create new restriction enzyme binding sites 
        # when the sequence is split, but we'll live with that for now.
        cc_result["SEARCH_SEQUENCE"] = cc_result["SEARCH_SEQUENCE"].upper()
        
        for rs in restriction_sequences:
            rs = rs.upper()    
            if rs in cc_result["SEARCH_SEQUENCE"]:
                if cut_position=="right":
                    cc_result["SEARCH_SEQUENCE"] = cc_result["SEARCH_SEQUENCE"].split( rs )[-1] # right-most piece of sequence, closest to cut
                else:
                    cc_result["SEARCH_SEQUENCE"] = cc_result["SEARCH_SEQUENCE"].split( rs )[0] # left-most piece of sequence, closest to cut
        
        
        # Write out input file for Primer3 and call Primer3
        fo = open( fn_primer3_input, 'w' )
        fo.write( "SEQUENCE_TEMPLATE=" + cc_result["SEARCH_SEQUENCE"] + '\n' )
        fo.write( "SEQUENCE_PRIMER_PAIR_OK_REGION_LIST=" + cc_result["PRIMER_LOCATION_RESTRICTION"] + '\n' )
        fo.write( "PRIMER_PRODUCT_SIZE_RANGE=" + product_size_bounds + '\n' )
        fo.write( "=\n" )
        fo.close()
        fo = open(fn_primer3_output, 'w')
        subprocess.call([fn_primer3, fn_primer3_input], stdout=fo)
        fo.close()
        
        cc_result["PRIMER3_RESULTS"] = parse_primer3_output( fn_primer3_output )

        for result in cc_result["PRIMER3_RESULTS"]:
            result_identifier = cc_result["GUIDE_IDX"] + "_" + result["PRIMER3_IDX"]
            f_isPcr.write( '\t'.join( [result_identifier, result["PRIMER_LEFT_SEQUENCE"], result["PRIMER_RIGHT_SEQUENCE"] ] ) + '\n')
        
        cc_results.append( cc_result )
        idx += 1

f.close()

f_isPcr.close()

print "# MESSAGE: Read " + str(len(cc_results)) + " guides"
print "# MESSAGE: in silico PCR input written to " + fn_isPcr_input
print "# MESSAGE: Calling in-silico PCR at " + fn_ispcr
print "# MESSAGE: Using reference: " + fn_ref_twobit

# call in-silico PCR using fn_isPcr_input, write BED-formatted output to fn_isPcr_output
# mark up cc_results with number of matches for each pair and then write final results
subprocess.call( [fn_ispcr, "-out=bed", fn_ref_twobit, fn_isPcr_input, fn_isPcr_output] )

print "# MESSAGE: Reading in-silico PCR output from:" + fn_isPcr_output
f = open(fn_isPcr_output)
for line in f:
    chrom, loc_begin, loc_end, result_identifier = line.split("\t")[0:4]
    guide_id, primer_id = result_identifier.split("_")
    location_length = chrom + ":" + loc_begin + "-" + loc_end + "_" + str(int(loc_end)-int(loc_begin))
    primer3_result = cc_results[int(guide_id)]["PRIMER3_RESULTS"][int(primer_id)]
    try:
        primer3_result["IN_SILICO_RESULTS"].append( location_length )
    except KeyError:
        primer3_result["IN_SILICO_RESULTS"] = [location_length]

# Write out final results

primer3_header = ["PRIMER_LEFT_SEQUENCE", "PRIMER_LEFT_TM", "PRIMER_LEFT", "PRIMER_RIGHT_SEQUENCE", "PRIMER_RIGHT_TM", "PRIMER_RIGHT", "PRIMER_PAIR_PRODUCT_SIZE"]

print "# MESSAGE: Wrote output to: " + fn_out
fo = open(fn_out, 'w')
for cc_result in cc_results:

    if cc_result["GUIDE_IDX"]=="0":
        fo.write( "\t".join( ["GUIDE_IDX", "GUIDE_SEQUENCE", "GUIDE_LOCATION"] ) )
        fo.write( '\t' + "\t".join(primer3_header) )
        fo.write( '\tNUMBER_IN_SILICO_MATCHES\tIN_SILICO_MATCH_DETAILS\n' )
    
    for primer3_result in cc_result["PRIMER3_RESULTS"]:
        fo.write( '\t'.join( [cc_result["GUIDE_IDX"], cc_result["GUIDE_SEQUENCE"], cc_result["GUIDE_LOCATION"] ]) )
        fo.write( '\t' + '\t'.join( [primer3_result[token] for token in primer3_header] )  )
        try:
            fo.write( '\t' + str( len(primer3_result["IN_SILICO_RESULTS"]) ) )
            fo.write( '\t' + ';'.join( primer3_result["IN_SILICO_RESULTS"] ) + '\n' )
        except KeyError:
            fo.write( '\t0\tno_PCR_results\n')

fo.close()

