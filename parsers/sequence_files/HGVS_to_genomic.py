import pyximport
pyximport.install()
import pyhgvs as hgvs
import pyhgvs.utils as hgvs_utils 
import pygr
from pygr.seqdb import SequenceFileDB  
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-i", "--input_string", dest="instr", help="format NM_007294.3:c.8T>G", default="")
(options, args) = parser.parse_args()
instr = options.instr

# Read genome sequence using pygr. 
fn_hg19='/raw/human_sequence_reference/UCSC/hg19.fa'
fn_refgene='/home/david/software/hgvs-master/pyhgvs/data/genes.refGene'
genome = SequenceFileDB(fn_hg19)
def get_transcript(name):
    return transcripts.get(name)

with open(fn_refgene) as infile:
    transcripts = hgvs_utils.read_transcripts(infile)

out = hgvs.parse_hgvs_name(instr,genome,get_transcript=get_transcript)  
print( "\t".join( [str(x) for x in out] ) )



