import vcf
from optparse import OptionParser
parser = OptionParser()
parser.add_option("-i", "--file_in", dest="fn_in", help="path to vcf", default="")
parser.add_option("-s", "--symbol_list", dest="symbol_restriction", help="restrict to these symbols, comma-delimited", default="")
parser.add_option("-f", "--filter_list", dest="filters", help="filters, comma-delimited, required to pass all", default="")
parser.add_option("-o", "--file_out", dest="fn_out", help="output vcf file to write", default="")
(options, args) = parser.parse_args()
fn_in = options.fn_in
fn_out = options.fn_out
filters = options.filters.split(",")
symbol_restriction=options.symbol_restriction

if symbol_restriction=="":
    symbol_restriction=[]
else:
    symbol_restriction = set( symbol_restriction.split(",") )

known_filters = ["", "cosmic", "nonsense_or_pathogenic_missense", "pass"]
for filter in filters:
    if not filter in known_filters:
        raise RuntimeError("unknown filter")

vcf_reader = vcf.Reader(open(fn_in, 'r'))
records = []
for record in vcf_reader:
    records.append( record )
    if len(records) % 10000 == 0:
        print( str(len(records)) )

def filter_nonsense(record):
    con = record.INFO["ExonicFunc.refGene"][0]
    if con is None:
        return False
    
    if not con is None:
        return( con == "stopgain" or con=="frameshift_insertion" or con=="frameshift_deletion" )

def filter_pass(record):
    if not record.FILTER is None:
        return len(record.FILTER)==0
    
    return False
    
def filter_missense(record):
    con = record.INFO["ExonicFunc.refGene"][0]
    if con is None:
        return False
    
    if not con is None:
        return( con=="nonsynonymous_SNV" )

def filter_cosmic(record):
    return not record.INFO["cosmic70"][0] is None

def filter_pathogenic(record):
    CLINSIG = record.INFO["CLINSIG"][0]
    if( CLINSIG is None ):
        return False
    
    return "pathogenic" in CLINSIG.lower()


vcf_writer = vcf.Writer(open(fn_out, 'w'), vcf_reader)
idx=1
for record in records:
    idx += 1
    if idx % 10000== 0:
        print idx

    keep = True
    if len(symbol_restriction)>0:
        symbol = record.INFO["Gene.refGene"][0]
        keep = not symbol is None and symbol in symbol_restriction
    
    for filter in filters:
        if keep and filter=="pass":
            keep = filter_pass(record)
        
        elif keep and filter=="cosmic":
            keep = filter_cosmic(record)
        
        elif keep and filter=="nonsense_or_pathogenic_missense":
            keep = filter_nonsense(record) or ( filter_missense(record) and filter_pathogenic(record) )
        
    if keep:
        vcf_writer.write_record( record )

vcf_writer.close()
