import argparse
parser = argparse.ArgumentParser(description='convert sanger copy number data to BED files')
parser.add_argument('dictionary', help='cell line dictionary file')
parser.add_argument('segments', help='input file from sanger')
parser.add_argument('dir_bed', help='output BED file directory')

args = parser.parse_args()

fn_dictionary = args.dictionary
fn_segments = args.segments
dir_bed = args.dir_bed

s2id = {}
f = open(fn_dictionary)
h = f.readline().strip('\r\n').split('\t')
idx_s = h.index("name_sanger")
for line in f:
    a = line.rstrip('\r\n').split('\t')
    s2id[ a[idx_s] ] = a[0]

f.close()
s2id["FADU"] ="FADU"

f = open(fn_segments)
#ID	ID_SAMPLE	ID_tumour	Primary site	Site subtype	Primary histology	Histology subtype	SAMPLE_NAME	TOTAL_CN	MINOR_ALLELE	MUT_TYPE	ID_STUDY	Chromosome:G_Start..G_Stop
h = f.readline().strip('\r\n').split('\t')

idx_seg = h.index("Chromosome:G_Start..G_Stop")
idx_min = h.index("MINOR_ALLELE")
idx_tot = h.index("TOTAL_CN")
idx_lin = h.index("SAMPLE_NAME")

current_cell_line = ""
for line in f:
    a = line.rstrip('\r\n').split('\t')
    seg, min, tot, cell_line = a[idx_seg], a[idx_min], a[idx_tot], a[idx_lin]
    cell_line = s2id[cell_line]
    if cell_line != current_cell_line:
        if current_cell_line != "":
            fo.close()
        
        current_cell_line = cell_line
        fo = open( dir_bed + '/' + cell_line + "___segments_PICNIC.BED", "w")
    
    chrom, be = seg.split(":")
    loc_start, loc_end = be.split("..")
    fo.write( chrom + '\t' + loc_start + '\t' + loc_end + '\t' +  'id:' + cell_line + ',minor:' + min + ',major:' + tot + '\n')

fo.close()
