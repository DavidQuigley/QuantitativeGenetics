# sample attributes must have a column named filename, and patient identifier is first column 
#
# parse_segments_to_BED --samples /datasets/human_esophageal_TCGA/CN/sample_attributes.txt \
#   --dir_segments /datasets/human_esophageal_TCGA/CN/segments \
#   --dir_beds  /datasets/human_esophageal_TCGA/CN/segment_beds

import argparse
parser = argparse.ArgumentParser(description='convert TCGA copy number data to BED files')
parser.add_argument('--samples', dest="fn_sa", help='path to sample attributes')
parser.add_argument('--dir_segments', dest="dir_segments", help='path to read segment files')
parser.add_argument('--dir_beds', dest="dir_beds", help='path to write bed files')

args = parser.parse_args()

dir_segments = args.dir_segments
dir_beds = args.dir_beds

f = open(args.fn_sa)
h = f.readline().rstrip('\r\n').split('\t')
idx_filename = h.index('filename')
segment_files = []
sample_names = []
for line in f:
    a = line.rstrip('\r\n').split('\t')
    segment_files.append( a[idx_filename] )
    sample_names.append( a[0] )

f.close()
for i, segment_file in enumerate(segment_files):
    f = open( dir_segments + '/' + segment_file )
    fo = open( dir_beds + '/' + sample_names[i] + '_segments.bed', 'w' )
    h = f.readline()
    for line in f:
        a = line.rstrip('\r\n').split('\t')
        fo.write( a[1] + '\t' + a[2] + '\t' + a[3] + '\t' + a[5] + '\n' )
    
    fo.close()
    f.close()
