# EGA delivers metadata as a pile of XML files. Connect the filename to the common name
# of the cell line
# EGA_link_files_to_common_names.py --dir_experiments /raw/human_lines_sanger/bam/EGAD00001001039/xmls/experiments \
#  --dir_samples /raw/human_lines_sanger/bam/EGAD00001001039/xmls/samples 
#  --fn_out /raw/human_lines_sanger/bam/EGA_file_lookup
# 
import os
import argparse
parser = argparse.ArgumentParser(description='connect EGA filename to common name')
parser.add_argument('--dir_experiments', dest="dir_experiments", help='path to xmls/experiments')
parser.add_argument('--dir_samples', dest="dir_samples", help='path to xmls/samples')
parser.add_argument('--fn_out', dest="fn_out", help='path to output file')
args = parser.parse_args()

EGA_alias2file = {}
EGA_alias2common = {}

dir_root_exp = args.dir_experiments
dir_root_sam = args.dir_samples
fn_out = args.fn_out

fn_out = "/raw/human_lines_sanger/bam/
import xml.etree.ElementTree as ET
for file in os.listdir(dir_root_exp):
    tree = ET.parse(dir_root_exp + '/' + file)
    root = tree.getroot()
    file_ID = root[0].attrib["alias"]
    sample_alias = root[0][2][1].attrib["refname"]
    EGA_alias2file[ sample_alias ] = file_ID

for file in os.listdir(dir_root_sam):
    tree = ET.parse(dir_root_sam + '/' + file)
    root = tree.getroot()
    common_name = root[0][2][4][1].text
    sample_alias=root[0].attrib["alias"]
    if not common_name is None:
        EGA_alias2common[ sample_alias ] = common_name

f = open(fn_out, 'w')
f.write("ALIAS\tfilename\tcommon_name\n")
for alias in EGA_alias2file:
    try:
        print alias + '\t' + EGA_alias2file[alias] + '\t' + EGA_alias2common[alias]
    except KeyError:
        print alias + '\t' + EGA_alias2file[alias] + '\tNA'

f.close()