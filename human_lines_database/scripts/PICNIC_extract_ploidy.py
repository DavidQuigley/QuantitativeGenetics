import os
import argparse

parser = argparse.ArgumentParser(description='extract ploidy')
parser.add_argument('dir_ploidy', help='input segment file directory')
parser.add_argument('fn_output', help='file to write')
args = parser.parse_args()

dir_o2 = args.dir_ploidy
fn_dirs = os.listdir(dir_o2)
fn_ploidy = args.fn_output 

if ".DS_Store" in fn_dirs:
    fn_dirs.remove(".DS_Store")
    
fn2p = {}
for fn_dir in fn_dirs:
    fn = fn_dir.split('/')[-1]
    fp = open( dir_o2 + '/' + fn )
    ploidy = fp.readline().split(',')[1]
    fn2p[fn] = ploidy

fo = open(fn_ploidy, 'w')
fo.write("IDENTIFIER\tploidy\n")
for fn in sorted(fn2p.keys()):
    cell_line_name = fn.replace("_feature.TXT","")
    cell_line_name = cell_line_name.replace("ploidy_", "")
    cell_line_name = cell_line_name.replace(".csv", "")
    fo.write( cell_line_name + '\t' + fn2p[fn] + '\n')        

fo.close()
