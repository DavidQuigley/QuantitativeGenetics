import os

dir_o2 = "/Volumes/2014_backup/datasets/human_lines_ccle/CN/output2"
fn_ploidy = '/datasets/human_lines_database/cell_line_attributes/CCLE_ploidy.txt'
fn_cells = '/datasets/human_lines_database/cell_line_attributes/cell_line_dictionary'
fn_dirs = os.listdir(dir_o2)
fn_dirs.remove(".DS_Store")
fn2p = {}
for fn_dir in fn_dirs:
    fn = fn_dir.split('/')[-1]
    fp = open( dir_o2 + '/' + fn_dir + '/ploidy_' + fn + '.csv')
    ploidy = fp.readline().split(',')[1]
    fn2p[fn] = ploidy

fo = open(fn_ploidy, 'w')
fo.write("IDENTIFIER\tploidy\n")
for fn in sorted(fn2p.keys()):
    fo.write( fn.replace("_feature.TXT","") + '\t' + fn2p[fn] + '\n')        

fo.close()
