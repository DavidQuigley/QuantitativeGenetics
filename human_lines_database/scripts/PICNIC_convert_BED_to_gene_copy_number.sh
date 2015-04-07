#!/bin/bash
for f in /datasets/human_lines_ccle/CN/segment_beds/*.BED; do
    intersectBed -wb -a /datasets/human_lines_database/annotations/refseq_genes_hg19.BED \
    -b $f | cut -f 4,8 > "$f"___gene_CN.txt
done
rm /datasets/human_lines_ccle/CN/gene_CN/*.txt
mv /datasets/human_lines_ccle/CN/segment_beds/*.txt  /datasets/human_lines_ccle/CN/gene_CN

for f in /datasets/human_lines_ccle/CN/gene_CN/*.txt; do
    cat $f | /datasets/human_lines_database/scripts/PICNIC_reformat_bed.awk > $f.expanded.txt
done

rm /datasets/human_lines_ccle/CN/gene_CN/*_CN.txt
