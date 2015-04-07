#python /datasets/human_lines_database/scripts/SANGER_convert_CN_to_BED.py \
# /datasets/human_lines_database/cell_line_attributes/cell_line_dictionary_2015_02_18.txt \
# /datasets/human_lines_sanger/CosmicCLPCompleteCNA_no_Y_chr.tsv \
# /datasets/human_lines_sanger/CN/segment_beds_extreme

DIR=/datasets/human_lines_sanger/CN

# ----------------------------------------------------------------------------------------
# Segment BED files I can make from sanger only contain extreme values (0 or >5). To 
# work sensibly with downstream analysis I need BED files that have complete coverage.
# Hopefully this code will be obsolete when I get the raw CEL files.
# ----------------------------------------------------------------------------------------

#ls $DIR/segment_beds_extreme/ | grep BED > $DIR/segment_beds/filenames 

#while read filename; do
#    cell_line=(${filename//___/ })
#    echo $cell_line
#    bedtools subtract -a /datasets/human_lines_database/annotations/hg19_chromosomes.BED \
#    -b $DIR/segment_beds_extreme/$filename | cat - $DIR/segment_beds_extreme/$filename | sort -n |  sed -e "s/hg19/$cell_line/g" > $DIR/segment_beds/$filename
#done < $DIR/segment_beds/filenames
 
    
ECHO "generated BED files"
# ----------------------------------------------------------------------------------------
# Create gene-level calls by Intersecting segment BED files with gene BED files
# Move files into CN/gene_CN
# reformat to expand out minor allele and total copy number
# remove BED file intermediates
# ----------------------------------------------------------------------------------------



ECHO "intersecting BED files with gene loci"
for f in $DIR/segment_beds/*.BED; do
    intersectBed -wb -a /datasets/human_lines_database/annotations/refseq_genes_hg19_2015_02_20.BED \
    -b $f | cut -f 4,8 > "$f"___gene_CN.txt
done
rm $DIR/gene_CN/*.txt
ECHO "reformatting gene CN files"
mv $DIR/segment_beds/*.txt  $DIR/gene_CN
for f in $DIR/gene_CN/*.txt; do
    cat $f | /datasets/human_lines_database/scripts/PICNIC_reformat_bed.awk > $f.expanded.txt
done
rm $DIR/gene_CN/*_CN.txt
    
# ----------------------------------------------------------------------------------------
# combine individual gene copy number files into a single matrix
# First do total copy, then minor copy
# ----------------------------------------------------------------------------------------
ECHO "generating consolidated CN file: total"
Rscript /datasets/human_lines_database/scripts/PICNIC_aggregate_CCLE_CN_files.R \
 /datasets/human_lines_database/annotations/refseq_genes_hg19_2015_02_20.BED \
 $DIR/gene_CN \
 /datasets/human_lines_database/copy_number/CN_total_sanger.txt \
 total

ECHO "generating consolidated CN file: minor"
 Rscript /datasets/human_lines_database/scripts/PICNIC_aggregate_CCLE_CN_files.R \
 /datasets/human_lines_database/annotations/refseq_genes_hg19_2015_02_20.BED \
 $DIR/gene_CN \
 /datasets/human_lines_database/copy_number/CN_minor_sanger.txt \
 minor
