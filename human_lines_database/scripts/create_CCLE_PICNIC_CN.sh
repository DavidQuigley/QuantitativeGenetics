#!/bin/bash


# ************************************************************************************
# ----------------------------------------------------------------------------------------
# RUN ON CLUSTER: generate jobs and run them
# ----------------------------------------------------------------------------------------

JAVA_BIN=/mnt/speed/quigleyd/software/PICNIC/Matlab_Compiler_Runtime/_jvm/bin/java
PICNIC_ROOT=/mnt/speed/quigleyd/software/PICNIC
OUT_ROOT=/mnt/speed/quigleyd/CN
DIR_HOME=/home/quigleyd
LOCAL_CN="/raw/human_lines_CCLE/CN"
DIR_SCRIPTS=/notebook/code/src/human_lines_database/scripts
DIR_CCLE="/datasets/human_lines_CCLE/CN"
DIR_BEDTOOLS=/opt/bedtools2/bin

# convert CEL files to text format

#$JAVA_BIN -Xmx2G \
#  -jar $PICNIC_ROOT/celConverter/CelFileConverter.jar \
#  -c $PICNIC_ROOT/cdf/GenomeWideSNP_6.Full.cdf \
#  -m $PICNIC_ROOT/celConverter/Snp6FeatureMappings.csv \
#  -s $OUT_ROOT/raw_data/CCLE \
#  -t $OUT_ROOT/PICNIC_outdir/raw

# rename feature intensity files; required because PICNIC expects this file format
# and mangles directory structure if not present

#cd $OUT_ROOT/PICNIC_outdir/raw
#for f in *.feature_intensity; do
#    mv $f CGP_$f 
#done

# python PICNIC_generate_jobs.py
# sh run_PICNIC_jobs.sh

# ************************************************************************************
# RUN LOCALLY
# ----------------------------------------------------------------------------------------
# download segment files by SCP
# ----------------------------------------------------------------------------------------
# sh PICNIC_download_segment_files.sh

# ----------------------------------------------------------------------------------------
# convert segments to BED format and lift over to hg19
# ----------------------------------------------------------------------------------------
python /datasets/human_lines_database/scripts/PICNIC_convert_segment_to_BED.py

# ----------------------------------------------------------------------------------------
# extract and store ploidy information
# ----------------------------------------------------------------------------------------
python $DIR_SCRIPTS/PICNIC_extract_ploidy.py $LOCAL_CN/ploidy $DIR_CCLE/ccle_ploidy.txt

# ----------------------------------------------------------------------------------------
# Create gene-level calls by Intersecting segment BED files with gene BED files
# Move files into CN/gene_CN
# reformat to expand out minor allele and total copy number
# remove BED file intermediates
# ----------------------------------------------------------------------------------------

DIR=/datasets/human_lines_ccle/CN

for f in $DIR/segment_beds/*.BED; do
    intersectBed -wb -a /datasets/human_lines_database/annotations/refseq_genes_hg19.BED \
    -b $f | cut -f 4,8 > "$f"___gene_CN.txt
done
rm $DIR/gene_CN/*.txt
mv $DIR/segment_beds/*.txt  $DIR/gene_CN
for f in $DIR/gene_CN/*.txt; do
    cat $f | /datasets/human_lines_database/scripts/PICNIC_reformat_bed.awk > $f.expanded.txt
done
rm $DIR/gene_CN/*_CN.txt

# ----------------------------------------------------------------------------------------
# combine individual gene copy number files into a single matrix
# First do total copy, then minor copy
# ----------------------------------------------------------------------------------------
Rscript /datasets/human_lines_database/scripts/PICNIC_aggregate_CCLE_CN_files.R \
/datasets/human_lines_database/annotations/refseq_genes_hg19_2015_02_20.BED \
$DIR/gene_CN \
/datasets/human_lines_database/copy_number/CN_total_CCLE.txt \
total

Rscript /datasets/human_lines_database/scripts/PICNIC_aggregate_CCLE_CN_files.R \
/datasets/human_lines_database/annotations/refseq_genes_hg19_2015_02_20.BED \
$DIR/gene_CN \
/datasets/human_lines_database/copy_number/CN_minor_CCLE.txt \
minor

# ----------------------------------------------------------------------------------------
# Update proberef to HG19 using liftover values, result is ProbeRef_hg19.csv
# ----------------------------------------------------------------------------------------
python /datasets/human_lines_database/scripts/PICNIC_lift_ProbeRef_to_HG19.py \
  $DIR/PICNIC/annotation/ProbeRef.csv \
  $DIR/PICNIC/annotation/ProbeRef_hg19.csv \
  $DIR/PICNIC/annotation/hg19_SNP6_liftover_results.bed \
  $DIR/PICNIC/annotation/hg19_SNP6_liftover_best_fits_for_deleted_loci.bed


# ----------------------------------------------------------------------------------------
# Combine all segment files into a single source
# ----------------------------------------------------------------------------------------
Rscript /datasets/human_lines_database/scripts/PICNIC_consolidate_segments.R