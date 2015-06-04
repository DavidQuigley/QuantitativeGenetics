#!/bin/bash

# ----------------------------------------------------------------------------------------
# RUN ON CLUSTER: generate jobs and run them
# ----------------------------------------------------------------------------------------

JAVA_BIN=/mnt/speed/quigleyd/software/PICNIC/Matlab_Compiler_Runtime/_jvm/bin/java
PICNIC_ROOT=/mnt/speed/quigleyd/software/PICNIC
OUT_ROOT=/mnt/speed/quigleyd/CN
DIR_HOME=/home/quigleyd
LOCAL_CN="/raw/human_lines_sanger/CN"
DIR_SCRIPTS=/notebook/code/src/human_lines_database/scripts
DIR_SANGER="/datasets/human_lines_sanger/CN"
DIR_BEDTOOLS=/opt/bedtools2/bin
# convert CEL files to text format

#$JAVA_BIN -Xmx2G \
#  -jar $PICNIC_ROOT/celConverter/CelFileConverter.jar \
#  -c $PICNIC_ROOT/cdf/GenomeWideSNP_6.Full.cdf \
#  -m $PICNIC_ROOT/celConverter/Snp6FeatureMappings.csv \
#  -s $OUT_ROOT/raw_data/sanger \
#  -t $OUT_ROOT/PICNIC_outdir/raw

# rename feature intensity files; required because PICNIC expects this file format
# and mangles directory structure if not present

#cd $OUT_ROOT/PICNIC_outdir/raw
#for f in *.feature_intensity; do
#    mv $f CGP_$f 
#done

# Generate commands to run PICNIC in job control system
python $PICNIC_ROOT/PICNIC_generate_jobs.py

# Submit jobs
sh $DIR_HOME/scripts/run_PICNIC_jobs.sh

ls $OUT_ROOT/PICNIC_outdir/output2 | sort > $DIR_HOME/output_dirlist

# ----------------------------------------------------------------------------------------
# RUN LOCALLY
# ----------------------------------------------------------------------------------------

# Pull down directory list
scp -P 2223 quigleyd@64.54.200.169:$DIR_HOME/output_dirlist $LOCAL_CN/output_dirlist

while read f; do
   scp -P 2223 quigleyd@64.54.200.169:$OUT_ROOT/PICNIC_outdir/output2/$f/segments_$f.csv $LOCAL_CN/segments/$f
done <$LOCAL_CN/output_dirlist

while read f; do
   scp -P 2223 quigleyd@64.54.200.169:$OUT_ROOT/PICNIC_outdir/output2/$f/ploidy_$f.csv $LOCAL_CN/ploidy/ploidy_$f.csv
done <$LOCAL_CN/output_dirlist

# sh PICNIC_download_segment_files.sh

# ----------------------------------------------------------------------------------------
# convert segments to BED format and lift over to hg19
# ----------------------------------------------------------------------------------------
python $DIR_SCRIPTS/PICNIC_convert_segment_to_BED.py \
       $LOCAL_CN/segments $DIR_SANGER/segment_beds name_sanger

# ----------------------------------------------------------------------------------------
# extract and store ploidy information
# ----------------------------------------------------------------------------------------
python $DIR_SCRIPTS/PICNIC_extract_ploidy.py $LOCAL_CN/ploidy $DIR_SANGER/sanger_ploidy.txt


# ----------------------------------------------------------------------------------------
# Create gene-level calls by Intersecting segment BED files with gene BED files
# Move files into CN/gene_CN
# reformat to expand out minor allele and total copy number
# remove BED file intermediates
# ----------------------------------------------------------------------------------------

for f in $DIR_SANGER/segment_beds/*.BED; do
    $DIR_BEDTOOLS/intersectBed -wb -a /datasets/human_lines_database/annotations/refseq_genes_hg19_2015_02_20.BED \
    -b $f | cut -f 4,8 > "$f"___gene_CN.txt
done
rm $DIR_SANGER/gene_CN/*.txt
mv $DIR_SANGER/segment_beds/*.txt  $DIR_SANGER/gene_CN
for f in $DIR_SANGER/gene_CN/*.txt; do
    cat $f | $DIR_SCRIPTS/PICNIC_reformat_bed.awk > $f.expanded.txt
done
rm $DIR_SANGER/gene_CN/*_CN.txt

# ----------------------------------------------------------------------------------------
# combine individual gene copy number files into a single matrix
# First do total copy, then minor copy
# ----------------------------------------------------------------------------------------
Rscript $DIR_SCRIPTS/PICNIC_aggregate_CCLE_CN_files.R \
/datasets/human_lines_database/annotations/refseq_genes_hg19_2015_02_20.BED \
$DIR_SANGER/gene_CN \
/datasets/human_lines_database/copy_number/CN_total_sanger.txt \
total
gzip /datasets/human_lines_database/copy_number/CN_total_sanger.txt

Rscript $DIR_SCRIPTS/PICNIC_aggregate_CCLE_CN_files.R \
/datasets/human_lines_database/annotations/refseq_genes_hg19_2015_02_20.BED \
$DIR_SANGER/gene_CN \
/datasets/human_lines_database/copy_number/CN_minor_sanger.txt \
minor
gzip /datasets/human_lines_database/copy_number/CN_minor_sanger.txt

# ----------------------------------------------------------------------------------------
# Update proberef to HG19 using liftover values, result is ProbeRef_hg19.csv
# ----------------------------------------------------------------------------------------
#python $DIR_SCRIPTS/PICNIC_lift_ProbeRef_to_HG19.py \
#  /datasets/human_lines_database/annotations/ProbeRef.csv \
#  /datasets/human_lines_database/annotations/ProbeRef_hg19.csv \
#  /datasets/human_lines_database/annotations/hg19_SNP6_liftover_results.bed \
#  /datasets/human_lines_database/annotations/hg19_SNP6_liftover_best_fits_for_deleted_loci.bed


# ----------------------------------------------------------------------------------------
# Combine all segment files into a single source
# ----------------------------------------------------------------------------------------
Rscript $DIR_SCRIPTS/PICNIC_consolidate_segments.R \
    $DIR_SANGER/segment_beds /datasets/human_lines_database/copy_number/sanger_segments

