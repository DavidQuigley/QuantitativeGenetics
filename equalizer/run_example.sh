#!/bin/bash
#BEDTOOLS_DIR="/opt/bedtools2/bin"

if [ -z ${BEDTOOLS_DIR+x} ]; then 
 {
   echo "ERROR: please edit this file and set the BEDTOOLS_DIR variable to the location of"
   echo "the bin subdirectory under which BEDTOOLS was installed, e.g. /opt/bedtools2/bin"
   exit 1
 }; 
 else echo "BEDTOOLS_DIR is set to '$BEDTOOLS_DIR'"; 
fi

python bin/equalizer.py \
-f gene \
-p mogene-1_1_nosprfvbSnps \
-c MoGene-1_1-st-v1 \
-v cdc26.vcf \
-a data/annotation/probe_bed/MoGene_probes_MM10.bed \
-b $BEDTOOLS_DIR \
-s na34 -g mm10 -d 38 -i 1.1 -w Mouse -r 'Mus musculus' \
-y data/annotation/MoGene-1_1/MoGene-1_1-st-v1.na34.mm10.probeset.csv \
-t data/annotation/MoGene-1_1/MoGene-1_1-st-v1.na34.mm10.transcript.csv \
-q data/annotation/MoGene-1_1/MoGene-1_1-st-v1.r4.pgf \
-m data/annotation/MoGene-1_1/MoGene-1_1-st-v1.r4.mps \
-k data/annotation/MoGene-1_1/MoGene-1_1-st-v1.r4.clf \
-u 'David Quigley' -l dquigley@cc.ucsf.edu \
-o results/MoGene_FVB_SPR

R --no-save  < results/MoGene_FVB_SPR/create_package.R
R --no-save  < normalize_example_microarrays.R

bin/eqtl \
-dresults/expression/mammary/MoGene_expr_above_bg_refseq.txt \
-fdata/sample_attributes_mammary.txt \
-gresults/expression/mammary/MoGene_gene_attributes_above_bg_refseq.txt \
-edata/genotypes/calls.txt \
-hdata/genotypes/sample_attributes.txt \
-sdata/genotypes/gene_attributes.txt \
-ysymbol -cmammary.id -n1000 \
-r10513583 \
-oresults/eqtl/eqtl_mammary_above_bg_refseq_1000.txt

bin/eqtl \
-dresults/expression/mammary/MoGene_nosprfvbsnps_expr_above_bg_refseq.txt \
-fdata/sample_attributes_mammary.txt \
-gresults/expression/mammary/MoGene_nosprfvbsnps_gene_attributes_above_bg_refseq.txt \
-edata/genotypes/calls.txt \
-hdata/genotypes/sample_attributes.txt \
-sdata/genotypes/gene_attributes.txt \
-ysymbol -cmammary.id -n1000 \
-r10513583 \
-oresults/eqtl/eqtl_mammary_CORRECTED_above_bg_refseq_1000.txt

more results/eqtl/eqtl_mammary_above_bg_refseq_1000.txt
more results/eqtl/eqtl_mammary_CORRECTED_above_bg_refseq_1000.txt
