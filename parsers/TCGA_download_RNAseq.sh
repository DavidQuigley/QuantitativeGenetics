BYY=$1
BMM=$2
BDD=$3
BTISSUE=$4
DATA_DIR=$5

wget http://gdac.broadinstitute.org/runs/stddata__${BYY}_${BMM}_${BDD}/data/${BTISSUE}/${BYY}${BMM}${BDD}/gdac.broadinstitute.org_${BTISSUE}.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.${BYY}${BMM}${BDD}00.0.0.tar.gz
gunzip gdac.broadinstitute.org_${BTISSUE}.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.${BYY}${BMM}${BDD}00.0.0.tar.gz
tar -xf gdac.broadinstitute.org_${BTISSUE}.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.${BYY}${BMM}${BDD}00.0.0.tar
mv gdac.broadinstitute.org_${BTISSUE}.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.${BYY}${BMM}${BDD}00.0.0/${BTISSUE}.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt rnaseq_${BTISSUE}_rsem_genes_normalized.txt
rm -r gdac.broadinstitute.org_${BTISSUE}.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.${BYY}${BMM}${BDD}00.0.0
rm gdac.broadinstitute.org_${BTISSUE}.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.${BYY}${BMM}${BDD}00.0.0.tar

wget http://gdac.broadinstitute.org/runs/stddata__${BYY}_${BMM}_${BDD}/data/${BTISSUE}/${BYY}${BMM}${BDD}/gdac.broadinstitute.org_${BTISSUE}.Merge_Clinical.Level_1.${BYY}${BMM}${BDD}00.0.0.tar.gz
gunzip gdac.broadinstitute.org_${BTISSUE}.Merge_Clinical.Level_1.${BYY}${BMM}${BDD}00.0.0.tar.gz
tar -xf gdac.broadinstitute.org_${BTISSUE}.Merge_Clinical.Level_1.${BYY}${BMM}${BDD}00.0.0.tar
rm gdac.broadinstitute.org_${BTISSUE}.Merge_Clinical.Level_1.${BYY}${BMM}${BDD}00.0.0.tar
mv gdac.broadinstitute.org_${BTISSUE}.Merge_Clinical.Level_1.${BYY}${BMM}${BDD}00.0.0/${BTISSUE}.clin.merged.txt sample_attributes_merged_${BTISSUE}.txt
rm -r gdac.broadinstitute.org_${BTISSUE}.Merge_Clinical.Level_1.${BYY}${BMM}${BDD}00.0.0

wget http://gdac.broadinstitute.org/runs/stddata__${BYY}_${BMM}_${BDD}/data/${BTISSUE}/${BYY}${BMM}${BDD}/gdac.broadinstitute.org_${BTISSUE}.Clinical_Pick_Tier1.Level_4.${BYY}${BMM}${BDD}00.1.0.tar.gz
gunzip gdac.broadinstitute.org_${BTISSUE}.Clinical_Pick_Tier1.Level_4.${BYY}${BMM}${BDD}00.1.0.tar.gz
tar -xf gdac.broadinstitute.org_${BTISSUE}.Clinical_Pick_Tier1.Level_4.${BYY}${BMM}${BDD}00.1.0.tar
rm gdac.broadinstitute.org_${BTISSUE}.Clinical_Pick_Tier1.Level_4.${BYY}${BMM}${BDD}00.1.0.tar
mv gdac.broadinstitute.org_${BTISSUE}.Clinical_Pick_Tier1.Level_4.${BYY}${BMM}${BDD}00.1.0/${BTISSUE}.clin.merged.picked.txt sample_attributes_${BTISSUE}.txt
rm -r gdac.broadinstitute.org_${BTISSUE}.Clinical_Pick_Tier1.Level_4.${BYY}${BMM}${BDD}00.1.0  

python /notebook/code/src/parsers/parse_TCGA_rnaseq_normalized.py \
--file_in $DATA_DIR/rnaseq_${BTISSUE}_rsem_genes_normalized.txt \
--file_sa_in $DATA_DIR/sample_attributes_${BTISSUE}.txt \
--file_genes_out $DATA_DIR/gene_attributes_rnaseq_${BTISSUE}_rsem_genes_normalized.txt \
--file_expr_out $DATA_DIR/expr_rnaseq2_${BTISSUE}_rsem_genes_normalized.txt \
--file_sa_out $DATA_DIR/sample_attributes_${BTISSUE}_rsem_genes_normalized.txt \
--require_symbol T
