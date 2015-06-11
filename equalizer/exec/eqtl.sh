bin/eqtl \
-d./results/expression/mammary/MoGene_expr_above_bg_refseq.txt \
-f./data/sample_attributes_mammary.txt \
-g./results/expression/mammary/MoGene_gene_attributes_above_bg_refseq.txt \
-e./data/genotypes/calls.txt \
-h./data/genotypes/sample_attributes.txt \
-s./data/genotypes/gene_attributes.txt \
-ysymbol \
-cmammary.id \
-n1000 \
-r10513583 \
-o./results/eqtl/eqtl_mammary_above_bg_refseq_1000.txt

bin/eqtl \
-d./results/expression/mammary/MoGene_nosprfvbsnps_expr_above_bg_refseq.txt \
-f./data/sample_attributes_mammary.txt \
-g./results/expression/mammary/MoGene_nosprfvbsnps_gene_attributes_above_bg_refseq.txt \
-e./data/genotypes/calls.txt \
-h./data/genotypes/sample_attributes.txt \
-s./data/genotypes/gene_attributes.txt \
-ysymbol \
-cmammary.id \
-n1000 \
-r10513583 \
-o./results/eqtl/eqtl_mammary_CORRECTED_above_bg_refseq_1000.txt

echo "uncorrected result:"
grep Cdc26 ./results/eqtl/eqtl_mammary_above_bg_refseq_1000.txt
echo "corrected result:"
grep Cdc26 ./results/eqtl/eqtl_mammary_CORRECTED_above_bg_refseq_1000.txt
