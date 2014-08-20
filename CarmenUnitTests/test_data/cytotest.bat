# Source CARMEN 1.0
# Date 2014_08_20_13_57_03
# Class_A  8 members.
# Class_B  0 members.
# Focus_Probe GENE1
# Min_var 0
# Delta_rho 0
# Abs_rho 0
# Abs_rho_Class_B 0
# Data_File ../../test_data/test_correlation_dataset.txt
# Sample_File ../../test_data/test_sample_attributes.txt
# Genes_File ../../test_data/test_correlation_gene_attributes.txt
# Max_pval_eqtl 1
# Rho_difference 
# Min_clique_size 0
# Min_percent_present 0
# Gene_name_column 
# Include_seed_neighbor_correlations True
# Limit_network_to_seeds False
# Min_clique_size 0
# N_perms 0
# Percent_required 0
# Limit_to_seeds False
# Seeds GENE1
"/cytoscape/cytoscape.exe" -N "../../test_data/cytotest.sif" -e "../../test_data/cytotest.eda" -e "../../test_data/cytotest_dir.eda" -V "/cytoscape/props.txt" -n "../../test_data/cytotest.noa" -e "../../test_data/cytotest_pval.eda"
