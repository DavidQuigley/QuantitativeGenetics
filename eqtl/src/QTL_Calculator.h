const int LBL_A = -1;
const int LBL_B = 1;
const int METHOD_REGRESSION=1;
const int METHOD_STUDENT_T_TEST=2;

// This is the primary data structure used to hold QTL results.
class rqtl{
public:
    rqtl();
    rqtl( int idx_genotype, int idx_probe, double pvalue, double min_permuted_pval, double t_stat );
	int idx_genotype;
	int idx_probe;
	double pvalue;
    double min_permuted_pval;
    double t_stat;
    double mean_a, mean_b, mean_c;
};


class qtl_by_chrom{
public:
    qtl_by_chrom(int idx_probe, int n_chrom);
    std::vector<double> T_obs;
    std::vector<int> snp_indexes;
    std::vector<bool> is_robust_adjustment;
    std::vector<double> pvalues_raw;
    std::vector<double> pvalues_perm;
    int idx_probe;
    
    double get_max_T_obs();
    void update_perm_pvalue(std::vector<double>& max_perm_tstats);
    void calculate_final_pvalues(int n_samples, int n_perms);
};


class QTL_Calculator{
public:
    QTL_Calculator(Dataset* data_expr, Attributes* sa_expr,	Attributes* ga_expr, Dataset* data_snps, Attributes* sa_snps, Attributes* ga_snps, std::string shared_g_p_col);
        
    void limit_ids_by_var(double min_var);
    void limit_ids_by_probe_id(std::string probe_id);
    void limit_ids_by_NA(double fraction_required);
    void limit_ids_by_NA_min_max(double fraction_required, double fraction_max);
    void limit_ids_by_subset(int this_fraction, int n_fractions);
    void restrict_probesets_and_snps( std::vector<std::string>& probesets, std::vector<std::string>& snps, std::vector<int> & idx_snps);
    
    void check_probes_in_KV( KeyValuesParser& KV);
    void get_matched_snp_sample_idx( std::vector<int>& idx );
    void get_matched_gene_sample_idx( std::vector<int>& idx );
    int get_min_obs_per_genotype();
    int get_n_matched_samples();
    int get_cis_interval();
    int get_size_gene_idx();
    void set_logfile( std::string fn_log );
    void set_analysis_method( int method );
    void set_min_obs_per_genotype(int min);
    void set_n_perms(int n_perms);
    void set_cis_interval(int cis_interval);
    void set_n_perms_max(int n_perms_max);    
    void set_robust(bool robust);
    void set_thread_count(int n_threads);
    void set_verbose(bool verbose);
    
    void calculate_specific_SNPs(Matrix<double>& results, std::vector<int>& idx_snps);
    void calculate_genome_wide( bool genomewide_by_chromosome);
    void calculate_permutations_for_gene_probe_pairs( KeyValuesParser& KV ); 
    void calculate_all_regressions(ClassMinerOptions* cmo_expr, ClassMinerOptions* cmo_snps, double fraction_required, double max_pvalue_to_report);

    void print_rqtl( ClassMinerOptions * cmo_expr, ClassMinerOptions * cmo_snps, double min_var, double fraction_required );
    void print_by_chromosome( ClassMinerOptions* cmo_expr, ClassMinerOptions* cmo_snps, double fraction_required);
    void print_all_regressions(ClassMinerOptions* cmo_expr, ClassMinerOptions* cmo_snps, double fraction_required );
    
    void reportNetResults(Matrix<double>& results, std::vector<int>& idx_snps, ClassMinerOptions* cmo_expr, ClassMinerOptions* cmo_snps );
    void reportFileList( std::string fn, std::vector<std::string>& keys, std::vector<std::string>& fn_out );
    
    void sort_QTLs();
    
private:
    
    // n_perms_max is used when we have a larger number of permutations to run if the gene is very significant after n_perms iterations
    //   The normal use case would be, "run 10 permutations, if P<=0.01, run 100"
    // expr_snps is the matched indexes between gene samples and SNP samples
    //
    
    boost::mutex io_mutex;
    boost::mutex results_mutex;
    boost::mutex thread_iter_mutex;
    
    std::vector<rqtl*> qtls; 
    std::vector<qtl_by_chrom*> qtls_by_chrom;
    std::vector<int> idx; // indexes to process;
    Dataset* data_expr;
	Attributes* sa_expr;
	Attributes* ga_expr;
	Dataset* data_snps;
	Attributes* sa_snps;
	Attributes* ga_snps;
    bool is_robust, verbose;
    int n_perms, n_perms_max;
    int n_threads;
    int cis_interval;
    int analysis_method;
    int index_for_next_thread;
    std::string fn_log;  
    std::ofstream fs_log;
    std::vector<std::pair<int, int>*> expr_snps;
    int min_obs_per_genotype;
    void match_expr_and_snps(std::string shared_geno_pheno_col); 
    
    void calculate_snp_idx2chr(boost::unordered_map<int, int>& snp_idx2chr, std::string& chr_colname, int& n_chr);
    void calculate_sample_indexes(std::vector<int>* true_sample_idx, std::vector< std::vector<int>* >& perm_sample_idx, int n_samples);
    
    double find_var(std::vector<double>* C);
    double find_variance(double mean, const vector<double>& acc);
    void find_mean(std::vector<double>* V, std::vector<int>* ab, double& mean_a, double& mean_b);
    void find_mean_three_way(std::vector<double>* V, std::vector<int>* ab, double& mean_a, double& mean_b, double& mean_c);
    void find_mean_three_way(std::vector<double>* V, std::vector<double>* ab, double& mean_a, double& mean_b, double& mean_c);
    void find_mean_sd_three_way(std::vector<double>* V, std::vector<double>* ab, double& mean_a, double& mean_b, double& mean_c, double& sd_a, double& sd_b, double& sd_c);

    void fill_XY( std::vector<double>& X, std::vector<double>& Y, int snps_g_idx, int expr_g_idx, std::vector<int>* sample_idx );
    void modify_XY_to_be_robust( std::vector<double>& X, std::vector<double>& Y);
    bool regress(vector<double>& X, vector<double>& Y, double& m, double& b, double& t_stat);
    bool regress_direct(std::vector<int>* sample_idx, int expr_g_idx, int snps_g_idx, double& m, double& b, double& t_stat);
    bool student_t_direct(std::vector<int>* sample_idx, int expr_g_idx, int snps_g_idx, double& t_stat);            
    
    bool calculation_wrapper( int snps_g_idx, int expr_g_idx, std::vector<int>*  sample_idx, double& t_stat, double value_to_exceed);

    int request_next_index_to_process();
    void report_thread_progress( int thread_id, int current, int total);
    void process_genomewide_in_thread(int thread_id, std::vector<int>* true_sample_idx, std::vector< std::vector<int>* >& perm_sample_idx);
    void process_by_chr_in_thread(int thread_id,  std::vector<int>* true_sample_idx, std::vector< std::vector<int>* >& perm_sample_idx);
    void process_cis_only_in_thread(int thread_id,  std::vector<int>* true_sample_idx, std::vector< std::vector<int>* >& perm_sample_idx, std::vector<int>* idx_all_snps, int cis_interval);
    void process_all_regressions_in_thread(int thread_id, double max_pvalue_to_report );
    
    void find_best_tstat_across_loci( int expr_g_idx, std::vector<int>* sample_idx, std::vector<int>* genotype_idx, int& idx_best_locus, double& t_obs );
    void find_best_tstat_by_chr( int expr_g_idx, std::vector<int>* sample_idx, std::vector<int>& idx_Tstat_maxes, std::vector<double>& Tstat_stopvals, std::vector<double>& Tstat_maxes, boost::unordered_map<int,int>& snp_idx2chr );    
    bool permutation_tstat_across_loci_exceeds( int expr_g_idx, std::vector<int>* expr_g_indexes, std::vector<int>* sample_idx, double tstat_observed, double& max_tstat_perm );
    double permute_by_chromosome( int expr_g_idx, double T_obs, int chr_idx, std::vector< std::vector<int>* >& perm_sample_idx, boost::unordered_map<int,int>& snp_idx2chr );
    
    std::string create_output_header(ClassMinerOptions* cmo_expr, ClassMinerOptions* cmo_snps, double max_var, double fraction_required);
    void spear_result_string(rqtl* q, std::string& spear_str);
    void spear_result_string(qtl_by_chrom* q, std::string& spear_str);
    void report_progress( int current, int total); // not thread-safe
    void write_to_log(std::string msg);
    void close_logfile_if_open();
};
