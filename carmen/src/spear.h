typedef boost::unordered_map<std::size_t, std::vector<int>* > HASH_SIZE_VEC;

class CacheHash{
public:
	CacheHash( Matrix<float>* data, std::vector<int> valid_cols );
	bool has( int idx, std::set<int>& vals );
	std::vector<int>* retreive( int idx, std::set<int>& vals );
	void store( int idx, std::set<int>& vals, std::vector<int>* ranks );
	void clear( int idx );
	std::set<int>* missing_by_idx( int idx );
	~CacheHash();
private:
	boost::unordered_map<int, HASH_SIZE_VEC*> hsh;
	boost::unordered_map< int, std::set<int>* > idx2missing;
};

class Spear{
public:
	Spear(double rho_a, double rho_b, int row1, int row2, double z_score);
    Spear(double rho_a, double rho_b, int row1, int row2, double z_score, double perm_pvalue);
	Spear(Spear* source);
	const double rho_a();
	const double rho_b();
	const int row_1();
	const int row_2();
	const double perm_p_value();
	const double delta();
	const double z_score();

    const int rhoa, rhob, row1, row2, zscore; // rhoa rhob and zscore are stored as ints to save space
    double perm_pvalue;
};

class RewiringResult{
public:
    RewiringResult(int row_1, double z_sum, double z_increase, double z_decrease, double p_perm, int n_z_lt5, int n_z_5_6, int n_z_6_7, int n_z_7_8, int n_z_8_9, int n_z_gt9 );
    
    int row_1();
    double z_sum();
    double z_increase();
    double z_decrease();
    double p_perm();
    int n_z_lt5();
    int n_z_5_6();
    int n_z_6_7();
    int n_z_7_8();
    int n_z_8_9();
    int n_z_gt9();
    void increment_p_perm(double value);
    double row1, rhoa, zsum, zinc, zdec, z_lt5, z_5_6, z_6_7, z_7_8, z_8_9, z_gt9, pperm;
};

class Sorter{
public:
	Sorter(double val, double index);
	double v[3];
};


class DifferentalCorrelationResult{
public:
	DifferentalCorrelationResult( std::string title, int N, double meanA, double meanB, double mean_diff_obs, double pval);
	std::string title;
	int N;
	double meanA, meanB, meanDiff, pval;
};

class ProbeSet{
public:
	ProbeSet(std::string title, std::vector<std::string>& probelist, std::vector<std::string>& probe_names);
	std::string title;
	std::vector<std::string> probelist;
	std::vector<std::string> probe_names;
};

class ProbeSets{
public:
	ProbeSets();
	void add_from_genelist( std::string set_title, Attributes* ga, std::vector<std::string> genelist, int min_probe_set_size );
	void add_from_probelist(std::string set_title, Attributes* ga, std::vector<std::string> probelist, int min_probe_set_size );
	void clear();
	ProbeSet* next_probeset();
	void reset_iterator();
	int size();
	~ProbeSets();
private:
	int itr_counter;
	std::vector<ProbeSet*> probesets;
};

class Rowpair_iterator{
public:
	Rowpair_iterator(int N, int initial_idx, int maximum_idx);
	Rowpair_iterator(std::vector< std::pair<int, int>* >& inbound_pairs);
	bool next_rowpair(int& row1, int& row2 );
	int get_ctr();
	~Rowpair_iterator();
protected:
	int ctr, i, j, N, maximum_idx;
	std::vector< std::pair<int, int>* > pairs;
};

class Spearman{
public:
	Spearman();
	Spearman(std::string fn_in, Attributes* ga);
	double correlation(Matrix<float>* A, Matrix<float>* B, std::string method);
	void calculate_rewiring_coefficient();
    void calculate_differential_correlation_by_probesets(ProbeSets& ps, std::vector<DifferentalCorrelationResult*>& results);
	double find_var(std::vector<double>* C);
	bool get_mean_difference(double& meanA, double& meanB);
	int get_n_probes_to_process();
	bool get_success();
	void load_data_for_run();
	void prune_data_by_seeds(); // called by run(), run_DC(); determines if idx should be modified
	void set_allow_uncorrelated_loci_with_eQTL(bool allow);
	void set_GWER_method(std::string method);
	void set_input_files(std::string fn_expr, std::string fn_sa, std::string fn_ga);
	void set_do_distribution(bool do_distribution);
	void set_limit_a(std::string class_a);
	void set_limit_b(std::string class_b);
	void set_method(std::string corr_type);
	void set_min_zscore(double min_zscore);
	void set_seeds(std::vector<std::string>& seeds);
	void set_min_var(double min_var);
	void set_max_eqtl_pval(double max_pval);
	void set_percent_required(double percent_required);
	void set_corr_diff(double corr_diff);
	void set_corr_abs_a(double corr_abs);
	void set_corr_abs_b(double corr_abs);
	void set_batch_extension(std::string ext);
	void set_fn_cytoscape(std::string fn);
	void set_fn_cytoscape_props(std::string fn);
	void set_fn_eqtl(std::string fn);
	void set_n_permutations(int n_perms);
	void set_gene_name_column(std::string gene_name_column );
	void set_column_extra_attributes(std::string fn_attributes);
	void set_include_seed_neighbor_correlations(bool include_seed_neighbor_correlations);
	void set_require_eqtl(bool req);
	void set_fn_out( std::string fn_out);
    void set_n_threads(int n);
	void set_verbose( bool is_verbose );
	void set_limit_network_to_seeds( bool limit_to_seeds );
	void set_min_clique_size( int min_clique );
	void run();
	void run_DC(Rowpair_iterator* RPI, std::vector<int>* pvalue_distribution, double early_stop_pvalue);
	void results(std::vector< Spear* >& spears);
	void print_spears();
	void write_spears();
    void print_rewiring_coefficient();
	void write_rewiring_coefficient();
	void print_distribution();
	void write_distribution();
	int n_spears();
	void write_to_cytoscape(double min_abs, std::string base_location, Attributes* ga);
	void write_to_cytoscape(double min_abs, std::string base_location, Attributes* ga, GOAnnotationParser* go, GeneAnnotationParser* gene_parser);
	void write_to_cytoscape(double min_abs, std::string fn_base, Attributes* ga, GOAnnotationParser* go, GeneAnnotationParser* gene_parser, std::string fn_eQTL, double max_perm_pval, bool require_eQTL);
	double find_spearman_from_ranks( Matrix<double>* ranks1, int row1, Matrix<double>* ranks2, int row2);
	double find_spearman_from_ranks( Matrix<double>* ranks, int row1, int row2, std::vector<int>& cols1, std::vector<int>& cols2 );
	void shuffle_ranks_between_labels(Matrix<double>* ranks_obs_a, Matrix<double>* ranks_obs_b, Matrix<double>* ranks_perm_a, Matrix<double>* ranks_perm_b );
	~Spearman();

private:
	Dataset* data;
	Attributes* sa;
	Attributes* ga;
	std::vector<Spear*> spears;
    std::vector<RewiringResult*> rewiring_results;
	std::vector<int> idx;
	std::vector<int> rho_distribution; // can be filled by find_DC_distribution()
	int n_before_NA, n_after_NA;
	bool do_distribution;
	std::string fn_expr, fn_sa, fn_ga, class_a, class_b, corr_type, fn_out, fn_eqtl, gene_name_column, GWER_method;
	std::string fn_cytoscape, fn_cytoscape_props, batch_extension, col_extra_attributes;
	std::vector<std::string> seeds;
	double min_var, percent_required, corr_diff, corr_abs_a, corr_abs_b, max_eqtl_pval, min_zscore;
	int n_perms, current_permutation, min_clique_size, n_threads;
	bool verbose, success, include_seed_neighbor_correlations, limit_network_to_seeds, require_eqtl, allow_uncorrelated_loci_with_eQTL;
    
    boost::mutex io_mutex;
    boost::mutex results_mutex;
    boost::mutex thread_iter_mutex;
    Matrix<int>* permutations_idx_a;
    Matrix<int>* permutations_idx_b;
    double z_score_sum_obs;
    double global_p_perm;
    
	double calculate_differential_correlations_probeset_pval(double mean_diff_obs);
	void calculate_differential_correlation_GWER();
	void calculate_ranks( Matrix<float>* raw_data, int row, std::vector<int>* idx, std::vector<int>& ranks);
	void convert_spears_to_graph(Graph* G, Attributes* ga, boost::unordered_map<int, std::vector<std::string>* >& g2p, double min_abs);
	void generate_result_header(std::string& header);
	void limit_ids_by_var();
    void limit_ids_by_NA();
    double fisher_zscore(double rhoA, double rhoB, int nA, int nB);
	void limit_ids_by_seeds(); // called within prune_data_by_seeds, does the actual modification of idx
	void load_data_for_differential_correlation();

	double mean(std::vector<double>& M);
	void find_ranks( Matrix<float>* raw_data, std::vector<int>& a_idx, Matrix<double>* ranks);
	double find_spearman( Matrix<float>* raw_data1, int row1, std::vector<int>* sample_idx1, Matrix<float>* raw_data2, int row2, std::vector<int>* sample_idx2);
	void find_DC_experimentwise();
	void find_DC_probewise_pvals();
	void find_DC_distribution();
	double find_correlation_r(Matrix<float>* raw_data1, int row1, std::vector<int>* idx_1, Matrix<float>* raw_data2, int row2, std::vector<int>* idx_2);
	void find_spears();
	void find_GWER();
	void find_DC_GWER();
	void permute_group_labels( std::vector<int>* a_idx, std::vector<int>* b_idx, std::vector<int>& a_idx_perm, std::vector<int>& b_idx_perm );
    void permute_group_labels_persistent( std::vector<int>* a_idx, std::vector<int>* b_idx, std::vector<int>& a_idx_perm, std::vector<int>& b_idx_perm, boost::mt19937& rng );
    
    int request_permutation_number();
    void process_rewiring_in_thread(int thread_id );

	void prepare_graphs(Graph* G, Graph* G_QTL, HASH_I_VECTOR_STR& g2p, double min_abs, Attributes* ga, std::string fn_eQTL, double max_perm_pval, bool require_eQTL);
	void remove_empty_seeds();
    void write_GO_terms_for_cytoscape( GeneAnnotationParser* gene_parser, GOAnnotationParser* go, std::string fn_bp, std::string fn_mf, std::string fn_cc, Graph* G );
    void write_extra_attribute_for_cytoscape( Attributes* ga, std::string fn_extra, Graph* G, HASH_I_VECTOR_STR& g2p );
    void get_gene_name(Graph* g, int idx, std::string& g1);

};

