class ClassMinerOptions{
public:
	ClassMinerOptions();
	std::string file_name_dis;	
	std::string file_name_sa;
	std::string file_name_ga;
	std::string file_name_out;
	std::string class_a;
	std::string class_b;	
	std::string discretization;
	std::string mine_type;
	std::string which_class;
	float disc_lower;
	float disc_upper;
	int min_sup;
	float min_conf;
	float min_imp;
	double max_chi2;
	std::string gene_limit;
	int max_depth;
	int timeout;
	int out_style;
	int n_tests;
	time_t start_time;
	bool verbose;
    std::string version;
};
