class Rule{
public:
	Rule();
	Rule(Rule* rule_to_clone);
	double chi2;
	double sup_a_p, sup_b_p;
	int sup_a_n, sup_b_n;
	double conf_p;
	double imp;
	double weight;
	double t_stat;
	int parent_idx;
	int rule_idx;
	int rule_class;
	std::vector<std::string> rules_g;
	std::vector<std::string> rules_i;
	std::vector<std::string> rules_r;

	static const int CLASS_A = -1;
	static const int CLASS_B = 1;
	
	void        clone(Rule* rule_to_clone);
	void        clone_out_probes( std::vector<Rule*>& rules_from_split);
	void        load_line(std::string line);
	void        load_line(Ant* ant, ClassifierDataset* data, Attributes* ga, ClassMinerOptions* cmo);
	std::string get_friendly_label();
	std::string get_identifier_label();
	std::string get_identifier_cutoff_label();
	void        evaluate( ClassifierDataset* data, bool disc_is_none );
};

class Edge{
public:
	Edge(std::string name1, std::string value1, std::string name2, std::string value2, int edge_class, double sup, double conf, double p_value);
	std::string name1;
	std::string value1;
	std::string name2;
	std::string value2;
	int edge_class;
	double sup;
	double conf;
	double p_value;
};

class RuleSet{
public:
	RuleSet();
	RuleSet(ClassMinerOptions* cmo, std::vector<Ant*>* ants, ClassifierDataset* data, Attributes* ga);
	RuleSet(std::string file_name, bool strip_edges = false);
	~RuleSet();
	std::vector<Rule*> rules;
	double minsup;
	double minconf;
	double minimp;
	std::string class_a;
	int class_a_size;
	std::string class_b;
	int class_b_size;
	std::string disc;
	double disc_lower;
	double disc_upper;
	double max_chi2;
	int max_depth;
	int n_tested;
	double bon_pval;
	std::string gene_limit;
	std::string file_name_dis;
	std::string file_name_ga;
	std::string file_name_sa;
	std::vector<Edge*>* edges;

	void build_rules_by_samples(Matrix<double>* r_by_s, ClassifierDataset* data);
	void find_edges();
	
	void write_crosstab(std::string file_name);
	void write_header(std::ofstream& f_out, int n_rules, double merge_threshold);
	void write_rule(std::ofstream& f_out, int rule_idx);
	void write_rules_by_samples(std::string filename, Matrix<double>* r_by_s, Attributes* sa);
	void write_network(std::string file_name);
	void write_expression(std::string file_name, bool verbose);
	void write_frequency_list(std::string fn_out, std::string fn_source, bool verbose);
	void write_gene_name_list(std::string file_name, bool verbose);
	void unit_tests( std::string basedir );
private:
	bool assertEqual(int a, int b);
	void report(std::string testname, int tst, bool passed );
};


class CrossValidation{
public:
	std::vector< RuleSet* >* rulesets;

	CrossValidation( std::string file_name_cv );
	~CrossValidation();
	void find_core_rules(RuleSet* core, Attributes* ga, double max_p);
};
