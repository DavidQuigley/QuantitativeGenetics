class Classifier{
public:
	std::vector<int> labels;
	std::vector<double> predictions, calls;
	std::vector<int> a_idx_train, b_idx_train, a_idx_test, b_idx_test;
	std::vector<Rule*> features_available;
    std::vector<Rule*> features_learned;
	bool verbose;
    std::string method;
    int n_rounds;
    std::string lbl_a_tr, lbl_b_tr, lbl_a_tst, lbl_b_tst;
    Evaluation eval;
    Matrix<double>* r_by_s;
    typedef std::pair<std::string, std::vector<int> > limit_pair;
    
    Classifier();

    void apply_to_dataset( std::string file_name, ClassifierDataset* ds );
    void classify_colSum(int top_t, bool optimize);
	void classify_adaboost(int rounds, bool robust);   
    void load_from_ruleset(RuleSet* rs, Matrix<double>* r_by_s, limit_pair& a_tr, limit_pair& b_tr, limit_pair& a_tst, limit_pair& b_tst);
	void print(ClassifierDataset* data, bool is_LOOCV, std::string fn_rules, std::string fn_classifier);
    void write(std::string fn_out, ClassifierDataset* data, bool is_LOOCV, std::string fn_rules, std::string fn_classifier);
    void print_labels(ClassifierDataset* data, Attributes* sa, std::string fn_classifier);
	void write_labels(std::string fn_out, ClassifierDataset* data, Attributes* sa, std::string fn_classifier);
    void write_rules_by_samples(std::string fn_out, Attributes* sa);
    void unit_tests(std::string dir_testdata);

    ~Classifier();

private:
    
	void evaluate();
    void process_rules();
    void emit(std::ofstream& f_out, std::string output);
    bool filter_out_idx(std::vector<int>* V, int target);
    bool assertEqual(bool a, bool b);
    bool assertEqual(int a, int b);
	void report(std::string testname, int tst, bool passed );
};
