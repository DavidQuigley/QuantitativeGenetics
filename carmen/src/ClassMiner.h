class ClassMiner{
public:
	void append_L1_L2_antecedents(ClassifierDataset * data, std::vector<Ant*>* ants, ClassMinerOptions* cmo, int& n_tests);
	void append_L_ants(ClassifierDataset * data, std::vector<Ant*>* ants, ClassMinerOptions* cmo, int L, int& n_tests);
	void assign_t_stats(std::vector<Ant*>* ants, ClassifierDataset* data);
	void filter_ants(std::vector<Ant*>* ants, ClassifierDataset * data, float min_conf, int min_sup, float min_imp, double max_chi2);
	void print_ants_stdout(std::vector<Ant*>* ants, ClassifierDataset* data, ClassMinerOptions* cmo);
	void print_ants_file(std::vector<Ant*>* ants, ClassifierDataset* data, ClassMinerOptions* cmo);
private:
	inline int round(float x);
	inline int round_up(float x);
};
