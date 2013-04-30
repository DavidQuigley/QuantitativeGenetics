class perm_result{
public:
	std::string name;
	std::string identifier;
	float var_a;
	float var_b;
	float mu_a;
	float mu_b;
	float t_val;
	float p_value;
	int idx;
};

class Perm{
public:
	Perm();

	int n_perm;
	int mean_trim;
	double max_p_value;
	int n_a;
	int n_b;
	std::string file_name_dis;
	std::string file_name_ga;
	std::string file_name_sa;
	std::string class_a;
	std::string class_b;

	void load_settings(std::string diff_file_name, double max_p_value);
	float find_p_value( float r, vector< vector<int>*> * perms_a, vector< vector<int>*> * perms_b, Matrix<float>* data, int row, int trim_percent);
	float find_t_stat( float mu_a, float var_a, int n_a, float mu_b, float var_b, int n_b);
    void permutations( Rawdata* rd, std::vector<int> & idx, vector<int>* idx_a, vector<int>* idx_b, Attributes* ga, bool verbose);
	void write( ClassMinerOptions* cmo );
	void write_as_ruleset( ClassMinerOptions* cmo, Rawdata* rd, vector<int>* idx_a, vector<int>* idx_b );
	void trimmed_stats(Matrix<float>* data, int row, vector<int>* cols, int trim_percent, float& mean, float& variance);
	void limit_ids_by_NA(std::vector<int> &idx, ClassifierDataset* data, double fraction_required);
	int rnd(int aRange);
	float find_variance(float mean, const vector<float>& acc, int idx_start, int idx_end);
	void shuffle(vector<int>* V);

private:
	
	std::vector<perm_result*> results;
	void permute( vector<int>* idx_a, vector<int>* idx_b, vector<int>* perm_a, vector<int>* perm_b);
	inline int round(float f);

};
