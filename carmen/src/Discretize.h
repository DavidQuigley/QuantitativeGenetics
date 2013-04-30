class Discretize{
public:
	Discretize(Rawdata* raw_data, std::vector<int>& idx, std::string method_str, float value_a, float value_b);
	~Discretize();
    std::string method;
    float mult_upper, mult_lower;
	Matrix<int>* dis;
	std::vector<float> cutoff_up;
	std::vector<float> cutoff_dn;

	float find_median(std::vector<float>* C);
	float find_mean(std::vector<float>* C);
	float find_stdev(float mean, std::vector<float>* C);

private:

	float find_mad(float median, std::vector<float>* C);
	float find_var(float mean, std::vector<float>* C);
	float find_t(float mu_a, float mu_b, float var_a, float var_b, float n_a, float n_b);
	void discretize_SD(float value_a, float value_b, std::vector<int>& idx, Rawdata* rd);
	void discretize_SD_Samples(float value_a, float value_b, std::vector<int>& idx, Rawdata* rd);
	void discretize_MAD(float value_a, float value_b, std::vector<int>& idx, Rawdata* rd);
	void discretize_PER(float value_a, float value_b, std::vector<int>& idx, Rawdata* rd);
	void discretize_ABS(float value_a, float value_b, std::vector<int>& idx, Rawdata* rd);
	void discretize_NONE(std::vector<int>& idx, Rawdata* rd);

	//void discretize_ttest(Rawdata* rd, std::vector<std::string>* ids_a, std::vector<std::string>* ids_b);
	int disc(float f, float dn, float up);
	inline int round(float x);
};
