class Rawdata{
public:
	Rawdata();
    std::string file_name;
	HASH_S_I identifier2idx;
	HASH_S_I sample_name2idx;
	std::vector<std::string> identifiers;
	std::vector<std::string> sample_names;
	Matrix<float>* data;

	void set_verbose(bool verbose);
	int index_of_identifier(std::string identifier);
	int index_of_sample(std::string sample_name);
	void load(std::string file_name);
	void unit_tests(std::string fn_data_dir);
	~Rawdata();

protected:
	void _load(std::string file_name, std::vector<float*>& rows, int& n_cols, bool &has_missing_values);
	float MISSING_DATA_FLOAT;
private:
	bool verbose;
};
