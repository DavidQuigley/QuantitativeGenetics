// Dataset keeps track of identifiers, sample_names separately from the raw_data 
// versions because we will probably restrict the set of identifiers or 
// samples under consideration.
class extract_sorter{
public:
	extract_sorter(double val, std::string s, int index);
	double val;
	int idx;
	std::string val_str;
};

class Dataset{
public:
	Dataset();

	std::vector<std::string> * identifiers;
	std::vector<std::string> * sample_names;
	std::vector<int>* a_idx;
	std::vector<int>* b_idx;
	bool a_is_target, b_is_target;
    std::string limits_a, limits_b;
	Rawdata* raw_data;
	Discretize* discrete;
	Attributes* sa;
	Attributes* ga;
	void load(Attributes* sa, Attributes* ga, ClassMinerOptions* cmo);
	void load(Attributes* sa, Attributes* ga, std::string fn_data, std::string limit_a, std::string limit_b);
	void extract(std::vector<std::vector<double>*>& values, std::vector<std::string>& samples, std::vector<std::string>& genes, std::vector<std::string> probes, std::string sample_limit, std::string group_by, bool sort_by_value, std::string label_attribute );
    void write_discretized(std::string fn, std::string separator, std::string discretization, float disc_lower, float disc_upper);
	~Dataset();

protected:
	void difference( std::vector<int>* source, std::vector<int>* diff);
	void intersection( std::vector<int>* source, std::vector<int>* diff);
	void check_sample_integrity(Rawdata* D, std::vector<std::string>* a, std::vector<std::string>* b);
	void check_gene_integrity(Rawdata* R, Attributes* ga );
	void restrict_identifiers_with_gene_disc(Rawdata* D, std::string class_limit, std::vector<std::string>* identifiers, Attributes* sa);
	void find_permitted_genes(HASH_S_I* permitted_genes, ClassMinerOptions* cmo, Attributes* ga);
	bool limits_contain_gene_restriction(std::string class_limit);
	void confirm_disjoint();
	void grouped_sort(std::vector<int>& idx_sorted, std::vector<std::string> samples, int idx_probe, std::string group_by, bool sort_within_groups);
	float MISSING_DATA_FLOAT;
};

// ClassifierDataset is a Numeric dataset with the addition of "feature" enumeration
// features is a binary matrix {feature x sample} with value 1 if feature is present, 
// 0 if features is absent.  We use translate to decode feature identity.
// Translate is a vector of two-element int arrays.  Element 0 is the row of the feature in 
// raw_data.  Element 1 is the discretized value of that row.  A single raw_data 
// row will typically have two rows in features and translate, corresponding to 
// the "up" and "down" states.
class ClassifierDataset : public Dataset{
public:
	ClassifierDataset();
	
	Matrix<int>* features; 
	std::vector<int*> * translate;
	
	void load(Attributes* sa, Attributes* ga, ClassMinerOptions* cmo); // overloaded
	void load_original(Attributes* sa, Attributes* ga, ClassMinerOptions* cmo); // deprecated
	void samples_with_gene_value(std::string identifier, int value, std::vector<int>* sample_idx);
	~ClassifierDataset();
};
