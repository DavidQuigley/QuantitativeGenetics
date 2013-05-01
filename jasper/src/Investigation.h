class Investigation{
	
public:
	typedef std::pair<std::string, std::string> prop_pair;

	Investigation();
	Investigation( std::string fn_properties );
	void clone_from( Investigation* source );
    void read_from_file( std::string fn_properties );
	void set_current( std::string investigation_name );
	void write_dataset( std::string fn_out, std::vector<std::string> probes, std::vector<std::string> samples, float percent_req );
	~Investigation();
	std::string investigation_name;
	std::string fn_properties;
	std::string fn_expr;
	std::string fn_ga;
	std::string fn_sa;
	std::string dir_results;
	std::string species;
	std::string data_type;
	std::string gene_name_column;
	bool is_mac; // true if loaded on a mac
	Attributes* sa;
	Attributes* ga;
	ConfigParser* cp;
	HASH_S_S nicknames;

private:
	void read_experiments_directory();
	void probes_with_percent_req(Rawdata& expr, boost::unordered_map<int,int>& probe_idx_present, float percent_req, std::vector<int>& idx_expr);
};
