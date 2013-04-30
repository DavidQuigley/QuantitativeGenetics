class GOAnnotation{
public:
	enum GO_branch{
		GO_BP=0,
		GO_MF,
		GO_CC
	};
	int idx;
	std::string go_id;
	GO_branch branch; // one of GO_BP, GO_MF, GO_CC
	std::string description;
	GOAnnotation();
	GOAnnotation( int idx, std::string go_id, GO_branch branch, std::string description );
};


class GOAnnotationParser{
	// Reads in a GO annotation file
	// Meant to be used as a lookup in conjunction with GeneAnnotationParser
public:
	GOAnnotationParser();
	GOAnnotationParser(std::string fn);
	void load(std::string fn);
	GOAnnotation* get_annotation(int idx);
	int get_idx(std::string go_id);
	void get_idx_matching_description(std::string searchterm, std::vector<int>& results);
	int size();
	~GOAnnotationParser();

private:
	std::vector<GOAnnotation*> GO;
	HASH_S_I go_id2idx;
};


class GeneAnnotation{
public:

	std::string symbol;
	std::string chr;
	std::string ucsc_chr;
	std::string full_name;
	std::string loc_start;
	std::vector< int > GO_annotations;
	// to reduce storage costs we load the annotations once and the refer to them by idx
	// these idx values are not persistent.

	GeneAnnotation(std::string symbol, std::string chr, std::string ucsc_chr, std::string full_name, std::string loc_start, std::vector<int>& GO );
};

class GeneAnnotationParser{
public:
	GeneAnnotationParser();
	GeneAnnotationParser(std::string fn, GOAnnotationParser& GOP);
	void load(std::string fn, GOAnnotationParser& GOP);
	GeneAnnotation* get_annotation(std::string symbol);
	void get_symbols_with_GO_idx(int target_go_idx, std::vector<std::string>& symbols);
private:
	boost::unordered_map<std::string, GeneAnnotation*> symbol2annotation;
};


class ConfigParser{
public:
	typedef std::pair<std::string, std::string> prop_pair;

	ConfigParser();
	ConfigParser(std::string fn_props);
	void add_section(std::string section_name);
	void clear();
	void create(std::string fn);
	void load(std::string fn);
	void set(std::string section, std::string A, std::string B);
	std::string get(std::string section, std::string A);
	void get_section(std::string section, std::vector< prop_pair* >*);
	void write();
	void write(std::string fn);
	void remove(std::string section, std::string A);
	~ConfigParser();
	boost::unordered_map<std::string, std::vector< prop_pair* >* > contents;
private:
	std::string fn_loaded;
	std::vector<std::string> sections;
};


class CarmenParser{
public:
	typedef boost::tokenizer<boost::char_separator<char> > tokenizer;

	CarmenParser();
	CarmenParser(std::string filename);
	
	std::string fn_loaded;

	int CountValueLines();
	bool get_value(std::string header, int idx, std::string& s);
	bool get_value_vector( std::string header_keyword, std::vector<std::string>* header_values );
    void get_column_names( std::vector<std::string>& column_names );
	void Load(std::string filename);
	void PrepareToReadValues();
	void ExtractProbes(std::vector<std::string>& probes);
	void ExtractGenes(std::vector<std::string>& genes);
	void ExtractGeneCounts(HASH_S_I& hsh_genes);
	void ExtractEQTL(std::vector<QTL*>& qtls, double max_eqtl_pval);
	bool ReadNextValue(std::vector<std::string>& values);
	void write_compressed_spear(std::string fn_out, HASH_S_I& p2id);

private:
	boost::unordered_map<std::string, std::vector<std::string>* > header2values;
	std::vector<std::string> column_names;
	bool file_has_values;
	bool has_CRLF;
	std::string extension;
	std::ifstream fs;
};

class KeyValuesParser{
// tab-delimited file where first element of each line is a key and following elements are values
public:
    KeyValuesParser(std::string fn_kv);
    void keys(std::vector<std::string>& k);
    void get(std::string key, std::vector<std::string>& values);
    bool has_key(std::string key);
private:
    std::vector<std::string> keys_in_order; // used to track order
    boost::unordered_map<std::string, std::vector<std::string>* > key2values;
    bool has_CRLF;
};
