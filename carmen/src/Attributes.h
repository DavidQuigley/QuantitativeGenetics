class Attributes{
public:
	std::vector<std::string> identifiers; // row labels
	std::vector<std::string> attrib;  // column labels
	HASH_S_I attribute2idx;
	HASH_S_I identifier2idx;
    std::string file_name;
	Attributes(std::string no_result);
	~Attributes();
	void find_identifiers_in_class(std::string limit, std::vector<std::string>* identifiers_in_class);
    void get_idx_by_genomic_range( std::string chromosome, int loc_start, int loc_end, std::vector<int>& matching_loci);
	std::string get_gene_name_column();
    std::string get_chr_column();
    bool get_chromosome_and_locus_by_identifier( std::string identifier, std::string& chromosome, int& locus);
	void insert(std::string attribute, int idx, std::string value);
	void indices_with_property(std::string attribute, std::string value, std::vector<int>* idx);
	void indices_with_property(std::string attribute, std::string value, std::vector<int>* idx, bool case_sensitive);
	void indices_without_property(std::string attribute, std::string value, std::vector<int>* idx);	
	void intersection( std::vector<int>* source, std::vector<int>* diff);
	void load(std::string file_name);
	std::string prop_for_identifier(std::string identifier, std::string attribute);
    void restrict_identifiers( std::vector<std::string> target );
	void set_gene_name_column(std::string col_name);
    void set_chromosome_and_locus( std::string col_chr, std::string col_locus );
	void write(std::string file_name);
	std::string no_result;
private:
	std::string col_gene_names, col_chr, col_locus;
	std::vector<std::string*> values;  // actual values of attributes
    // Use chromosome and locus to store chr:loc_start. 
    // Since not all rows of a gene attributes file will have a valid chromosome and locus,
    // length(chromosome) may not be equal to length(identifiers); store the actual index 
    // in idx_of_locus
    std::vector<std::string> chromosome;
    std::vector<int> locus;
    std::vector<int> idx_of_locus;
    
};
