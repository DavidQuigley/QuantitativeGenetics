class GridCorrelation : public wxGrid{
public:
	GridCorrelation(wxWindow* parent, Investigation* investigation, std::string filename, wxWindowID id, wxSize size);
	void Update(std::string gene_name, double min_abs_corr);
	bool write_to_file(std::string filepath);
    std::string current_probe_id;

private:
	wxWindow* m_parent;
	std::string filename;
	Graph G;
	Investigation* investigation;
	boost::unordered_map<int, std::string> id2p;
	HASH_S_S p2g;
	
	void Clear();
	void read_file_into_graph(std::string filename);
    
	static const int BORDER_PXL = 5;
	static const int FLAGS = wxEXPAND | wxALL | wxALIGN_LEFT | wxALIGN_CENTER_VERTICAL;
};
