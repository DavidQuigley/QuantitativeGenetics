enum {
	ID_INVEST_PROPERTIES_EXPR = 10000,
	ID_INVEST_PROPERTIES_GA,
	ID_INVEST_PROPERTIES_SA,
	ID_INVEST_PROPERTIES_RESULTS,
	ID_INVEST_PROPERTIES_SPECIES,
	ID_INVEST_PROPERTIES_DATA_TYPE,
	ID_INVEST_PROPERTIES_GENENAME
};

class InvestigationPropertiesPanel: public wxPanel{

public:
	InvestigationPropertiesPanel(wxWindow* parent, int id, Investigation* investigation);
	bool CheckConsistency(std::string & err);
	std::string get_fn_expr();
	std::string get_fn_sa();
	std::string get_fn_ga();
	std::string get_dir_results();
	std::string get_species();
	std::string get_data_type();
	std::string get_gene_name_column();
	void set_fn_expr(std::string fn);
	void set_fn_sa(std::string fn);
	void set_fn_ga(std::string fn);
	void set_dir_results(std::string fn);
	void set_species(std::string fn);
	void set_data_type(std::string fn);
	void set_enabled(bool is_enabled);
	void set_gene_name_column(std::string gnc);

private:
	wxWindow* parent;
	Investigation* investigation;
	wxFlexGridSizer* sizer_top;
	wxStaticText* lbl_expr_value;
	wxStaticText* lbl_ga_value;
	wxStaticText* lbl_sa_value;
	wxStaticText* lbl_results_value;
	wxStaticText* lbl_gene_name_column;
	wxButton* btn_raw_data;
	wxButton* btn_ga;
	wxButton* btn_sa;
	wxButton* btn_results;
	wxChoice* cho_gene_name_column;
	wxChoice* cho_data_type;
	wxChoice* cho_attributes;
	wxChoice* cho_annotation;
	std::vector<std::string> gene_name_column;
	wxArrayString ars_annotations;
	wxArrayString ars_data_types;
	std::string most_recent_folder;
	void redraw();
	void OnClickRawData( wxCommandEvent& event );
	void OnClickSA( wxCommandEvent& event );
	void OnClickGA( wxCommandEvent& event );
	void OnClickResults( wxCommandEvent& event );
	static const int BORDER_PXL = 1;
	static const int FLAGS = wxALL | wxALIGN_LEFT | wxALIGN_CENTER_VERTICAL;
	static const int FLAGS_RIGHT = wxALL | wxALIGN_RIGHT | wxALIGN_CENTER_VERTICAL;
};
