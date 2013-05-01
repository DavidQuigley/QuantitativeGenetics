
enum {
	ID_TXT_FILE_NAME= 10000,
	ID_TXT_SUPPORT,
	ID_TXT_CONFIDENCE,
	ID_TXT_IMPROVEMENT,
	ID_TXT_MAX_CHI,
	ID_TXT_MAX_DEPTH,
	ID_CHO_SHOW_RESULTS,
	ID_DISCRETIZATION,
	ID_LBL_GENE_LIST,
	ID_BTN_GENE_LIST,
	ID_PANEL_A, 
	ID_PANEL_B
};

class RulesetDialog: public wxDialog{
	DECLARE_CLASS( RulesetDialog )

public:
	RulesetDialog( Investigation* investigation, std::string filename );
	void set_controls_from_file( std::string filename );
	void CreateControls();
	void OnClickOk( wxCommandEvent& event );
	void OnChooseDisc( wxCommandEvent& event );
	Investigation* investigation;

	wxTextCtrl* txt_file_name;
	wxTextCtrl* txt_support;
	wxTextCtrl* txt_confidence;
	wxTextCtrl* txt_improvement;
	wxTextCtrl* txt_max_chi;
	wxTextCtrl* txt_max_depth;
	wxChoice* cho_show_results;
    std::string new_filename;
	PanelDiscretization* discretization;
	
	PanelLimits* limit_A;
	PanelLimits* limit_B;
	std::vector< std::pair<std::string, std::string> > bounds;
	std::vector< std::pair<std::string, std::string> > bounds_per;

private:
	static const int BORDER_PXL = 2;
	static const int SMALL_BORDER=1;
	static const int FLAGS = wxALL | wxALIGN_LEFT | wxALIGN_CENTER_VERTICAL;
	static const int FLAGS_RIGHT = wxALL | wxALIGN_RIGHT | wxALIGN_CENTER_VERTICAL;
};
