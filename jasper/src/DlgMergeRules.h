enum{
	ID_TXT_THRESHOLD = 1000,
	ID_TXT_MERGE_RULE_TREE,
	ID_BTN_REFRESH,
	ID_BTN_SAVE_MERGED,
	ID_LBL_RESULTS
};


class MergeRulesDialog: public wxDialog{
	DECLARE_CLASS( MergeRulesDialog )
public:
	MergeRulesDialog( Investigation* investigation, std::string filename );
	bool cancel_load;
	bool saved_new_file;
    std::string new_filename;

private:
	Matrix<double>* r_by_s;
	RuleSet* ruleset;
	std::vector< std::vector<int>* > rules;
	Investigation* investigation;
	double threshold;	
	void merge_rules( double threshold );
	void merge_rules_old( double threshold );
	bool CreateControls(std::string filename);
	void OnClickOk( wxCommandEvent& WXUNUSED(event) );
	void OnClickSaveMerged( wxCommandEvent& WXUNUSED(event) );
	void OnClickRefresh( wxCommandEvent& WXUNUSED(event) );
	void redraw();
	std::string filename;
	wxImageList imagelist;
	wxTreeItemId root;
	wxTextCtrl* txt_threshold;
	wxTreeCtrl* tree;
	wxStaticText* lbl_results;
	int img_folder;
	int img_file_open;
	int img_file_normal;
	
	static const int BORDER_PXL = 2;
	static const int SMALL_BORDER=1;
	static const int FLAGS = wxALL | wxALIGN_LEFT | wxALIGN_CENTER_VERTICAL;
	static const int FLAGS_RIGHT = wxALL | wxALIGN_RIGHT | wxALIGN_CENTER_VERTICAL;
};
