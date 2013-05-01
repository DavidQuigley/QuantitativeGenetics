enum{
	CTL_BTN_REMOVE = 10010,
	CTL_BTN_ADD_SAMPLE,
	CTL_BTN_ADD_GENE,
	CTL_TXT_LIMITS,
	ID_CHO_ATTRIBUTES,
	ID_CHO_IS,
	ID_CHO_VALUES,
	
	ID_TXT_GENE_NAME,
	ID_LBX_IDENTIFIERS,
	ID_CHO_IS_ISNOT,
	ID_CHO_UP_DOWN,
	ID_BTN_SEARCH,
	ID_GOFILTER_OK,

	ID_GOFILTER_SEARCHTERMS,
	ID_GOFILTER_RESULTS,
	ID_GOFILTER_ANNOTATION,
	ID_GOFILTER_GENES,
	ID_GOFILTER_PROBES,
	ID_GOFILTER_BTNSEARCH
};


class GOFilterDialog: public wxDialog{
public:
	DECLARE_CLASS(GOFilterDialog)
	GOFilterDialog(Investigation* investigation);
	std::string get_probe_list();
	~GOFilterDialog();
	std::vector<std::string> genes, probes;
private:
	Investigation* investigation;
	wxTextCtrl* txt_searchterms;
	wxListBox* lbx_results;
	wxChoice* cho_annotation;
	wxTextCtrl* txt_genes;
	wxTextCtrl* txt_probes;
	wxButton* btn_search;
	wxButton* btn_cancel;
	wxButton* btn_ok;
	std::string loaded_annotation;
	wxStaticText* lbl_probes;
	wxStaticText* lbl_genes;
	GOAnnotationParser* GO_parser;
	GeneAnnotationParser* gene_parser;
	boost::unordered_map <string, std::vector<int>* > gene2probe_idx;
	wxArrayString ars_annotations;
	std::vector<std::string> fn_annotations;
	std::vector<int> result_GO_idx;

	void CreateControls();
	void OnClickSearch(wxCommandEvent& event );
	void OnClickOk(wxCommandEvent& event );
	void OnClickResults(wxCommandEvent& event );

	static const int BORDER_PXL = 4;
	static const int FLAGS = wxALL | wxALIGN_LEFT | wxALIGN_CENTER_VERTICAL;
	static const int FLAGS_TOP = wxALL | wxALIGN_LEFT | wxALIGN_TOP;
	static const int FLAGS_RIGHT = wxALL | wxALIGN_RIGHT | wxALIGN_CENTER_VERTICAL;
};

class GeneFilterDialog: public wxDialog
{
public:
	DECLARE_CLASS(GeneFilterDialog)
	GeneFilterDialog(Investigation* investigation);
	std::string limit;
private:
	Investigation* investigation;
	wxTextCtrl* txt_gene_name;
	wxChoice* cho_is;
	wxChoice* cho_up_down;
	wxListBox* lbx_identifiers;
	std::vector<std::pair<std::string, std::string>*> g_p;
	void CreateControls();
	void OnClickSearch(wxCommandEvent& event );
	void OnClickOk(wxCommandEvent& event );
	void OnCloseWindow(wxCommandEvent& event);

	static const int BORDER_PXL = 4;
	static const int FLAGS = wxALL | wxALIGN_LEFT | wxALIGN_CENTER_VERTICAL;
	static const int FLAGS_TOP = wxALL | wxALIGN_LEFT | wxALIGN_TOP;
	static const int FLAGS_RIGHT = wxALL | wxALIGN_RIGHT | wxALIGN_CENTER_VERTICAL;
};


class SampleFilterDialog: public wxDialog
{
public:
	DECLARE_CLASS(SampleFilterDialog)
	SampleFilterDialog(Investigation* investigation);
	std::string limit;
private:
	Investigation* investigation;
	wxChoice* cho_attributes;
	wxChoice* cho_is;
	wxChoice* cho_values;
	std::vector<std::string> values;
	void CreateControls();
	void OnChangeAttribute(wxCommandEvent& event );
	void OnClickOk(wxCommandEvent& event );
	void OnClickCancel(wxCommandEvent& event );
	void redraw();
	static const int BORDER_PXL = 4;
	static const int FLAGS = wxALL | wxALIGN_LEFT | wxALIGN_CENTER_VERTICAL;
	static const int FLAGS_RIGHT = wxALL | wxALIGN_RIGHT | wxALIGN_CENTER_VERTICAL;
};


class PanelLimits : public wxPanel{
public:
	PanelLimits(wxWindow* parent, int id, Investigation* investigation);
	void OnSelected( wxCommandEvent& event );
	void OnRemove( wxCommandEvent& event );
	void OnAddSample( wxCommandEvent& event );
	void OnAddGene( wxCommandEvent& event );
    void disable_gene_limits();
    void Clear();
	void set_limits(std::string limits);
    void set_investigation(std::string investigation_name);
	std::string get_limits();
    Investigation* investigation;

private:
	wxWindow* parent;
	wxListBox* lbx_limits;
	wxButton* btn_remove;
	wxButton* btn_add_sample;
	wxButton* btn_add_gene;
	wxBoxSizer* sizer_top;
	void redraw();
	static const int BORDER_PXL = 1;
	static const int FLAGS = wxALL | wxALIGN_LEFT | wxALIGN_CENTER_VERTICAL;
	static const int FLAGS_RIGHT = wxALL | wxALIGN_RIGHT | wxALIGN_CENTER_VERTICAL;
};
