enum{
	ID_CORR_NET_TXT_FILE_NAME = 1000,
	ID_CORR_NET_TXT_MIN_CORR,
	ID_CORR_NET_TXT_SEEDS,
    ID_CORR_NET_TXT_GENE,
    ID_CORR_NET_TXT_MIN_PRESENT,
	ID_CORR_NEW_FILE,
	ID_CORR_NET_PANEL_A,
    ID_CORR_NET_BTN_ADD_GENE,
    ID_CORR_NET_BTN_GO,
	ID_CORR_NET_CHK_LIMIT_SEEDS,
	ID_CORR_NET_CHO_MIN_CLIQUE,
	ID_CORR_NET_CHO_ATTRIBUTES,
	ID_CORR_NET_CHO_ANNOTATION,
	ID_CORR_NET_TXT_MAXIMUM_PVALUE,
	ID_CORR_NET_BTN_EQTL,
	ID_CORR_NET_CHK_REQUIRE_EQTL,
	ID_CORR_NET_CHK_ALLOW_UNCORRELATED
};


class CorrelationNetworkDialog: public wxDialog
{
	DECLARE_CLASS( CorrelationNetworkDialog )

public:
    CorrelationNetworkDialog( Investigation* investigation, std::string fn_spear );
    std::string batch_filename;

private:
	void OnClickOk( wxCommandEvent& evt );
    void OnClickAddGene( wxCommandEvent& evt );
    void OnClickGO( wxCommandEvent& evt );
    void OnClickEQTL( wxCommandEvent& evt );
    std::vector<std::pair<std::string, std::string>*> g_p;
    std::string fn_spear, fn_attributes, fn_go, fn_eqtl;
    std::vector<std::string> fn_annotations;
	Investigation* investigation;
	void CreateControls();
	wxBoxSizer* sizer_top;
	wxTextCtrl* txt_file_name, *txt_gene, *txt_maximum_pvalue;
    wxStaticText* lbl_eqtl_file_name_value, * lbl_eqtl_file_name, *lbl_maximum_pvalue, *lbl_require_eqtl, *lbl_allow_uncorrelated;
    wxButton* btn_add_gene, *btn_eqtl, *btn_GO;
    wxChoice* cho_attributes, *cho_annotation, *cho_new_file, *cho_min_clique;
    wxTextCtrl* txt_min_corr, *txt_seeds, *txt_min_present;
	PanelLimits* limit_A;
	wxCheckBox* chk_limit_seeds, *chk_require_eqtl, *chk_allow_uncorrelated;

	static const int BORDER_PXL = 2;
	static const int SMALL_BORDER=1;
	static const int FLAGS = wxALL | wxALIGN_LEFT | wxALIGN_CENTER_VERTICAL;
	static const int FLAGS_RIGHT = wxALL | wxALIGN_RIGHT | wxALIGN_CENTER_VERTICAL;
	static const int FLAGS_TOP_LEFT = wxALL | wxALIGN_LEFT | wxALIGN_TOP;
};
