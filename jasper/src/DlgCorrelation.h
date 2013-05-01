enum{
	ID_CORR_TXT_FILE_NAME = 1000,
	ID_CORR_TXT_MIN_VAR,
	ID_CORR_TXT_MIN_CORR,
	ID_CORR_TXT_DELTA_CORR,
	ID_CORR_TXT_MIN_PRESENT,
	ID_CORR_CHK_DIFFERENTIAL,
	ID_CORR_CHO_METHOD,
	ID_CORR_BTN_GO,
	ID_CORR_PROBELIST,
	ID_CORR_PANEL_A,
	ID_CORR_PANEL_B,
	ID_RADIO_SINGLE,
	ID_PROBE_ID
};


class CorrelationDialog: public wxDialog
{
	DECLARE_CLASS( CorrelationDialog )

public:
	CorrelationDialog( Investigation* investigation );
    std::string new_filename;

private:
	void OnClickOk( wxCommandEvent& evt );
	void OnClickGO( wxCommandEvent& evt );
	void OnUpdateLabel(wxCommandEvent& evt );
	void OnClickSingle(wxCommandEvent& evt );
	void OnCheckDifferential(wxCommandEvent& evt );
	void OnChar( wxKeyEvent& event );
	Investigation* investigation;
	void CreateControls();
	void redraw();
	wxTextCtrl* txt_file_name;
	wxTextCtrl* txt_min_var;
	wxTextCtrl* txt_min_corr;
	wxTextCtrl* txt_min_delta_corr;
	wxChoice* cho_method;
	wxStaticText* lbl_min_delta_corr;
	wxStaticText* lbl_probe_id;
	wxStaticText* lbl_condition1;
	wxStaticText* lbl_condition2;
	wxStaticText* lbl_probe_list;
	wxTextCtrl* txt_probe_list;
	wxButton* btn_GO;

	wxTextCtrl* txt_min_present;
	wxCheckBox* chk_differential;
	wxRadioBox* rdo_single;
	wxTextCtrl* txt_probe_id;
	PanelLimits* limit_A;
	PanelLimits* limit_B;
	wxBoxSizer* sizer_top;
	static const int BORDER_PXL = 2;
	static const int SMALL_BORDER=1;
	static const int FLAGS = wxALL | wxALIGN_LEFT | wxALIGN_CENTER_VERTICAL;
	static const int FLAGS_RIGHT = wxALL | wxALIGN_RIGHT | wxALIGN_CENTER_VERTICAL;
	static const int FLAGS_TOP_LEFT = wxALL | wxALIGN_LEFT | wxALIGN_TOP;
};
