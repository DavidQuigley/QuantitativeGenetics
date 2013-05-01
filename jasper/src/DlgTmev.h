enum{
	ID_TMEV_TXT_FILE_NAME = 1000,
	ID_TMEV_TXT_MIN_PERCENT_PRESENT,
	ID_TMEV_BTN_GO,
	ID_TMEV_RDO_ALL,
    ID_TMEV_RDO_FILE,
    ID_TMEV_RDO_PASTE,
    ID_TMEV_TXT_PROBES,
	ID_TMEV_LIMIT_A
};


class TmevDialog: public wxDialog
{
	DECLARE_CLASS( TmevDialog )

public:
	TmevDialog( Investigation* investigation, std::string filename );
    std::string new_filename;
private:
	void OnClickOk( wxCommandEvent& evt );
	void OnClickGO(wxCommandEvent& evt );
	void OnUpdateLabel(wxCommandEvent& evt );
    void OnChangeAttribute(wxCommandEvent& evt);
    Investigation* investigation;
	void CreateControls(std::string filename);
	void redraw();
    std::string rs_filename;
	//wxTextCtrl* txt_file_name;
	wxTextCtrl* txt_min_percent_present;
    //wxTextCtrl* txt_perms;
    wxTextCtrl* txt_probes;
    wxRadioButton* rdo_all;
    wxButton* btn_GO;
    wxRadioButton* rdo_file;
    wxRadioButton* rdo_paste;
	PanelLimits* limit_A;
	wxBoxSizer* sizer_top;
	static const int BORDER_PXL = 2;
	static const int SMALL_BORDER=1;
	static const int FLAGS = wxALL | wxALIGN_LEFT | wxALIGN_CENTER_VERTICAL;
	static const int FLAGS_RIGHT = wxALL | wxALIGN_RIGHT | wxALIGN_CENTER_VERTICAL;
	static const int FLAGS_TOP_LEFT = wxALL | wxALIGN_LEFT | wxALIGN_TOP;
};
