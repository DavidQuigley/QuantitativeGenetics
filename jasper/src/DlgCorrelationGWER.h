enum{
	ID_CORR_GWER_TXT_FILE_NAME = 1000,
	ID_CORR_GWER_N_PERMS,
	ID_CORR_GWER_PANEL_A
};

class CorrelationGWERDialog: public wxDialog
{
	DECLARE_CLASS( CorrelationGWERDialog )

public:
	CorrelationGWERDialog( Investigation* investigation );
    std::string new_filename;
private:
    Investigation* investigation;
    wxArrayString ars_perms;
    wxTextCtrl* txt_file_name;
    wxChoice * cho_n_perms;
	PanelLimits* limit_A;
	static const int BORDER_PXL = 2;
	static const int SMALL_BORDER=1;
	static const int FLAGS = wxALL | wxALIGN_LEFT | wxALIGN_CENTER_VERTICAL;
	void OnClickOk( wxCommandEvent& evt );
	void CreateControls();
};
