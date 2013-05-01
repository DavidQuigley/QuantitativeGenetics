enum{
	ID_CORR1P_TXT_FILE_NAME = 1000,
	ID_CORR1P_PROBE_ID,
    ID_CORR1P_TXT_MIN_CORR,
	ID_CORR1P_PANEL_A,
	ID_CORR1P_RADIO_SINGLE
};

class CorrelationOneProbeDialog: public wxDialog
{
	DECLARE_CLASS( CorrelationOneProbeDialog )

public:
	CorrelationOneProbeDialog( Investigation* investigation );
    std::string new_filename;

private:
	void OnClickOk( wxCommandEvent& evt );
	void OnUpdateLabel(wxCommandEvent& evt );
	
	Investigation* investigation;
	void CreateControls();
	wxTextCtrl* txt_file_name;
	wxTextCtrl* txt_probe_id;
	wxTextCtrl* txt_min_corr;
	static const int BORDER_PXL = 2;
	static const int SMALL_BORDER=1;
	static const int FLAGS = wxALL | wxALIGN_LEFT | wxALIGN_CENTER_VERTICAL;
	static const int FLAGS_RIGHT = wxALL | wxALIGN_RIGHT | wxALIGN_CENTER_VERTICAL;
};
