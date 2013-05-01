enum{
	ID_DIFF_TXT_FILE_NAME = 1000,
	ID_DIFF_TXT_MAX_P,
	ID_DIFF_LIMIT_A,
	ID_DIFF_LIMIT_B,
	ID_DIFF_TRIM,
	ID_DIFF_CHK_DOPERMS,
	ID_DISCRETIZATION_PANEL
};


class DifferenceDialog: public wxDialog
{
	DECLARE_CLASS( DifferenceDialog )

public:
	DifferenceDialog( Investigation* investigation, std::string filename );
    std::string new_filename;
private:
	void OnClickOk( wxCommandEvent& evt );
	void OnCheckDoPerms( wxCommandEvent& evt);
	void OnUpdateLabel(wxCommandEvent& evt );
	Investigation* investigation;
	void CreateControls(std::string filename);
	void redraw();
	wxTextCtrl* txt_file_name;
	wxTextCtrl* txt_max_p;
    wxTextCtrl* txt_perms;
    wxTextCtrl* txt_trim;
    wxCheckBox* chk_doperms;
	PanelDiscretization* discretization;
	PanelLimits* limit_A;
	PanelLimits* limit_B;
	static const int BORDER_PXL = 2;
	static const int SMALL_BORDER=1;
	static const int FLAGS = wxALL | wxALIGN_LEFT | wxALIGN_CENTER_VERTICAL;
	static const int FLAGS_RIGHT = wxALL | wxALIGN_RIGHT | wxALIGN_CENTER_VERTICAL;
};
