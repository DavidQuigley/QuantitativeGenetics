enum{
	ID_CLASS_A_TXT_FILE_NAME= 1000,
	ID_CLASS_A_CHO_INVESTIGATIONS,
	ID_CLASS_A_LIMIT_A,
	ID_CLASS_A_LIMIT_B,
    ID_CLASS_A_DISC
};

class ClassificationApplyDialog: public wxDialog{
	DECLARE_CLASS( ClassificationApplyDialog )

public:
	ClassificationApplyDialog( std::string filename, Investigation* investigation );
    std::string new_filename;

private:
    std::string filename;
    Investigation* investigation;
    std::vector<std::string> investigations;
    void reload_investigation();
    void OnClickOk(wxCommandEvent& WXUNUSED(event));
    void OnChooseInvestigation( wxCommandEvent& event );
    void CreateControls(std::string filename);
    wxTextCtrl* txt_file_name;
    wxChoice* cho_investigations;
    PanelDiscretization* discretization;
    PanelLimits* limit_A;
    PanelLimits* limit_B;
    static const int BORDER_PXL = 1;
	static const int FLAGS = wxALL | wxALIGN_LEFT | wxALIGN_CENTER_VERTICAL;
	static const int FLAGS_RIGHT = wxALL | wxALIGN_RIGHT | wxALIGN_CENTER_VERTICAL;
};
