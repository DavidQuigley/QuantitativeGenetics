enum{
	ID_CLASS_TXT_FILE_NAME = 1000,
	ID_CLASS_CHO_METHOD,
	ID_CLASS_TXT_ROUNDS,
	ID_CLASS_TXT_MAX_RULES,
	ID_CLASS_CHK_LOOCV
};


class ClassifierDialog: public wxDialog
{
	DECLARE_CLASS( ClassifierDialog )

public:
	ClassifierDialog( Investigation* investigation, std::string filename );
    std::string new_filename;
private:
	std::string filename;
	void OnClickOk( wxCommandEvent& evt );
	void OnClickMethod(wxCommandEvent& evt );
	Investigation* investigation;
	void CreateControls(std::string filename);
	void redraw();
	wxStaticText* lbl_rounds;
	wxStaticText* lbl_max_rules;
	wxTextCtrl* txt_file_name;
	wxTextCtrl* txt_rounds;
	wxTextCtrl* txt_max_rules;
	wxChoice* cho_method;
	wxCheckBox* chk_loocv;
	static const int BORDER_PXL = 2;
	static const int SMALL_BORDER=1;
	static const int FLAGS = wxALL | wxALIGN_LEFT | wxALIGN_CENTER_VERTICAL;
	static const int FLAGS_RIGHT = wxALL | wxALIGN_RIGHT | wxALIGN_CENTER_VERTICAL;
};
