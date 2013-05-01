enum{
	ID_ANNOTATION_FILE_NAME = 1000,
	ID_ANNOTATION_ANNOTATION,
    ID_ANNOTATION_RDO_ALL,
    ID_ANNOTATION_RDO_FILE,
    ID_ANNOTATION_RDO_PASTE,
    ID_ANNOTATION_PROBES,

};


class AnnotationDialog: public wxDialog
{
	DECLARE_CLASS( AnnotationDialog )

public:
	AnnotationDialog( Investigation* investigation, std::string filename );
    std::string new_filename;
private:
	void OnClickOk( wxCommandEvent& evt );
	void OnClickCancel( wxCommandEvent& evt );
    void OnChangeAttribute(wxCommandEvent& evt);
    Investigation* investigation;
	void CreateControls(std::string filename);
	std::string rs_filename;
	//wxTextCtrl* txt_file_name;
	wxChoice* cho_annotation;
	wxArrayString ars_annotations;
	std::vector<std::string> fn_annotations;
    wxTextCtrl* txt_probes;
    wxRadioButton* rdo_file;
    wxRadioButton* rdo_paste;
	static const int BORDER_PXL = 2;
	static const int SMALL_BORDER=1;
	static const int FLAGS = wxALL | wxALIGN_LEFT | wxALIGN_CENTER_VERTICAL;
	static const int FLAGS_RIGHT = wxALL | wxALIGN_RIGHT | wxALIGN_CENTER_VERTICAL;
};
