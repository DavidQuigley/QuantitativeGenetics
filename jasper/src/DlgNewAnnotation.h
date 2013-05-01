//#include "wx/dialog.h"

enum{
	ID_NEW_ANNOTATION_SET=10000,
	ID_NEW_ANNOTATION_NAME
};

class NewAnnotationDialog: public wxDialog{
public:
	DECLARE_CLASS( NewAnnotationDialog )
	NewAnnotationDialog(Investigation* investigation);
	std::string fn_annotation;
	std::string name_annotation;
	void OnClickOk( wxCommandEvent& event );
	void OnClickSet( wxCommandEvent& event );
	void CreateControls();

private:
	Investigation* investigation;
	wxBoxSizer* sizer_top;
	wxTextCtrl* txt_annotation_name;
	wxStaticText* lbl_fn_annotation_value;
	wxButton* btn_set;
	static const int BORDER_PXL = 4;
	static const int FLAGS = wxALL | wxALIGN_LEFT | wxALIGN_CENTER_VERTICAL;
	static const int FLAGS_RIGHT = wxALL | wxALIGN_RIGHT | wxALIGN_CENTER_VERTICAL;
};
