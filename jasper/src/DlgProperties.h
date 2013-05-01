
enum {
	ID_PROPERTIES_CYTOSCAPE = 10000,
	ID_PROPERTIES_VIZMAP,
	ID_PROPERTIES_EXCEL,
	ID_PROPERTIES_TEXT_EDITOR,
	ID_PROPERTIES_INVESTIGAION,
	ID_PROPERTIES_GO,
	ID_PROPERTIES_ANNOTATION_CHO,
	ID_PROPERTIES_ANNOTATION_BTN

};

class PropertyDialog: public wxDialog
{
public:
	DECLARE_CLASS( PropertyDialog )
	PropertyDialog( Investigation* investigation);
	void CreateControls();
	std::vector<std::string> new_annotation_files;
	std::vector<std::string> new_annotation_names;

private:
	Investigation* investigation;
	wxArrayString ars_annot;
	std::vector<std::string> fn_annotations;
	wxStaticText* lbl_cytoscape_value;
	wxStaticText* lbl_vizmap_value;
	wxStaticText* lbl_excel_value;
	wxStaticText* lbl_text_editor_value;
	wxStaticText* lbl_go_value;
	wxStaticText* lbl_annotation_value;
	wxChoice* cho_annotation;

	void OnClickCytoscape(wxCommandEvent& event );
	void OnClickVizmap(wxCommandEvent& event );
	void OnClickExcel(wxCommandEvent& event );
	void OnClickTextEditor(wxCommandEvent& event );
	void OnClickGO(wxCommandEvent& event );
	void OnChangeAnnotation( wxCommandEvent& event );
	void OnNewAnnotation(wxCommandEvent& event );
	static const int BORDER_PXL = 5;
	static const int FLAGS = wxALL | wxALIGN_LEFT | wxALIGN_CENTER_VERTICAL;
	static const int FLAGS_RIGHT = wxALL | wxALIGN_RIGHT | wxALIGN_CENTER_VERTICAL;
};
