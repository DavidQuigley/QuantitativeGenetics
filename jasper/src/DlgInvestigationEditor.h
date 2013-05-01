enum{
	ID_INVESTIGATION_EDITOR_INVESTIGAION = 10000
};

class InvestigationEditorDialog: public wxDialog
{
public:
	DECLARE_CLASS( InvestigationEditorDialog )
	InvestigationEditorDialog( Investigation* investigation);
	void CreateControls();
	std::string fn_expr, fn_ga, fn_sa, dir_results;
	std::string species, data_type, investigation_name, gene_name_column;

private:
	Investigation* investigation;
	InvestigationPropertiesPanel* pnl_investigation;
	std::vector<std::string> values;
	void OnClickOK(wxCommandEvent& event );
	static const int BORDER_PXL = 5;
	static const int FLAGS = wxALL | wxALIGN_LEFT | wxALIGN_CENTER_VERTICAL;
	static const int FLAGS_RIGHT = wxALL | wxALIGN_RIGHT | wxALIGN_CENTER_VERTICAL;
};
