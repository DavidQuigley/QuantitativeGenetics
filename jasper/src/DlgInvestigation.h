#include "wx/dialog.h"

// Control identifiers

enum {
	ID_INVESTIGATION = 10000
};

class InvestigationDialog: public wxDialog{
public:
	DECLARE_CLASS( InvestigationDialog )
	InvestigationDialog(Investigation* investigation);
	std::string fn_expr;
	std::string fn_ga;
	std::string fn_sa;
	std::string dir_results;
	std::string investigation_name;
	std::string data_type;
	std::string species;
	std::string gene_name_column;
	void OnClickOk( wxCommandEvent& event );
	void CreateControls();

private:
	wxTextCtrl* txt_investigations;
	Investigation* investigation;
	InvestigationPropertiesPanel* pnl_investigation;
	static const int BORDER_PXL = 4;
	static const int FLAGS = wxALL | wxALIGN_LEFT | wxALIGN_CENTER_VERTICAL;
	static const int FLAGS_RIGHT = wxALL | wxALIGN_RIGHT | wxALIGN_CENTER_VERTICAL;
};
