enum {
	DISC_DLG_DISC_PANEL = 10000,
    DISC_DLG_DISC_LIMIT_A,
    DISC_DLG_FILE_NAME
};

class DiscretizationDialog: public wxDialog
{
	DECLARE_CLASS( DiscretizationDialog )

public:
    DiscretizationDialog(wxWindow* parent, Investigation* investigation );
	std::string disc_method, disc_upper, disc_lower, fn_discretization;
    std::string limits;
private:
    Investigation* investigation;
    void OnClickOk( wxCommandEvent& evt );
	wxTextCtrl* txt_file_name;
    PanelDiscretization* panel_disc;
    PanelLimits* limit_A;

	static const int BORDER_PXL = 5;
	static const int FLAGS = wxALL | wxALIGN_LEFT | wxALIGN_CENTER_VERTICAL;
    static const int FLAGS_RIGHT = wxALL | wxALIGN_RIGHT | wxALIGN_CENTER_VERTICAL;
};
