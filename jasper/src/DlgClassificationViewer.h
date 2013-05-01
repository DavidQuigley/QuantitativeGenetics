enum{
	ID_CLASS_V_GRID_CORR = 1000,
	ID_CLASS_V_BTN_EXPORT_SA
};

class ClassificationViewerDialog: public wxDialog
{
	DECLARE_CLASS( ClassificationViewerDialog )

public:
	ClassificationViewerDialog( std::string filename, Investigation* investigation );
private:
    Investigation* investigation;
    wxButton* btn_close;
    GridClassifierResults* grid_classification;
    void CreateControls(std::string filename);
    
static const int BORDER_PXL = 2;
	static const int SMALL_BORDER=1;
	static const int FLAGS = wxALL | wxALIGN_LEFT | wxALIGN_CENTER_VERTICAL;
	static const int FLAGS_RIGHT = wxALL | wxALIGN_RIGHT | wxALIGN_CENTER_VERTICAL;
};
