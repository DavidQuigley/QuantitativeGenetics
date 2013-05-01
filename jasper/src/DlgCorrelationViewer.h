enum{
	ID_COR_V_TXT_GENE_NAME = 1000,
	ID_COR_V_TXT_MIN_ABS_CORR,
	ID_COR_V_GRID_CORR,
	ID_COR_V_BTN_GO,
    ID_COR_V_BTN_SAVE
};

class CorrelationViewerDialog: public wxDialog
{
	DECLARE_CLASS( CorrelationViewerDialog )

public:
	CorrelationViewerDialog( std::string filename, Investigation* investigation );
	
private:
	wxTextCtrl* txt_gene_name;
	wxTextCtrl* txt_min_abs_corr;
	wxStaticText* lbl_current_probe_id_value;
	Investigation* investigation;
	wxButton* btn_go;
    wxButton* btn_save;
	wxButton* btn_close;
	GridCorrelation* grid_correlation;
	void CreateControls(std::string filename);
	void OnClickGo( wxCommandEvent& evt );
    void OnClickSave( wxCommandEvent& evt );
	//void OnClickClose( wxCommandEvent& evt );
	static const int BORDER_PXL = 2;
	static const int SMALL_BORDER=1;
	static const int FLAGS = wxALL | wxALIGN_LEFT | wxALIGN_CENTER_VERTICAL;
	static const int FLAGS_RIGHT = wxALL | wxALIGN_RIGHT | wxALIGN_CENTER_VERTICAL;
};
