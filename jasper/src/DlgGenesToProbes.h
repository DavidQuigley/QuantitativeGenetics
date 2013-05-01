enum{
	ID_G2P_BTN_CONVERT= 1000,
	ID_G2P_TXT_PROBES,
	ID_G2P_TXT_GENES
};

class GenesToProbesDialog: public wxDialog{
	DECLARE_CLASS( GenesToProbesDialog )
public:
	GenesToProbesDialog( Investigation* investigation);
	void CreateControls();

private:
	Investigation* investigation;
	wxTextCtrl* txt_genes;
	wxTextCtrl* txt_probes;
	wxButton* btn_convert;
	void OnClickConvert(wxCommandEvent& event );
	static const int BORDER_PXL = 5;
	static const int FLAGS = wxALL | wxALIGN_LEFT | wxALIGN_CENTER_VERTICAL;
	static const int FLAGS_RIGHT = wxALL | wxALIGN_RIGHT | wxALIGN_CENTER_VERTICAL;
	static const int FLAGS_CENTER = wxALL | wxALIGN_CENTER_VERTICAL;
};
