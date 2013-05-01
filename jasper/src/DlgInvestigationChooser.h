
enum {
	ID_CHO_INVESTIGATIONS = 10000,
    ID_CHO_SORTBY = 10001

};

class InvestigationChooserDialog: public wxDialog
{
	DECLARE_CLASS( InvestigationChooserDialog )

public:
	InvestigationChooserDialog( );
	InvestigationChooserDialog( wxWindow* parent, 
		ConfigParser* cp,
		wxWindowID id = wxID_ANY,
		const wxString& caption = wxT("Choose an Investigation to open"),
		const wxPoint& pos = wxDefaultPosition,
		const wxSize& size = wxDefaultSize,
		long style = wxCAPTION|wxRESIZE_BORDER|wxSYSTEM_MENU 
	);
	std::string investigation;

private:
	void Init(ConfigParser* cp);
	
	bool Create( wxWindow* parent,
		wxWindowID id = wxID_ANY,
		const wxString& caption = wxT("Investigation Browser"),
		const wxPoint& pos = wxDefaultPosition,
		const wxSize& size = wxDefaultSize,
		long style = wxCAPTION|wxRESIZE_BORDER|wxSYSTEM_MENU 
	);
	void CreateControls();
	wxButton* btn_Ok;
	wxButton* btn_Cancel;
	wxChoice* cho_investigations;
	wxChoice* cho_sortby;
    void OnChangeInvestigation( wxCommandEvent& event );
	void OnChangeSortBy( wxCommandEvent& event );
    void OnChar( wxKeyEvent& event );
	std::vector<std::string> investigations;
    std::vector<std::string> investigations_current_sorted_order;

	static const int BORDER_PXL = 5;
	static const int FLAGS = wxALL | wxALIGN_LEFT | wxALIGN_CENTER_VERTICAL;
	static const int FLAGS_RIGHT = wxALL | wxALIGN_RIGHT | wxALIGN_CENTER_VERTICAL;
};
