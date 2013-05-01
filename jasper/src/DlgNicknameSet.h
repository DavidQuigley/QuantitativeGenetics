enum {
	ID_TXT_NICKNAME = 10000
};

class NicknameSetDialog: public wxDialog
{
	DECLARE_CLASS( NicknameSetDialog )

public:
	NicknameSetDialog( );
	NicknameSetDialog( wxWindow* parent, 
		wxWindowID id = wxID_ANY,
		const wxString& caption = wxT("Set a nickname for this query"),
		const wxPoint& pos = wxDefaultPosition,
		const wxSize& size = wxDefaultSize,
		long style = wxCAPTION|wxRESIZE_BORDER|wxSYSTEM_MENU 
	);	
	bool Create( wxWindow* parent,
		wxWindowID id = wxID_ANY,
		const wxString& caption = wxT("Investigation Browser"),
		const wxPoint& pos = wxDefaultPosition,
		const wxSize& size = wxDefaultSize,
		long style = wxCAPTION|wxRESIZE_BORDER|wxSYSTEM_MENU 
	);
	
	void CreateControls();
	void OnClickOk( wxCommandEvent& event );
	std::string new_nickname;
	static const int BORDER_PXL = 5;
	static const int FLAGS = wxALL | wxALIGN_LEFT | wxALIGN_CENTER_VERTICAL;
};
