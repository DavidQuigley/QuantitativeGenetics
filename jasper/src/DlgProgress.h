enum{
	ID_LBX_PROGRESS = 10000,
	ID_BTN_CANCEL,
	ID_PROCESS,
	ID_TIMER
};

#define POPEN2_READ 0
#define POPEN2_WRITE 1

class ProgressDialog: public wxDialog
{
	DECLARE_CLASS( ProgressDialog )

public:
	ProgressDialog( wxWindow* parent, std::vector<std::string>& cmd, Investigation* investigation, bool show_progress );
	void OnProcessEnded( wxProcessEvent& event );
	void OnIdle( wxIdleEvent& event );
	void OnTimer(wxTimerEvent& event);
	void OnClickCancel( wxCommandEvent& event );
private:
	wxListBox* lbx_progress;
	bool cancelled, show_progress;
	std::string error;
	wxProcess* process;
	std::string call;
	wxTimer* timer;
	int file_descriptor;
	long pid;
	FILE *stream;
	Investigation* investigation;
	void tick(wxTimerEvent& evt);
	long popen2(std::string cmd, int &file_descriptor);
	bool ReadFromProcess();
	void CloseDialog();
	static const int BORDER_PXL = 2;
	static const int FLAGS = wxALL | wxALIGN_LEFT | wxALIGN_CENTER_VERTICAL;
	static const int FLAGS_RIGHT = wxALL | wxALIGN_RIGHT | wxALIGN_CENTER_VERTICAL;
};
