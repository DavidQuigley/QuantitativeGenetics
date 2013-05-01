enum{
	ID_CHOICE_DISC = 1000,
	ID_CHOICE_DISC_BOUNDS
};
enum{
	DISC_NONE=0,
	DISC_SD,
	DISC_MAD,
	DISC_ABS,
	DISC_PER,
	DISC_SD_SAMPLES
};
class PanelDiscretization: public wxPanel{
public:
	PanelDiscretization(wxWindow* parent, int id, Investigation* investigation);
	int get_discretization_code();
    void get_discretization_method_name(std::string& method);
	void get_discretization_bounds(std::string& lower, std::string& upper);
	void set_discretization(std::string& disc, std::string& lower, std::string& upper);

private:
	void redraw();
	void OnChooseDisc(wxCommandEvent& evt);
	wxWindow* parent;
	wxChoice* cho_disc;
	wxChoice* cho_disc_bounds;
    std::vector<std::string> method_abbreviations;
	std::vector< std::pair<std::string, std::string> > bounds;
	std::vector< std::pair<std::string, std::string> > bounds_per;
	static const int BORDER_PXL = 2;
	static const int SMALL_BORDER=1;
	static const int FLAGS = wxALL | wxALIGN_LEFT | wxALIGN_CENTER_VERTICAL;
	static const int FLAGS_RIGHT = wxALL | wxALIGN_RIGHT | wxALIGN_CENTER_VERTICAL;
};
