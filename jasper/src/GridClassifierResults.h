class GridClassifierResults : public wxGrid{
public:
	GridClassifierResults(wxWindow* parent, std::string filename, wxWindowID id, wxSize size);

private:
	wxWindow* m_parent;
	std::string filename;
    void set_row_color(int row, int r, int g, int b);
	static const int BORDER_PXL = 5;
	static const int FLAGS = wxEXPAND | wxALL | wxALIGN_LEFT | wxALIGN_CENTER_VERTICAL;
};
