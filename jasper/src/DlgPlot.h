enum{
	ID_PLOT_TXT_PROBES = 1000,
	ID_PLOT_TXT_GENE,
    ID_PLOT_BTN_ADD,
    ID_PLOT_CHO_GROUP_BY,
    ID_PLOT_CHO_LIMIT_BY,
    ID_PLOT_CHO_IS_OR_IS_NOT,
    ID_PLOT_CHO_LIMIT_TO_VALUE,
    ID_PLOT_CHK_SORT,
    ID_PLOT_CHK_LINES,
    ID_PLOT_BTN_PLOT,
    ID_PLOT_BTN_EXPORT,
    ID_PLOT_BTN_CLOSE,
    ID_PLOT_BMP,
    ID_PLOT_CHO_LABEL_ATTRIBUTE,
    ID_PLOT_KEY_LOCATION,
    ID_PLOT_CHO_WIDTH,
    ID_PLOT_CHO_HEIGHT
};

const int IMG_HEIGHT = 600;
const int IMG_WIDTH = 800;

class PlotDialog: public wxDialog
{
	DECLARE_CLASS( PlotDialog )

public:
	PlotDialog( Investigation* investigation, Dataset* dataset );

private:
    Dataset* dataset;
    wxTextCtrl* txt_probes;
	wxTextCtrl* txt_gene;
    Investigation* investigation;
    wxButton* btn_add_gene;
    wxChoice* cho_group_by;
    wxChoice* cho_limit_to;
    wxArrayString ars_limit_to;
    wxArrayString ars_limit_to_value;
    wxArrayString ars_key_locations;
    wxArrayString ars_height;
    wxArrayString ars_width;
    wxChoice* cho_is_or_is_not;
    wxChoice* cho_limit_to_value;
    wxChoice* cho_key_location;
    wxChoice* cho_label_attribute;
    wxChoice* cho_height;
    wxChoice* cho_width;
    wxCheckBox* chk_sort;
    wxCheckBox* chk_lines;
    wxButton* btn_export;
    wxButton* btn_plot;
    wxButton* btn_close;
    Plotter* plotter;
    std::vector<std::pair<std::string, std::string>*> g_p;
    std::vector<std::string> limits;
    std::vector<std::string> values;
    int plotter_width, total_width, plotter_height, total_height;
    void OnClickAdd( wxCommandEvent& evt );
    void OnClickPlot( wxCommandEvent& evt );
    void OnClickExport( wxCommandEvent& evt );
    void OnClickClose( wxCommandEvent& evt );
    void OnChangeAttribute( wxCommandEvent& evt );
    void OnChar(wxKeyEvent& event);
    void ShowTextFile(std::string fn);
    static const int BORDER_PXL = 2;
    static const int BORDER_NONE = 0;
	static const int SMALL_BORDER=1;
	static const int FLAGS = wxALL | wxALIGN_LEFT | wxALIGN_CENTER_VERTICAL;
	static const int FLAGS_RIGHT = wxALL | wxALIGN_RIGHT | wxALIGN_CENTER_VERTICAL;
	static const int FLAGS_TOP = wxALL | wxALIGN_RIGHT | wxALIGN_TOP;
};
