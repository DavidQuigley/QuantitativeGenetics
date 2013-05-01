class Plotter : public wxPanel{
public:
    Plotter(wxWindow* parent, int id, int width, int height);
    ~Plotter();
    void reset();
    void add_valueset( std::string name, std::vector<double>* values);
    void set_plot_values(std::vector<std::string> sample_labels, std::vector<std::string> gene_labels, std::vector<std::vector<double>*>& values);
    void render(wxDC& dc);
    void set_key_location(int keyloc);
    void SetDimensions(int height, int weight);
    void DrawLines(bool drawlines);
    void Repaint();
    DECLARE_EVENT_TABLE()
    enum KEY_LOCATIONS{ KEY_UPPER_RIGHT=0, KEY_UPPER_LEFT=1, KEY_LOWER_RIGHT=2, KEY_LOWER_LEFT=3 };

private:
    wxWindow* parent;
    int panel_width, panel_height;
    int scale_v, scale_h;                   // pixels per unit, scaled to window size.
    int margin_x, margin_y_upper, margin_y_lower;
    int y_orig;                             // useful because mac won't reverse Y axis
    std::vector<std::vector<double>*> values; // values to plot
    double min_y, max_y;                    // global highest and lowest points on plot
    int pen_width;
    int key_location;
    bool draw_lines_on_plot;
    std::vector<wxColour*> colors;

    std::vector<std::string> sample_labels;
    std::vector<std::string> value_labels;
    double NA_VALUE;
    const int DEFAULT_MARGIN_Y_LOWER;
    const int DEFAULT_MARGIN_Y_UPPER;
    const int DEFAULT_MARGIN_X;
    void draw_sample_labels(wxDC& dc);
    void draw_axis(wxDC& dc);
    void draw_plot(wxDC& dc);
    void set_scale();
    void set_sample_labels();
    void draw_key(wxDC& dc);
    void paintEvent(wxPaintEvent & evt);
};
