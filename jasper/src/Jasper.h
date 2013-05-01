
#ifndef __WXMSW__
	#include "sample.xpm";
#endif

class JasperGrid : public wxGrid{
public:
	JasperGrid(wxWindow* parent, wxWindowID id, wxSize size, bool is_mac);
	void Update(std::string filename);
	void Clear();
private:
	wxWindow* m_parent;
	static const int BORDER_PXL = 5;
	static const int FLAGS = wxEXPAND | wxALL | wxALIGN_LEFT | wxALIGN_CENTER_VERTICAL;
};


class TreePanel : public wxPanel{
public:
	TreePanel(wxWindow*  parent, Investigation* investigation, wxSize size);
	void OnSelectionChanged(wxTreeEvent& evt);
	JasperTree* tree;
private:
	Investigation* investigation;
	wxWindow* m_parent;
	static const int BORDER_PXL = 5;
	static const int FLAGS = wxEXPAND | wxALL | wxALIGN_LEFT | wxALIGN_CENTER_VERTICAL;
};


class JasperApp : public wxApp
{
public:
    virtual bool OnInit(); 
};

DECLARE_APP(JasperApp)

class JasperFrame : public wxFrame
{
	//DECLARE_CLASS(JasperFrame);
public:
	
    JasperFrame(const wxString& title);
	void react_to_tree_selection(int leaf_type);	
	Investigation* investigation;
	TreePanel* m_treepanel;
	JasperGrid* grid;

private:
	void ShowTextFile(std::string fn);
	void OnNewInvestigation( wxCommandEvent& event );
	void OnOpenInvestigation( wxCommandEvent& event );
	void OnFileRemove( wxCommandEvent& event );
	void OnFileShowProperties( wxCommandEvent& event );
	void OnFileEditProperties( wxCommandEvent& event );
	void OnFileEditInvestigation( wxCommandEvent& event );
	void OnFileShowSA( wxCommandEvent& event );
	void OnFileShowGA( wxCommandEvent& event );
	void OnFileShowRaw( wxCommandEvent& event );
    void OnFileShowRawDisc( wxCommandEvent& event );
	void OnQuit( wxCommandEvent& event );
	void OnAbout( wxCommandEvent& event );
	
    void OnMiningSetNickname( wxCommandEvent& event );
	void OnMiningClearNickname( wxCommandEvent& event );
	void OnMiningRulesetNew( wxCommandEvent& event );
	void OnMiningDelete( wxCommandEvent& event );
	void OnMiningMergeRules( wxCommandEvent& event );
	void OnMiningCorrelationNew( wxCommandEvent& event );
	void OnMiningCoreNew( wxCommandEvent& event );
	void OnMiningDifferenceNew( wxCommandEvent& event );
	
    void OnClassificationClassifierNew(wxCommandEvent& event);
    void OnClassificationApply(wxCommandEvent& event);
    void OnClassificationViewer(wxCommandEvent& event);

    void OnCorrelationOneProbe(wxCommandEvent& event);
    void OnCorrelationGWER(wxCommandEvent& event);
    void OnAnalysisCorrelationViewer(wxCommandEvent& event);
    
    void OnAnnotationAnnotate(wxCommandEvent& event);
    void OnAnnotationGOBrowser(wxCommandEvent& event);
    void OnAnnotationExport(wxCommandEvent& event);
    void OnAnnotationGenesToProbes(wxCommandEvent& event);

	void OnExportNotepad( wxCommandEvent& event );
	void OnExportExcel( wxCommandEvent& event );
	void OnExportCytoscape( wxCommandEvent& event );
	void OnExportFrequency( wxCommandEvent& event );
	void OnExportCrosstab( wxCommandEvent& event );
	void OnExportRulesBySamples( wxCommandEvent& event );
    void OnExportPlot( wxCommandEvent& event);
    void OnExportTmev( wxCommandEvent& event);
    void OnExportCorrelationNetwork( wxCommandEvent& event);
    void OnDocumentation( wxCommandEvent& event );

	void EnableInvestigationMenuItems();
    void fork_shell_command(std::string cmd);
    void initialize_properties(std::string prop_path, std::string homedir, bool is_mac);
    void SetUpMenu();
	void get_file_name_without_extension(std::string& fn);
	bool loaded_properties;
	wxMenuBar* menubar;
	wxSplitterWindow * m_root;
	wxBoxSizer* m_sizer_top;
};
