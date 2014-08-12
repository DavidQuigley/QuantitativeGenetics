#define _SCL_SECURE_NO_DEPRECATE 1
#define _CRT_SECURE_NO_DEPRECATE 1
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <set>
#include <sys/types.h>
#ifndef WIN32
#include <pwd.h>
#endif

using namespace std;
#include "wx/wx.h"
#include <wx/grid.h>
#include <wx/splitter.h>
#include <wx/process.h>
#include <wx/imaglist.h>
#include <wx/treectrl.h>
#include <wx/mimetype.h>
#include <wx/gdicmn.h>
#include <wx/font.h>
#include <wx/mimetype.h>
#include <wx/stdpaths.h>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>
#include <boost/random.hpp>
#include <boost/thread/mutex.hpp>
namespace alg = boost::algorithm;
namespace fs = boost::filesystem;
using namespace boost;
#include "control_ids.h"
#include "carmen/src/DataStructures.h"
#include "carmen/src/Attributes.h"
#include "carmen/src/ClassMinerOptions.h"
#include "carmen/src/Rawdata.h"
#include "carmen/src/Discretize.h"
#include "carmen/src/Dataset.h"
#include "carmen/src/chi2.h"
#include "carmen/src/Rule.h"
#include "carmen/src/graph.h"
#include "carmen/src/Parser.h"
#include "carmen/src/spear.h"
#include "Investigation.h"
#include "JasperTree.h"
#include "PanelInvestigationProperties.h"
#include "DlgInvestigation.h"
#include "DlgInvestigationChooser.h"
#include "DlgNicknameSet.h"
#include "PanelLimits.h"
#include "PanelDiscretization.h"
#include "DlgProgress.h"
#include "DlgRuleset.h"
#include "DlgDifference.h"
#include "DlgCorrelation.h"
#include "DlgMergeRules.h"
#include "DlgClassifier.h"
#include "GridCorrelation.h"
#include "DlgCorrelationViewer.h"
#include "GridClassifierResults.h"
#include "DlgClassificationViewer.h"
#include "DlgClassificationApply.h"
#include "DlgDiscretization.h"
#include "DlgCorrelationOneProbe.h"
#include "DlgCorrelationGWER.h"
#include "DlgCorrelationNetwork.h"
#include "PanelPlotter.h"
#include "DlgPlot.h"
#include "DlgTmev.h"
#include "DlgNewAnnotation.h"
#include "DlgAnnotate.h"
#include "DlgProperties.h"
#include "DlgInvestigationEditor.h"
#include "DlgGenesToProbes.h"
#include "Jasper.h"

IMPLEMENT_APP(JasperApp)

void JasperFrame::fork_shell_command(std::string cmd){
	// This code is only called from the macintosh; the PC side uses wxExecute
	// This is because wxExecute seems only to know how to call .app packages
	// on the OS X side.
#ifdef WIN32
	return;  // This code would never be called from windows
#else
	int	writepipe[2] = {-1,-1},	/* parent -> child */
		readpipe [2] = {-1,-1};	/* child -> parent */
	long childpid;
	writepipe[0] = -1;
	if (pipe(readpipe) < 0 || pipe(writepipe) < 0)
	    return;
	#define	PARENT_READ readpipe[0]
	#define	CHILD_WRITE readpipe[1]
	#define CHILD_READ writepipe[0]
	#define PARENT_WRITE writepipe[1]

	childpid = fork();
	if (childpid < 0)
	    return;
	else if (childpid == 0){
		// in child
	    close(PARENT_WRITE);
	    close(PARENT_READ);
	    dup2(CHILD_READ,  0);
	    close(CHILD_READ);
	    dup2(CHILD_WRITE, 1);
	    close(CHILD_WRITE);
	    execl("/bin/sh", "sh", "-c", cmd.c_str(), NULL);
	    // should not get to this point, as execl only returns on fatal error
	    std::stringstream ss;
	    ss << "Fatal error in command: " << cmd.c_str();
	    wxMessageBox(wxString::FromAscii(ss.str().c_str()));
	    perror("execl");
	}
	else{
		// in parent
		close(CHILD_READ);
		close(CHILD_WRITE);
		//fd = PARENT_READ;
		// normally would return fd, but we don't care about it here.
	}
#endif
}

JasperGrid::JasperGrid(wxWindow* parent, wxWindowID id, wxSize size, bool is_mac) : 
	wxGrid(parent, id, wxPoint(-1,-1), size){
	this->m_parent = parent;
	this->CreateGrid( 0, 0 );
	this->SetRowLabelSize(0);
    this->SetColLabelSize(15);
    this->EnableEditing(false);
	if(is_mac){
	    // DEBUG
		wxFont font(this->GetDefaultCellFont());
		font.SetPointSize(wxSMALL_FONT->GetPointSize());
		this->SetDefaultCellFont(font);
	}
	this->Fit();
}


void JasperGrid::Clear(){
	int n_rows = this->GetNumberRows();
	int n_cols = this->GetNumberCols();
	if( n_rows>0 )
		this->DeleteRows(0, n_rows);
	if( n_cols>0 )
		this->DeleteCols(0, n_cols);
}

void JasperGrid::Update(std::string filename){

	std::string n_rules_str;
	std::string lblA, lblB;
	std::vector<std::string> values;
	int i=0;
	float v;
	char vs[10];
	std::vector<std::string> file_name_parts;
	wxBeginBusyCursor();
	boost::algorithm::split(file_name_parts, filename, boost::algorithm::is_any_of(".") );
	std::string extension = file_name_parts.at( file_name_parts.size()-1 );	
	
	this->Clear();

	CarmenParser cp(filename);
	this->SetColLabelSize(20);

	if( extension.compare("ruleset")==0){
		cp.get_value("Class_A", 0, lblA);  
		cp.get_value("Class_B", 0, lblB);
		if( lblA.size() > 0 && lblB.size() > 0 ){
			this->AppendCols(6);
			this->SetColLabelValue(0, wxString::FromAscii("Rule"));
			this->SetColLabelValue(1, wxString::FromAscii("Chi2"));
			std::stringstream ss, ssb;
			ss << "Support " << lblA ;
			this->SetColLabelValue(2, wxString::FromAscii( ss.str().c_str()));
			this->AutoSizeColumn(2);

			ssb << "Support " << lblB ;
			this->SetColLabelValue(3, wxString::FromAscii( ssb.str().c_str()));
			this->AutoSizeColumn(3);
			this->SetColLabelValue(4, wxString::FromAscii("Confidence"));
			this->SetColLabelValue(5, wxString::FromAscii("Improvement"));
		}
		else{
			this->AppendCols(5);
			this->SetColLabelValue(0, wxString::FromAscii("Rule"));
			this->SetColLabelValue(1, wxString::FromAscii("Chi2"));
			std::stringstream ss;
			ss << "Support " << lblA ;
			this->SetColLabelValue(2, wxString::FromAscii( ss.str().c_str()));
			this->AutoSizeColumn(2);
			this->SetColLabelValue(3, wxString::FromAscii("Confidence"));
			this->SetColLabelValue(4, wxString::FromAscii("Improvement"));
		}
		cp.get_value("Generated", 0, n_rules_str); 
		int n_rules = boost::lexical_cast<int>( n_rules_str );
		if( n_rules>0 ){
			this->AppendRows(n_rules);	
			cp.PrepareToReadValues();
			while( cp.ReadNextValue(values) ){
				this->SetCellValue(i, 0, wxString::FromAscii( values.at(7).c_str() ) );
				this->SetCellValue(i, 1, wxString::FromAscii( values.at(0).c_str() ) );
				this->SetCellValue(i, 2, wxString::FromAscii( values.at(3).c_str() ) );
				this->SetCellValue(i, 3, wxString::FromAscii( values.at(5).c_str() ) );
				this->SetCellValue(i, 4, wxString::FromAscii( values.at(1).c_str() ) );
				v = boost::lexical_cast<float>( values.at(6) );
				if( v < 1 )
					v = v * 100; // XXX earlier versions of CARMEN wrote improvement/100 to file.
				//sprintf_s(vs, 10  "%2.3f", v);
				sprintf(vs, "%2.3f", v);
				this->SetCellValue(i, 5, wxString::FromAscii( vs ) );
				++i;
			}
		}
		this->AutoSizeColumn(0);
		if(this->GetColSize(0)>300)
			this->SetColSize(0, 300);
	}
	else if(extension.compare("classifier")==0){
		this->AppendCols(3);
		this->SetColLabelValue(0, wxString::FromAscii("Rule"));
		cp.get_value("Class_A_train", 0, lblA);
		cp.get_value("Class_B_train", 0, lblB);
		this->SetColLabelValue(1, wxString::FromAscii(lblA.c_str()));
		this->SetColLabelValue(2, wxString::FromAscii(lblB.c_str()));
		cp.PrepareToReadValues();
		float f;
		cp.get_value("n_rules", 0, n_rules_str); 
		int n_rules = boost::lexical_cast<int>( n_rules_str );
		this->AppendRows(n_rules);
		while( cp.ReadNextValue(values) ){
			this->SetCellValue(i, 0, wxString::FromAscii( values.at(0).c_str() ) );
			f = boost::lexical_cast<float>( values.at(3) );
			if(f>0){
				this->SetCellValue(i,1, wxString::FromAscii(values.at(3).c_str()));
				this->SetCellBackgroundColour(i,2,wxColor(255, 255, 255));
				this->SetCellBackgroundColour(i,1,wxColor(153, 204, 255));
                this->SetCellValue(i, 2, wxString::FromAscii("0"));
			}
			else{
				this->SetCellValue(i,2, wxString::FromAscii(values.at(3).c_str()));
				this->SetCellBackgroundColour(i,1,wxColor(255, 255, 255));
				this->SetCellBackgroundColour(i,2,wxColor(153, 204, 255));
                this->SetCellValue(i, 1, wxString::FromAscii("0"));
			}
			++i;
		}
		this->AutoSizeColumn(0);
		if(this->GetColSize(0)>300)
			this->SetColSize(0, 300);
		this->AutoSizeColumn(1);
		this->AutoSizeColumn(2);
	}
    else if(extension.compare("labels")==0){
        /// XXX so far, only bothered to code for single limiting factor
        this->AppendCols(4);
        this->SetColLabelValue(0, wxString::FromAscii("Sample"));
        this->SetColLabelValue(1, wxString::FromAscii("Called from training"));
        this->SetColLabelValue(2, wxString::FromAscii("Label in test"));
        this->SetColLabelValue(3, wxString::FromAscii("Weight"));
        this->AppendRows(cp.CountValueLines()-1); // minus one for header
        cp.PrepareToReadValues();
        std::string class_a;
        std::vector<std::string> limit_parts;
        cp.get_value("Class_A_test", 0, class_a);
        alg::split(limit_parts, class_a, boost::algorithm::is_any_of("=") );
        std::vector<std::string> header;
        cp.ReadNextValue(header);
        int idx_limit;
        for(idx_limit=0; idx_limit<(int)header.size();idx_limit++){
            if( header.at(idx_limit).compare(limit_parts.at(0))==0)
                break;
        }
        while(cp.ReadNextValue(values)){
            this->SetCellValue(i,0, wxString::FromAscii(values.at(0).c_str()));
            this->SetCellValue(i,1, wxString::FromAscii(values.at(1).c_str()));
            if( idx_limit < (int)header.size() )
                this->SetCellValue(i,2, wxString::FromAscii(values.at(idx_limit).c_str()));
            else
                this->SetCellValue(i,2, wxString::FromAscii("?"));
            this->SetCellValue(i,3, wxString::FromAscii(values.at(2).c_str()));
            ++i;
        }
    }
	else if(extension.compare("difference")==0){
		cp.get_value("Class_A", 0, lblA);
		cp.get_value("Class_B", 0, lblB);
        std::string n_perms;
        cp.get_value("N_Permutations", 0, n_perms);
        bool did_perms = n_perms.compare(std::string("1")) != 0;
        if(did_perms)
            this->AppendCols(6);
        else {
            this->AppendCols(5);
        }
		this->SetColLabelValue(0, wxString::FromAscii("symbol") );
		this->SetColLabelValue(1, wxString::FromAscii("probe") );
		this->SetColLabelValue(2, wxString::FromAscii(lblA.c_str()));
		this->SetColLabelValue(3, wxString::FromAscii(lblB.c_str()));
		this->SetColLabelValue(4, wxString::FromAscii("t stat"));
		if(did_perms){
            this->SetColLabelValue(5, wxString::FromAscii("p-value"));
        }
		int n_rules = cp.CountValueLines();
		cp.PrepareToReadValues();
		this->AppendRows(n_rules);
		while( cp.ReadNextValue(values) ){
			this->SetCellValue(i, 0, wxString::FromAscii( values.at(1).c_str() ) );
            this->SetCellValue(i, 1, wxString::FromAscii( values.at(0).c_str() ) );
			this->SetCellValue(i, 2, wxString::FromAscii( values.at(2).c_str() ) );
			this->SetCellValue(i, 3, wxString::FromAscii( values.at(4).c_str() ) );
			this->SetCellValue(i, 4, wxString::FromAscii( values.at(6).c_str() ) );
			if(did_perms){
                this->SetCellValue(i, 5, wxString::FromAscii( values.at(7).c_str() ) );
            }
			++i;
		}
	}
	else if(extension.compare("spear")==0){
		cp.get_value("Class_A", 0, lblA);  
		cp.get_value("Class_B", 0, lblB);
		bool show_full = false;
		if( lblA.size() > 0 && lblB.size() > 0 ){
			this->AppendCols(8);
			this->SetColLabelValue(0, wxString::FromAscii("symbol 1"));
			this->SetColLabelValue(1, wxString::FromAscii("probe 1"));
			this->SetColLabelValue(2, wxString::FromAscii("symbol 2"));
			this->SetColLabelValue(3, wxString::FromAscii("probe 2"));
			std::stringstream ss;
			ss << "r (" << lblA << ")";
			this->SetColLabelValue(4, wxString::FromAscii(ss.str().c_str()));
			std::stringstream ss2;
			ss2 << "r (" << lblB << ")";
			this->SetColLabelValue(5, wxString::FromAscii(ss2.str().c_str()));
			this->SetColLabelValue(6, wxString::FromAscii("r (change)"));
            this->SetColLabelValue(7, wxString::FromAscii("Z score"));
			show_full = true;
		}
		else{
			this->AppendCols(6);
			std::stringstream ss;
			this->SetColLabelValue(0, wxString::FromAscii("symbol 1"));
			this->SetColLabelValue(1, wxString::FromAscii("probe 1"));
			this->SetColLabelValue(2, wxString::FromAscii("symbol 2"));
			this->SetColLabelValue(3, wxString::FromAscii("probe 2"));
			ss << "r (" << lblA << ")";
			this->SetColLabelValue(4, wxString::FromAscii(ss.str().c_str()));
            this->SetColLabelValue(5, wxString::FromAscii("Z score"));
		}
		int i=0;
		int n_spear = cp.CountValueLines();
		int truncate = wxNO;
		if( n_spear > 100000 ){
			truncate = wxMessageBox(wxString::FromAscii("File has more than 100,000 entries.  Truncate to 100,000?"), wxString::FromAscii("WARNING: File is very large."), wxYES_NO|wxICON_QUESTION, this);
		}
		cp.PrepareToReadValues();
		this->AppendRows(n_spear);
		while( cp.ReadNextValue(values) ){
			this->SetCellValue(i, 0, wxString::FromAscii( values.at(0).c_str() ) );
            this->SetCellValue(i, 1, wxString::FromAscii( values.at(1).c_str() ) );
			this->SetCellValue(i, 2, wxString::FromAscii( values.at(2).c_str() ) );
			this->SetCellValue(i, 3, wxString::FromAscii( values.at(3).c_str() ) );
			this->SetCellValue(i, 4, wxString::FromAscii( values.at(4).c_str() ) );
			if( show_full ){
				this->SetCellValue(i, 5, wxString::FromAscii( values.at(5).c_str() ) );
				this->SetCellValue(i, 6, wxString::FromAscii( values.at(6).c_str() ) );
                this->SetCellValue(i, 7, wxString::FromAscii( values.at(8).c_str() ) ); // values.at(7) is perm_p
			}
            else{
                this->SetCellValue(i, 5, wxString::FromAscii( values.at(8).c_str() ) ); // values.at(7) is perm_p
            }
			++i;
			if( truncate==wxYES && i>100000 )
				break;
		}
	}
	wxEndBusyCursor();
	this->Layout();
	this->m_parent->FitInside();
}



TreePanel::TreePanel(wxWindow* parent, Investigation* investigation, wxSize size)
	: wxPanel(parent, ID_Ctrl_TreePanel, wxPoint(-1, -1), size, wxBORDER_DOUBLE){
	this->m_parent = parent;
	this->investigation = investigation;
	this->tree = new JasperTree(this, ID_Ctrl_Tree, investigation, size);
	Connect(ID_Ctrl_Tree, wxEVT_COMMAND_TREE_SEL_CHANGED, wxTreeEventHandler(TreePanel::OnSelectionChanged));
	wxFlexGridSizer* sizer_top = new wxFlexGridSizer(1,1,0,0);
	sizer_top->AddGrowableRow(0,1);
	sizer_top->AddGrowableCol(0,1);
	sizer_top->Add(this->tree, 1, wxGROW & FLAGS, BORDER_PXL);
	this->SetSizerAndFit(sizer_top);
}


void TreePanel::OnSelectionChanged(wxTreeEvent & evt){
	TreeData* d = (TreeData*) this->tree->GetItemData(evt.GetItem());
	JasperFrame* p = (JasperFrame*)this->m_parent;	
	if( d->leaf_type==LEAF_FILE && d->str.size()>0){
		p->grid->Update(d->str);
	}
	((JasperFrame*)this->m_parent)->react_to_tree_selection(d->leaf_type);
}


bool JasperApp::OnInit(){
	if ( !wxApp::OnInit() ) // parse a few common command-line options
        return false;
	JasperFrame *frame = new JasperFrame(wxT("CARMEN: Exploring systems genetics"));
	frame->Show(true);
	SetTopWindow( frame );
	return true;
}


void JasperFrame::initialize_properties(std::string prop_path, std::string homedir, bool is_mac){
	// create a new properties file from scratch
	ConfigParser cp;
	cp.create(prop_path);
	cp.set(std::string("History"), std::string("most_recent_fileset"), std::string("") );
	cp.add_section("Internal");
	cp.set(std::string("Internal"), std::string("homedir"), homedir );
	cp.add_section("External");
	if( !is_mac )
		cp.set(std::string("External"), std::string("text_editor"), "notepad.exe");
	cp.add_section(std::string("Nicknames"));
	cp.add_section(std::string("Filesets"));
	cp.add_section(std::string("Annotations"));
	std::stringstream ss_go, ss_mouse, ss_human;
	if( is_mac ){
		ss_go << homedir << "/Carmen.app/Contents/MacOS/gene_ontology.obo.compressed";
		ss_mouse << homedir << "/Carmen.app/Contents/MacOS/mouse_annotation.txt";
		ss_human << homedir << "/Carmen.app/Contents/MacOS/human_annotation.txt";
	}
	else{
		ss_go << homedir << "\\gene_ontology.obo.compressed";
		ss_mouse << homedir << "\\mouse_annotation.txt";
		ss_human << homedir << "\\human_annotation.txt";
	}
	cp.set("Annotations", "go", ss_go.str() );
	cp.set("Annotations", "mouse", ss_mouse.str() );
	cp.set("Annotations", "human", ss_human.str() );
	cp.write();
}


JasperFrame::JasperFrame(const wxString& title): wxFrame(NULL, wxID_ANY, title){

	SetIcon(wxICON(sample));
	std::string carmen_base_dir( boost::filesystem::current_path().string() );
	std::string prop_path, prop_path_old_location;
#ifdef WIN32
	bool is_mac=false;
	prop_path = carmen_base_dir + "\\.properties";
#else
	bool is_mac=true;
	struct passwd *pwd;
	pwd = getpwuid(getuid());
	std::string user_home_dir = pwd->pw_dir;
	prop_path = user_home_dir + "/.carmen_properties";
	prop_path_old_location = carmen_base_dir + "/Carmen.app/Contents/MacOS/.properties";
	bool properties_in_app = boost::filesystem::exists( boost::filesystem::path(prop_path_old_location) );
	bool properties_in_home = boost::filesystem::exists( boost::filesystem::path(prop_path) );
	// On the Mac, old location for .properties was in App bundle
	// Going forward, this file lives in the user's home folder
	if( properties_in_app && !properties_in_home ){
		std::cout << "Writing .carmen_properties to " << user_home_dir << "\n";
		ConfigParser cp(prop_path_old_location);
		cp.write(prop_path);
	}
#endif
	// NOTE: the call to investigaion->set_current is in JasperTree.cpp,
	// JasperTree::refresh_from_directory(). This is eventually called during
	// initialization and any other time investigation details may have changed.
	try{
		this->investigation = new Investigation(prop_path);
		this->loaded_properties = true;
	}
	catch( std::string msg ){
		std::stringstream ss;
		std::string fn_carmen;
		this->loaded_properties = false;
		ss << "Error attempting to load property file at\n" << prop_path << "\nGenerate new property file?";
		int yn = wxMessageBox(wxString::FromAscii(ss.str().c_str()), wxString::FromAscii("Error attempting to load property file"), wxYES|wxNO, this);
		if( yn==wxYES){
			try{
				this->initialize_properties(prop_path, carmen_base_dir, is_mac);
				this->investigation = new Investigation(prop_path);
				this->loaded_properties = true;
			}
			catch( std::string err){
				wxMessageBox(wxString::FromAscii(err.c_str()), wxString::FromAscii("Error writing .properties file"), wxOK , this);
				this->Close(true);
				return;
			}
		}
		else{
			wxMessageBox(wxString::FromAscii("This file is required; exiting CARMEN."));
			this->Close(true);
			return;
		}
	}
	if( this->investigation->cp->get(std::string("Internal"), std::string("debug")).compare(std::string("true"))==0 )
		std::cout << "Properties loaded from: " << prop_path << "\n";
	// Re-set homedir to our current homedir if we've moved since the last time we were run.
	// REMOVED THIS CODE. If you run from the command line in another directory, this blows up our ability to
	// find software that the Mac verison expects to find inside of the app bundle.
	//if( this->investigation->cp->get(std::string("Internal"), std::string("homedir") ).compare(carmen_base_dir) != 0 ){
	//	this->investigation->cp->set(std::string("Internal"), std::string("homedir"), carmen_base_dir );
	//	this->investigation->cp->write();
	//}
	this->SetUpMenu();
    CreateStatusBar(2);
	SetStatusText(_T("CARMEN 1.2"));
	wxFlexGridSizer *sizer_top = new wxFlexGridSizer(2,1,0,0);
	sizer_top->AddGrowableRow(1,1);
	sizer_top->AddGrowableCol(0,1);
	this->m_treepanel = new TreePanel(this, this->investigation, wxSize(600,300));
	this->grid = new JasperGrid(this, ID_Ctrl_Grid, wxSize(600,300), this->investigation->is_mac); 
	sizer_top->Add(this->m_treepanel, 1, wxGROW);
	sizer_top->Add(this->grid, 1, wxGROW);
	this->SetSizerAndFit(sizer_top);
	this->react_to_tree_selection(-1); // force a redraw

}


void JasperFrame::SetUpMenu(){
	this->menubar = new wxMenuBar();
    wxMenu* fileMenu = new wxMenu();
    menubar->Append(fileMenu, _T("&Investigation"));
	fileMenu->Append(Menu_File_New, _T("&New"), _T("Create a new Investigation"));
	fileMenu->Append(Menu_File_Open, _T("&Open"), _T("Open an Investigation"));
	fileMenu->Append(Menu_File_Edit_Investigation, _T("&Edit Current Investigation"), _T("Edit Current Investigation"));
	fileMenu->AppendSeparator();
#ifdef WIN32
	fileMenu->Append(wxID_PREFERENCES, _T("Edit CARMEN &Preferences"), _T("Edit CARMEN Preferences"));
#endif
	fileMenu->Append(Menu_File_Show_Properties, _T("Show CARMEN Proper&ties file"), _T("Show CARMEN properties file"));
	fileMenu->Append(Menu_File_Show_GA, _T("Show &Gene Attributes"), _T("Show gene attributes file"));
	fileMenu->Append(Menu_File_Show_SA, _T("&Show Sample Attributes"), _T("Show Sample Attributes file"));
	fileMenu->Append(Menu_File_Show_raw, _T("Show &Raw Data"), _T("Show Raw Data"));
    fileMenu->AppendSeparator();
    fileMenu->Append( Menu_Mining_Delete,  _T("&Delete Result File"), _T("Delete result file"));
    fileMenu->AppendSeparator();
	fileMenu->Append(Menu_File_Remove, _T("&Remove"), _T("Remove this Investigation"));
	fileMenu->AppendSeparator();
	fileMenu->Append(wxID_EXIT, _T("E&xit\tAlt-X"), _T("Quit this program"));
	
	wxMenu* miningMenu = new wxMenu;
	menubar->Append(miningMenu, _T("&Mining"));
	miningMenu->Append(Menu_Mining_Ruleset_New, _T("&New Ruleset"), _T("Create a new ruleset"));
	miningMenu->Append(Menu_Mining_Core_New, _T("New C&ore Ruleset"), _T("Create a new ruleset"));
	miningMenu->Append(Menu_Mining_Merge_Rules, _T("Merge Redundant Rules"), _T("Merge redundant rules"));
	miningMenu->AppendSeparator();
	miningMenu->Append(Menu_Mining_Set_Nickname, _T("&Set Nickname\tF2"), _T("Set Nickname for query"));
	miningMenu->Append(Menu_Mining_Clear_Nickname, _T("&Clear Nickname"), _T("Clear Nickname for query"));

	//wxMenu* classificationMenu = new wxMenu;
	//menubar->Append(classificationMenu, _T("&Classification"));
	//classificationMenu->Append(Menu_Classification_Classifier_New, _T("&New Classifier"), _T("Create a new classifier"));
    //classificationMenu->Append(Menu_Classification_Apply, _T("&Apply"), _T("Apply classifier to another dataset"));
    //classificationMenu->Append(Menu_Classification_Viewer, _T("&View Classifier Results"), _T("View classifier results"));

	wxMenu* correlationMenu = new wxMenu;
	menubar->Append(correlationMenu, _T("&Statistics"));
    correlationMenu->Append(Menu_Correlation_New, _T("&Correlation"), _T("Create a new correlation measurement"));
    correlationMenu->Append(Menu_Correlation_New_One_Probe, _T("Correlation with One &Probe"), _T("Create a new correlation measurement, one probe vs. all"));
	correlationMenu->Append(Menu_Correlation_Viewer,  _T("Correlation Vie&wer"), _T("Explore Correlation"));
	correlationMenu->Append(Menu_Correlation_GWER,  _T("Calculate &GWER"), _T("Calculate GWER"));
	correlationMenu->Append(Menu_Export_Correlation_Network, _T("Correlation &Network"), _T("Export Correlation network to Cytoscape"));
    correlationMenu->AppendSeparator();
    correlationMenu->Append(Menu_Mining_Difference_New, wxString::FromAscii("Di&fferential Expression"), wxString::FromAscii("Test for differential expression"));

    wxMenu* annotationMenu  = new wxMenu;
    menubar->Append(annotationMenu, _T("&Annotation"));
    annotationMenu->Append(Menu_Annotation_Annotate, _T("&Annotate genes"), _T("Annotate gene list"));
    annotationMenu->Append(Menu_Annotation_GO_Browser, _T("&Find genes by GO"), _T("Browse Gene Ontology"));
    annotationMenu->Append(Menu_Annotation_Export_GO, _T("Show GO &file"), _T("Show Gene Ontology file"));
    annotationMenu->AppendSeparator();
    annotationMenu->Append(Menu_Annotation_Genes_To_Probes, _T("&Convert genes to probes"), _T("Convert genes to probes"));

	wxMenu* exportMenu = new wxMenu;
	menubar->Append(exportMenu, _T("&Export"));
	exportMenu->Append(Menu_Export_Notepad, _T("To &Text Editor"), _T("Export file to Text Editor"));
	exportMenu->Append(Menu_Export_Excel, _T("To &Spreadsheet"), _T("Export file to Spreadsheet"));
    exportMenu->Append(Menu_Export_Tmev, _T("To &Datasheet"), _T("Export genes, samples, and data in one file"));
    exportMenu->Append(Menu_File_Show_raw_disc, _T("Raw &Data Discretized"), _T("Export Raw Data Discretized"));
    exportMenu->AppendSeparator();
	exportMenu->Append(Menu_Export_Cytoscape, _T("To &Cytoscape"), _T("Export file to Cytoscape"));
	exportMenu->Append(Menu_Export_Crosstab, _T("C&rosstab of rules"), _T("Export crosstab of rules by samples"));
	exportMenu->Append(Menu_Export_Frequency, _T("Gene &Frequency List"), _T("Export ranked list of gene frequencies"));
	exportMenu->Append(Menu_Export_Rules_by_Samples, _T("R&ules by Samples"), _T("Export matrix of {rules x samples}"));
    exportMenu->AppendSeparator();
    exportMenu->Append(Menu_Export_Plot, _T("&Plot"), _T("Plot gene expression and export values"));

	wxMenu* helpMenu = new wxMenu;
	menubar->Append(helpMenu, _T("&Help"));
	helpMenu->Append(Menu_Help_Documentation, _T("&Documentation"), _T("Show Documentation"));
	helpMenu->Append(wxID_ABOUT, _T("&About\tF1"), _T("Show about dialog"));
	Connect(wxID_ABOUT, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(JasperFrame::OnAbout));
	Connect(wxID_EXIT, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(JasperFrame::OnQuit));
	Connect(wxID_NEW, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(JasperFrame::OnNewInvestigation));
	Connect(wxID_OPEN, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(JasperFrame::OnOpenInvestigation));
    Connect(Menu_File_Edit_Investigation, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(JasperFrame::OnFileEditInvestigation));
	Connect(wxID_PREFERENCES, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(JasperFrame::OnFileEditProperties));
	Connect(Menu_File_Show_Properties, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(JasperFrame::OnFileShowProperties));
	Connect(Menu_File_Remove, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(JasperFrame::OnFileRemove));
	Connect(Menu_File_Show_GA, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(JasperFrame::OnFileShowGA));
	Connect(Menu_File_Show_SA, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(JasperFrame::OnFileShowSA));
	Connect(Menu_File_Show_raw, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(JasperFrame::OnFileShowRaw));

	Connect(Menu_Mining_Delete, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(JasperFrame::OnMiningDelete));

    Connect(Menu_Mining_Merge_Rules, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(JasperFrame::OnMiningMergeRules));
	Connect(Menu_Mining_Set_Nickname, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(JasperFrame::OnMiningSetNickname));
	Connect(Menu_Mining_Clear_Nickname, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(JasperFrame::OnMiningClearNickname));
	Connect(Menu_Mining_Ruleset_New, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(JasperFrame::OnMiningRulesetNew));
	Connect(Menu_Mining_Core_New, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(JasperFrame::OnMiningCoreNew));
	Connect(Menu_Mining_Difference_New, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(JasperFrame::OnMiningDifferenceNew));

    Connect(Menu_Correlation_New, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(JasperFrame::OnMiningCorrelationNew));
	Connect(Menu_Correlation_New_One_Probe, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(JasperFrame::OnCorrelationOneProbe));
    Connect(Menu_Correlation_Viewer, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(JasperFrame::OnAnalysisCorrelationViewer));
    Connect(Menu_Correlation_GWER, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(JasperFrame::OnCorrelationGWER));
    //Connect(Menu_Classification_Classifier_New, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(JasperFrame::OnClassificationClassifierNew));
    //Connect(Menu_Classification_Apply, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(JasperFrame::OnClassificationApply));
    //Connect(Menu_Classification_Viewer, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(JasperFrame::OnClassificationViewer));

    Connect(Menu_Annotation_Annotate, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(JasperFrame::OnAnnotationAnnotate));
    Connect(Menu_Annotation_GO_Browser, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(JasperFrame::OnAnnotationGOBrowser));
    Connect(Menu_Annotation_Export_GO, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(JasperFrame::OnAnnotationExport));
    Connect(Menu_Annotation_Genes_To_Probes, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(JasperFrame::OnAnnotationGenesToProbes));

    Connect(Menu_Export_Notepad, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(JasperFrame::OnExportNotepad));
	Connect(Menu_Export_Excel, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(JasperFrame::OnExportExcel));
	Connect(Menu_File_Show_raw_disc, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(JasperFrame::OnFileShowRawDisc));
	Connect(Menu_Export_Tmev, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(JasperFrame::OnExportTmev));
    Connect(Menu_Export_Cytoscape, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(JasperFrame::OnExportCytoscape));
	Connect(Menu_Export_Crosstab, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(JasperFrame::OnExportCrosstab));
	Connect(Menu_Export_Frequency, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(JasperFrame::OnExportFrequency));
	Connect(Menu_Export_Rules_by_Samples, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(JasperFrame::OnExportRulesBySamples));
    Connect(Menu_Export_Plot, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(JasperFrame::OnExportPlot));
    Connect(Menu_Export_Correlation_Network, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(JasperFrame::OnExportCorrelationNetwork));
    Connect(Menu_Help_Documentation, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(JasperFrame::OnDocumentation));
	SetMenuBar(menubar);
}

void JasperFrame::OnQuit(wxCommandEvent& WXUNUSED(event)){
    Close(true);  
}

void JasperFrame::OnAbout(wxCommandEvent& WXUNUSED(event)){
	wxMessageBox(wxString::FromAscii("CARMEN 1.3\nDavid Quigley\nBalmain Lab, UCSF\n2014"), _T("About Carmen"), wxOK | wxICON_INFORMATION, this);
}


void JasperFrame::OnFileRemove(wxCommandEvent& WXUNUSED(event)){
	if ( wxMessageBox(wxString::FromAscii("Remove this investigation? (Underlying files will not be touched)"), wxString::FromAscii("Remove?"), wxYES_NO | wxICON_QUESTION, this) ==wxYES ){
		std::string i_name(this->investigation->investigation_name );
		alg::to_lower(i_name);
		this->investigation->cp->remove( "Filesets", i_name+ "<>dir" );
		this->investigation->cp->remove( "Filesets", i_name + "<>raw" );
		this->investigation->cp->remove( "Filesets", i_name + "<>ga" );
		this->investigation->cp->remove( "Filesets", i_name + "<>sa" );
		std::string names = this->investigation->cp->get("Filesets", "names");
		alg::to_lower(names);
		std::vector<std::string> arr_names;
		alg::split(arr_names, names, alg::is_any_of(","));
		std::stringstream ss;
		std::string set_to="";
		for(int i=0; i<(int)arr_names.size(); i++){
			if( arr_names.at(i)!=i_name ){
				if(arr_names.at(i).size()>0){
					set_to = arr_names.at(i);
					ss << arr_names.at(i);
					if(i<(int)arr_names.size()-1)
						ss << ",";
				}
			}
		}
		this->investigation->cp->set ( "Filesets", "names", ss.str() );
		try{
			if( set_to.size()>0 ){ // one for removed, one for some other investigation
				this->investigation->cp->set("History", "most_recent_fileset", set_to);
				this->investigation->set_current( set_to );
				
			}
			else{
				this->investigation->cp->set("History", "most_recent_fileset", "");
				this->investigation->investigation_name="";
			}
		}
		catch(std::string err){
			wxMessageBox(wxString::FromAscii(err.c_str() ), wxString::FromAscii("ERROR"), wxOK | wxICON_EXCLAMATION, this);
		}
		this->investigation->cp->write();
		this->m_treepanel->tree->refresh_from_directory();
		react_to_tree_selection(LEAF_ROOT);
	}
}

void JasperFrame::OnOpenInvestigation(wxCommandEvent & WXUNUSED(event)){
	InvestigationChooserDialog dialog( NULL, this->investigation->cp );
	std::string current_investigation( this->investigation->investigation_name ); // may be empty

	if( dialog.ShowModal() == wxID_OK ){
		if( dialog.investigation.size() > 0 ){
			alg::to_lower(dialog.investigation);
			this->investigation->cp->set("History", "most_recent_fileset", dialog.investigation);
			this->investigation->cp->write();
			this->grid->Clear();
			((JasperTree*)this->m_treepanel->tree)->refresh_from_directory();
			if( this->investigation->investigation_name.size()==0 ){
				// load of new investigation failed for some reason
				// attempt to roll back to previously loaded one.
				this->investigation->cp->set("History", "most_recent_fileset", current_investigation);
				this->investigation->investigation_name = current_investigation;
				((JasperTree*)this->m_treepanel->tree)->refresh_from_directory();
			}
		}
	}
}


void JasperFrame::OnNewInvestigation(wxCommandEvent & WXUNUSED(event)){
	InvestigationDialog dialog(this->investigation);
	if (dialog.ShowModal() == wxID_OK){
		std::string fs("Filesets");
		std::string i_name = dialog.investigation_name;
		alg::to_lower(i_name);
		// check to see if we are overwriting an existing fileset
		std::string names = this->investigation->cp->get(fs, "names");
		alg::to_lower(names);
		std::vector<std::string> arr_names;
		alg::split(arr_names, names, alg::is_any_of(","));
		bool exists=false;
		for(int i=0; i<(int)arr_names.size(); i++)
			if( arr_names.at(i)==i_name )
				exists = true;
		if(exists){
			int result = wxMessageBox(wxString::FromAscii("There is already an investigation with that name.  Overwrite?"), wxString::FromAscii("Overwrite?"), wxYES_NO | wxICON_QUESTION, this);
			if( result!=wxYES )
				return;
		}
		std::string sa = i_name + "<>sa";
		std::string ga = i_name + "<>ga";
		std::string raw = i_name + "<>raw";
		std::string dir_results = i_name + "<>dir";
		std::string species = i_name + "<>species";
		std::string data_type = i_name + "<>data_type";
		std::string gene_name_column = i_name + "<>gene_name_column";
		this->investigation->cp->set(fs, sa, dialog.fn_sa );
		this->investigation->cp->set(fs, ga, dialog.fn_ga );
		this->investigation->cp->set(fs, raw, dialog.fn_expr );
		this->investigation->cp->set(fs, dir_results, dialog.dir_results);
		this->investigation->cp->set(fs, species, dialog.species);
		this->investigation->cp->set(fs, data_type, dialog.data_type);
		this->investigation->cp->set(fs, gene_name_column, dialog.gene_name_column );
		this->investigation->cp->set(std::string("History"), std::string("most_recent_fileset"), i_name );
		if( !exists )
			this->investigation->cp->set(fs, std::string("names"), names + "," + i_name);
		this->investigation->cp->write();
		this->investigation->set_current(i_name);
		this->m_treepanel->tree->refresh_from_directory();

		this->EnableInvestigationMenuItems(); // for special case when this is the first investigation
	}
}


void JasperFrame::OnMiningCoreNew( wxCommandEvent& evt ){
	std::string filepath;
	std::string file_noext;
	this->get_file_name_without_extension(file_noext);
	this->m_treepanel->tree->get_selected_filename(filepath);
	std::stringstream ss;
	if(this->investigation->is_mac)
		ss << this->investigation->dir_results << '/' << file_noext << ".core.ruleset";
	else
		ss << this->investigation->dir_results << '\\' << file_noext << ".core.ruleset";
	std::string fn_core(ss.str()); 

	std::vector<std::string> cmd;
	cmd.push_back( "coreminer" );
    cmd.push_back("-r" + filepath);
	cmd.push_back("-t0"); // haven't implemented UI for tolorance
    cmd.push_back("-o" + fn_core );
	cmd.push_back("-vT");

    ProgressDialog dialog( this, cmd, this->investigation, true );
	if( dialog.ShowModal() == wxID_OK){   
		this->m_treepanel->tree->refresh_from_directory();
		this->m_treepanel->tree->open_tree_to_file( fn_core );
	}
}


void JasperFrame::OnMiningSetNickname(wxCommandEvent& WXUNUSED(event)){
	if( this->m_treepanel->tree->get_selected_leaftype() == LEAF_QUERY){
		NicknameSetDialog dialog(NULL, wxID_ANY, wxT("Set Nickname for query"));
		if (dialog.ShowModal() == wxID_OK){
			std::string orig_label;
			std::string colon(":");
			std::string colon_word("<colon>");
			this->m_treepanel->tree->get_selected_original_label(orig_label);
			std::string name_lower = this->investigation->investigation_name;
			boost::algorithm::to_lower(name_lower);
			std::string key = "\"" + name_lower + "<>" + orig_label + "\"";
			std::string value = "\"" + dialog.new_nickname + "\"";
			boost::algorithm::replace_all(key, colon, colon_word);
			boost::algorithm::replace_all(value, colon, colon_word);
			boost::algorithm::to_lower(key);
			boost::algorithm::to_lower(value);
			this->investigation->cp->set("Nicknames", key, value);
			this->investigation->cp->write();
			this->investigation->read_from_file(this->investigation->fn_properties);
			this->investigation->set_current(this->investigation->investigation_name);
			this->m_treepanel->tree->refresh_from_directory();
		}
	}
}

void JasperFrame::OnMiningClearNickname(wxCommandEvent& WXUNUSED(event)){
	if( wxMessageBox(_T("CONFIRM: clear nickname?"), _T("Clear?"), wxYES_NO|wxICON_QUESTION)==wxYES){
		std::string orig_label;
		std::string colon(":");
		std::string colon_word("<colon>");

		this->m_treepanel->tree->get_selected_original_label(orig_label);
		std::string name_lower = this->investigation->investigation_name;
		boost::algorithm::to_lower(name_lower);
		std::string key = "\"" + name_lower + "<>" + orig_label + "\"";
		boost::algorithm::replace_all(key, colon, colon_word);
		boost::algorithm::to_lower(key);

		this->investigation->cp->remove("Nicknames", key);
		this->investigation->cp->write();
		this->m_treepanel->tree->refresh_from_directory();
	}
}


void JasperFrame::OnMiningRulesetNew(wxCommandEvent& WXUNUSED(event)){
	std::string filename;
	std::string ruleset(".ruleset");
	this->m_treepanel->tree->get_selected_filename(filename);
	if( !alg::find_first( filename, ruleset) )
		filename = "";

	RulesetDialog dialog(this->investigation, filename);
	if (dialog.ShowModal() == wxID_OK){
		this->m_treepanel->tree->refresh_from_directory();
        this->m_treepanel->tree->open_tree_to_file( dialog.new_filename );
	}
}


void JasperFrame::OnMiningMergeRules(wxCommandEvent& WXUNUSED(event)){
	std::string fn;
	this->m_treepanel->tree->get_selected_filename(fn);
	MergeRulesDialog dialog( this->investigation, fn );
	if( !dialog.cancel_load ){
		if( dialog.ShowModal() == wxID_OK){
            if( dialog.saved_new_file ){
				this->m_treepanel->tree->refresh_from_directory();
                this->m_treepanel->tree->open_tree_to_file( dialog.new_filename );
            }
		}
	}
}


void JasperFrame::OnMiningCorrelationNew(wxCommandEvent& WXUNUSED(event)){
	CorrelationDialog dialog(this->investigation);
	if (dialog.ShowModal() == wxID_OK){
		this->m_treepanel->tree->refresh_from_directory();
        this->m_treepanel->tree->open_tree_to_file( dialog.new_filename );
	}
}


void JasperFrame::OnCorrelationOneProbe(wxCommandEvent& WXUNUSED(event)){
	CorrelationOneProbeDialog dialog(this->investigation);
	if (dialog.ShowModal() == wxID_OK){
		this->m_treepanel->tree->refresh_from_directory();
        this->m_treepanel->tree->open_tree_to_file( dialog.new_filename );
	}
}

void JasperFrame::OnCorrelationGWER(wxCommandEvent& WXUNUSED(event)){
	CorrelationGWERDialog dialog(this->investigation);
	if (dialog.ShowModal() == wxID_OK){
		this->m_treepanel->tree->refresh_from_directory();
        this->m_treepanel->tree->open_tree_to_file( dialog.new_filename );
	}
}

void JasperFrame::OnMiningDifferenceNew(wxCommandEvent& WXUNUSED(event)){
	std::string filename;
	this->m_treepanel->tree->get_selected_filename(filename);
	DifferenceDialog dialog(this->investigation, filename);
	if (dialog.ShowModal() == wxID_OK){
		this->m_treepanel->tree->refresh_from_directory();
        this->m_treepanel->tree->open_tree_to_file( dialog.new_filename );
	}
}

void JasperFrame::OnMiningDelete(wxCommandEvent& WXUNUSED(event)){
	if( wxMessageBox(_T("CONFIRM: Delete result file?"), _T("Delete?"), wxYES_NO|wxICON_QUESTION)==wxYES){
		std::string fs;
		this->m_treepanel->tree->get_selected_filename( fs );
		boost::filesystem::path rs(fs);
		boost::filesystem::remove( rs );
		this->m_treepanel->tree->refresh_from_directory();
	}
}


void JasperFrame::OnClassificationClassifierNew(wxCommandEvent& WXUNUSED(event)){
	std::string filename;
	this->m_treepanel->tree->get_selected_filename(filename);
	ClassifierDialog dialog(this->investigation, filename);
	if (dialog.ShowModal() == wxID_OK){
		this->m_treepanel->tree->refresh_from_directory();
        this->m_treepanel->tree->open_tree_to_file( dialog.new_filename );
	}
}

void JasperFrame::OnClassificationApply(wxCommandEvent& WXUNUSED(event)){
	std::string filename;
	this->m_treepanel->tree->get_selected_filename(filename);
	ClassificationApplyDialog dialog(filename, this->investigation);
	if (dialog.ShowModal() == wxID_OK){
		this->m_treepanel->tree->refresh_from_directory();
        this->m_treepanel->tree->open_tree_to_file( dialog.new_filename );
	}
}


void JasperFrame::OnClassificationViewer(wxCommandEvent& event){
    std::string filename;
	this->m_treepanel->tree->get_selected_filename(filename);
	ClassificationViewerDialog dialog(filename, this->investigation);
	dialog.ShowModal();
}

void JasperFrame::ShowTextFile(std::string fn){
	std::stringstream ss;
	if(this->investigation->is_mac){
		ss << "open -t \"" << fn << "\""; // -t flag opens with default text editor
	}
	else{
		std::string loc_text = this->investigation->cp->get(std::string("External"), std::string("text_editor"));;
		ss << loc_text << " \"" << fn << "\"";
	}
	wxExecute(wxString::FromAscii(ss.str().c_str()));
}

void JasperFrame::OnDocumentation(wxCommandEvent& WXUNUSED(event)){
	std::stringstream ss;
	std::string homedir = this->investigation->cp->get(std::string("Internal"), std::string("homedir"));
	std::string helpfile = this->investigation->cp->get(std::string("Internal"), std::string("helpfile"));
	if( helpfile.size()==0 )
		helpfile = std::string("CARMEN_documentation.pdf");
	if( this->investigation->is_mac ){
		ss << "open \"" << homedir << "/Carmen.app/Contents/MacOS/" << helpfile << "\"";
        std::cout << ss.str();
	}
	else{
		ss << "\"" << homedir << "\\" << helpfile << "\"";
	}
	wxExecute( wxString::FromAscii(ss.str().c_str()) );
}

void JasperFrame::OnFileEditInvestigation(wxCommandEvent& WXUNUSED(event)){
	InvestigationEditorDialog dialog(this->investigation);
	if( dialog.ShowModal() == wxID_OK){
		this->investigation->cp->set( "Filesets", dialog.investigation_name+ "<>dir", dialog.dir_results);
		this->investigation->cp->set( "Filesets", dialog.investigation_name + "<>raw", dialog.fn_expr  );
		this->investigation->cp->set( "Filesets", dialog.investigation_name + "<>ga", dialog.fn_ga );
		this->investigation->cp->set( "Filesets", dialog.investigation_name + "<>sa", dialog.fn_sa);
		this->investigation->cp->set( "Filesets", dialog.investigation_name + "<>species", dialog.species);
		this->investigation->cp->set( "Filesets", dialog.investigation_name + "<>data_type", dialog.data_type);
		this->investigation->cp->set( "Filesets", dialog.investigation_name + "<>gene_name_column", dialog.gene_name_column);
		this->investigation->fn_expr = dialog.fn_expr;
		this->investigation->fn_ga = dialog.fn_ga;
		this->investigation->fn_sa = dialog.fn_sa;
		this->investigation->dir_results = dialog.dir_results;
		this->investigation->species = dialog.species;
		this->investigation->gene_name_column = dialog.gene_name_column;
		this->investigation->data_type = dialog.data_type;
		delete this->investigation->ga;
		delete this->investigation->sa;
		this->investigation->ga = new Attributes("NA");
		this->investigation->sa = new Attributes("NA");
		this->investigation->ga->load(dialog.fn_ga);
		this->investigation->sa->load(dialog.fn_sa);
		((JasperTree*)this->m_treepanel->tree)->refresh_from_directory();
	}
}

void JasperFrame::OnFileEditProperties(wxCommandEvent& WXUNUSED(event)){
	PropertyDialog dialog(this->investigation);
	if( dialog.ShowModal() == wxOK ){
		if( dialog.new_annotation_files.size()>0 ){
			for(int i=0; i<int(dialog.new_annotation_files.size() ); i++ ){
				this->investigation->cp->set(std::string("Annotations"), dialog.new_annotation_names.at(i), dialog.new_annotation_files.at(i) );
			}
		}
		this->investigation->cp->write();
	}
}

void JasperFrame::OnFileShowProperties(wxCommandEvent& WXUNUSED(event)){
	this->ShowTextFile( this->investigation->fn_properties );
}

void JasperFrame::OnFileShowSA(wxCommandEvent& WXUNUSED(event)){
	this->ShowTextFile( this->investigation->fn_sa );
}

void JasperFrame::OnFileShowGA(wxCommandEvent& WXUNUSED(event)){
	this->ShowTextFile(this->investigation->fn_ga);
}

void JasperFrame::OnFileShowRaw(wxCommandEvent& WXUNUSED(event)){
	this->ShowTextFile( this->investigation->fn_expr );
}

void JasperFrame::OnFileShowRawDisc(wxCommandEvent& WXUNUSED(event)){
	
	DiscretizationDialog dialog(this, this->investigation);
	if( dialog.ShowModal() != wxID_OK)
		return;
	std::string limit_a;
	
	if( dialog.limits.size() >0 )
		limit_a = std::string(dialog.limits);
	else
		limit_a = std::string("IDENTIFIER!NULL");
	ClassMinerOptions * cmo = new ClassMinerOptions();
	cmo->file_name_dis = this->investigation->fn_expr;
	cmo->file_name_sa = this->investigation->fn_sa;
	cmo->file_name_ga = this->investigation->fn_ga;
	cmo->class_a = limit_a;
	try{
		Dataset* data = new Dataset();
		Attributes* sa = new Attributes("NA");
		Attributes* ga = new Attributes("NA");
		sa->load(cmo->file_name_sa);
		ga->load(cmo->file_name_ga);
		float disc_lower, disc_upper;
		try{
			disc_upper = boost::lexical_cast<float>( dialog.disc_upper );
			disc_lower = boost::lexical_cast<float>( dialog.disc_lower );
		}
		catch( boost::bad_lexical_cast blc){
			disc_upper = 0;
			disc_lower = 0;
		}
		cmo->discretization = dialog.disc_method;
		cmo->disc_lower = disc_lower;
		cmo->disc_upper = disc_upper;
		data->load(sa, ga, cmo);
		data->write_discretized(dialog.fn_discretization, std::string("\t"), dialog.disc_method, disc_lower, disc_upper);
		this->ShowTextFile( dialog.fn_discretization );
		delete data;
		delete ga;
		delete sa;
	}
	catch(std::string err){
		wxMessageBox(wxString::FromAscii(err.c_str()));
	}
	delete cmo;
}


void JasperFrame::OnAnalysisCorrelationViewer(wxCommandEvent& event){
	std::string fn_spear;
	this->m_treepanel->tree->get_selected_filename(fn_spear);
	std::string fn_compressed( fn_spear );
	fn_compressed += ".compressed";
	std::ifstream f( fn_compressed.c_str() );
	wxBeginBusyCursor();
	CarmenParser cp(fn_spear);
	cp.write_compressed_spear(fn_compressed, this->investigation->ga->identifier2idx );
	wxEndBusyCursor();
	
	CorrelationViewerDialog dialog(fn_compressed, this->investigation);
	dialog.ShowModal();
}


void JasperFrame::OnAnnotationAnnotate(wxCommandEvent& WXUNUSED(event)){
	std::string filename;
	this->m_treepanel->tree->get_selected_filename(filename);
	AnnotationDialog dialog(this->investigation, filename);
	if( dialog.ShowModal() == wxID_OK)
		this->ShowTextFile( dialog.new_filename );
}

void JasperFrame::OnAnnotationGOBrowser(wxCommandEvent& WXUNUSED(event)){
	GOFilterDialog dialog(this->investigation);
	dialog.ShowModal();
}

void JasperFrame::OnAnnotationExport(wxCommandEvent& WXUNUSED(event)){
	this->ShowTextFile( this->investigation->cp->get(std::string("Annotations"), std::string("go")) );
}

void JasperFrame::OnAnnotationGenesToProbes(wxCommandEvent& WXUNUSED(event)){
	GenesToProbesDialog dialog(this->investigation);
	dialog.ShowModal();
}

void JasperFrame::OnExportNotepad(wxCommandEvent& WXUNUSED(event)){
	std::string fn;
	this->m_treepanel->tree->get_selected_filename(fn);
	this->ShowTextFile( fn );
}


void JasperFrame::get_file_name_without_extension(std::string& fn){
	
	std::stringstream file_name;
	std::string filepath;
	this->m_treepanel->tree->get_selected_filename(filepath);
	std::vector<std::string> path_parts, file_name_parts;
	if( filepath.find_first_of("/") == std::string::npos )
		alg::split( path_parts, filepath, alg::is_any_of("\\") );
	else
		alg::split( path_parts, filepath, alg::is_any_of("/") );
	alg::split( file_name_parts, path_parts.at(path_parts.size()-1), alg::is_any_of(".") );
	for(int i=0; i<(int)file_name_parts.size()-1; i++){
		file_name << file_name_parts.at(i);
		if(i<(int)file_name_parts.size()-2)
			file_name << ".";
	}
	fn = file_name.str();
}


void JasperFrame::OnExportCytoscape(wxCommandEvent& WXUNUSED(event)){

    // is cytoscape present?
	std::string forward("/");
	std::string back("\\");
	std::string slash;
    std::string loc_cytoscape = this->investigation->cp->get(std::string("External"), std::string("cytoscape"));
    if( loc_cytoscape.size()==0 ){
		wxMessageBox(wxString::FromAscii("Cytoscape location not set; please select Cytoscape application."), wxString::FromAscii("Cytoscape location not set"), wxOK|wxICON_EXCLAMATION );
		wxString cloc = wxFileSelector(wxString::FromAscii("Select file"), wxString::FromAscii(""), wxString::FromAscii(""), wxString::FromAscii(""), wxString::FromAscii(""), wxFD_OPEN, this);
		if( cloc.Len() > 0 ){
			loc_cytoscape = std::string(cloc.ToAscii());
			this->investigation->cp->set(std::string("External"), std::string("cytoscape"), loc_cytoscape);
			this->investigation->cp->write();
		}
		else{
			return;
		}
	}
    if( this->investigation->is_mac )
    	alg::replace_all(loc_cytoscape, std::string("Cytoscape.app"), std::string("Cytoscape.sh") );
    if( this->investigation->cp->get(std::string("Internal"), std::string("debug")).compare(std::string("true"))==0 )
    	std::cout << loc_cytoscape << "\n";
	if(this->investigation->is_mac){
		slash = forward;
		alg::replace_all(loc_cytoscape, back, forward);
	}
	else{
		slash = back;
		alg::replace_all(loc_cytoscape, forward, back);
	}
	if( loc_cytoscape.at(0) == '"' && loc_cytoscape.at( loc_cytoscape.size()-1 )=='"')
        loc_cytoscape = loc_cytoscape.substr(1, loc_cytoscape.size()-2 );

    // where will we write the files?
    std::string fn;
	get_file_name_without_extension(fn);
	wxString default_location = wxString::FromAscii(this->investigation->dir_results.c_str());
	std::stringstream ss;
	ss << fn + "_cyto";
	wxString default_name = wxString::FromAscii(ss.str().c_str());
	wxString default_extension = wxString::FromAscii("");
	wxString filepath = wxFileSelector(wxString::FromAscii("Base filename for Cytoscape"), default_location, default_name, default_extension, wxString::FromAscii("*.*"), wxFD_SAVE, this);
	std::string location(filepath.ToAscii() );
	if(location.size()>0){
        std::string filepath;
    	this->m_treepanel->tree->get_selected_filename(filepath);	
        std::vector<std::string> file_name_parts;
	    boost::algorithm::split(file_name_parts, filepath, boost::algorithm::is_any_of(".") );
	    std::string extension = file_name_parts.at( file_name_parts.size()-1 );
		std::string home_dir = this->investigation->cp->get(std::string("Internal"), std::string("homedir"));
		std::string cyto_viz = this->investigation->cp->get(std::string("Internal"), std::string("cytoscape_vizmap"));
        std::vector<std::string> cmd;
		if( extension.compare("spear")==0 || extension.compare("SPEAR")==0 ){
    		Spearman sp(filepath,  this->investigation->ga);
    		sp.set_fn_cytoscape(loc_cytoscape);
    		sp.set_fn_cytoscape_props(cyto_viz);
    		std::vector<std::string> cmd;

    		if( this->investigation->is_mac){
    			sp.set_batch_extension(std::string("sh"));
    			cmd.push_back(location + ".sh");
    		}
    		else{
    			sp.set_batch_extension(std::string("bat"));
    			cmd.push_back(location + ".bat");
    		}
    		sp.write_to_cytoscape(0, location, this->investigation->ga);
    		std::stringstream batch;
			if(this->investigation->is_mac){
				batch << "sh " << location << ".sh";
				fork_shell_command( batch.str() );
			}
			else{
    			batch << location << ".bat";
				wxExecute(wxString::FromAscii(batch.str().c_str()));
			}
        }
        else{
            // from rules file
        	try{

				RuleSet ruleset(filepath);
				ruleset.write_network(location);
				std::stringstream batch;
				batch << loc_cytoscape;
				batch << " -N \"" << location << ".sif\"";
				batch << " -n \"" << location << ".noa\"";
				batch << " -n \"" << location << "_names.noa\"";
				batch << " -V \"" << cyto_viz << "\"";
				batch << " -e \"" << location << "_class.eda\"";
				batch << " -e \"" << location << "_conf.eda\"";
				batch << " -e \"" << location << "_sup.eda\"";
				batch << " -e \"" << location << "_pval.eda\"";
				if(this->investigation->is_mac)
					fork_shell_command( batch.str() );
				else
					wxExecute(wxString::FromAscii(batch.str().c_str()));
			}
        	catch(std::string err){
        		wxMessageBox(wxString::FromAscii(err.c_str()));
        	}
        	catch(std::string* err){
        		wxMessageBox(wxString::FromAscii(err->c_str()));
        	}
        }
	}
}


bool comp_gene_counts(const std::pair<std::string, int> *a, const std::pair<std::string, int> *b){
	return a->second > b->second;
}

void JasperFrame::OnExportFrequency(wxCommandEvent& WXUNUSED(event)){
	std::string fn;
	get_file_name_without_extension(fn);
	std::string filepath;
	this->m_treepanel->tree->get_selected_filename(filepath);
	std::stringstream ss;
	ss << fn << ".frequency.txt";
	wxString default_location = wxString::FromAscii(this->investigation->dir_results.c_str());
	wxString default_name = wxString::FromAscii(ss.str().c_str());
	wxString default_extension = wxString::FromAscii("txt");
	wxString fileloc = wxFileSelector(wxString::FromAscii("Save To"), default_location, default_name, default_extension, wxString::FromAscii("*.*"), wxFD_SAVE, this);
	if( fileloc.Len()>0 ){
		CarmenParser cp(filepath);
		cp.PrepareToReadValues();
		HASH_S_I hsh_genes;
		cp.ExtractGeneCounts(hsh_genes);
		std::vector< std::pair<std::string, int>* >  counts;
		for( HASH_S_I::iterator itr=hsh_genes.begin(); itr != hsh_genes.end(); itr++ ){
			std::pair<std::string, int>* p = new std::pair<std::string, int>();
			p->first = (*itr).first;
			p->second = (*itr).second;
			counts.push_back(p);
		}
		std::sort(counts.begin(), counts.end(), comp_gene_counts);
		std::string fn_out(fileloc.ToAscii() );
		std::ofstream f_out(fn_out.c_str());
		if( !f_out.is_open() ){
			std::stringstream ss;
			ss << "Unable to open for writing: " << fn_out;
			wxMessageBox(wxString::FromAscii(ss.str().c_str()), wxString::FromAscii("Unable to open"), wxOK|wxICON_EXCLAMATION, this);
			return;
		}
		f_out << "# Gene frequency list from source file " << filepath << "\n";
		for(int i=0; i<(int)counts.size(); i++){
			f_out << counts.at(i)->first << "\t" << counts.at(i)->second << "\n";
			delete counts.at(i);
		}
		f_out.close();
		this->ShowTextFile( std::string( fileloc.ToAscii()) );
	}
}


void JasperFrame::OnExportTmev(wxCommandEvent& WXUNUSED(event)){
    std::string filename;
	this->m_treepanel->tree->get_selected_filename(filename);
	TmevDialog dialog(this->investigation,  filename);
	if( dialog.ShowModal() == wxID_OK)
		this->ShowTextFile( dialog.new_filename );
}


void JasperFrame::OnExportCorrelationNetwork(wxCommandEvent& WXUNUSED(event)){
    // is cytoscape present?
	std::string forward("/");
	std::string back("\\");
    std::string loc_cytoscape = this->investigation->cp->get(std::string("External"), std::string("cytoscape"));
	// if not present, let user tell us where it is
    if( loc_cytoscape.size()==0 ){
		wxMessageBox(wxString::FromAscii("Cytoscape location not set; please select Cytoscape application."), wxString::FromAscii("Cytoscape location not set"), wxOK|wxICON_EXCLAMATION );
		wxString cloc = wxFileSelector(wxString::FromAscii("Select file"), wxString::FromAscii(""), wxString::FromAscii(""), wxString::FromAscii(""), wxString::FromAscii(""), wxFD_OPEN, this);
		if( cloc.Len() > 0 ){
			loc_cytoscape = std::string(cloc.ToAscii());
			this->investigation->cp->set(std::string("External"), std::string("cytoscape"), loc_cytoscape);
			this->investigation->cp->write();
		}
		else{
			return;
		}
	}

    // Used to use cytoscape.sh; this is wrong. If found, replace with .app
    if( this->investigation->is_mac ){
    	alg::replace_all(loc_cytoscape, std::string("Cytoscape.app"), std::string("Cytoscape.sh") );
    	std::cout << "replacing with " << loc_cytoscape << "\n";
    	this->investigation->cp->set(std::string("External"), std::string("cytoscape"), loc_cytoscape);
    	this->investigation->cp->write();
    }
    std::string loc_properties = this->investigation->cp->get(std::string("Internal"), std::string("cytoscape_vizmap"));
    if( loc_properties.size()==0 ){
    	std::stringstream ss;
    	if( this->investigation->is_mac )
    		ss << this->investigation->cp->get(std::string("Internal"), std::string("homedir")) << "/Carmen.app/Contents/MacOS/ClassAnalysis.props";
    	else
    		ss << this->investigation->cp->get(std::string("Internal"), std::string("homedir")) << "\\ClassAnalysis.props";
    	this->investigation->cp->set(std::string("Internal"), std::string("cytoscape_vizmap"), ss.str());
    	this->investigation->cp->write();
    }

	if( loc_cytoscape.at(0) == '"' && loc_cytoscape.at( loc_cytoscape.size()-1 )=='"')
		loc_cytoscape = loc_cytoscape.substr(1, loc_cytoscape.size()-2 );
	if(!this->investigation->is_mac)
		alg::replace_all(loc_cytoscape, forward, back);

	bool exists_as_file = wxFileExists( wxString::FromAscii( loc_cytoscape.c_str() ) );
	bool exists_as_dir = wxDirExists( wxString::FromAscii( loc_cytoscape.c_str() ) );
	if( loc_cytoscape.size()==0 || (!exists_as_file && !exists_as_dir )  ){
		std::stringstream ss;
		ss << "ERROR: Cytoscape not found at location: " << loc_cytoscape;
		wxMessageBox(wxString::FromAscii(ss.str().c_str()), wxString::FromAscii("Cannot find Cytoscape"), wxOK|wxICON_EXCLAMATION, this);
		return;
	}

	std::string filename;
	this->m_treepanel->tree->get_selected_filename(filename);
	CorrelationNetworkDialog dialog(this->investigation,  filename);
	if( dialog.ShowModal() == wxID_OK ){
		boost::filesystem::path sh_file(dialog.batch_filename.c_str());
		if( boost::filesystem::exists( sh_file ) ){
			if(this->investigation->is_mac){
				std::stringstream ss;
				ss << "sh " << dialog.batch_filename;
				fork_shell_command( ss.str() );
			}
			else{
				wxExecute(wxString::FromAscii(dialog.batch_filename.c_str()));
			}
		}
		else{
			std::stringstream ss;
			ss << "ERROR: File not created: " << dialog.batch_filename;
			wxMessageBox(wxString::FromAscii(ss.str().c_str()), wxString::FromAscii("File not created correctly"), wxOK|wxICON_EXCLAMATION, this);
		}
	}
}


void JasperFrame::OnExportPlot(wxCommandEvent& WXUNUSED(event)){
	Dataset* dataset = new Dataset();
	try{ dataset->load(this->investigation->sa, this->investigation->ga, this->investigation->fn_expr, std::string("IDENTIFIER!NULL"), std::string("")); }
	catch(std::string e){
		wxMessageBox(wxString::FromAscii(e.c_str()));
		return;
	}
	PlotDialog dialog( this->investigation, dataset );
    dialog.ShowModal();
    delete dataset;
}


void JasperFrame::OnExportCrosstab(wxCommandEvent& WXUNUSED(event)){
	std::string fn;
	get_file_name_without_extension(fn);
	std::string filepath, slash;
	this->m_treepanel->tree->get_selected_filename(filepath);
	if(this->investigation->is_mac)
		slash = std::string("/");
	else
		slash = std::string("\\");
	std::stringstream ss;
	ss << fn << ".crosstab.txt";
	wxString default_location = wxString::FromAscii(this->investigation->dir_results.c_str());
	wxString default_name = wxString::FromAscii(ss.str().c_str());
	wxString default_extension = wxString::FromAscii("txt");
	wxString fileloc = wxFileSelector(wxString::FromAscii("Save To"), default_location, default_name, default_extension, wxString::FromAscii("*.*"), wxFD_SAVE, this);
	if( fileloc.size() > 0 ){
		try{
			RuleSet ruleset(filepath);
			ruleset.write_crosstab(std::string(fileloc.ToAscii()));
			this->ShowTextFile(std::string(fileloc.ToAscii()) );
		}
		catch(std::string err){
			wxMessageBox(wxString::FromAscii(err.c_str()));
		}
	}
}


void JasperFrame::OnExportExcel(wxCommandEvent& WXUNUSED(event)){
	std::string fn;
	std::string forward("/");
	std::string back("\\");
	this->m_treepanel->tree->get_selected_filename(fn);
	std::string loc_excel = this->investigation->cp->get(std::string("External"), std::string("excel"));
	if( this->investigation->is_mac )
		alg::replace_all(loc_excel, back, forward );
	else
		alg::replace_all(loc_excel, forward, back);
	if( loc_excel.size()==0){
		wxMessageBox(wxString::FromAscii("Excel location not set; please select Excel application."), wxString::FromAscii("Excel location not set"), wxOK|wxICON_EXCLAMATION );
		wxString cloc = wxFileSelector(wxString::FromAscii("Select file"), wxString::FromAscii(""), wxString::FromAscii(""), wxString::FromAscii(""), wxString::FromAscii(""), wxFD_OPEN, this);
		if( cloc.Len() > 0 ){
			loc_excel = std::string(cloc.ToAscii());
			this->investigation->cp->set(std::string("External"), std::string("excel"), loc_excel);
			this->investigation->cp->write();
		}
		else{
			return;
		}
	}
	std::stringstream ss;
	std::string homedir = this->investigation->cp->get(std::string("Internal"), std::string("homedir"));
	if(this->investigation->is_mac){
		if(loc_excel.at(0)=='\"' )
			ss << "open -a " << loc_excel << " \"" << fn << "\"";
		else
			ss << "open -a \"" << loc_excel << "\" \"" << fn << "\"";
		this->fork_shell_command(ss.str());
	}
	else{
		ss << loc_excel << " \"" << fn << "\"";
		std::string cmd = ss.str();
		alg::replace_all(cmd, forward, back);
		wxExecute(wxString::FromAscii(cmd.c_str()));
	}
	
}

void JasperFrame::OnExportRulesBySamples(wxCommandEvent& WXUNUSED(event)){
	std::string fn;
	get_file_name_without_extension(fn);
	std::string filepath;
	this->m_treepanel->tree->get_selected_filename(filepath);
	std::stringstream ss;
	ss << fn << ".txt" ;
	wxString default_location = wxString::FromAscii(this->investigation->dir_results.c_str());
	wxString default_name = wxString::FromAscii(ss.str().c_str());
	wxString default_extension = wxString::FromAscii("txt");
	wxString fileloc = wxFileSelector(wxString::FromAscii("Save To"), default_location, default_name, default_extension, wxString::FromAscii("*.*"), wxFD_SAVE, this);
	if( fileloc.Len()>0 ){
		try{
			std::string fn(fileloc.ToAscii());
			RuleSet* ruleset = new RuleSet(filepath);
			Matrix<double>* r_by_s = new Matrix<double>();
			ClassifierDataset* data = new ClassifierDataset();
			Attributes* sa = new Attributes("NA");
			Attributes* ga = new Attributes("NA");
			sa->load(ruleset->file_name_sa);
			ga->load(ruleset->file_name_ga);
			ClassMinerOptions * cmo = new ClassMinerOptions();
			cmo->file_name_dis = ruleset->file_name_dis;
			cmo->file_name_sa = ruleset->file_name_sa;
			cmo->file_name_ga = ruleset->file_name_ga;
			cmo->discretization = ruleset->disc;
			cmo->disc_lower = (float)ruleset->disc_lower;
			cmo->disc_upper = (float)ruleset->disc_upper;
			cmo->gene_limit = ruleset->gene_limit;
			cmo->class_a = ruleset->class_a;
			cmo->class_b = ruleset->class_b;
			data->load(sa, ga, cmo);
			ruleset->build_rules_by_samples(r_by_s, data);
			ruleset->write_rules_by_samples(fn, r_by_s, sa);
			delete cmo;
			delete data;
			delete sa;
			delete ga;
			this->ShowTextFile(fn);
		}
		catch(std::string err){
			wxMessageBox(wxString::FromAscii(err.c_str()));
		}
	}
}


void JasperFrame::EnableInvestigationMenuItems(){
	// Enable various menu items if there is an investigation open
	this->menubar->Enable(Menu_File_Edit_Investigation, true);
	this->menubar->Enable(Menu_File_Remove, true);
	this->menubar->Enable(Menu_File_Open, true);
	this->menubar->Enable(Menu_File_Show_GA, true);
	this->menubar->Enable(Menu_File_Show_SA, true);
	this->menubar->Enable(Menu_File_Show_raw, true);
    this->menubar->Enable(Menu_File_Show_raw_disc, true);
    this->menubar->Enable(Menu_Mining_Ruleset_New, true);
	this->menubar->Enable(Menu_Correlation_New, true);
    this->menubar->Enable(Menu_Correlation_New_One_Probe, true);
	this->menubar->Enable(Menu_Export_Correlation_Network, true);
	this->menubar->Enable(Menu_Mining_Difference_New, true);
    this->menubar->Enable(Menu_Export_Tmev, true);
    this->menubar->Enable(Menu_Export_Plot, true);
    this->menubar->Enable(Menu_Correlation_GWER, true);
    this->menubar->Enable(Menu_Annotation_Genes_To_Probes, true);
}

void JasperFrame::react_to_tree_selection(int leaf_type){
	// Enable/disable various menu items depending on what type of file is selected

	this->menubar->Enable(Menu_File_Show_Properties, this->loaded_properties);
	this->menubar->Enable(Menu_File_Remove, false);
	this->menubar->Enable(Menu_File_Edit_Investigation, false);
	this->menubar->Enable(Menu_Mining_Ruleset_New, false);
	this->menubar->Enable(Menu_File_Show_GA, false);
	this->menubar->Enable(Menu_File_Show_SA, false);
	this->menubar->Enable(Menu_File_Show_raw, false);
    this->menubar->Enable(Menu_File_Show_raw_disc, false);
    this->menubar->Enable(Menu_Mining_Delete, false);
	this->menubar->Enable(Menu_Mining_Core_New, false);
	this->menubar->Enable(Menu_Mining_Difference_New, false);

	this->menubar->Enable(Menu_Mining_Merge_Rules, false);
	this->menubar->Enable(Menu_Mining_Set_Nickname, false);
	this->menubar->Enable(Menu_Mining_Clear_Nickname, false);

    this->menubar->Enable(Menu_Correlation_New, false);
    this->menubar->Enable(Menu_Correlation_New_One_Probe, false);
    this->menubar->Enable(Menu_Correlation_Viewer, false);
    this->menubar->Enable(Menu_Correlation_GWER, false);

	this->menubar->Enable(Menu_Export_Notepad, false);
	this->menubar->Enable(Menu_Export_Excel, false);
    this->menubar->Enable(Menu_Export_Tmev, false);
	this->menubar->Enable(Menu_Export_Cytoscape, false);
	this->menubar->Enable(Menu_Export_Frequency, false);
	this->menubar->Enable(Menu_Export_Crosstab, false);
	this->menubar->Enable(Menu_Export_Rules_by_Samples, false);
	this->menubar->Enable(Menu_Export_Cytoscape, false);
    this->menubar->Enable(Menu_Export_Plot, false);
    this->menubar->Enable(Menu_Export_Correlation_Network, false);

    this->menubar->Enable(Menu_Annotation_Annotate, false);
    this->menubar->Enable(Menu_Annotation_GO_Browser, false);
    this->menubar->Enable(Menu_Annotation_Export_GO, false);
    this->menubar->Enable(Menu_Annotation_Genes_To_Probes, false);

    if( this->investigation->cp->get(std::string("Annotations"), std::string("go")).size() > 0 ){
    	// has Gene Ontology definition file
    	this->menubar->Enable(Menu_Annotation_GO_Browser, true);
    	this->menubar->Enable(Menu_Annotation_Export_GO, true);
    	typedef std::pair<std::string, std::string> prop_pair;
    	std::vector< prop_pair* >* annotations = new std::vector<prop_pair*>();
    	this->investigation->cp->get_section(std::string("Annotations"), annotations);
    	for(int i=0; i<(int)annotations->size(); i++){
    		if( annotations->at(i)->first.compare("go") != 0 ){
    			this->menubar->Enable(Menu_Annotation_Annotate, true); // at least one species annotation
    		}
    		delete annotations->at(i);
    	}
    }

	if( this->investigation->investigation_name.size()>0){
		this->EnableInvestigationMenuItems();
    }

	if( this->investigation->cp->get("Filesets", "names").size() == 0 ){
		this->menubar->Enable(Menu_File_Open, false);
	}

	if( leaf_type == LEAF_QUERY ){
		this->menubar->Enable(Menu_Mining_Set_Nickname, true);
		this->menubar->Enable(Menu_Mining_Clear_Nickname, true);
	}
	else if( leaf_type == LEAF_FILE ){
		this->menubar->Enable(Menu_Mining_Delete, true);
		this->menubar->Enable(Menu_Export_Notepad, true);
		this->menubar->Enable(Menu_Export_Excel, true);
		std::string filename;
		this->m_treepanel->tree->get_selected_filename(filename);
		if( alg::find_first( filename, ".ruleset" ) ){
			this->menubar->Enable(Menu_Mining_Core_New, true);
			this->menubar->Enable(Menu_Export_Cytoscape, true);
			this->menubar->Enable(Menu_Export_Crosstab, true);
			this->menubar->Enable(Menu_Export_Frequency, true);
			this->menubar->Enable(Menu_Mining_Merge_Rules, true);
			this->menubar->Enable(Menu_Export_Rules_by_Samples, true);
		}
		else if( alg::find_first(filename, ".spear" )){
            this->menubar->Enable(Menu_Correlation_Viewer, true);
            this->menubar->Enable(Menu_Export_Frequency, true);
            this->menubar->Enable(Menu_Export_Cytoscape, true);
		}
        //else if( alg::find_first(filename, ".classifier" ) || alg::find_first(filename, ".labels" )){
        //    this->menubar->Enable(Menu_Classification_Viewer, true);
        //    this->menubar->Enable(Menu_Classification_Apply, true);
        //}
	}
	else{
		this->grid->Clear();
	}
}
