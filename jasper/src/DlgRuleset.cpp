#define _SCL_SECURE_NO_DEPRECATE 1
#define _CRT_SECURE_NO_DEPRECATE 1
#include <time.h>
#include <wx/dialog.h>
#include <wx/sizer.h>
#include <wx/stattext.h>
#include <wx/choice.h>
#include <wx/button.h>
#include <wx/textctrl.h>
#include <wx/arrstr.h>
#include <wx/panel.h>
#include <wx/listbox.h>
#include <wx/wx.h>
#include <wx/statbox.h>
#include <wx/process.h>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/random.hpp>
namespace alg = boost::algorithm;
#include <fstream>
#include <string>
#include <vector>
#include "carmen/src/DataStructures.h"
#include "carmen/src/Attributes.h"
#include "carmen/src/Rawdata.h"
#include "carmen/src/Parser.h"
#include "Investigation.h"
using namespace boost;
using namespace std;
#include "PanelLimits.h"
#include "PanelDiscretization.h"
#include "DlgProgress.h"
#include "DlgRuleset.h"

IMPLEMENT_CLASS( RulesetDialog, wxDialog )

RulesetDialog::RulesetDialog( Investigation* investigation, std::string filename ){
	this->investigation = investigation;
	SetExtraStyle(wxWS_EX_BLOCK_EVENTS|wxDIALOG_EX_CONTEXTHELP);
	wxDialog::Create( NULL, wxID_ANY, wxT("Create new Ruleset"), wxDefaultPosition, wxDefaultSize, wxCAPTION|wxRESIZE_BORDER|wxSYSTEM_MENU );
	CreateControls();
	GetSizer()->Fit(this); // This fits the dialog to the minimum size dictated by the sizers
	GetSizer()->SetSizeHints(this); // This ensures that the dialog cannot be sized smaller than the minimum size
	Centre();
	if(filename.size()>0)
		set_controls_from_file( filename );
}


void RulesetDialog::set_controls_from_file(std::string filename ){
	CarmenParser cp(filename);
	std::string class_a, class_b, sci_str, max_chi, max_depth, disc, disc_mult, which_class;
	cp.get_value(std::string("Class_A"), 0, class_a);
	cp.get_value(std::string("Class_B"), 0, class_b);
	cp.get_value(std::string("S_C_I"), 0, sci_str);
	std::vector<std::string> sci;
	alg::split(sci, sci_str, alg::is_any_of(",") );
	cp.get_value(std::string("Maximum_Chi2"), 0, max_chi);
	cp.get_value(std::string("Max_Depth"), 0, max_depth);
	cp.get_value(std::string("Discretization"), 0, disc);
	cp.get_value(std::string("Which_Class"), 0, which_class);
	if( which_class.size()==0 )
		which_class="both";
	else if( which_class.compare("a")==0 || which_class.compare("A")==0 )
		this->cho_show_results->SetSelection(1);
	else if( which_class.compare("b")==0 || which_class.compare("B")==0 )
		this->cho_show_results->SetSelection(2);
	this->txt_support->SetValue( wxString::FromAscii(sci.at(0).c_str() ) );
	this->txt_confidence->SetValue( wxString::FromAscii(sci.at(1).c_str() ) );
	this->txt_improvement->SetValue( wxString::FromAscii(sci.at(2).c_str() ) );
	this->txt_max_depth->SetValue( wxString::FromAscii(max_depth.c_str()) );
	this->txt_max_chi->SetValue( wxString::FromAscii(max_chi.c_str())  );
	this->limit_A->set_limits( class_a );
	this->limit_B->set_limits( class_b );
	
	if( disc.compare("none")==0 ){
		std::string none("none");
		std::string empty_string("");
		this->discretization->set_discretization(none, empty_string, empty_string);
	}
	else{
		cp.get_value(std::string("Discretization"), 3, disc_mult);
		std::vector<std::string> mu;
		alg::split(mu, disc_mult, alg::is_any_of(",") );
		std::string lower = mu.at(0);
		std::string upper = mu.at(1);
		this->discretization->set_discretization(disc, lower, upper);
	}
}

void RulesetDialog::CreateControls()
{
	char timebuf[80];
	struct tm* newtime;
	time_t long_time;
	time( &long_time );
	newtime = localtime( &long_time );
	strftime(timebuf, 80, "%Y_%m_%d_%H_%M_%S", newtime);
	wxStaticText* lbl_file_name = new wxStaticText( this, wxID_ANY, wxString::FromAscii("File Name:") );
	this->txt_file_name = new wxTextCtrl( this, ID_TXT_FILE_NAME, wxString::FromAscii(timebuf), wxDefaultPosition, wxSize(200, 22) );
	this->txt_file_name->ChangeValue(wxString::FromAscii(timebuf));
	wxStaticText* lbl_support = new wxStaticText( this, wxID_ANY, wxString::FromAscii("Support:"));
	this->txt_support = new wxTextCtrl( this, ID_TXT_SUPPORT, wxString::FromAscii("100"), wxDefaultPosition, wxSize(70, 22) );
	this->txt_support->SetValue( wxT("100") );
	wxStaticText* lbl_confidence = new wxStaticText( this, wxID_ANY, wxString::FromAscii("Confidence:"));
	this->txt_confidence = new wxTextCtrl( this, ID_TXT_CONFIDENCE, wxString::FromAscii("100"), wxDefaultPosition, wxSize(70, 22) );
	this->txt_confidence->SetValue( wxT("100") );
	wxStaticText* lbl_improvement = new wxStaticText( this, wxID_ANY, wxString::FromAscii("Improvement:"));
	this->txt_improvement = new wxTextCtrl( this, ID_TXT_IMPROVEMENT,wxString::FromAscii("1"), wxDefaultPosition, wxSize(70, 22) );
	this->txt_improvement->SetValue( wxT("1") );
	wxStaticText* lbl_max_chi = new wxStaticText( this, wxID_ANY, wxString::FromAscii("Max Chi-sq P value:"));
	this->txt_max_chi = new wxTextCtrl( this, ID_TXT_MAX_CHI, wxString::FromAscii("1"), wxDefaultPosition, wxSize(70, 22) );
	this->txt_max_chi->SetValue( wxT("1") );
	wxStaticText* lbl_max_depth = new wxStaticText( this, wxID_ANY, wxString::FromAscii("Max Depth:"));
	this->txt_max_depth = new wxTextCtrl( this, ID_TXT_MAX_DEPTH, wxString::FromAscii("3"), wxDefaultPosition, wxSize(70, 22) );
	this->txt_max_depth->SetValue( wxT("3") );
	wxStaticText* lbl_show_results = new wxStaticText( this, wxID_ANY, wxString::FromAscii("Show Results For:") );
	wxArrayString ars_results;
	ars_results.Add(wxString::FromAscii("Both classes"));
	ars_results.Add(wxString::FromAscii("Only Class A"));
	ars_results.Add(wxString::FromAscii("Only Class B"));
	this->cho_show_results = new wxChoice( this, ID_CHO_SHOW_RESULTS, wxDefaultPosition, wxDefaultSize, ars_results );
	this->cho_show_results->SetSelection(0);
	wxStaticText* lbl_disc = new wxStaticText( this, wxID_ANY, wxString::FromAscii("Discretization:") );
	this->discretization = new PanelDiscretization(this, ID_DISCRETIZATION, this->investigation);

	wxStaticBox* static_a = new wxStaticBox(this, wxID_ANY, wxString::FromAscii("Class A Definition"));
    wxStaticBoxSizer* static_class_a_sizer = new wxStaticBoxSizer(static_a, wxHORIZONTAL);
	wxStaticBox* static_b = new wxStaticBox(this, wxID_ANY, wxString::FromAscii("Class B Definition"));
    wxStaticBoxSizer* static_class_b_sizer = new wxStaticBoxSizer(static_b, wxHORIZONTAL);
	if(investigation->is_mac){
		wxFont font_small(lbl_file_name->GetFont());
		font_small.SetPointSize(wxSMALL_FONT->GetPointSize());
		lbl_file_name->SetFont( font_small );
		lbl_support->SetFont( font_small );
		lbl_confidence->SetFont( font_small );
		lbl_improvement->SetFont( font_small );
		lbl_max_chi->SetFont( font_small );
		lbl_max_depth->SetFont( font_small );
		lbl_show_results->SetFont( font_small );
		cho_show_results->SetFont( font_small );
		lbl_disc->SetFont( font_small );
//		wxFont font(wxNORMAL_FONT->GetPointSize(), wxFONTFAMILY_SWISS, wxNORMAL, wxNORMAL);
//		static_a->SetFont( font );
//		static_b->SetFont( font );
	}
	this->limit_A = new PanelLimits(this, ID_PANEL_A, this->investigation);
	this->limit_B = new PanelLimits(this, ID_PANEL_B, this->investigation);

	wxSizer* sizer_ok = CreateButtonSizer(wxOK|wxCANCEL);

	Connect(wxID_OK, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(RulesetDialog::OnClickOk));

	// POSITIONING
	
	wxFlexGridSizer* sizer_upper = new wxFlexGridSizer(8, 2, 0, 0);
	sizer_upper->Add(lbl_file_name, 0, FLAGS, BORDER_PXL);
	sizer_upper->Add(this->txt_file_name, 0, FLAGS, BORDER_PXL);
	sizer_upper->Add(lbl_support, 0, FLAGS, BORDER_PXL);
	sizer_upper->Add(this->txt_support, 0, FLAGS, BORDER_PXL);
	sizer_upper->Add(lbl_confidence, 0, FLAGS, BORDER_PXL);
	sizer_upper->Add(this->txt_confidence, 0, FLAGS, BORDER_PXL);
	sizer_upper->Add(lbl_improvement, 0, FLAGS, BORDER_PXL);
	sizer_upper->Add(this->txt_improvement, 0, FLAGS, BORDER_PXL);
	sizer_upper->Add(lbl_max_chi, 0, FLAGS, BORDER_PXL);
	sizer_upper->Add(this->txt_max_chi, 0, FLAGS, BORDER_PXL);
	sizer_upper->Add(lbl_max_depth, 0, FLAGS, BORDER_PXL);
	sizer_upper->Add(this->txt_max_depth, 0, FLAGS, BORDER_PXL);
	sizer_upper->Add(lbl_show_results, 0, FLAGS, BORDER_PXL);
	sizer_upper->Add(this->cho_show_results, 0, FLAGS, BORDER_PXL);
	sizer_upper->Add(lbl_disc, 0, FLAGS, BORDER_PXL);
	sizer_upper->Add(this->discretization, 0, wxRIGHT|wxTOP|wxBOTTOM|wxALIGN_LEFT|wxALIGN_CENTER_VERTICAL, BORDER_PXL);
	wxBoxSizer* sizer_top = new wxBoxSizer(wxVERTICAL);
	sizer_top->Add(sizer_upper, 0, FLAGS, BORDER_PXL);

	wxBoxSizer* sizer_limits = new wxBoxSizer(wxHORIZONTAL);
	static_class_a_sizer->Add(this->limit_A, 0, FLAGS, BORDER_PXL);
	static_class_b_sizer->Add(this->limit_B, 0, FLAGS, BORDER_PXL);
	sizer_limits->Add(static_class_a_sizer, 0, FLAGS, BORDER_PXL);
	sizer_limits->Add(static_class_b_sizer, 0, FLAGS, BORDER_PXL);
	sizer_top->Add(sizer_limits, 0, FLAGS_RIGHT, BORDER_PXL);
	sizer_top->Add(sizer_ok, 0, FLAGS_RIGHT, BORDER_PXL);
	
	this->SetSizer(sizer_top);
}


void RulesetDialog::OnClickOk( wxCommandEvent& WXUNUSED(event) ){
	
	std::string d( this->investigation->fn_expr );
	std::string f( this->investigation->fn_sa );
	std::string g( this->investigation->fn_ga );

	std::string a = this->limit_A->get_limits();
	std::string b = this->limit_B->get_limits();
	std::string r(""); // Not currently implemented
	std::string s( this->txt_support->GetValue().ToAscii() );
    std::string c( this->txt_confidence->GetValue().ToAscii() );
	std::string i( this->txt_improvement->GetValue().ToAscii() );
	std::string n( this->txt_max_depth->GetValue().ToAscii() );
	std::string o( this->txt_file_name->GetValue().ToAscii() );
	std::string p( this->txt_max_chi->GetValue().ToAscii() );
	std::string w;
	if(this->cho_show_results->GetCurrentSelection()==0)
		w= "both";
	else if( this->cho_show_results->GetCurrentSelection()==1 )
		w = "a";
	else if( this->cho_show_results->GetCurrentSelection()==2 )
		w = "b";

	if(a.size()==0 && b.size()==0){
		wxMessageBox(wxString::FromAscii("Please select a limit for Class A or Class B."), _T("Not Ready"), wxOK | wxICON_EXCLAMATION, this);
		return;
	}
	if( o.size()==0 ){
		wxMessageBox(wxString::FromAscii("Please enter a name for the ruleset file that will be generated"), _T("Not Ready"), wxOK | wxICON_EXCLAMATION, this);
		return;
	}
	if( ! alg::find_first( o, ".ruleset") )
		o = o + ".ruleset";
	
	std::stringstream ss;
	if(this->investigation->is_mac)
		ss << this->investigation->dir_results + '/' + o;
	else
		ss << this->investigation->dir_results + '\\' + o;
	std::string attempted_filename = ss.str();
	int s_int, c_int, i_int;
	try{
		s_int = boost::lexical_cast<int>(s);
		c_int = boost::lexical_cast<int>(c);
		i_int = boost::lexical_cast<int>(i);
	}
	catch( ... ){
		wxMessageBox(wxString::FromAscii("Support, Confidence, and Improvement must be integers"), _T("Not Ready"), wxOK | wxICON_EXCLAMATION, this);
		return;
	}
	if( s_int >100 || s_int < 0 ){
		wxMessageBox(wxString::FromAscii("Support must be between 0 and 100"), _T("Not Ready"), wxOK | wxICON_EXCLAMATION, this);
		return;
	}
	if( c_int >100 || s_int < 0 ){
		wxMessageBox(wxString::FromAscii("Confidence must be between 0 and 100"), _T("Not Ready"), wxOK | wxICON_EXCLAMATION, this);
		return;
	}
	if( i_int >100 || s_int < 0 ){
		wxMessageBox(wxString::FromAscii("Improvement must be between 0 and 100"), _T("Not Ready"), wxOK | wxICON_EXCLAMATION, this);
		return;
	}

	try{
		float p_float = boost::lexical_cast<float>(p);
		p_float=0; // eliminate unused variable warning
	}
	catch( ... ){
		wxMessageBox(wxString::FromAscii("Max Chi2 must be a real-valued number."), _T("Not Ready"), wxOK | wxICON_EXCLAMATION, this);
		return;
	}
	try{
		int n_int = boost::lexical_cast<int>(n);
		n_int = 0; // eliminate unused variable warning
	}
	catch( ... ){
		wxMessageBox(wxString::FromAscii("Maximum depth must be an integer (including zero)."), _T("Not Ready"), wxOK | wxICON_EXCLAMATION, this);
		return;
	}

	std::vector<std::string> cmd;
	cmd.push_back( "miner" );
	cmd.push_back("-d" + d );
	cmd.push_back("-f" + f );
	cmd.push_back("-g" + g );
    cmd.push_back("-a" + a );
    cmd.push_back("-b" + b );
	cmd.push_back("-s" + s );
	cmd.push_back("-c" + c );
	cmd.push_back("-i" + i );
	cmd.push_back("-n" + n );
	cmd.push_back("-p" + p );
	cmd.push_back("-w" + w );
	cmd.push_back("-vT");
    cmd.push_back("-o" + attempted_filename );
    
	//if( r.compare("(none)")==0)  GENE LIST CONSTRAINT NOT IMPLEMENTED
	//	cmd.push_back("-r" + r );
	//cmd.push_back("-k" + k );   TIMEOUT NOT IMPLEMENTED
    std::string disc_method;
    this->discretization->get_discretization_method_name(disc_method);
    //int d_idx = this->discretization->get_discretization_code();
    cmd.push_back( std::string("-m") + disc_method);
    if( disc_method.compare("none")!=0 ){
		std::string lower, upper;
		this->discretization->get_discretization_bounds(lower, upper);
		cmd.push_back("-l" + lower );
        cmd.push_back("-u" + upper );
	}
	ProgressDialog dialog(this, cmd, this->investigation, true);
	if (dialog.ShowModal() == wxID_OK){
        this->new_filename = attempted_filename;
		this->AcceptAndClose();
	}
}
