#define _CRT_SECURE_NO_DEPRECATE 1
#include <time.h>
#include <wx/dialog.h>
#include <wx/sizer.h>
#include <wx/stattext.h>
#include <wx/button.h>
#include <wx/textctrl.h>
#include <wx/radiobox.h>
#include <wx/panel.h>
#include <wx/statbox.h>
#include <wx/wx.h>
#include <wx/process.h>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/random.hpp>
namespace alg = boost::algorithm;
#include <fstream>
#include <set>
#include <string>
#include <vector>
#include "carmen/src/DataStructures.h"
#include "carmen/src/Attributes.h"
#include "carmen/src/Rawdata.h"
#include "carmen/src/ClassMinerOptions.h"
#include "carmen/src/Discretize.h"
#include "carmen/src/Dataset.h"
#include "carmen/src/Parser.h"
#include "carmen/src/graph.h"
#include "carmen/src/spear.h"
#include "Investigation.h"
using namespace boost;
using namespace std;
#include "PanelLimits.h"
#include "DlgProgress.h"
#include "DlgCorrelation.h"

IMPLEMENT_CLASS( CorrelationDialog, wxDialog )

CorrelationDialog::CorrelationDialog( Investigation* investigation ){
	this->investigation = investigation;
	SetExtraStyle(wxWANTS_CHARS);
	wxDialog::Create( NULL, wxID_ANY, wxT("Calculate Correlation"), wxDefaultPosition, wxDefaultSize, wxCAPTION|wxRESIZE_BORDER|wxSYSTEM_MENU );
	CreateControls();
	GetSizer()->Fit(this); 
	GetSizer()->SetSizeHints(this); 
	Centre();
	redraw();
}

void CorrelationDialog::CreateControls(){
	char timebuf[80];
	struct tm* newtime;
	time_t long_time;
	time( &long_time );
	newtime = localtime( &long_time );
	strftime(timebuf, 80, "%Y_%m_%d_%H_%M_%S", newtime);

	wxStaticBoxSizer* static_sample_sizer = new wxStaticBoxSizer(wxVERTICAL, this, wxString::FromAscii("Sample Limits") );
	wxStaticBoxSizer* static_correlation_sizer = new wxStaticBoxSizer(wxVERTICAL, this, wxString::FromAscii("Correlation Settings") );
	wxStaticBoxSizer* static_probe_sizer = new wxStaticBoxSizer(wxVERTICAL, this, wxString::FromAscii("Probe Limits") );

	wxStaticText* lbl_file_name = new wxStaticText( this, wxID_ANY, wxT("File Name:") );
	if(investigation->is_mac){
    	wxFont font(lbl_file_name->GetFont());
    	font.SetPointSize(wxNORMAL_FONT->GetPointSize());
    	font.SetFamily(wxFONTFAMILY_SWISS);
    	font.SetWeight(wxFONTWEIGHT_BOLD);
		static_correlation_sizer->GetStaticBox()->SetFont( font );
		static_sample_sizer->GetStaticBox()->SetFont( font );
		static_correlation_sizer->GetStaticBox()->SetFont( font );
		static_probe_sizer->GetStaticBox()->SetFont( font );
	}
	this->txt_file_name = new wxTextCtrl( this, ID_CORR_TXT_FILE_NAME, wxString::FromAscii(timebuf), wxDefaultPosition, wxSize(200, 20) );
	this->txt_file_name->ChangeValue(wxString::FromAscii(timebuf));
	wxStaticText* lbl_min_corr = new wxStaticText( this, wxID_ANY, wxT("Minimum Correlation: +/-"));
	this->txt_min_corr = new wxTextCtrl( this, ID_CORR_TXT_MIN_CORR, wxT("0.7"), wxDefaultPosition, wxSize(50, 20) );
	this->txt_min_corr->SetValue( wxT("0.7") );
	wxArrayString methods;
	methods.Add(wxString::FromAscii("Spearman rank correlation") );
	methods.Add(wxString::FromAscii("Pearson correlation") );
	this->cho_method = new wxChoice( this, ID_CORR_CHO_METHOD, wxDefaultPosition, wxDefaultSize, methods);
	this->cho_method->SetSelection(0);

	wxString choices[2];
	choices[0] = wxString::FromAscii("All pairwise correlations");
	choices[1] = wxString::FromAscii("One probe vs. all others");
	this->rdo_single = new wxRadioBox(this, ID_RADIO_SINGLE, wxString::FromAscii(""), wxDefaultPosition, wxSize(200, 70), 2, choices, 2, wxRA_SPECIFY_ROWS);
	this->lbl_probe_id = new wxStaticText( this, wxID_ANY, wxT("Seed:"));
	this->txt_probe_id= new wxTextCtrl( this, ID_PROBE_ID, wxT(""), wxDefaultPosition, wxSize(200, 20) );
	this->lbl_probe_list = new wxStaticText( this, wxID_ANY, wxT("Limit to probes:"));
	this->txt_probe_list = new wxTextCtrl( this, ID_CORR_PROBELIST, wxT(""), wxDefaultPosition, wxSize(150, 100), wxTE_MULTILINE);

	this->chk_differential= new wxCheckBox(this, ID_CORR_CHK_DIFFERENTIAL, wxString::FromAscii(""));
	wxStaticText* lbl_differential = new wxStaticText( this, wxID_ANY, wxT("Calculate differential correlation"));

	this->lbl_min_delta_corr = new wxStaticText( this, wxID_ANY, wxT("Minimum Change in Correlation:"));
	this->txt_min_delta_corr = new wxTextCtrl( this, ID_CORR_TXT_DELTA_CORR, wxT("0.0"), wxDefaultPosition, wxSize(50, 20) );
	this->lbl_condition1 = new wxStaticText( this, wxID_ANY, wxT("Sample group 1"));
	this->lbl_condition2 = new wxStaticText( this, wxID_ANY, wxT("Sample group 2"));
	this->limit_A = new PanelLimits(this, ID_CORR_PANEL_A, investigation);
	this->limit_B = new PanelLimits(this, ID_CORR_PANEL_B, investigation);
    this->limit_A->disable_gene_limits();
    this->limit_B->disable_gene_limits();

    wxStaticText* lbl_min_present = new wxStaticText( this, wxID_ANY, wxT("Minimum % Present:"));
	this->txt_min_present = new wxTextCtrl( this, ID_CORR_TXT_MIN_PRESENT, wxT("90"), wxDefaultPosition, wxSize(50, 20) );
	this->txt_min_present->SetValue( wxT("90") );
	wxStaticText* lbl_min_var = new wxStaticText( this, wxID_ANY, wxT("Minimum Variance:"));
	this->txt_min_var = new wxTextCtrl( this, ID_CORR_TXT_MIN_VAR, wxT("0.0"), wxDefaultPosition, wxSize(50, 20) );
	this->txt_min_var->SetValue( wxT("0.0") );
	this->btn_GO = new wxButton( this, ID_CORR_BTN_GO, wxT("&Add by Ontology"), wxDefaultPosition, wxSize(120,20));

	wxSizer* sizer_ok = CreateButtonSizer(wxOK|wxCANCEL);

	if(investigation->is_mac){
		wxFont font_small(lbl_file_name->GetFont());
		font_small.SetPointSize(wxSMALL_FONT->GetPointSize());
		lbl_file_name->SetFont( font_small );
		cho_method->SetFont( font_small );
		//rdo_single->SetFont( font_small );
		lbl_min_corr->SetFont( font_small );
		lbl_probe_id->SetFont( font_small );
		lbl_differential->SetFont( font_small );
		lbl_min_delta_corr->SetFont( font_small );
		lbl_condition1->SetFont( font_small );
		lbl_condition2->SetFont( font_small );
		lbl_min_present->SetFont( font_small );
		lbl_min_var->SetFont( font_small );
		btn_GO->SetFont( font_small );
	}

	Connect(wxID_OK, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(CorrelationDialog::OnClickOk));
	Connect(ID_CORR_PANEL_A, wxEVT_COMMAND_TEXT_UPDATED, wxCommandEventHandler(CorrelationDialog::OnUpdateLabel));
	Connect(ID_CORR_PANEL_B, wxEVT_COMMAND_TEXT_UPDATED, wxCommandEventHandler(CorrelationDialog::OnUpdateLabel));
	Connect(ID_RADIO_SINGLE, wxEVT_COMMAND_RADIOBOX_SELECTED, wxCommandEventHandler(CorrelationDialog::OnClickSingle));
	Connect(ID_CORR_CHK_DIFFERENTIAL, wxEVT_COMMAND_CHECKBOX_CLICKED, wxCommandEventHandler(CorrelationDialog::OnCheckDifferential));
	Connect(ID_CORR_BTN_GO, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(CorrelationDialog::OnClickGO));
	Connect(wxEVT_KEY_DOWN, wxKeyEventHandler(CorrelationDialog::OnChar));

	this->sizer_top = new wxBoxSizer(wxVERTICAL);
	wxBoxSizer* sizer_filename = new wxBoxSizer(wxHORIZONTAL);
	sizer_filename->Add(lbl_file_name, 0, FLAGS, BORDER_PXL);
	sizer_filename->Add(this->txt_file_name, 0, FLAGS, BORDER_PXL);
	this->sizer_top->Add( sizer_filename, 0, FLAGS, BORDER_PXL);


	// CORRELATION SETTINGS
	wxBoxSizer* sizer_min_corr = new wxBoxSizer(wxHORIZONTAL);
	sizer_min_corr->Add(lbl_min_corr, 0, FLAGS, BORDER_PXL);
	sizer_min_corr->Add(this->txt_min_corr, 0, FLAGS, 0);
	static_correlation_sizer->Add(sizer_min_corr, 0, FLAGS, BORDER_PXL);
	static_correlation_sizer->Add(this->cho_method, 0, FLAGS, BORDER_PXL);

	// radio, seed probe, probe list
	static_correlation_sizer->Add(this->rdo_single, 0, FLAGS, BORDER_PXL);
	wxBoxSizer* sizer_seed = new wxBoxSizer(wxHORIZONTAL);
	sizer_seed->Add(this->lbl_probe_id, 0, FLAGS, BORDER_PXL);
	sizer_seed->Add(this->txt_probe_id, 0, FLAGS, BORDER_PXL);
	rdo_single->SetSelection(1);
	static_correlation_sizer->Add(sizer_seed, 0, FLAGS, BORDER_PXL);
	static_correlation_sizer->Add(this->lbl_probe_list, 0, FLAGS_TOP_LEFT, BORDER_PXL);
	wxBoxSizer* sizer_probe_list= new wxBoxSizer(wxHORIZONTAL);
	sizer_probe_list->Add(this->txt_probe_list, 0, FLAGS, BORDER_PXL);
	sizer_probe_list->Add(this->btn_GO, 0, FLAGS_TOP_LEFT, BORDER_PXL);
	static_correlation_sizer->Add(sizer_probe_list, 0, FLAGS, BORDER_PXL);

	wxBoxSizer* sizer_differential= new wxBoxSizer(wxHORIZONTAL);
	sizer_differential->Add(this->chk_differential, 0, FLAGS, BORDER_PXL);
	sizer_differential->Add(lbl_differential, 0, FLAGS, BORDER_PXL);
	static_correlation_sizer->Add(sizer_differential, 0, FLAGS, BORDER_PXL);

	wxBoxSizer* sizer_delta= new wxBoxSizer(wxHORIZONTAL);
	sizer_delta->Add(lbl_min_delta_corr, 0, FLAGS, BORDER_PXL);
	sizer_delta->Add(this->txt_min_delta_corr, 0, FLAGS, BORDER_PXL);
	static_correlation_sizer->Add(sizer_delta, 0, FLAGS, BORDER_PXL);

	this->sizer_top->Add( static_correlation_sizer, 0, FLAGS, BORDER_PXL);

	// SAMPLE LIMITS
	wxFlexGridSizer* sizer_limits = new wxFlexGridSizer(2, 2, 0, 0);
	sizer_limits->Add(this->lbl_condition1, 0, FLAGS, BORDER_PXL);
	sizer_limits->Add(this->limit_A, 0, FLAGS, BORDER_PXL);
	sizer_limits->Add(this->lbl_condition2, 0, FLAGS, BORDER_PXL);
	sizer_limits->Add(this->limit_B, 0, FLAGS, BORDER_PXL);
	static_sample_sizer->Add(sizer_limits, 0, FLAGS, BORDER_PXL);
	this->sizer_top->Add( static_sample_sizer, 0, FLAGS, BORDER_PXL);

	// PROBE LIMITS
	wxFlexGridSizer* sizer_probe_limits= new wxFlexGridSizer(3, 2, 0, 0);
	sizer_probe_limits->Add(lbl_min_present, 0, FLAGS, BORDER_PXL);
	sizer_probe_limits->Add(this->txt_min_present, 0, FLAGS, BORDER_PXL);
	sizer_probe_limits->Add(lbl_min_var, 0, FLAGS, BORDER_PXL);
	sizer_probe_limits->Add(this->txt_min_var, 0, FLAGS, BORDER_PXL);

	static_probe_sizer->Add(sizer_probe_limits, 0, FLAGS, BORDER_PXL);
	this->sizer_top->Add( static_probe_sizer, 0, FLAGS, BORDER_PXL);
	this->sizer_top->Add( sizer_ok, 0, FLAGS_RIGHT, BORDER_PXL);
	this->sizer_top->Layout();
	this->SetSizer(this->sizer_top);

	this->redraw();
}

void CorrelationDialog::OnChar(wxKeyEvent& event){
	if (event.GetKeyCode() == WXK_ESCAPE){
		this->EndModal(wxCANCEL);
	}
	event.Skip();
};


void CorrelationDialog::redraw(){
	this->txt_min_delta_corr->Show(false);
	this->lbl_min_delta_corr->Show(false);
	this->lbl_probe_id ->Show(false);
	this->txt_probe_id->Show(false);
	this->limit_B->Show(false);
	this->lbl_condition1->Show(false);
	this->lbl_condition2->Show(false);
	this->lbl_probe_list->Show(false);
	this->txt_probe_list->Show(false);
	this->btn_GO->Show(false);
	if(this->rdo_single->GetSelection()==0){
		this->txt_probe_id->SetValue(wxString::FromAscii(""));
		this->lbl_probe_list->Show(true);
		this->txt_probe_list->Show(true);
		this->btn_GO->Show(true);
	}
	else{
		this->txt_probe_id->Show(true);
		this->lbl_probe_id->Show(true);
	}
	if( this->chk_differential->IsChecked() ){
		this->limit_B->Show(true);
		this->txt_min_delta_corr->Show(true);
		this->lbl_min_delta_corr->Show(true);
		this->lbl_condition1->Show(true);
		this->lbl_condition2->Show(true);
	}
	else{
		this->limit_B->Clear();
		this->txt_min_delta_corr->SetValue(wxString::FromAscii("0.0") );
	}
	this->sizer_top->Fit(this);
}

void CorrelationDialog::OnClickGO(wxCommandEvent& event){
	GOFilterDialog dialog(this->investigation);
	if (dialog.ShowModal() == wxID_OK){
		this->txt_probe_list->AppendText( wxString::FromAscii(dialog.get_probe_list().c_str()) );
	}
	this->redraw();
}

void CorrelationDialog::OnClickOk(wxCommandEvent& event){
	std::string d( this->investigation->fn_expr );
	std::string f( this->investigation->fn_sa );
	std::string g( this->investigation->fn_ga );
	std::string a = this->limit_A->get_limits();
	std::string b = this->limit_B->get_limits();

	if( this->chk_differential->IsChecked() ){
		if( a.size()==0 || b.size()==0 ){
			wxMessageBox(wxString::FromAscii("For differential correlation you must place limits on both sample groups."), _T("Not Ready"), wxOK | wxICON_EXCLAMATION, this);
			return;
		}
	}
	else if( a.size()==0 ){
		a = std::string("IDENTIFIER*NULL");
	}
	std::string o( this->txt_file_name->GetValue().ToAscii() );
	std::string m( this->txt_min_var->GetValue().ToAscii() );
	std::string s( this->txt_min_corr->GetValue().ToAscii() );
	std::string c( this->txt_min_delta_corr->GetValue().ToAscii() );
	std::string present( this->txt_min_present->GetValue().ToAscii() );
	std::string p("");
	std::string r("F");
	std::string l("F");
	std::string y( this->investigation->ga->get_gene_name_column() );
	std::string t("spearman");
	if( this->cho_method->GetSelection()==1 )
		t = std::string("pearson");
	wxArrayString aProbes;
	if(this->rdo_single->GetSelection()==0){
		std::vector<std::string> probes_to_use_v;
		std::string probes_to_use( this->txt_probe_list->GetValue().ToAscii() );
		if(probes_to_use.size()>0){
			alg::replace_all(probes_to_use, std::string("\n"), std::string(","));
			p = probes_to_use;
		}
		if(p.size()>0){
			r = std::string("T");
			l = std::string("T");
		}
	}
	else{
		p = this->txt_probe_id->GetValue().ToAscii();
		std::string searchterm(p);
		std::string searchterm_LC(searchterm);
		std::string searchterm_UC(searchterm);
		boost::algorithm::to_lower(searchterm_LC);
		boost::algorithm::to_upper(searchterm_UC);
		if( this->investigation->ga->identifier2idx.find(searchterm) != this->investigation->ga->identifier2idx.end() ){
			aProbes.Add(wxString::FromAscii(p.c_str()));
		}
		else if( this->investigation->ga->identifier2idx.find(searchterm_LC) != this->investigation->ga->identifier2idx.end() ){
			aProbes.Add(wxString::FromAscii(p.c_str()));  // check lower case version
		}
		else if( this->investigation->ga->identifier2idx.find(searchterm_UC) != this->investigation->ga->identifier2idx.end() ){
			aProbes.Add(wxString::FromAscii(p.c_str()));  // check upper case version
		}
		else{
			std::string gene_name_column = this->investigation->ga->get_gene_name_column();
			for(int i=0; i<(int)this->investigation->ga->identifiers.size(); i++){
				std::string identifier = this->investigation->ga->identifiers.at(i);
				std::string gene = this->investigation->ga->prop_for_identifier(identifier, gene_name_column );
				boost::algorithm::to_lower(gene);
				if( gene.compare(searchterm_LC)==0 ){
					aProbes.Add( wxString::FromAscii(identifier.c_str()) );
				}
			}
		}
		if( aProbes.size()==0 ){
			wxMessageBox(wxString::FromAscii("Specified gene was not found in data set."), _T("Not Ready"), wxOK | wxICON_EXCLAMATION, this);
			return;
		}

		if( aProbes.size()==1 ){
			p = aProbes.Item(0).ToAscii();
		}
		else{
			//wxString chosen_probe = wxGetSingleChoice(wxString::FromAscii("Choose a probe"), wxString::FromAscii("Choose a probe"), aProbes).ToAscii();
            wxString chosen_probe = wxT(wxGetSingleChoice(wxString(wxT("Choose a probe")), wxString(wxT("Choose a probe")), aProbes));

			if( chosen_probe.Len()==0 ){
				return;
			}
			else
				p = chosen_probe.ToAscii();
		}
	}
	if( ! alg::find_first( o, ".spear") )
		o = o + ".spear";
	std::string attempted_filename;
	std::stringstream ss;
	if(this->investigation->is_mac)
		ss << this->investigation->dir_results << "/" << o;
	else
		ss << this->investigation->dir_results << "\\" << o;
	attempted_filename = ss.str();

	float m_f, c_f, s_f, present_f;
	try{
		m_f = boost::lexical_cast<float>(m);
		c_f = boost::lexical_cast<float>(c);
		s_f = boost::lexical_cast<float>(s);
		present_f = boost::lexical_cast<float>(present);
	}
	catch( ... ){
		wxMessageBox(wxString::FromAscii("Min Var, Min Corr, Min Corr Change, and Precent Present must all be real-valued numbers."), _T("Not Ready"), wxOK | wxICON_EXCLAMATION, this);
		return;
	}
	if( m_f >1 || m_f < 0 ){
		wxMessageBox(wxString::FromAscii("Minimum variance must be between 0 and 1 (inclusive)."), _T("Not Ready"), wxOK | wxICON_EXCLAMATION, this);
		return;
	}
	if( s_f >1 || s_f < 0 ){
		wxMessageBox(wxString::FromAscii("Minimum absolute correlation must be between 0 and 1 (inclusive)."), _T("Not Ready"), wxOK | wxICON_EXCLAMATION, this);
		return;
	}
	if( c_f >2 || c_f < -2 ){
		wxMessageBox(wxString::FromAscii("Minimum delta correlation must be between -2 and 2 (inclusive)."), _T("Not Ready"), wxOK | wxICON_EXCLAMATION, this);
		return;
	}
	if( present_f >100 || present_f < 0 ){
		wxMessageBox(wxString::FromAscii("Minimum percent present must be between 0 and 100 (inclusive)."), _T("Not Ready"), wxOK | wxICON_EXCLAMATION, this);
		return;
	}
	std::vector<std::string> cmd;
	cmd.push_back( "spear" );
	cmd.push_back("-d" + d );
	cmd.push_back("-f" + f );
	cmd.push_back("-g" + g );
	cmd.push_back("-a" + a );
    cmd.push_back("-b" + b );
	cmd.push_back("-m" + m );
	cmd.push_back("-r" + r );
	cmd.push_back("-l" + l );
	present = boost::lexical_cast<std::string>(present_f / 100);
	cmd.push_back("-n" + present);
	cmd.push_back("-c" + c );
	cmd.push_back("-s" + s );
	cmd.push_back("-p" + p );
	cmd.push_back("-y" + y );
	cmd.push_back("-t" + t );
	cmd.push_back("-vT");
	cmd.push_back("-o" + attempted_filename );
	ProgressDialog dialog(this, cmd, this->investigation, true);
	if (dialog.ShowModal() == wxID_OK){
        this->new_filename = attempted_filename;
        this->AcceptAndClose();
	}
}

void CorrelationDialog::OnUpdateLabel(wxCommandEvent& evt){
	redraw();
}

void CorrelationDialog::OnClickSingle(wxCommandEvent& evt){
	redraw();
}

void CorrelationDialog::OnCheckDifferential(wxCommandEvent& evt){
	redraw();
}
