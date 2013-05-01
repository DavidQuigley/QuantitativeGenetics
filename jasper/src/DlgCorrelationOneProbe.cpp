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
#include <string>
#include <vector>
#include <set>
#include "carmen/src/DataStructures.h"
#include "carmen/src/Attributes.h"
#include "carmen/src/Rawdata.h"
#include "carmen/src/ClassMinerOptions.h"
#include "carmen/src/Discretize.h"
#include "carmen/src/Parser.h"
#include "carmen/src/Dataset.h"
#include "carmen/src/graph.h"
#include "carmen/src/spear.h"
#include "Investigation.h"
using namespace boost;
using namespace std;
#include "DlgProgress.h"
#include "DlgCorrelationOneProbe.h"

IMPLEMENT_CLASS( CorrelationOneProbeDialog, wxDialog )

CorrelationOneProbeDialog::CorrelationOneProbeDialog( Investigation* investigation ){
	this->investigation = investigation;
	SetExtraStyle(wxWS_EX_BLOCK_EVENTS|wxDIALOG_EX_CONTEXTHELP);
	wxDialog::Create( NULL, wxID_ANY, wxT("Calculate Correlation, One Probe"), wxDefaultPosition, wxDefaultSize, wxCAPTION|wxRESIZE_BORDER|wxSYSTEM_MENU );
	CreateControls();
	GetSizer()->Fit(this); 
	GetSizer()->SetSizeHints(this); 
	Centre();
}

void CorrelationOneProbeDialog::CreateControls(){
	char timebuf[80];
	struct tm* newtime;
	time_t long_time;
	time( &long_time );
	newtime = localtime( &long_time );
	strftime(timebuf, 80, "%Y_%m_%d_%H_%M_%S", newtime);

	wxStaticText* lbl_file_name = new wxStaticText( this, wxID_ANY, wxT("File Name:") );
	this->txt_file_name = new wxTextCtrl( this, ID_CORR1P_TXT_FILE_NAME, wxString::FromAscii(timebuf), wxDefaultPosition, wxSize(200, 20) );
	this->txt_file_name->ChangeValue(wxString::FromAscii(timebuf));
    wxStaticText* lbl_probe_id = new wxStaticText( this, wxID_ANY, wxT("Probe or Gene:") );
    this->txt_probe_id= new wxTextCtrl( this, ID_CORR1P_PROBE_ID, wxT(""));
    wxStaticText* lbl_min_corr = new wxStaticText( this, wxID_ANY, wxT("Minimum Correlation: +/-") );
    this->txt_min_corr= new wxTextCtrl( this, ID_CORR1P_TXT_MIN_CORR, wxT("0.5"));
    this->txt_min_corr->SetValue( wxT("0.5") );
	wxSizer* sizer_ok = CreateButtonSizer(wxOK|wxCANCEL);

	if( this->investigation->is_mac ){
		wxFont font(lbl_file_name->GetFont());
		font.SetPointSize(wxSMALL_FONT->GetPointSize());
		lbl_file_name->SetFont(font);
		lbl_probe_id->SetFont(font);
		lbl_min_corr->SetFont(font);
	}

	Connect(wxID_OK, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(CorrelationOneProbeDialog::OnClickOk));

	wxBoxSizer* sizer_top = new wxBoxSizer(wxVERTICAL);
	wxFlexGridSizer* sizer_upper = new wxFlexGridSizer(4, 2, 0, 0);
	sizer_upper->Add(lbl_file_name, 0, FLAGS, BORDER_PXL);
	sizer_upper->Add(this->txt_file_name, 0, FLAGS, BORDER_PXL);
	sizer_upper->Add(lbl_probe_id, 0, FLAGS, BORDER_PXL);
	sizer_upper->Add(this->txt_probe_id, 0, FLAGS, BORDER_PXL);
	sizer_upper->Add(lbl_min_corr, 0, FLAGS, BORDER_PXL);
	sizer_upper->Add(this->txt_min_corr, 0, FLAGS, BORDER_PXL);
	
	sizer_top->Add( sizer_upper, 0, FLAGS, BORDER_PXL);
	sizer_top->Add( sizer_ok, 0, FLAGS_RIGHT, BORDER_PXL);
	SetSizer(sizer_top);
}

void CorrelationOneProbeDialog::OnClickOk(wxCommandEvent& event){

	std::string p( this->txt_probe_id->GetValue().ToAscii() );
	std::string searchterm(p);
	std::string o( this->txt_file_name->GetValue().ToAscii() );
	std::string s( this->txt_min_corr->GetValue().ToAscii() );
	wxArrayString aProbes;
	std::string gene_name_col = this->investigation->ga->get_gene_name_column();
	if( this->investigation->ga->identifier2idx.find(searchterm) == this->investigation->ga->identifier2idx.end() ){
		boost::algorithm::to_lower(searchterm);
		for(int i=0; i<(int)this->investigation->ga->identifiers.size(); i++){
			std::string identifier = this->investigation->ga->identifiers.at(i);
			std::string gene = this->investigation->ga->prop_for_identifier(identifier, gene_name_col);
			boost::algorithm::to_lower(gene);
			if( gene.compare(searchterm)==0 ){
				aProbes.Add( wxString::FromAscii(identifier.c_str()) );
			}
		}
		if( aProbes.size()==0 ){
			wxMessageBox(wxString::FromAscii("Specified gene was not found in data set."), _T("Not Ready"), wxOK | wxICON_EXCLAMATION, this);
			return;
		}
	}
	else
		aProbes.Add( wxString::FromAscii(searchterm.c_str()) );
	if( aProbes.size()==1 ){
        p = aProbes.Item(0).ToAscii();
    }
	else{
		p = wxGetSingleChoice(wxString::FromAscii("Choose a probe"), wxString::FromAscii("Choose a probe"), aProbes).ToAscii();
	}
	float s_f;
	try{
		s_f = boost::lexical_cast<float>(s);
	}
	catch( ... ){
		wxMessageBox(wxString::FromAscii("Minimum correlationmust be a real-valued number."), _T("Not Ready"), wxOK | wxICON_EXCLAMATION, this);
		return;
	}
	if( s_f >1 || s_f < 0 ){
		wxMessageBox(wxString::FromAscii("Minimum absolute correlation must be between 0 and 1 (inclusive)."), _T("Not Ready"), wxOK | wxICON_EXCLAMATION, this);
		return;
	}
	if( ! alg::find_first( o, ".spear") )
		o = o + ".spear";

	std::stringstream new_file;
	if(this->investigation->is_mac)
		new_file << this->investigation->dir_results << "/" << o ;
	else
		new_file << this->investigation->dir_results << "\\" << o ;

	std::vector<std::string> cmd;
	cmd.push_back( "spear" );
	cmd.push_back("-d" + this->investigation->fn_expr );
	cmd.push_back("-f" + this->investigation->fn_sa );
	cmd.push_back("-g" + this->investigation->fn_ga );
	cmd.push_back("-aIDENTIFIER*NULL" );
	cmd.push_back("-s" + s );
	cmd.push_back("-p" + p );
	cmd.push_back("-rF" );
	cmd.push_back("-y" + this->investigation->ga->get_gene_name_column() );
	cmd.push_back("-vT");
	cmd.push_back("-o" + new_file.str() );

	ProgressDialog dialog(this, cmd, this->investigation, true);
	if (dialog.ShowModal() == wxID_OK){
        this->new_filename = new_file.str();
        this->AcceptAndClose();
	}
}
