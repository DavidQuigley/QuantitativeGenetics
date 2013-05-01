#define _CRT_SECURE_NO_DEPRECATE 1
#include <wx/wx.h>
#include <time.h>
#include <wx/dialog.h>
#include <wx/sizer.h>
#include <wx/stattext.h>
#include <wx/button.h>
#include <wx/textctrl.h>
#include <wx/radiobox.h>
#include <wx/panel.h>
#include <wx/arrstr.h>
#include <wx/choice.h>
#include <wx/listbox.h>
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
#include <set>
#include <vector>
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
#include "PanelLimits.h"
#include "DlgProgress.h"
#include "DlgCorrelationGWER.h"

IMPLEMENT_CLASS( CorrelationGWERDialog, wxDialog )

CorrelationGWERDialog::CorrelationGWERDialog( Investigation* investigation ){
	this->investigation = investigation;
	SetExtraStyle(wxWS_EX_BLOCK_EVENTS|wxDIALOG_EX_CONTEXTHELP);
	wxDialog::Create( NULL, wxID_ANY, wxT("Calculate Correlation GWER"), wxDefaultPosition, wxDefaultSize, wxCAPTION|wxRESIZE_BORDER|wxSYSTEM_MENU );
	CreateControls();
	GetSizer()->Fit(this);
	GetSizer()->SetSizeHints(this);
	Centre();
}

void CorrelationGWERDialog::CreateControls(){
	char timebuf[80];
	struct tm* newtime;
	time_t long_time;
	time( &long_time );
	newtime = localtime( &long_time );
	strftime(timebuf, 80, "%Y_%m_%d_%H_%M_%S", newtime);

	wxStaticText* lbl_file_name = new wxStaticText( this, wxID_ANY, wxT("File Name:") );
	this->txt_file_name = new wxTextCtrl( this, ID_CORR_GWER_TXT_FILE_NAME, wxString::FromAscii(timebuf), wxDefaultPosition, wxSize(200, 20) );
	this->txt_file_name->ChangeValue(wxString::FromAscii(timebuf));
	wxStaticText* lbl_n_perms = new wxStaticText( this, wxID_ANY, wxT("Permutations:") );

	this->ars_perms.Add(wxString::FromAscii("10"));
	this->ars_perms.Add(wxString::FromAscii("100"));
	this->ars_perms.Add(wxString::FromAscii("250"));
	this->ars_perms.Add(wxString::FromAscii("1000"));
	this->ars_perms.Add(wxString::FromAscii("10000"));
	this->cho_n_perms = new wxChoice( this, ID_CORR_GWER_N_PERMS, wxDefaultPosition, wxDefaultSize, this->ars_perms);
	this->cho_n_perms->SetSelection(3);

	wxStaticText* lbl_limit_a = new wxStaticText( this, wxID_ANY, wxT("Sample limits:"));
	this->limit_A = new PanelLimits(this, ID_CORR_GWER_PANEL_A, this->investigation);
	this->limit_A->disable_gene_limits();

	wxFlexGridSizer* sizer_top = new wxFlexGridSizer(4, 2, 0, 0);
	sizer_top->Add(lbl_file_name, 0, FLAGS, BORDER_PXL);
	sizer_top->Add(this->txt_file_name, 0, FLAGS, BORDER_PXL);

	sizer_top->Add(lbl_n_perms, 0, FLAGS, BORDER_PXL);
	sizer_top->Add(this->cho_n_perms, 0, FLAGS, BORDER_PXL);

	sizer_top->Add(lbl_limit_a, 0, FLAGS, BORDER_PXL);
	sizer_top->Add(this->limit_A, 0, FLAGS, BORDER_PXL);

	sizer_top->Add( new wxStaticText(this, wxID_ANY, wxT("")), 0, FLAGS, BORDER_PXL);
	wxSizer* sizer_ok = CreateButtonSizer(wxOK|wxCANCEL);
	sizer_top->Add(sizer_ok, 0, FLAGS, BORDER_PXL);

	if(this->investigation->is_mac){
		wxFont font(lbl_file_name->GetFont());
		font.SetPointSize(wxSMALL_FONT->GetPointSize());
		lbl_file_name->SetFont(font);
		lbl_n_perms->SetFont(font);
		cho_n_perms->SetFont(font);
	}
	Connect(wxID_OK, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(CorrelationGWERDialog::OnClickOk));
	this->SetSizer(sizer_top);

}

void CorrelationGWERDialog::OnClickOk(wxCommandEvent& event){

	std::string z( this->ars_perms.Item(  this->cho_n_perms->GetCurrentSelection() ).ToAscii() );
	std::string o( this->txt_file_name->GetValue().ToAscii() );
	if( ! boost::algorithm::find_first( o, ".spear") )
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
	cmd.push_back("-a" + this->limit_A->get_limits() );
	cmd.push_back("-z" + z);
	cmd.push_back("-vT");
	cmd.push_back("-o" + new_file.str() );

	ProgressDialog dialog(this, cmd, this->investigation, true);
	if (dialog.ShowModal() == wxID_OK){
		this->new_filename = new_file.str();
		this->AcceptAndClose();
	}
}
