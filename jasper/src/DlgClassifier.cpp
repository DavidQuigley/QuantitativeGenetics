#define _CRT_SECURE_NO_DEPRECATE 1
#define _SCL_SECURE_NO_DEPRECATE 1
#include <wx/dialog.h>
#include <wx/sizer.h>
#include <wx/stattext.h>
#include <wx/button.h>
#include <wx/textctrl.h>
#include <wx/choice.h>
#include <wx/checkbox.h>
#include <wx/wx.h>
#include <wx/arrstr.h>
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
#include "DlgProgress.h"
#include "DlgClassifier.h"

IMPLEMENT_CLASS( ClassifierDialog, wxDialog )

ClassifierDialog::ClassifierDialog( Investigation* investigation, std::string filename ){
	this->investigation = investigation;
	this->filename = filename;
	SetExtraStyle(wxWS_EX_BLOCK_EVENTS|wxDIALOG_EX_CONTEXTHELP);
	wxDialog::Create( NULL, wxID_ANY, wxT("Create Classifier"), wxDefaultPosition, wxDefaultSize, wxCAPTION|wxRESIZE_BORDER|wxSYSTEM_MENU );
	CreateControls(filename);
	GetSizer()->Fit(this); 
	GetSizer()->SetSizeHints(this); 
	Centre();
	redraw();
}

void ClassifierDialog::CreateControls(std::string filename){
	
	std::vector<std::string> path_parts;
	alg::split( path_parts, filename, alg::is_any_of("\\") );
	std::string fn_ruleset = path_parts.at( path_parts.size()-1 );
	alg::replace_all(fn_ruleset, ".ruleset", "");
	
	wxStaticText* lbl_file_name = new wxStaticText( this, wxID_ANY, wxT("File Name:") );
	this->txt_file_name = new wxTextCtrl( this, ID_CLASS_TXT_FILE_NAME, wxString::FromAscii(fn_ruleset.c_str()), wxDefaultPosition, wxSize(200, 20) );
	
	wxStaticText* lbl_method = new wxStaticText( this, wxID_ANY, wxT("Method:") );
	wxArrayString ars_method;
	ars_method.Add(wxT("Boosting"));
	ars_method.Add(wxT("Sum of rules"));
	this->cho_method = new wxChoice( this, ID_CLASS_CHO_METHOD, wxDefaultPosition, wxDefaultSize, ars_method );
	this->cho_method->SetSelection(0);
	
	this->lbl_rounds = new wxStaticText( this, wxID_ANY, wxT("Boosting rounds:") );
	this->txt_rounds = new wxTextCtrl( this, ID_CLASS_TXT_ROUNDS, wxT("20"), wxDefaultPosition, wxSize(70, 20) );
	
	this->lbl_max_rules = new wxStaticText( this, wxID_ANY, wxT("Max rules per class:") );
	this->txt_max_rules = new wxTextCtrl( this, ID_CLASS_TXT_MAX_RULES, wxT("100"), wxDefaultPosition, wxSize(70, 20) );
	
	wxStaticText* lbl_loocv = new wxStaticText( this, wxID_ANY, wxT("Leave One Out CV?") );
	this->chk_loocv = new wxCheckBox(this, ID_CLASS_CHK_LOOCV, wxString::FromAscii(""));

	wxSizer* sizer_ok = CreateButtonSizer( wxOK | wxCANCEL );

	Connect(wxID_OK, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(ClassifierDialog::OnClickOk));
	Connect(ID_CLASS_CHO_METHOD, wxEVT_COMMAND_CHOICE_SELECTED, wxCommandEventHandler(ClassifierDialog::OnClickMethod));

	wxBoxSizer* sizer_top = new wxBoxSizer(wxVERTICAL);
	
	wxFlexGridSizer* sizer_upper = new wxFlexGridSizer(5, 2, 0, 0);
	sizer_upper->Add(lbl_file_name, 0, FLAGS, BORDER_PXL);
	sizer_upper->Add(this->txt_file_name, 0, FLAGS, BORDER_PXL);

	sizer_upper->Add(lbl_method, 0, FLAGS, BORDER_PXL);
	sizer_upper->Add(this->cho_method, 0, FLAGS, BORDER_PXL);

	sizer_upper->Add(lbl_rounds, 0, FLAGS, BORDER_PXL);
	sizer_upper->Add(this->txt_rounds, 0, FLAGS, BORDER_PXL);

	sizer_upper->Add(lbl_max_rules, 0, FLAGS, BORDER_PXL);
	sizer_upper->Add(this->txt_max_rules, 0, FLAGS, BORDER_PXL);

	sizer_upper->Add(lbl_loocv, 0, FLAGS, BORDER_PXL);
	sizer_upper->Add(this->chk_loocv, 0, FLAGS, BORDER_PXL);
	
	sizer_top->Add(sizer_upper, 0, FLAGS, BORDER_PXL);
	sizer_top->Add(sizer_ok, 0, FLAGS_RIGHT, BORDER_PXL);

	SetSizer(sizer_top);
}

void ClassifierDialog::OnClickOk(wxCommandEvent& event){
	std::string d( this->investigation->fn_expr );
	std::string f( this->investigation->fn_sa );
	std::string g( this->investigation->fn_ga );
	std::string rounds( this->txt_rounds->GetValue().ToAscii() );
	std::string max_rules( this->txt_max_rules->GetValue().ToAscii() );
	std::string o( this->txt_file_name->GetValue().ToAscii() );

	std::vector<std::string> cmd;
	cmd.push_back( "buildclassifier" );
	cmd.push_back("-d" + d );
	cmd.push_back("-f" + f );
	cmd.push_back("-g" + g );
	cmd.push_back("-r" + this->filename);
	if( this->chk_loocv->IsChecked())
		cmd.push_back("-xT");
	if( this->cho_method->GetCurrentSelection()==0){
		cmd.push_back("-mboosting");
		int r_int;
		try{
			r_int = boost::lexical_cast<int>(rounds);
		}
		catch( ... ){
			wxMessageBox(wxString::FromAscii("Number of rounds must be a number greater than zero."), _T("Not Ready"), wxOK | wxICON_EXCLAMATION, this);
			return;
		}
		if( r_int < 0 ){
			wxMessageBox(wxString::FromAscii("Number of rounds must be a number greater than zero."), _T("Not Ready"), wxOK | wxICON_EXCLAMATION, this);
			return;
		}
		cmd.push_back("-n" + rounds);
	}
	else{
		cmd.push_back("-mcolsums");
		int r_int;
		try{
			r_int = boost::lexical_cast<int>(max_rules);
		}
		catch( ... ){
			wxMessageBox(wxString::FromAscii("Maximum number of rules per category must be a number greater than zero."), _T("Not Ready"), wxOK | wxICON_EXCLAMATION, this);
			return;
		}
		if( r_int < 0 ){
			wxMessageBox(wxString::FromAscii("Maximum number of rules per category must be a number greater than zero."), _T("Not Ready"), wxOK | wxICON_EXCLAMATION, this);
			return;
		}
		cmd.push_back("-t" + max_rules);
	}
	
	if( o.size()==0 ){
		wxMessageBox(wxString::FromAscii("Please enter a name for the classifier file that will be generated"), _T("Not Ready"), wxOK | wxICON_EXCLAMATION, this);
		return;
	}
	if( ! alg::find_first( o, ".classifier") )
		o = o + ".classifier";
	std::stringstream ss;
	if( this->investigation->is_mac)
		ss << this->investigation->dir_results << '/' << o;
	else
		ss << this->investigation->dir_results << '\\' << o;
	std::string attempted_filename = ss.str();
	cmd.push_back("-vT");
    cmd.push_back("-o" + attempted_filename );
	ProgressDialog dialog(this, cmd, this->investigation, false);
	if (dialog.ShowModal() == wxID_OK){
        this->new_filename = attempted_filename;
		this->AcceptAndClose();
	}
}

void ClassifierDialog::redraw(){
	this->txt_max_rules->Enable(false);
	this->txt_rounds->Enable(false);
	this->lbl_max_rules->Enable(false);
	this->lbl_rounds->Enable(false);


	if( this->cho_method->GetCurrentSelection() == 0 ){
		this->txt_rounds->Enable(true);
		this->lbl_rounds->Enable(true);
	}
	else{
		this->txt_max_rules->Enable(true);
		this->lbl_max_rules->Enable(true);
	}
}

void ClassifierDialog::OnClickMethod(wxCommandEvent& evt){
	this->redraw();
}
