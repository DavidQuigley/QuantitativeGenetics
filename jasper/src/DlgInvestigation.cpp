#define _CRT_SECURE_NO_DEPRECATE 1
#include <string>
#include <vector>
#include "wx/wx.h"
#include "wx/sizer.h"
#include "wx/stattext.h"
#include "wx/textctrl.h"
#include "wx/button.h"
using namespace std;
#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/random.hpp>
namespace alg = boost::algorithm;
#include "carmen/src/DataStructures.h"
#include "carmen/src/Attributes.h"
#include "carmen/src/Rawdata.h"
#include "carmen/src/Parser.h"
#include "Investigation.h"
#include "PanelInvestigationProperties.h"
#include "DlgInvestigation.h"

IMPLEMENT_CLASS( InvestigationDialog, wxDialog )

InvestigationDialog::InvestigationDialog( Investigation* investigation){
	this->investigation = investigation;
	SetExtraStyle(wxWS_EX_BLOCK_EVENTS|wxDIALOG_EX_CONTEXTHELP);
	wxDialog::Create( NULL, wxID_ANY, wxT("Create New Investigation"), wxDefaultPosition, wxDefaultSize, wxCAPTION|wxRESIZE_BORDER|wxSYSTEM_MENU );
	this->fn_expr = std::string("");
	this->fn_sa = std::string("");
	this->fn_ga = std::string("");
	this->dir_results = std::string("");
	this->species = std::string("");
	this->data_type = std::string("");
	this->gene_name_column = std::string("");
	CreateControls();
	GetSizer()->Fit(this);
	GetSizer()->SetSizeHints(this);
	Centre();
}


void InvestigationDialog::CreateControls()
{
	// CONTROLS
	wxStaticText* lbl_investigations = new wxStaticText( this, wxID_STATIC, wxT("Investigation Name:"), wxDefaultPosition, wxDefaultSize, 0 );
	this->txt_investigations = new wxTextCtrl ( this, ID_INVESTIGATION, wxEmptyString, wxDefaultPosition, wxSize(200,22), 0 );
	this->pnl_investigation = new InvestigationPropertiesPanel(this, wxID_ANY, this->investigation);
	wxSizer* sizer_ok = CreateButtonSizer(wxOK|wxCANCEL);
	if(investigation->is_mac){
		// mac uses giant fonts by default; PC is better.
		wxFont font(this->txt_investigations->GetFont());
		font.SetPointSize(wxSMALL_FONT->GetPointSize());
		this->txt_investigations->SetFont(font);
	}
	Connect(wxID_OK, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(InvestigationDialog::OnClickOk));

	// POSITIONING
	wxBoxSizer* sizer_top = new wxBoxSizer(wxVERTICAL);
	wxBoxSizer* sizer_upper = new wxBoxSizer(wxHORIZONTAL);
	sizer_upper->Add(lbl_investigations, 0, FLAGS, BORDER_PXL);
	sizer_upper->Add(txt_investigations, 0, FLAGS, BORDER_PXL);
	sizer_top->Add(sizer_upper, 0, FLAGS, BORDER_PXL);
	sizer_top->Add(this->pnl_investigation, 0, FLAGS, BORDER_PXL);
	sizer_top->Add(sizer_ok, 0, FLAGS_RIGHT, BORDER_PXL);
	this->SetSizer(sizer_top);
}


void InvestigationDialog::OnClickOk( wxCommandEvent& event ){

	std::string err;
	if( this->txt_investigations->GetValue().Len() == 0 ){
		wxMessageBox(wxString::FromAscii("Please choose a name for the investigation"), wxString::FromAscii("ERROR: missing investigation name"), wxOK | wxICON_EXCLAMATION, this);
		return;
	}
	if( this->pnl_investigation->CheckConsistency(err) ){
		this->fn_expr = this->pnl_investigation->get_fn_expr();
		this->fn_ga = this->pnl_investigation->get_fn_ga();
		this->fn_sa = this->pnl_investigation->get_fn_sa();
		this->dir_results = this->pnl_investigation->get_dir_results();
		this->species = this->pnl_investigation->get_species();
		this->data_type = this->pnl_investigation->get_data_type();
		this->gene_name_column = this->pnl_investigation->get_gene_name_column();
		this->investigation_name = txt_investigations->GetValue().ToAscii();
		this->AcceptAndClose();
	}
	else{
		wxMessageBox(wxString::FromAscii(err.c_str() ), wxString::FromAscii(err.c_str()), wxOK | wxICON_EXCLAMATION, this);
	}
}


