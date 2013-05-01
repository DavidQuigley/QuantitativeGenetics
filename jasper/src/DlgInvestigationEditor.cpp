#define _SCL_SECURE_NO_DEPRECATE 1
#include <wx/wx.h>
#include <wx/dialog.h>
#include <wx/sizer.h>
#include <wx/stattext.h>
#include <wx/textctrl.h>
#include <wx/choice.h>
#include <wx/button.h>
#include <wx/arrstr.h>
#include <boost/algorithm/string.hpp>
namespace alg = boost::algorithm;
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/random.hpp>
using namespace boost;
using namespace std;
#include "carmen/src/DataStructures.h"
#include "carmen/src/Attributes.h"
#include "carmen/src/Rawdata.h"
#include "carmen/src/Parser.h"
#include "Investigation.h"
#include "PanelInvestigationProperties.h"
#include "DlgInvestigationEditor.h"

IMPLEMENT_CLASS( InvestigationEditorDialog, wxDialog )

InvestigationEditorDialog::InvestigationEditorDialog( Investigation* investigation ){
	this->investigation = investigation;
	this->investigation_name = this->investigation->investigation_name;
	SetExtraStyle(wxWS_EX_BLOCK_EVENTS|wxDIALOG_EX_CONTEXTHELP);
	wxDialog::Create( NULL, wxID_ANY, wxT("Edit Investigation"), wxDefaultPosition, wxDefaultSize, wxCAPTION|wxRESIZE_BORDER|wxSYSTEM_MENU );
	CreateControls();
	GetSizer()->Fit(this);
	GetSizer()->SetSizeHints(this);
	Centre();
}

void InvestigationEditorDialog::CreateControls(){
	wxStaticBoxSizer* static_investigation_sizer = new wxStaticBoxSizer(wxVERTICAL, this, wxString::FromAscii("Current Investigation") );
	wxStaticText* lbl_investigation = new wxStaticText( this, wxID_STATIC, wxT("Investigation Name:"), wxDefaultPosition, wxDefaultSize, 0 );
	wxStaticText* lbl_investigation_value = new wxStaticText( this, wxID_STATIC, wxString::FromAscii(this->investigation->investigation_name.c_str()), wxDefaultPosition, wxDefaultSize, 0 );
	this->pnl_investigation = new InvestigationPropertiesPanel(this, wxID_ANY, this->investigation);
	this->pnl_investigation->set_fn_expr( this->investigation->fn_expr );
	this->pnl_investigation->set_fn_sa( this->investigation->fn_sa );
	this->pnl_investigation->set_fn_ga( this->investigation->fn_ga );
	this->pnl_investigation->set_dir_results( this->investigation->dir_results );
	this->pnl_investigation->set_species( this->investigation->species );
	this->pnl_investigation->set_data_type( this->investigation->data_type );
	this->pnl_investigation->set_gene_name_column( this->investigation->gene_name_column );

	wxSizer* sizer_ok = CreateButtonSizer(wxOK|wxCANCEL);
	Connect(wxID_OK, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(InvestigationEditorDialog::OnClickOK));

	// POSITIONING
	wxBoxSizer* sizer_top = new wxBoxSizer(wxVERTICAL);
	wxBoxSizer* sizer_investigations = new wxBoxSizer(wxHORIZONTAL);
	sizer_investigations->Add(lbl_investigation, 0, FLAGS, BORDER_PXL);
	sizer_investigations->Add(lbl_investigation_value, 0, FLAGS, BORDER_PXL);
	static_investigation_sizer->Add(sizer_investigations, 0, FLAGS, BORDER_PXL);
	static_investigation_sizer->Add(this->pnl_investigation, 0, FLAGS, BORDER_PXL);
	sizer_top->Add(static_investigation_sizer, 0, FLAGS, BORDER_PXL);
	sizer_top->Add(sizer_ok, 0, FLAGS_RIGHT, BORDER_PXL);
	this->SetSizer(sizer_top);
}

void InvestigationEditorDialog::OnClickOK(wxCommandEvent& event ){
	std::string err;
	if( this->pnl_investigation->CheckConsistency(err) ){
		this->fn_expr = this->pnl_investigation->get_fn_expr();
		this->fn_ga = this->pnl_investigation->get_fn_ga();
		this->fn_sa = this->pnl_investigation->get_fn_sa();
		this->dir_results = this->pnl_investigation->get_dir_results();
		this->species = this->pnl_investigation->get_species();
		this->data_type = this->pnl_investigation->get_data_type();
		this->gene_name_column = this->pnl_investigation->get_gene_name_column();
		this->AcceptAndClose();
	}
	else{
		wxMessageBox(wxString::FromAscii(err.c_str()));
	}
}
