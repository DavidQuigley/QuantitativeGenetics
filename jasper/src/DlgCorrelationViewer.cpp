#define _CRT_SECURE_NO_DEPRECATE 1
#include <stdio.h>
#include <string>
#include <algorithm>
#include <iostream>
#include <math.h>
#include <fstream>
#include <vector>
#include <set>
using namespace std;
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/random.hpp>
using namespace boost;
namespace alg = boost::algorithm;
#include <wx/wx.h>
#include <wx/grid.h>
#include <wx/choicdlg.h>
#include <wx/process.h>
#include <wx/mimetype.h>
#include "carmen/src/DataStructures.h"
#include "carmen/src/Attributes.h"
#include "carmen/src/Rawdata.h"
#include "carmen/src/graph.h"
#include "carmen/src/Parser.h"
#include "Investigation.h"
#include "GridCorrelation.h"
#include "DlgCorrelationViewer.h"

IMPLEMENT_CLASS( CorrelationViewerDialog, wxDialog )


CorrelationViewerDialog::CorrelationViewerDialog( std::string filename, Investigation* investigation ){
	this->investigation = investigation;
	SetExtraStyle(wxWS_EX_BLOCK_EVENTS|wxDIALOG_EX_CONTEXTHELP);
	wxDialog::Create( NULL, wxID_ANY, wxT("Correlation Explorer"), wxDefaultPosition, wxSize(500,500), wxCAPTION|wxRESIZE_BORDER|wxSYSTEM_MENU );
	CreateControls(filename);
	GetSizer()->Fit(this); 
	GetSizer()->SetSizeHints(this); 
	Centre();
}

void CorrelationViewerDialog::CreateControls(std::string filename){
	wxStaticText* lbl_gene_name = new wxStaticText( this, wxID_ANY, wxT("Gene Name:") );
	this->txt_gene_name = new wxTextCtrl( this, ID_COR_V_TXT_GENE_NAME, wxT(""), wxDefaultPosition, wxSize(100, 22) );
	
	wxStaticText* lbl_min_abs_corr = new wxStaticText( this, wxID_ANY, wxT("Minimum Correlation: +/-") );
	this->txt_min_abs_corr = new wxTextCtrl( this, ID_COR_V_TXT_MIN_ABS_CORR, wxT("0.7"), wxDefaultPosition, wxSize(70, 22) );
	
    this->btn_save = new wxButton( this, ID_COR_V_BTN_SAVE, wxT("Save") );
	this->btn_go = new wxButton( this, ID_COR_V_BTN_GO, wxT("Go") );
	this->btn_go->SetDefault();

	wxStaticText* lbl_current_probe_id = new wxStaticText( this, wxID_ANY, wxT("Probe ID:") );
	this->lbl_current_probe_id_value = new wxStaticText( this, wxID_ANY, wxT("(none)") );

	this->grid_correlation = new GridCorrelation(this, this->investigation, filename, ID_COR_V_GRID_CORR, wxSize(400,500));
	this->btn_close = new wxButton( this, wxID_CANCEL, wxT("Close") );

	Connect(ID_COR_V_BTN_GO, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(CorrelationViewerDialog::OnClickGo));
    Connect(ID_COR_V_BTN_SAVE, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(CorrelationViewerDialog::OnClickSave));

	wxBoxSizer* sizer_top = new wxBoxSizer(wxVERTICAL);
	wxFlexGridSizer* sizer_upper = new wxFlexGridSizer(4, 2, 0, 0);
	sizer_upper->Add(lbl_gene_name, 0, FLAGS, BORDER_PXL);
	sizer_upper->Add(this->txt_gene_name, 0, FLAGS, BORDER_PXL);

	sizer_upper->Add(lbl_min_abs_corr, 0, FLAGS, BORDER_PXL);
	sizer_upper->Add(this->txt_min_abs_corr, 0, FLAGS, BORDER_PXL);
	
	sizer_upper->Add(lbl_current_probe_id, 0, FLAGS, BORDER_PXL);
	sizer_upper->Add(this->lbl_current_probe_id_value, 0, FLAGS, BORDER_PXL);

	sizer_upper->Add(new wxStaticText( this, wxID_ANY, wxT("") ), 0, FLAGS, BORDER_PXL);
	sizer_upper->Add(this->btn_go, 0, FLAGS, BORDER_PXL);

	sizer_top->Add(sizer_upper, 0, FLAGS, BORDER_PXL);
	sizer_top->Add(this->grid_correlation, 0, FLAGS, BORDER_PXL);
	
    wxBoxSizer* sizer_bot = new wxBoxSizer(wxHORIZONTAL);

    sizer_bot->Add(this->btn_save, 0, FLAGS, BORDER_PXL);
    sizer_bot->Add(this->btn_close, 0, FLAGS_RIGHT, BORDER_PXL);
    sizer_top->Add(sizer_bot, 0, FLAGS, BORDER_PXL);
	SetSizer(sizer_top);
}


void CorrelationViewerDialog::OnClickSave( wxCommandEvent& evt){

	std::stringstream ss;
	if(this->investigation->is_mac)
		ss << this->investigation->dir_results << "/" << "output.txt";
	else
		ss << this->investigation->dir_results << "\\" << "output.txt";
    wxString wxlocation = wxGetTextFromUser(wxString::FromAscii("File name for correlations"), wxString::FromAscii("File name for correlations"), wxString::FromAscii(ss.str().c_str()), this);
	std::string location( wxlocation.ToAscii() );
    if(wxlocation.size()>0){
        bool success = this->grid_correlation->write_to_file(location);
        if( success ){
        	std::stringstream ssfn;
        	if(this->investigation->is_mac){
        		ssfn << "open -t " << location; // -t flag opens with default text editor
        		wxExecute(wxString::FromAscii(ssfn.str().c_str()));
        	}
        	else{
        		std::string loc_text = this->investigation->cp->get(std::string("External"), std::string("text_editor"));;
        		ssfn << loc_text << " " << location;
        		wxExecute(wxString::FromAscii(ssfn.str().c_str()));
        	}
        }
    }
}


void CorrelationViewerDialog::OnClickGo( wxCommandEvent& evt ){
	std::string gene_name( this->txt_gene_name->GetValue().ToAscii() );
	if( gene_name.size()==0 ){
		wxMessageBox(wxString::FromAscii("ERROR: You must enter a gene symbol or probe id"), _T("Not Ready"), wxOK | wxICON_EXCLAMATION, this);
		return;
	}
	double rho;
	try{
		std::string s( this->txt_min_abs_corr->GetValue().ToAscii() );
		rho = boost::lexical_cast<double>( s );
	}
	catch( ... ){
		wxMessageBox(wxString::FromAscii("ERROR: maximum absolute correlation must be a number"), _T("Not Ready"), wxOK | wxICON_EXCLAMATION, this);
		return;
	}
	
	if( rho < 0 || rho > 1 ){
		wxMessageBox(wxString::FromAscii("ERROR: maximum absolute correlation must be a number between 0 and 1"), _T("Not Ready"), wxOK | wxICON_EXCLAMATION, this);
		return;
	}
	this->grid_correlation->Update(gene_name, rho);
	this->lbl_current_probe_id_value->SetLabel( wxString::FromAscii( this->grid_correlation->current_probe_id.c_str() ));
}
