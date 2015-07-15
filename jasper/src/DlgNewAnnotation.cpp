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
#include "DlgNewAnnotation.h"

IMPLEMENT_CLASS( NewAnnotationDialog, wxDialog )

NewAnnotationDialog::NewAnnotationDialog(Investigation* investigation){
	this->investigation = investigation;
	SetExtraStyle(wxWS_EX_BLOCK_EVENTS|wxDIALOG_EX_CONTEXTHELP);
	wxDialog::Create( NULL, wxID_ANY, wxT("Create New Investigation"), wxDefaultPosition, wxDefaultSize, wxCAPTION|wxRESIZE_BORDER|wxSYSTEM_MENU );
	this->fn_annotation = std::string("");
	this->name_annotation = std::string("");
	CreateControls();
	GetSizer()->Fit(this);
	GetSizer()->SetSizeHints(this);
	Centre();
}


void NewAnnotationDialog::CreateControls()
{
	// CONTROLS
	wxStaticText* lbl_annotation_name= new wxStaticText( this, wxID_STATIC, wxT("Annotation Name:"), wxDefaultPosition, wxDefaultSize, 0 );
	this->txt_annotation_name = new wxTextCtrl ( this, ID_NEW_ANNOTATION_NAME, wxEmptyString, wxDefaultPosition, wxSize(200,22), 0 );
	wxStaticText* lbl_fn_annotation = new wxStaticText( this, wxID_STATIC, wxT("Annotation Name:"), wxDefaultPosition, wxDefaultSize, 0 );
	wxButton* btn_annotation = new wxButton( this, ID_NEW_ANNOTATION_SET, wxT("Browse"), wxDefaultPosition, wxDefaultSize, 0 );
	this->lbl_fn_annotation_value = new wxStaticText( this, wxID_STATIC, wxT("(not set)"), wxDefaultPosition, wxDefaultSize, 0 );

	wxSizer* sizer_ok = CreateButtonSizer(wxOK|wxCANCEL);

	Connect(wxID_OK, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(NewAnnotationDialog::OnClickOk));
	Connect(ID_NEW_ANNOTATION_SET, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(NewAnnotationDialog::OnClickSet));

	// POSITIONING
	this->sizer_top = new wxBoxSizer(wxVERTICAL);

	wxBoxSizer* sizer_name = new wxBoxSizer(wxHORIZONTAL);
	sizer_name->Add(lbl_annotation_name, 0, FLAGS, BORDER_PXL);
	sizer_name->Add(this->txt_annotation_name, 0, FLAGS, BORDER_PXL);
	this->sizer_top->Add(sizer_name, 0, FLAGS, BORDER_PXL);

	wxBoxSizer* sizer_fn = new wxBoxSizer(wxHORIZONTAL);
	sizer_fn->Add(lbl_fn_annotation, 0, FLAGS, BORDER_PXL);
	sizer_fn->Add(btn_annotation, 0, FLAGS, BORDER_PXL);
	sizer_fn->Add(lbl_fn_annotation_value, 0, FLAGS, BORDER_PXL);
	this->sizer_top->Add(sizer_fn, 0, FLAGS, BORDER_PXL);

	this->sizer_top->Add(sizer_ok, 0, FLAGS_RIGHT, BORDER_PXL);
	this->SetSizer(this->sizer_top);
}


void NewAnnotationDialog::OnClickSet( wxCommandEvent& event ){
	wxString cloc = wxFileSelector(wxString::FromAscii("Select file"), wxString::FromAscii(""), wxString::FromAscii(""), wxString::FromAscii(""), wxString::FromAscii(""), wxFD_OPEN, this);
	if( cloc.Len() > 0 ){
		this->lbl_fn_annotation_value->SetLabel(cloc);
		this->Fit();
	}
}


void NewAnnotationDialog::OnClickOk( wxCommandEvent& event ){

	// check to make sure name doesn't already exist
	this->name_annotation = this->txt_annotation_name->GetValue().ToAscii();
	this->fn_annotation = this->lbl_fn_annotation_value->GetLabel().ToAscii();

	if( this->name_annotation.size()==0 ){
		wxMessageBox(wxString::FromAscii("New annotation name must not be empty"), wxString::FromAscii("ERROR: no name provided"), wxOK | wxICON_EXCLAMATION, this);
		return;
	}

	typedef std::pair<std::string, std::string> prop_pair;
	std::vector< prop_pair* >* annotations = new std::vector<prop_pair*>();
	this->investigation->cp->get_section(std::string("Annotations"), annotations);
	bool exists = false;
	std::string fn_go;
	for(int i=0; i<(int)annotations->size(); i++){
		if( annotations->at(i)->first.compare(name_annotation) == 0 ){
			exists=true;
		}
		else if(annotations->at(i)->first.compare(std::string("go")) == 0 ){
			fn_go = annotations->at(i)->second;
		}
		delete annotations->at(i);
	}
	delete annotations;
	if( exists ){
		wxMessageBox(wxString::FromAscii("There is already an annotation with that name"), wxString::FromAscii("ERROR: There is already an annotation with that name"), wxOK | wxICON_EXCLAMATION, this);
		return;
	}

	GOAnnotationParser go(fn_go);
	try{
		GeneAnnotationParser gene_parser(fn_annotation, go);
	}
	catch( std::string err){
		std::stringstream ss;
		ss << "Error loading requested annotation:\n" << err;
		wxMessageBox(wxString::FromAscii(ss.str().c_str()), wxString::FromAscii("ERROR loading annotation file"), wxOK | wxICON_EXCLAMATION, this);
		return;
	}
	this->AcceptAndClose();
}


