#define _CRT_SECURE_NO_DEPRECATE 1
#define _SCL_SECURE_NO_DEPRECATE 1
#include <wx/dialog.h>
#include <wx/sizer.h>
#include <wx/stattext.h>
#include <wx/button.h>
#include <wx/textctrl.h>
#include <wx/radiobox.h>
#include <wx/wx.h>
#include <wx/process.h>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/random.hpp>
namespace alg = boost::algorithm;
#include <iostream>
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
#include "DlgAnnotate.h"

IMPLEMENT_CLASS( AnnotationDialog, wxDialog )

AnnotationDialog::AnnotationDialog( Investigation* investigation, std::string filename ){
	this->investigation = investigation;

	SetExtraStyle(wxWS_EX_BLOCK_EVENTS|wxDIALOG_EX_CONTEXTHELP);
	wxDialog::Create( NULL, wxID_ANY, wxT("Annotate gene list"), wxDefaultPosition, wxDefaultSize, wxCAPTION|wxRESIZE_BORDER|wxSYSTEM_MENU );
	CreateControls(filename);
	GetSizer()->Fit(this);
	GetSizer()->SetSizeHints(this);
	Centre();
}

void AnnotationDialog::CreateControls(std::string rs_filename)
{
	// CONTROLS
	//wxStaticText* lbl_file_name = new wxStaticText( this, wxID_ANY, wxT("File Name:") );
	//this->txt_file_name = new wxTextCtrl( this, ID_ANNOTATION_FILE_NAME, wxString::FromAscii("annotation.txt"), wxDefaultPosition, wxSize(200, 20) );
	//this->txt_file_name->ChangeValue(wxString::FromAscii("annotation.txt"));
	wxStaticText* lbl_annotation = new wxStaticText( this, wxID_ANY, wxT("Annotation:"));
	typedef std::pair<std::string, std::string> prop_pair;
	std::vector< prop_pair* >* annotations = new std::vector<prop_pair*>();
	this->investigation->cp->get_section(std::string("Annotations"), annotations);
	int species_idx = 0;
	for(int i=0; i<(int)annotations->size(); i++){
	   	if( annotations->at(i)->first.compare("go") != 0 ){
	   		this->ars_annotations.Add(wxString::FromAscii(annotations->at(i)->first.c_str()));
	   		this->fn_annotations.push_back( annotations->at(i)->second );
	   		if( annotations->at(i)->first.compare( this->investigation->species ) == 0 )
	   			species_idx = this->ars_annotations.Count()-1;
	   	}
	   	delete annotations->at(i);
	}

	this->cho_annotation = new wxChoice( this, ID_ANNOTATION_ANNOTATION, wxDefaultPosition, wxDefaultSize, this->ars_annotations);
	this->cho_annotation->SetSelection(species_idx);

    this->rdo_file = new wxRadioButton(this, ID_ANNOTATION_RDO_FILE, wxString::FromAscii("Genes in currently selected file"));
    this->rdo_paste = new wxRadioButton(this, ID_ANNOTATION_RDO_PASTE, wxString::FromAscii("Genes pasted into box below"));
    this->rs_filename = rs_filename;
    if( rs_filename.size()==0 ){
        this->rdo_file->Enable(false);
    }
    this->rdo_paste->SetValue(true);
    this->txt_probes = new wxTextCtrl( this, ID_ANNOTATION_PROBES, wxT(""), wxDefaultPosition, wxSize(150, 150), wxTE_MULTILINE );


	wxSizer* sizer_ok = CreateButtonSizer(wxOK|wxCANCEL);
	Connect(wxID_OK, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(AnnotationDialog::OnClickOk));
	Connect(wxID_CANCEL, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(AnnotationDialog::OnClickCancel));
	Connect(ID_ANNOTATION_RDO_FILE, wxEVT_COMMAND_RADIOBUTTON_SELECTED, wxCommandEventHandler(AnnotationDialog::OnChangeAttribute));
	Connect(ID_ANNOTATION_RDO_PASTE, wxEVT_COMMAND_RADIOBUTTON_SELECTED, wxCommandEventHandler(AnnotationDialog::OnChangeAttribute));
	SetEscapeId(wxID_CANCEL);

    wxBoxSizer* sizer_top = new wxBoxSizer(wxVERTICAL);

	wxFlexGridSizer* sizer_upper = new wxFlexGridSizer(1, 2, 0, 0);
	//sizer_upper->Add(lbl_file_name, 0, FLAGS, BORDER_PXL);
	//sizer_upper->Add(this->txt_file_name, 0, FLAGS, BORDER_PXL);
	sizer_upper->Add(lbl_annotation , 0, FLAGS, BORDER_PXL);
	sizer_upper->Add(this->cho_annotation, 0, FLAGS, BORDER_PXL);
    sizer_top->Add(sizer_upper, 0, FLAGS, BORDER_PXL);
    sizer_top->Add(this->rdo_file, 0, FLAGS, BORDER_PXL);
    sizer_top->Add(this->rdo_paste, 0, FLAGS, BORDER_PXL);
	sizer_top->Add(this->txt_probes, 0, FLAGS, BORDER_PXL);
    sizer_top->Add(sizer_ok, 0, FLAGS_RIGHT, BORDER_PXL);
    if(this->investigation->is_mac){
    	// mac uses giant fonts by default; PC is better.
    	wxFont font(this->GetFont());
    	font.SetPointSize(wxSMALL_FONT->GetPointSize());
    	this->SetFont(font);
        this->cho_annotation->SetFont(font);
        this->rdo_file->SetFont( font );
        this->rdo_paste->SetFont( font );
        //lbl_file_name->SetFont( font );
        //this->txt_file_name->SetFont( font );
        //lbl_annotation->SetFont(font);
   	}
	this->SetSizer(sizer_top);
}

void AnnotationDialog::OnChangeAttribute(wxCommandEvent& WXUNUSED(event)){
	if( this->rdo_paste->GetValue() ){
		this->txt_probes->Enable(true);
	}
	else{
		this->txt_probes->Enable(false);
		this->txt_probes->Clear();
	}
}

void AnnotationDialog::OnClickCancel( wxCommandEvent& WXUNUSED(event) ){
	EndModal(wxID_CANCEL);
}

void AnnotationDialog::OnClickOk( wxCommandEvent& WXUNUSED(event) ){
	//std::string o( this->txt_file_name->GetValue().ToAscii() );
	wxString default_location = wxString::FromAscii(this->investigation->dir_results.c_str());
	wxString default_name = wxString::FromAscii("annotation.txt");
	wxString default_extension = wxString::FromAscii("txt");
	wxString filename = wxFileSelector(wxString::FromAscii("Save To"), default_location, default_name, default_extension, wxString::FromAscii("*.*"), wxFD_SAVE, this);
	if( filename.size()==0 ){
		return;
	}
	std::string fn_out(filename.ToAscii());
	std::ofstream f_out(fn_out.c_str());
	if( !f_out.is_open() ){
		wxMessageBox(wxString::FromAscii("Unable to open file for writing"));
		return;
	}
	GOAnnotationParser gop;
	GeneAnnotationParser genep;
	std::vector<std::string> genes;
	std::string f;
	try{
		gop.load(this->investigation->cp->get(std::string("Annotations"), std::string("go")));
		genep.load(this->fn_annotations.at( this->cho_annotation->GetSelection() ), gop);
		if(this->rdo_file->GetValue()){
			f = this->rs_filename;
			CarmenParser cp;
			cp.Load(f);
			cp.PrepareToReadValues();
			cp.ExtractGenes(genes);
		}
		else if(this->rdo_paste->GetValue() ){
			wxString items = this->txt_probes->GetValue();
		    std::string genes_str( items.ToAscii() );
		    std::vector<std::string> genes_raw;
		    alg::split(genes_raw, genes_str, alg::is_any_of("\n") );
		    for(int i=0; i<(int)genes_raw.size(); i++){
		    	if( genes_raw.at(i).compare( std::string("") )  != 0 ){
		    		genes.push_back(genes_raw.at(i));
		    	}
		    }
		}
	}
	catch(std::string err){
		f_out.close();
		wxMessageBox(wxString::FromAscii(err.c_str()));
		return;
	}
	GeneAnnotation* gene=NULL;
	GOAnnotation* go=NULL;
	std::vector<std::string> genes_not_found;
	f_out << "Symbol\tName\tChr\tloc_start\tGO_BP\tGO_MF\tGO_CC\n";
	std::vector< std::string > GO_BP, GO_MF, GO_CC;
	std::string join_str;
	for(int i=0; i<(int)genes.size(); i++){
		try{
			gene=genep.get_annotation(genes.at(i));
			GO_MF.clear();
			GO_BP.clear();
			GO_CC.clear();
			f_out << gene->symbol << "\t" << gene->full_name << "\t" << gene->chr << "\t" << gene->loc_start;
			for(int j=0; j<(int)gene->GO_annotations.size(); j++){
				go = gop.get_annotation( gene->GO_annotations.at(j) );
				if(go->branch == GOAnnotation::GO_MF )
					GO_MF.push_back(go->description);
				else if(go->branch == GOAnnotation::GO_BP )
					GO_BP.push_back(go->description);
				else
					GO_CC.push_back(go->description);
			}

			if(GO_BP.size()==0)
				f_out << "\tNA";
			else{
				f_out << "\t";
				stringhasher::join(join_str, &GO_BP, std::string(","));
				f_out << join_str;
			}
			if(GO_MF.size()==0)
				f_out << "\tNA";
			else{
				f_out << "\t";
				stringhasher::join(join_str, &GO_MF, std::string(","));
				f_out << join_str;
			}
			if(GO_CC.size()==0)
				f_out << "\tNA";
			else{
				f_out << "\t";
				stringhasher::join(join_str, &GO_CC, std::string(","));
				f_out << join_str;
			}

			f_out << "\n";
		}
		catch(std::string err){
			f_out << genes.at(i) << "\tNA\tNA\tNA\tNA\tNA\tNA\n";
			genes_not_found.push_back(genes.at(i));
		}
	}
	f_out.close();
	if(genes_not_found.size()>0){
		std::stringstream ss;
		ss << "Unable to identify the following genes:\n";
		for(int i=0; i<(int)genes_not_found.size(); i++){
			ss << " " << genes_not_found.at(i);
			if(i%5==0 && i>0)
				ss << "\n";
		}
		wxMessageBox(wxString::FromAscii(ss.str().c_str()));
	}
	if(genes.size()>genes_not_found.size()){
		this->new_filename = fn_out;
		this->AcceptAndClose();
	}
}
