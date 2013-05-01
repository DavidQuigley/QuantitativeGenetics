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
#include "DlgNewAnnotation.h"
#include "DlgProperties.h"

IMPLEMENT_CLASS( PropertyDialog, wxDialog )

PropertyDialog::PropertyDialog( Investigation* investigation ){
	this->investigation = investigation;
	SetExtraStyle(wxWS_EX_BLOCK_EVENTS|wxDIALOG_EX_CONTEXTHELP);
	wxDialog::Create( NULL, wxID_ANY, wxT("Edit Carmen Properties"), wxDefaultPosition, wxDefaultSize, wxCAPTION|wxRESIZE_BORDER|wxSYSTEM_MENU );
	CreateControls();
	GetSizer()->Fit(this);
	GetSizer()->SetSizeHints(this);
	Centre();
}


void PropertyDialog::CreateControls()
{
	wxStaticBoxSizer* static_applications_sizer = new wxStaticBoxSizer(wxVERTICAL, this, wxString::FromAscii("Application Locations") );
	wxStaticBoxSizer* static_annotations_sizer = new wxStaticBoxSizer(wxVERTICAL, this, wxString::FromAscii("Annotation Properties") );

	std::string cytoscape( this->investigation->cp->get(std::string("External"), std::string("cytoscape") ) );
	std::string excel( this->investigation->cp->get(std::string("External"), std::string("excel") ) );
	std::string vizmap( this->investigation->cp->get(std::string("Internal"), std::string("cytoscape_vizmap") ) );
	std::string text_editor( this->investigation->cp->get(std::string("External"), std::string("text_editor") ) );
	if( cytoscape.size()==0 )
		cytoscape = std::string("(not set)");
	if( excel.size()==0 )
		excel = std::string("(not set)");
	if( vizmap.size()==0 )
		vizmap = std::string("(not set)");
	if( text_editor.size()==0 )
		text_editor = std::string("(not set)");

	wxStaticText* lbl_cytoscape = new wxStaticText( this, wxID_STATIC, wxT("Cytoscape:"), wxDefaultPosition, wxDefaultSize, 0 );
	this->lbl_cytoscape_value = new wxStaticText( this, wxID_STATIC, wxString::FromAscii(cytoscape.c_str()), wxDefaultPosition, wxDefaultSize, 0 );
	wxButton* btn_cytoscape = new wxButton( this, ID_PROPERTIES_CYTOSCAPE, wxT("Browse"), wxDefaultPosition, wxDefaultSize, 0 );

	wxStaticText* lbl_vizmap = new wxStaticText( this, wxID_STATIC, wxT("Cytoscape vizmap:"), wxDefaultPosition, wxDefaultSize, 0 );
	this->lbl_vizmap_value = new wxStaticText( this, wxID_STATIC, wxString::FromAscii(vizmap.c_str()), wxDefaultPosition, wxDefaultSize, 0 );
	wxButton* btn_vizmap = new wxButton( this, ID_PROPERTIES_VIZMAP, wxT("Browse"), wxDefaultPosition, wxDefaultSize, 0 );

	wxStaticText* lbl_excel = new wxStaticText( this, wxID_STATIC, wxT("Spreadsheet:"), wxDefaultPosition, wxDefaultSize, 0 );
	this->lbl_excel_value = new wxStaticText( this, wxID_STATIC, wxString::FromAscii(excel.c_str()), wxDefaultPosition, wxDefaultSize, 0 );
	wxButton* btn_excel = new wxButton( this, ID_PROPERTIES_EXCEL, wxT("Browse"), wxDefaultPosition, wxDefaultSize, 0 );

	wxStaticText* lbl_text_editor = new wxStaticText( this, wxID_STATIC, wxT("Text editor:"), wxDefaultPosition, wxDefaultSize, 0 );
	this->lbl_text_editor_value = new wxStaticText( this, wxID_STATIC, wxString::FromAscii(text_editor.c_str()), wxDefaultPosition, wxDefaultSize, 0 );
	wxButton* btn_text_editor = new wxButton( this, ID_PROPERTIES_TEXT_EDITOR, wxT("Browse"), wxDefaultPosition, wxDefaultSize, 0 );

	std::string go = this->investigation->cp->get(std::string("Annotations"), std::string("go") );
	wxStaticText* lbl_go = new wxStaticText( this, wxID_STATIC, wxT("GO Annotation file:"), wxDefaultPosition, wxDefaultSize, 0 );
	this->lbl_go_value = new wxStaticText( this, wxID_STATIC, wxString::FromAscii(go.c_str()), wxDefaultPosition, wxDefaultSize, 0 );
	wxButton* btn_go = new wxButton( this, ID_PROPERTIES_GO, wxT("Browse"), wxDefaultPosition, wxDefaultSize, 0 );

	wxStaticText* lbl_annotation = new wxStaticText( this, wxID_ANY, wxT("Annotations"));
	typedef std::pair<std::string, std::string> prop_pair;
	std::vector< prop_pair* >* annotations = new std::vector<prop_pair*>();
	this->investigation->cp->get_section(std::string("Annotations"), annotations);
	for(int i=0; i<(int)annotations->size(); i++){
		if( annotations->at(i)->first.compare("go") != 0 ){
			ars_annot.Add(wxString::FromAscii(annotations->at(i)->first.c_str()));
			this->fn_annotations.push_back( annotations->at(i)->second );
			delete annotations->at(i);
		}
	}
	this->cho_annotation = new wxChoice( this, ID_PROPERTIES_ANNOTATION_CHO, wxDefaultPosition, wxDefaultSize, ars_annot );
	std::string fn_annotation("");
	if(this->fn_annotations.size()>0){
		this->cho_annotation->SetSelection(0);
		fn_annotation = this->fn_annotations.at(0);
	}
	this->lbl_annotation_value = new wxStaticText( this, wxID_STATIC, wxString::FromAscii(fn_annotation.c_str()), wxDefaultPosition, wxDefaultSize, 0 );
	wxButton* btn_annotation = new wxButton( this, ID_PROPERTIES_ANNOTATION_BTN, wxT("New Annotation"), wxDefaultPosition, wxDefaultSize, 0 );

	wxSizer* sizer_ok = CreateButtonSizer(wxOK|wxCANCEL);

	if( this->investigation->is_mac ){
		lbl_text_editor->Show(false);
		lbl_text_editor_value->Show(false);
		btn_text_editor->Show(false);
		// mac uses giant fonts by default; PC is better.
		wxFont font(btn_cytoscape->GetFont());
		font.SetPointSize(wxSMALL_FONT->GetPointSize());
		lbl_cytoscape->SetFont(font);
		lbl_vizmap->SetFont(font);
		lbl_excel->SetFont(font);
		lbl_text_editor->SetFont(font);
		lbl_go->SetFont(font);
		lbl_annotation->SetFont(font);
		btn_annotation->SetFont(font);
		btn_go->SetFont(font);
		btn_cytoscape->SetFont(font);
		btn_vizmap->SetFont(font);
		btn_excel->SetFont(font);
		lbl_cytoscape_value->SetFont(font);
		lbl_vizmap_value->SetFont(font);
		lbl_excel_value->SetFont(font);
		this->lbl_go_value->SetFont(font);
		this->lbl_annotation_value->SetFont(font);
	}

	Connect(ID_PROPERTIES_CYTOSCAPE, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(PropertyDialog::OnClickCytoscape));
	Connect(ID_PROPERTIES_VIZMAP, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(PropertyDialog::OnClickVizmap));
	Connect(ID_PROPERTIES_EXCEL, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(PropertyDialog::OnClickExcel));
	Connect(ID_PROPERTIES_TEXT_EDITOR, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(PropertyDialog::OnClickTextEditor));
	Connect(ID_PROPERTIES_GO, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(PropertyDialog::OnClickGO));
	Connect(ID_PROPERTIES_ANNOTATION_CHO, wxEVT_COMMAND_CHOICE_SELECTED, wxCommandEventHandler(PropertyDialog::OnChangeAnnotation));
	Connect(ID_PROPERTIES_ANNOTATION_BTN, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(PropertyDialog::OnNewAnnotation));

	// POSITIONING
	wxBoxSizer* sizer_top = new wxBoxSizer(wxVERTICAL);

	wxFlexGridSizer* sizer_app = new wxFlexGridSizer(4, 3, 0, 0);
	sizer_app->Add(lbl_cytoscape, 0, FLAGS, BORDER_PXL);
	sizer_app->Add(btn_cytoscape, 0, FLAGS, BORDER_PXL);
	sizer_app->Add(lbl_cytoscape_value, 0, FLAGS, BORDER_PXL);

	sizer_app->Add(lbl_vizmap, 0, FLAGS, BORDER_PXL);
	sizer_app->Add(btn_vizmap, 0, FLAGS, BORDER_PXL);
	sizer_app->Add(lbl_vizmap_value, 0, FLAGS, BORDER_PXL);

	sizer_app->Add(lbl_excel, 0, FLAGS, BORDER_PXL);
	sizer_app->Add(btn_excel, 0, FLAGS, BORDER_PXL);
	sizer_app->Add(lbl_excel_value, 0, FLAGS, BORDER_PXL);

	sizer_app->Add(lbl_text_editor, 0, FLAGS, BORDER_PXL);
	sizer_app->Add(btn_text_editor, 0, FLAGS, BORDER_PXL);
	sizer_app->Add(lbl_text_editor_value, 0, FLAGS, BORDER_PXL);

	static_applications_sizer->Add(sizer_app, 0, FLAGS, BORDER_PXL);
	sizer_top->Add(static_applications_sizer, 0, FLAGS, BORDER_PXL);

	wxFlexGridSizer* sizer_annot = new wxFlexGridSizer(2, 3, 0, 0);
	sizer_annot->Add(lbl_go, 0, FLAGS, BORDER_PXL);
	sizer_annot->Add(btn_go, 0, FLAGS, BORDER_PXL);
	sizer_annot->Add(lbl_go_value, 0, FLAGS, BORDER_PXL);
	sizer_annot->Add(lbl_annotation, 0, FLAGS, BORDER_PXL);
	sizer_annot->Add(cho_annotation, 0, FLAGS, BORDER_PXL);
	sizer_annot->Add(lbl_annotation_value, 0, FLAGS, BORDER_PXL);
	static_annotations_sizer->Add(sizer_annot,0,FLAGS, BORDER_PXL);
	static_annotations_sizer->Add(btn_annotation,0,FLAGS, BORDER_PXL);
	sizer_top->Add(static_annotations_sizer, 0, FLAGS, BORDER_PXL);

	sizer_top->Add(sizer_ok, 0, FLAGS_RIGHT, BORDER_PXL);

	this->SetSizer(sizer_top);
}

void PropertyDialog::OnClickCytoscape(wxCommandEvent& event ){
	std::string app("app");
	std::string sh("sh");
	std::string loc_cytoscape;
	wxString cloc = wxFileSelector(wxString::FromAscii("Select file"), wxString::FromAscii(""), wxString::FromAscii(""), wxString::FromAscii(""), wxString::FromAscii("*.sh|*.app|*.exe"), wxFD_OPEN, this);
	if( cloc.Len() > 0 ){
		loc_cytoscape = std::string(cloc.ToAscii());
		alg::replace_all(loc_cytoscape, app, sh);
		this->investigation->cp->set(std::string("External"), std::string("cytoscape"), loc_cytoscape);
		this->investigation->cp->write();
		this->lbl_cytoscape_value->SetLabel( cloc );
	}
}

void PropertyDialog::OnClickVizmap(wxCommandEvent& event ){
	wxString cloc = wxFileSelector(wxString::FromAscii("Select file"), wxString::FromAscii(""), wxString::FromAscii(""), wxString::FromAscii(""), wxString::FromAscii("*.props"), wxFD_OPEN, this);
	if( cloc.Len() > 0 ){
		this->investigation->cp->set(std::string("Internal"), std::string("cytoscape_vizmap"), std::string(cloc.ToAscii()));
		this->investigation->cp->write();
		this->lbl_vizmap_value->SetLabel( cloc );
	}
}

void PropertyDialog::OnClickExcel(wxCommandEvent& event ){
	wxString cloc;
	if( this->investigation->is_mac)
		cloc = wxFileSelector(wxString::FromAscii("Select file"), wxString::FromAscii(""), wxString::FromAscii(""), wxString::FromAscii(""), wxString::FromAscii("App files (*.app)|*.app"), wxFD_OPEN, this);
	else
		cloc = wxFileSelector(wxString::FromAscii("Select file"), wxString::FromAscii(""), wxString::FromAscii(""), wxString::FromAscii(""), wxString::FromAscii("Exe files (*.exe)|*.exe"), wxFD_OPEN, this);
	if( cloc.Len() > 0 ){
		this->investigation->cp->set(std::string("External"), std::string("excel"), std::string(cloc.ToAscii()));
		this->investigation->cp->write();
		this->lbl_excel_value->SetLabel( cloc );
	}
}

void PropertyDialog::OnClickTextEditor(wxCommandEvent& event ){
	wxString cloc = wxFileSelector(wxString::FromAscii("Select file"), wxString::FromAscii(""), wxString::FromAscii(""), wxString::FromAscii(""), wxString::FromAscii("*.exe"), wxFD_OPEN, this);
	if( cloc.Len() > 0 ){
		this->investigation->cp->set(std::string("External"), std::string("text_editor"), std::string(cloc.ToAscii()));
		this->investigation->cp->write();
		this->lbl_text_editor_value->SetLabel( cloc );
	}
}


void PropertyDialog::OnClickGO(wxCommandEvent& event ){
	wxString cloc = wxFileSelector(wxString::FromAscii("Select file"), wxString::FromAscii(""), wxString::FromAscii(""), wxString::FromAscii(""), wxString::FromAscii("*.txt"), wxFD_OPEN, this);
	if( cloc.Len() > 0 ){
		this->investigation->cp->set(std::string("Annotations"), std::string("go"), std::string(cloc.ToAscii()));
		this->investigation->cp->write();
		this->lbl_go_value->SetLabel( cloc );
	}
}


void PropertyDialog::OnChangeAnnotation( wxCommandEvent& event ){
	int cur_sel = this->cho_annotation->GetCurrentSelection();
	this->lbl_annotation_value->SetLabel( wxString::FromAscii(this->fn_annotations.at(cur_sel).c_str() ) );
	this->Fit();
}


void PropertyDialog::OnNewAnnotation(wxCommandEvent& event ){
	NewAnnotationDialog dialog( this->investigation );
	if( dialog.ShowModal() == wxID_OK){
		this->new_annotation_files.push_back( dialog.fn_annotation );
		this->new_annotation_names.push_back( dialog.name_annotation );
		this->ars_annot.Add( wxString::FromAscii(dialog.name_annotation.c_str()) );
		this->fn_annotations.push_back(dialog.fn_annotation);
		this->cho_annotation->SetSelection( this->cho_annotation->GetCount()-1 );
		this->Fit();
	}
}
