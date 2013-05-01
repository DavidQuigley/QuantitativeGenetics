#define _SCL_SECURE_NO_DEPRECATE 1
#include <wx/sizer.h>
#include <wx/stattext.h>
#include <wx/button.h>
#include <wx/textctrl.h>
#include <wx/panel.h>
#include <wx/listbox.h>
#include <wx/dialog.h>
#include <wx/choice.h>
#include <wx/dynarray.h>
#include <wx/wx.h>
#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/random.hpp>
#include <boost/filesystem/path.hpp>
#include "carmen/src/DataStructures.h"
#include "carmen/src/Attributes.h"
#include "carmen/src/ClassMinerOptions.h"
#include "carmen/src/Rawdata.h"
#include "carmen/src/Discretize.h"
#include "carmen/src/Dataset.h"
#include "carmen/src/Parser.h"
#include "Investigation.h"
#include "PanelInvestigationProperties.h"

InvestigationPropertiesPanel::InvestigationPropertiesPanel(wxWindow* parent, int id, Investigation* investigation) : wxPanel(parent, id){
	this->parent = parent;
	this->investigation = investigation;
	this->most_recent_folder = std::string("");
	wxStaticText* lbl_raw_data = new wxStaticText( this, wxID_STATIC, wxT("Raw Data:"), wxDefaultPosition, wxDefaultSize, 0 );
	this->lbl_expr_value = new wxStaticText( this, wxID_STATIC, wxT("(not set)"), wxDefaultPosition, wxDefaultSize, 0 );
	this->btn_raw_data = new wxButton( this, ID_INVEST_PROPERTIES_EXPR, wxT("Browse"), wxDefaultPosition, wxDefaultSize, 0 );

	wxStaticText* lbl_sa = new wxStaticText( this, wxID_STATIC, wxT("Sample Attributes:"), wxDefaultPosition, wxDefaultSize, 0 );
	this->lbl_sa_value = new wxStaticText( this, wxID_STATIC, wxT("(not set)"), wxDefaultPosition, wxDefaultSize, 0 );
	this->btn_sa = new wxButton( this, ID_INVEST_PROPERTIES_SA, wxT("Browse"), wxDefaultPosition, wxDefaultSize, 0 );

	wxStaticText* lbl_ga = new wxStaticText( this, wxID_STATIC, wxT("Probe Attributes:"), wxDefaultPosition, wxDefaultSize, 0 );
	this->lbl_ga_value = new wxStaticText( this, wxID_STATIC, wxT("(not set)"), wxDefaultPosition, wxDefaultSize, 0 );
	this->btn_ga = new wxButton( this, ID_INVEST_PROPERTIES_GA, wxT("Browse"), wxDefaultPosition, wxDefaultSize, 0 );

	wxStaticText* lbl_results = new wxStaticText( this, wxID_STATIC, wxT("Results Directory:"), wxDefaultPosition, wxDefaultSize, 0 );
	this->lbl_results_value = new wxStaticText( this, wxID_STATIC, wxT("(not set)"), wxDefaultPosition, wxDefaultSize, 0 );
	this->btn_results = new wxButton( this, ID_INVEST_PROPERTIES_RESULTS, wxT("Browse"), wxDefaultPosition, wxDefaultSize, 0 );

	wxStaticText* lbl_annotation = new wxStaticText( this, wxID_ANY, wxT("Species:"));
	typedef std::pair<std::string, std::string> prop_pair;
	std::vector< prop_pair* >* annotations = new std::vector<prop_pair*>();
	this->investigation->cp->get_section(std::string("Annotations"), annotations);
	for(int i=0; i<(int)annotations->size(); i++){
	   	if( annotations->at(i)->first.compare("go") != 0 ){
	   		this->ars_annotations.Add(wxString::FromAscii(annotations->at(i)->first.c_str()));
	   	}
	   	delete annotations->at(i);
	}
	this->cho_annotation = new wxChoice( this, ID_INVEST_PROPERTIES_SPECIES, wxDefaultPosition, wxDefaultSize, this->ars_annotations);
	this->cho_annotation->SetSelection(0);

	this->ars_data_types.Add(wxString::FromAscii("Gene Expression"));
	this->ars_data_types.Add(wxString::FromAscii("Genotype"));
	this->ars_data_types.Add(wxString::FromAscii("DNA Copy Number"));
	wxStaticText* lbl_data_type = new wxStaticText( this, wxID_ANY, wxT("Data Type:"));
	this->cho_data_type= new wxChoice( this, ID_INVEST_PROPERTIES_DATA_TYPE, wxDefaultPosition, wxDefaultSize, this->ars_data_types);
	this->cho_data_type->SetSelection(0);

	this->lbl_gene_name_column = new wxStaticText( this, wxID_STATIC, wxT("Gene Symbol Column:"), wxDefaultPosition, wxDefaultSize, 0 );
	this->cho_gene_name_column = new wxChoice( this, ID_INVEST_PROPERTIES_GENENAME, wxDefaultPosition, wxDefaultSize);

	this->gene_name_column.clear();
	this->cho_gene_name_column->Enable(false);

	if(investigation->is_mac){
	    // mac uses giant fonts by default; PC is better.
		wxFont font(this->btn_raw_data->GetFont());
		font.SetPointSize(wxSMALL_FONT->GetPointSize());
	    this->btn_raw_data->SetFont(font);
	    lbl_raw_data->SetFont(font);
	    lbl_sa->SetFont(font);
	    lbl_ga->SetFont(font);
	    lbl_data_type->SetFont(font);
	    lbl_results->SetFont(font);
	    lbl_annotation->SetFont(font);
	    cho_annotation->SetFont(font);
	    cho_data_type->SetFont(font);
	    cho_gene_name_column->SetFont(font);
	    this->btn_sa->SetFont(font);
	    this->btn_ga->SetFont(font);
	    this->btn_results->SetFont(font);
	    this->lbl_expr_value->SetFont(font);
	    this->lbl_sa_value->SetFont(font);
	    this->lbl_ga_value->SetFont(font);
	    this->lbl_results_value->SetFont(font);
	    this->lbl_gene_name_column->SetFont(font);
	}

	Connect(ID_INVEST_PROPERTIES_EXPR, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(InvestigationPropertiesPanel::OnClickRawData));
	Connect(ID_INVEST_PROPERTIES_SA, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(InvestigationPropertiesPanel::OnClickSA));
	Connect(ID_INVEST_PROPERTIES_GA, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(InvestigationPropertiesPanel::OnClickGA));
	Connect(ID_INVEST_PROPERTIES_RESULTS, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(InvestigationPropertiesPanel::OnClickResults));

	this->sizer_top =  new wxFlexGridSizer(7, 3, 0, 0);

	sizer_top->Add(lbl_raw_data, 0, FLAGS, BORDER_PXL);
	sizer_top->Add(this->btn_raw_data, 0, FLAGS, BORDER_PXL);
	sizer_top->Add(this->lbl_expr_value, 0, FLAGS, BORDER_PXL);

	sizer_top->Add(lbl_sa, 0, FLAGS, BORDER_PXL);
	sizer_top->Add(this->btn_sa, 0, FLAGS, BORDER_PXL);
	sizer_top->Add(this->lbl_sa_value, 0, FLAGS, BORDER_PXL);

	sizer_top->Add(lbl_ga, 0, FLAGS, BORDER_PXL);
	sizer_top->Add(this->btn_ga, 0, FLAGS, BORDER_PXL);
	sizer_top->Add(this->lbl_ga_value, 0, FLAGS, BORDER_PXL);

	sizer_top->Add(lbl_results, 0, FLAGS, BORDER_PXL);
	sizer_top->Add(this->btn_results, 0, FLAGS, BORDER_PXL);
	sizer_top->Add(this->lbl_results_value, 0, FLAGS, BORDER_PXL);

	sizer_top->Add(lbl_annotation, 0, FLAGS, BORDER_PXL);
	sizer_top->Add(cho_annotation, 0, FLAGS, BORDER_PXL);
	sizer_top->Add(new wxStaticText( this, wxID_STATIC, wxT("") ), 0, FLAGS, BORDER_PXL);

	sizer_top->Add(lbl_data_type, 0, FLAGS, BORDER_PXL);
	sizer_top->Add(cho_data_type, 0, FLAGS, BORDER_PXL);
	sizer_top->Add(new wxStaticText( this, wxID_STATIC, wxT("") ), 0, FLAGS, BORDER_PXL);

	sizer_top->Add(lbl_gene_name_column, 0, FLAGS, BORDER_PXL);
	sizer_top->Add(cho_gene_name_column, 0, FLAGS, BORDER_PXL);
	sizer_top->Add(new wxStaticText( this, wxID_STATIC, wxT("") ), 0, FLAGS, BORDER_PXL);

	this->SetSizer(this->sizer_top);
	this->sizer_top->Fit(this);
}

std::string InvestigationPropertiesPanel::get_fn_expr(){
	if( strcmp( this->lbl_expr_value->GetLabel().ToAscii(), "(not set)" )==0 )
		return std::string("");
	return std::string( this->lbl_expr_value->GetLabel().ToAscii() );
}

std::string InvestigationPropertiesPanel::get_fn_ga(){
	if( strcmp( this->lbl_ga_value->GetLabel().ToAscii(), "(not set)" )==0 )
		return std::string("");
	return std::string( this->lbl_ga_value->GetLabel().ToAscii() );
}

std::string InvestigationPropertiesPanel::get_fn_sa(){
	if( strcmp( this->lbl_sa_value->GetLabel().ToAscii(), "(not set)" )==0 )
		return std::string("");
	return std::string( this->lbl_sa_value->GetLabel().ToAscii() );
}

std::string InvestigationPropertiesPanel::get_dir_results(){
	if( strcmp( this->lbl_results_value->GetLabel().ToAscii(), "(not set)" )==0 )
		return std::string("");
	return std::string( this->lbl_results_value->GetLabel().ToAscii() );
}

std::string InvestigationPropertiesPanel::get_species(){
	return std::string( this->ars_annotations.Item( this->cho_annotation->GetSelection()).ToAscii() );
}

std::string InvestigationPropertiesPanel::get_data_type(){
	return std::string( this->ars_data_types.Item( this->cho_data_type->GetSelection()).ToAscii() );
}

std::string InvestigationPropertiesPanel::get_gene_name_column(){
	return this->gene_name_column.at( this->cho_gene_name_column->GetSelection() );
}

void InvestigationPropertiesPanel::set_fn_expr(std::string fn){
	this->lbl_expr_value->SetLabel(wxString::FromAscii(fn.c_str() ) );
	this->redraw();
}

void InvestigationPropertiesPanel::set_fn_ga(std::string fn){
	try{
		Attributes* ga = new Attributes(std::string("NA"));
		ga->load(fn);
		int idx=0;
		this->gene_name_column.clear();
		this->cho_gene_name_column->Clear();
		for(int i=0; i<int(ga->attrib.size() ); i++){
			this->gene_name_column.push_back(ga->attrib.at(i));
			this->cho_gene_name_column->Append(wxString::FromAscii(ga->attrib.at(i).c_str()));
			if( ga->attrib.at(i).compare("Gene Name")==0 ) // for legacy reasons
				idx = i;
			else if( ga->attrib.at(i).compare("Symbol")==0 )
				idx = i;
			else if( ga->attrib.at(i).compare("symbol")==0 )
				idx = i;
		}
		this->cho_gene_name_column->SetSelection(idx);
		this->cho_gene_name_column->Enable(true);
		this->lbl_ga_value->SetLabel(wxString::FromAscii(fn.c_str() ) );
		this->redraw();
	}
	catch( std::string err){
		std::stringstream ss;
		ss << "ERROR loading attributes:\n" << err;
		wxMessageBox(wxString::FromAscii(ss.str().c_str()), wxString::FromAscii("ERROR loading file"), wxOK | wxICON_EXCLAMATION, this);
		return;
	}
}

void InvestigationPropertiesPanel::set_fn_sa(std::string fn){
	this->lbl_sa_value->SetLabel(wxString::FromAscii(fn.c_str() ) ) ;
	this->redraw();
}

void InvestigationPropertiesPanel::set_dir_results(std::string fn){
	this->lbl_results_value->SetLabel(wxString::FromAscii(fn.c_str() ) );
	this->redraw();
}

void InvestigationPropertiesPanel::set_gene_name_column(std::string column_name){
	int idx=-1;
	for(int i=0; i < int(this->gene_name_column.size()); i++){
		if( this->gene_name_column.at(i).compare( column_name ) == 0 )
			idx = i;
	}
	if(idx>-1)
		this->cho_gene_name_column->SetSelection(idx);
	else
		this->cho_gene_name_column->SetSelection(0);
}


void InvestigationPropertiesPanel::set_species(std::string fn){
	for(int i=0; i<int(this->cho_annotation->GetCount()); i++){
		if( this->cho_annotation->GetString(i).Cmp(wxString::FromAscii(fn.c_str()))==0 )
			this->cho_annotation->SetSelection(i);
	}
}

void InvestigationPropertiesPanel::set_data_type(std::string fn){
	for(int i=0; i<int(this->cho_data_type->GetCount()); i++){
		if( this->cho_data_type->GetString(i).Cmp(wxString::FromAscii(fn.c_str()))==0 )
			this->cho_data_type->SetSelection(i);
	}
}

bool InvestigationPropertiesPanel::CheckConsistency(std::string& err){
	std::string fn_raw_data( this->get_fn_expr() );
	std::string fn_ga( this->get_fn_ga() );
	std::string fn_sa(this->get_fn_sa() );
	std::string fn_results( this->get_dir_results() );
	if(fn_raw_data.size()==0){
		err = std::string("Raw data file not set.");
		return false;
	}
	if(fn_ga.size()==0){
		err = std::string("Gene Attribute file not set.");
		return false;
	}
	if(fn_sa.size()==0){
		err = std::string("Sample Attribute file not set.");
		return false;
	}
	if( fn_results.size()==0 ){
		err = std::string("Results directory not set.");
		return false;
	}
	Attributes ga("NA");
	Attributes sa("NA");
	// First check: do three data files load?
	try{ ga.load( fn_ga ); }
	catch( std::string err_msg ){
		err = err_msg;
		return false;
	}
	try{ sa.load( fn_sa ); }
	catch( std::string err_msg ){
		err = err_msg;
		return false;
	}
	// Next check: do identifiers match each other?
	Dataset* dataset = new Dataset();
	try{ dataset->load(&sa, &ga, fn_raw_data, std::string("IDENTIFIER!NULL"), std::string("")); }
	catch(std::string err_msg){
	  	err = err_msg;
	  	delete dataset;
	   	return false;
	}
	delete dataset;
	return true;
}


void InvestigationPropertiesPanel::OnClickRawData( wxCommandEvent& event ){

	wxString rd = wxFileSelector(wxString::FromAscii("Select file"), wxString::FromAscii(this->most_recent_folder.c_str()), wxString::FromAscii(""), wxString::FromAscii("txt"), wxString::FromAscii("*.txt"), wxFD_OPEN, this);
	if( rd.Len() > 0 ){
		std::string path( rd.ToAscii() );
		Rawdata raw;
		try{
			raw.load( path );
		}
		catch( std::string err){
			std::stringstream ss;
			ss << "ERROR loading data:\n" << err;
			wxMessageBox(wxString::FromAscii(ss.str().c_str()), wxString::FromAscii("ERROR loading file"), wxOK | wxICON_EXCLAMATION, this);
			return;
		}
		boost::filesystem::path rd_path(path);
		this->most_recent_folder = rd_path.parent_path().string();
		this->lbl_expr_value->SetLabel( rd );
		redraw();
	}
}


void InvestigationPropertiesPanel::OnClickSA( wxCommandEvent& event ){
	wxString sa = wxFileSelector(wxString::FromAscii("Select file"), wxString::FromAscii(this->most_recent_folder.c_str()), wxString::FromAscii(""), wxString::FromAscii("txt"), wxString::FromAscii("*.txt"), wxFD_OPEN, this);
	if( sa.Len() > 0 ){
		Attributes a("NA");
		try{
			a.load(std::string(sa.ToAscii()));
		}
		catch( std::string err){
			std::stringstream ss;
			ss << "ERROR loading attributes:\n" << err;
			wxMessageBox(wxString::FromAscii(ss.str().c_str()), wxString::FromAscii("ERROR loading file"), wxOK | wxICON_EXCLAMATION, this);
			return;
		}
		this->lbl_sa_value->SetLabel( sa );
		redraw();
	}
}


void InvestigationPropertiesPanel::OnClickGA( wxCommandEvent& event ){
	wxString fn_ga = wxFileSelector(wxString::FromAscii("Select file"), wxString::FromAscii(this->most_recent_folder.c_str()), wxString::FromAscii(""), wxString::FromAscii("txt"), wxString::FromAscii("*.txt"), wxFD_OPEN, this);
	if( fn_ga.Len() > 0 ){
		this->set_fn_ga(std::string(fn_ga.ToAscii()));
	}
}


void InvestigationPropertiesPanel::OnClickResults( wxCommandEvent& event ){
	wxString r = wxDirSelector(wxString::FromAscii("Select experimental results directory"), wxString::FromAscii(this->most_recent_folder.c_str()), wxDD_NEW_DIR_BUTTON, wxDefaultPosition, this);
	if( r.Len() > 0 ){
		this->lbl_results_value->SetLabel( r );
		redraw();
	}
}


void InvestigationPropertiesPanel::redraw(){
	this->Fit();
	this->parent->Fit();
}

void InvestigationPropertiesPanel::set_enabled(bool is_enabled){
	this->btn_raw_data->Enable(is_enabled);
	this->btn_ga->Enable(is_enabled);
	this->btn_sa->Enable(is_enabled);
	this->btn_results->Enable(is_enabled);
	this->redraw();
}
