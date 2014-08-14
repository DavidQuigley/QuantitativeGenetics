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
using namespace boost;
namespace alg = boost::algorithm;
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
using namespace std;
#include "carmen/src/DataStructures.h"
#include "carmen/src/Attributes.h"
#include "carmen/src/Rawdata.h"
#include "carmen/src/Parser.h"
#include "Investigation.h"
#include "PanelLimits.h"

IMPLEMENT_CLASS( SampleFilterDialog, wxDialog )
IMPLEMENT_CLASS( GeneFilterDialog, wxDialog )
IMPLEMENT_CLASS( GOFilterDialog, wxDialog )

GOFilterDialog::GOFilterDialog(Investigation* investigation){

	this->investigation= investigation;
	SetExtraStyle(wxWS_EX_BLOCK_EVENTS|wxDIALOG_EX_CONTEXTHELP);
	wxDialog::Create( NULL, wxID_ANY, wxString::FromAscii("Search by GO terms"), wxDefaultPosition, wxDefaultSize, wxCAPTION);

	this->GO_parser = new GOAnnotationParser();
	this->gene_parser = new GeneAnnotationParser();
	try{
		this->GO_parser->load(this->investigation->cp->get(std::string("Annotations"), std::string("go")));
	}
	catch(std::string err){
		wxMessageBox(wxString::FromAscii(err.c_str()), _T("ERROR loading GO annotation"), wxOK | wxICON_EXCLAMATION, this);
		return;
	}
	std::vector<int>* v;
	std::string gene_name_column = this->investigation->ga->get_gene_name_column();
	for(int i=0; i<(int)this->investigation->ga->identifiers.size(); i++){
		std::string identifier = this->investigation->ga->identifiers.at(i);
		std::string gene = this->investigation->ga->prop_for_identifier(identifier, gene_name_column);
		boost::algorithm::to_lower(gene);
		if( this->gene2probe_idx.find(gene)==this->gene2probe_idx.end() ){
			v = new std::vector<int>();
			gene2probe_idx[gene] = v;
		}
		else{
			v = gene2probe_idx[gene];
		}
		v->push_back(i);
	}

	CreateControls();
	GetSizer()->Fit(this); // This fits the dialog to the minimum size dictated by the sizers
	GetSizer()->SetSizeHints(this); // This ensures that the dialog cannot be sized smaller than the minimum size
	Centre();
}


void GOFilterDialog::CreateControls(){

	typedef std::pair<std::string, std::string> prop_pair;
	std::vector< prop_pair* >* annotations = new std::vector<prop_pair*>();
	this->investigation->cp->get_section(std::string("Annotations"), annotations);
	int species_idx=0;
	for(int i=0; i<(int)annotations->size(); i++){
	   	if( annotations->at(i)->first.compare("go") != 0 ){
	   		this->ars_annotations.Add(wxString::FromAscii(annotations->at(i)->first.c_str()));
	   		this->fn_annotations.push_back( annotations->at(i)->second );
	   		if( annotations->at(i)->first.compare( this->investigation->species ) == 0 )
	   			species_idx = this->ars_annotations.Count()-1;
	   	}
	   	delete annotations->at(i);
	}
	this->gene_parser->load( this->fn_annotations.at(0), *(this->GO_parser));
	this->loaded_annotation = this->fn_annotations.at(0);

	wxStaticText* lbl_go_terms = new wxStaticText( this, wxID_ANY, wxString::FromAscii("GO Term:"));
	this->txt_searchterms = new wxTextCtrl( this, ID_GOFILTER_SEARCHTERMS, wxString::FromAscii(""), wxDefaultPosition, wxSize(200, 20) );
	this->cho_annotation = new wxChoice( this, ID_GOFILTER_ANNOTATION, wxDefaultPosition, wxDefaultSize, this->ars_annotations);
	this->cho_annotation->SetSelection(species_idx);
	this->btn_search = new wxButton( this, ID_GOFILTER_BTNSEARCH, wxString::FromAscii("&Search"));

	wxStaticText* lbl_results = new wxStaticText( this, wxID_ANY, wxString::FromAscii("Matching Gene Ontology terms:"));
	this->lbx_results = new wxListBox(this, ID_GOFILTER_RESULTS, wxDefaultPosition, wxSize(450, 200), 0, NULL, wxLB_MULTIPLE);
	wxStaticText* lbl_mac = new wxStaticText( this, wxID_ANY, wxString::FromAscii("Use command-click to unselect a result"));
	this->lbl_genes = new wxStaticText( this, wxID_ANY, wxString::FromAscii("Genes (0):"));
	this->txt_genes= new wxTextCtrl( this, ID_GOFILTER_GENES, wxT(""), wxDefaultPosition, wxSize(150, 100), wxTE_MULTILINE);
	this->lbl_probes= new wxStaticText( this, wxID_ANY, wxString::FromAscii("Probes (0):"));
	this->txt_probes = new wxTextCtrl( this, ID_GOFILTER_PROBES, wxT(""), wxDefaultPosition, wxSize(150, 100), wxTE_MULTILINE);

	this->btn_cancel = new wxButton(this, wxID_CANCEL, wxString::FromAscii("&Cancel"));
	this->btn_ok = new wxButton(this, ID_GOFILTER_OK, wxString::FromAscii("&OK"));

	Connect(ID_GOFILTER_RESULTS, wxEVT_COMMAND_LISTBOX_SELECTED, wxCommandEventHandler(GOFilterDialog::OnClickResults));
	Connect(ID_GOFILTER_BTNSEARCH, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(GOFilterDialog::OnClickSearch));
	Connect(ID_GOFILTER_OK, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(GOFilterDialog::OnClickOk));
	wxBoxSizer* sizer_top = new wxBoxSizer(wxVERTICAL);

	wxBoxSizer* sizer_search = new wxBoxSizer(wxHORIZONTAL);
	sizer_search->Add(lbl_go_terms, 0, FLAGS_TOP, BORDER_PXL);
	sizer_search->Add(this->txt_searchterms, 0, FLAGS_TOP, BORDER_PXL);
	sizer_search->Add(this->cho_annotation, 0, FLAGS_TOP, BORDER_PXL);
	sizer_search->Add(this->btn_search, 0, FLAGS_TOP, BORDER_PXL);
	sizer_top->Add(sizer_search, 0, FLAGS, BORDER_PXL);

	sizer_top->Add(lbl_results, 0, FLAGS_TOP, BORDER_PXL);
	sizer_top->Add(this->lbx_results, 0, FLAGS_TOP, BORDER_PXL);
	sizer_top->Add(lbl_mac, 0, FLAGS_TOP, BORDER_PXL);
	if( !this->investigation->is_mac )
		lbl_mac->Show(false);

	wxFlexGridSizer* sizer_gp = new wxFlexGridSizer(2, 2, 0, 0);
	sizer_gp->Add(lbl_genes, 0, FLAGS_TOP, BORDER_PXL);
	sizer_gp->Add(lbl_probes, 0, FLAGS_TOP, BORDER_PXL);
	sizer_gp->Add(this->txt_genes, 0, FLAGS_TOP, BORDER_PXL);
	sizer_gp->Add(this->txt_probes, 0, FLAGS_TOP, BORDER_PXL);

	sizer_top->Add(sizer_gp, 0, FLAGS, BORDER_PXL);
	wxBoxSizer* sizer_btn = new wxBoxSizer(wxHORIZONTAL);
	sizer_btn->Add(this->btn_cancel, 0, FLAGS_TOP, BORDER_PXL);
	sizer_btn->Add(this->btn_ok, 0, FLAGS_TOP, BORDER_PXL);

	sizer_top->Add(sizer_btn, 0, FLAGS_RIGHT, BORDER_PXL);
	this->SetSizer(sizer_top);
}

void GOFilterDialog::OnClickResults(wxCommandEvent& WXUNUSED(event)){
	// read which symbols are selected
	std::vector<int> sel;
	for( int i=0; i<int(this->lbx_results->GetCount()); i++){
		if( this->lbx_results->IsSelected(i) ){
			sel.push_back(i);
		}
	}


	this->genes.clear();
	this->probes.clear();

	if( this->fn_annotations.at( this->cho_annotation->GetSelection() ).compare( this->loaded_annotation) != 0 ){
		this->loaded_annotation = this->fn_annotations.at(this->cho_annotation->GetSelection());
		this->gene_parser->load( this->fn_annotations.at(this->cho_annotation->GetSelection()) , *(this->GO_parser));
	}
	this->lbx_results->Refresh();
	this->lbx_results->Update();
	HASH_S_I seen;
	std::vector<std::string> search_results;
	for( int i=0; i<int(sel.size()); i++){
		int idx = this->result_GO_idx.at(sel.at(i));
		this->gene_parser->get_symbols_with_GO_idx(idx, search_results);
		for(int j=0; j<(int)search_results.size(); j++){
			if( seen.find(search_results.at(j)) == seen.end() ){
				seen[search_results.at(j)]=1;
				this->genes.push_back( search_results.at(j) );
			}
		}
	}

	std::sort(this->genes.begin(), this->genes.end() );
	std::vector<int>* v;
	std::string probe;
	wxString ssp, ssg;
	for(int i=0; i<(int)this->genes.size(); i++){
		if( this->gene2probe_idx.find( this->genes.at(i) ) != this->gene2probe_idx.end() ){
			v = this->gene2probe_idx[this->genes.at(i)];
			for(int j=0; j<(int)v->size(); j++){
				probe = this->investigation->ga->identifiers.at(v->at(j));
				ssp << wxString::FromAscii(probe.c_str() ) << wxString::FromAscii("\n");
				this->probes.push_back(probe);
			}
		}
		boost::algorithm::to_upper(genes.at(i));
		ssg <<  wxString::FromAscii(this->genes.at(i).c_str()) << wxString::FromAscii("\n");
	}
	this->txt_genes->SetValue( ssg );
	this->txt_probes->SetValue( ssp );
	std::stringstream ngenes, nprobes;
	ngenes << "Genes (" << this->genes.size() << "):";
	nprobes << "Probes (" << this->probes.size() << "):";
	this->lbl_genes->SetLabel( wxString::FromAscii(ngenes.str().c_str() ) );
	this->lbl_probes->SetLabel( wxString::FromAscii(nprobes.str().c_str() ) );
}


void GOFilterDialog::OnClickSearch(wxCommandEvent& WXUNUSED(event)){
	this->txt_genes->Clear();
	this->txt_probes->Clear();
	this->lbx_results->Clear();
	this->result_GO_idx.clear();

	this->gene_parser->load(this->fn_annotations.at( this->cho_annotation->GetSelection() ), *this->GO_parser);
	// DO SEARCH; simple version for now assumes one term
	std::string searchterm( this->txt_searchterms->GetValue().ToAscii() );
	this->GO_parser->get_idx_matching_description(searchterm, this->result_GO_idx);
	GOAnnotation* go;
	std::vector<std::string> results;
	for(int i=0; i<(int)this->result_GO_idx.size(); i++){
		go = this->GO_parser->get_annotation(this->result_GO_idx.at(i)  );
		results.push_back( go->description );
	}
	for( int i=0; i<(int)results.size(); i++){
		this->lbx_results->Append( wxString::FromAscii( results.at(i).c_str() ) );
	}
}


GOFilterDialog::~GOFilterDialog(){

}

std::string GOFilterDialog::get_probe_list(){
	std::stringstream ss;
	for(int i=0; i<(int)this->probes.size(); i++){
		ss << probes.at(i) << "\n";
	}
	return ss.str();
}

void GOFilterDialog::OnClickOk(wxCommandEvent& WXUNUSED(event)){
	for (boost::unordered_map <string,std::vector<int>*>::iterator i=this->gene2probe_idx.begin(); i!=gene2probe_idx.end(); ++i){
		delete i->second;
	}
	this->AcceptAndClose();
}


GeneFilterDialog::GeneFilterDialog(Investigation* investigation){
	this->investigation = investigation;
	SetExtraStyle(wxWS_EX_BLOCK_EVENTS|wxDIALOG_EX_CONTEXTHELP);
	wxDialog::Create( NULL, wxID_ANY, wxString::FromAscii("Choose Sample Filter"), wxDefaultPosition, wxDefaultSize, wxCAPTION);
	CreateControls();
	GetSizer()->Fit(this); // This fits the dialog to the minimum size dictated by the sizers
	GetSizer()->SetSizeHints(this); // This ensures that the dialog cannot be sized smaller than the minimum size
	Centre();
	wxBeginBusyCursor();
	std::string gene_name_column = this->investigation->ga->get_gene_name_column();
	for(int i=0; i<(int)this->investigation->ga->identifiers.size(); i++){
		std::string identifier = this->investigation->ga->identifiers.at(i);
		std::string gene = this->investigation->ga->prop_for_identifier(identifier, gene_name_column);
		boost::algorithm::to_lower(gene);
		this->g_p.push_back( new std::pair<std::string, std::string>(gene, identifier) );
	}
	wxEndBusyCursor();
}


void GeneFilterDialog::CreateControls(){

	wxStaticText* lbl_find = new wxStaticText( this, wxID_ANY, wxString::FromAscii("Find identifiers for gene:"));
	this->txt_gene_name = new wxTextCtrl( this, ID_TXT_GENE_NAME, wxString::FromAscii(""), wxDefaultPosition, wxSize(70, 20) );
	wxButton* btn_search = new wxButton( this, ID_BTN_SEARCH, wxString::FromAscii("&Search"));

	this->lbx_identifiers = new wxListBox(this, ID_LBX_IDENTIFIERS, wxDefaultPosition, wxSize(130, 70), 0, NULL, wxLB_SINGLE );

	wxArrayString ars_is;
	ars_is.Add(wxString::FromAscii( "is") );
	ars_is.Add(wxString::FromAscii( "is not") );
	this->cho_is = new wxChoice( this, ID_CHO_IS_ISNOT, wxDefaultPosition, wxDefaultSize, ars_is );
	this->cho_is->SetSelection(0);

	wxArrayString ars_up_down;
	ars_up_down.Add(wxString::FromAscii( "Up") );
	ars_up_down.Add(wxString::FromAscii( "Down") );
	this->cho_up_down = new wxChoice( this, ID_CHO_UP_DOWN, wxDefaultPosition, wxDefaultSize, ars_up_down );
	this->cho_up_down->SetSelection(0);

	wxSizer* sizer_lower = CreateButtonSizer(wxOK|wxCANCEL);

	Connect(ID_BTN_SEARCH, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(GeneFilterDialog::OnClickSearch));
	Connect(wxID_OK, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(GeneFilterDialog::OnClickOk));

	wxBoxSizer* sizer_top = new wxBoxSizer(wxVERTICAL);
	wxBoxSizer* sizer_upper = new wxBoxSizer(wxHORIZONTAL);
	wxBoxSizer* sizer_middle = new wxBoxSizer(wxHORIZONTAL);

	sizer_upper->Add(lbl_find, 0, FLAGS, BORDER_PXL);
	sizer_upper->Add(this->txt_gene_name, 0, FLAGS, BORDER_PXL);
	sizer_upper->Add(btn_search, 0, FLAGS, BORDER_PXL);
	
	sizer_middle->Add(this->lbx_identifiers, 0, FLAGS_TOP, BORDER_PXL);
	sizer_middle->Add(this->cho_is, 0, FLAGS_TOP, BORDER_PXL);
	sizer_middle->Add(this->cho_up_down, 0, FLAGS_TOP, BORDER_PXL);

	sizer_top->Add(sizer_upper, 0, FLAGS, BORDER_PXL);
	sizer_top->Add(sizer_middle, 0, FLAGS, BORDER_PXL);
	sizer_top->Add(sizer_lower, 0, FLAGS_RIGHT, BORDER_PXL);
	this->SetSizer(sizer_top);
}


void GeneFilterDialog::OnCloseWindow(wxCommandEvent& WXUNUSED(event)){
	for(int i=0; i<(int)this->g_p.size(); i++)
		delete g_p.at(i);
}

void GeneFilterDialog::OnClickOk(wxCommandEvent& WXUNUSED(event)){
	wxArrayString items = this->lbx_identifiers->GetStrings();
	wxArrayInt sel;
	wxString id;
	if( this->lbx_identifiers->GetSelections(sel)==0 ){
		wxMessageBox(wxString::FromAscii("Please select an identifier in the text box."), _T("ERROR: no identifier selected"), wxOK | wxICON_EXCLAMATION, this);
		return;
	}
	for(int i=0; i<(int)items.size(); i++){
		if( i == sel[0] )
			id = items[i];
	}
	this->limit = "gene:" + std::string( id.ToAscii() );
	if( this->cho_is->GetSelection()==0 )
		this->limit += "=";
	else
		this->limit += "!";
	if( this->cho_up_down->GetSelection()==0)
		this->limit += "1";
	else
		this->limit += "-1";
	this->AcceptAndClose();
}


void GeneFilterDialog::OnClickSearch(wxCommandEvent& WXUNUSED(event)){
	std::string searchterm( this->txt_gene_name->GetValue().ToAscii() );
	if( searchterm.size() == 0 ){
		wxMessageBox(wxString::FromAscii("Please enter a search term in the text box."), _T("ERROR: no search term"), wxOK | wxICON_EXCLAMATION, this);
		return;
	}
	this->lbx_identifiers->Clear();
	boost::algorithm::to_lower(searchterm);
	for(int i=0; i<int(this->g_p.size()); i++){
		if( this->g_p.at(i)->first.compare(searchterm)==0 ){
			this->lbx_identifiers->Append( wxString::FromAscii(this->g_p.at(i)->second.c_str() ));
		}
	}
	if( this->lbx_identifiers->GetCount() == 1 ){
		this->lbx_identifiers->SetSelection(0);
	}
}


SampleFilterDialog::SampleFilterDialog(Investigation* investigation){
	this->investigation = investigation;
	SetExtraStyle(wxWS_EX_BLOCK_EVENTS|wxDIALOG_EX_CONTEXTHELP);
	wxDialog::Create( NULL, wxID_ANY, wxString::FromAscii("Choose Sample Filter"), wxDefaultPosition, wxDefaultSize, wxCAPTION);
	CreateControls();
	GetSizer()->Fit(this); 
	GetSizer()->SetSizeHints(this);
	Centre(); 
}


void SampleFilterDialog::CreateControls(){
	wxArrayString ars_attributes;
	for(int i=0; i<(int)this->investigation->sa->attrib.size(); i++)
		ars_attributes.Add(wxString::FromAscii( this->investigation->sa->attrib.at(i).c_str() ));
	this->cho_attributes = new wxChoice( this, ID_CHO_ATTRIBUTES, wxDefaultPosition, wxDefaultSize, ars_attributes );
	this->cho_attributes->SetSelection(0);

	wxArrayString ars_is;
	ars_is.Add(wxString::FromAscii( "is") );
	ars_is.Add(wxString::FromAscii( "is not") );
	this->cho_is = new wxChoice( this, ID_CHO_IS, wxDefaultPosition, wxDefaultSize, ars_is );
	this->cho_is->SetSelection(0);

	this->cho_values = new wxChoice( this, ID_CHO_VALUES);
	
	//wxButton* btn_ok = new wxButton( this, ID_BTN_OK, wxString::FromAscii("&Ok"));
	//wxButton* btn_cancel = new wxButton( this, wxID_CANCEL, wxString::FromAscii("&Cancel") );
	//btn_ok->SetDefault();
	//this->SetEscapeId(wxID_CANCEL);
	//wxStdDialogButtonSizer* sizer_btn = CreateStdDialogButtonSizer(wxOK | wxCLOSE);
	wxSizer* sizer_ok = CreateButtonSizer(wxOK|wxCANCEL);

	Connect(ID_CHO_ATTRIBUTES, wxEVT_COMMAND_CHOICE_SELECTED, wxCommandEventHandler(SampleFilterDialog::OnChangeAttribute));
	Connect(wxID_OK, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(SampleFilterDialog::OnClickOk));
	Connect(wxID_CANCEL, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(SampleFilterDialog::OnClickCancel));

	wxBoxSizer* sizer_top = new wxBoxSizer(wxVERTICAL);
	wxBoxSizer* sizer_cbo = new wxBoxSizer(wxHORIZONTAL);
	sizer_cbo->Add(this->cho_attributes, 0, FLAGS, BORDER_PXL);
	sizer_cbo->Add(this->cho_is, 0, FLAGS, BORDER_PXL);
	sizer_cbo->Add(this->cho_values, 0, FLAGS, BORDER_PXL);
	
	//wxBoxSizer* sizer_btn = new wxBoxSizer(wxHORIZONTAL);
	//sizer_btn->Add(btn_ok, 0, FLAGS, BORDER_PXL);
	//sizer_btn->Add(btn_cancel, 0, FLAGS, BORDER_PXL);
	
	sizer_top->Add(sizer_cbo, 0, FLAGS, BORDER_PXL);
	sizer_top->Add(sizer_ok, 0, FLAGS_RIGHT, BORDER_PXL);
	this->SetSizer(sizer_top);
	redraw();
}

void SampleFilterDialog::OnChangeAttribute(wxCommandEvent& WXUNUSED(event)){
	this->redraw();
}

void SampleFilterDialog::OnClickOk( wxCommandEvent& event ){
	this->limit = this->investigation->sa->attrib.at( this->cho_attributes->GetSelection() );
	if( this->cho_is->GetSelection()==0 )
		this->limit += "=";
	else
		this->limit += "!";
	this->limit += this->values.at( this->cho_values->GetSelection() );
	this->AcceptAndClose();
}

void SampleFilterDialog::OnClickCancel( wxCommandEvent& event ){
	this->EndModal(wxID_CANCEL);
}

void SampleFilterDialog::redraw(){
	this->cho_values->Clear();
	this->values.clear();
	std::string attribute = this->investigation->sa->attrib.at( this->cho_attributes->GetSelection() );
	HASH_S_I seen;
	std::string value;
	
	for(int i=0; i<(int)this->investigation->sa->identifiers.size(); i++){
		value = this->investigation->sa->prop_for_identifier(this->investigation->sa->identifiers.at(i),attribute);
		if(seen.find( value ) == seen.end()){
			seen[value] = 1;
			this->values.push_back(value);
		}
	}
	std::sort(this->values.begin(), this->values.end());
	for(int i=0; i<(int)this->values.size(); i++)
		this->cho_values->Append(wxString::FromAscii(this->values.at(i).c_str()));
	this->cho_values->SetSelection(0);
}


PanelLimits::PanelLimits(wxWindow* parent, int id, Investigation* investigation) : wxPanel(parent, id){
	this->parent = parent;
    
	this->investigation = new Investigation();
    this->investigation->clone_from(investigation);
	this->lbx_limits = new wxListBox(this, CTL_TXT_LIMITS, wxDefaultPosition, wxSize(200, 50), 0, NULL, wxLB_SINGLE );
	int btn_size = 160;
	if( this->investigation->is_mac){
	    btn_size=120;
	}
	this->btn_remove = new wxButton( this, CTL_BTN_REMOVE, wxString::FromAscii("Remove"), wxDefaultPosition, wxSize(btn_size, 20));
	this->btn_add_sample = new wxButton( this, CTL_BTN_ADD_SAMPLE, wxString::FromAscii("Add Sample Limit"), wxDefaultPosition, wxSize(btn_size, 20));
	this->btn_add_gene = new wxButton( this, CTL_BTN_ADD_GENE, wxString::FromAscii("Add Gene Limit"), wxDefaultPosition, wxSize(btn_size, 20));
	this->btn_remove->Enable(false);
	wxFont font(this->btn_add_gene->GetFont());
	font.SetPointSize(wxSMALL_FONT->GetPointSize());
	this->btn_add_gene->SetFont( font );
	this->btn_add_sample->SetFont( font );
	this->btn_remove->SetFont( font );
	this->lbx_limits->SetFont( font );
	this->sizer_top = new wxBoxSizer(wxHORIZONTAL);
	this->sizer_top->Add(this->lbx_limits, 0, FLAGS, BORDER_PXL);
	
	wxBoxSizer* sizer_right = new wxBoxSizer(wxVERTICAL);
	sizer_right->Add(this->btn_remove, 0, FLAGS, BORDER_PXL);
	sizer_right->Add(this->btn_add_sample, 0, FLAGS, BORDER_PXL);
	sizer_right->Add(this->btn_add_gene, 0, FLAGS, BORDER_PXL);
	this->sizer_top->Add(sizer_right, 0, FLAGS, BORDER_PXL);
    
	Connect(CTL_TXT_LIMITS, wxEVT_COMMAND_LISTBOX_SELECTED, wxCommandEventHandler(PanelLimits::OnSelected));
	Connect(CTL_BTN_REMOVE, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(PanelLimits::OnRemove));
	Connect(CTL_BTN_ADD_SAMPLE, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(PanelLimits::OnAddSample));
	Connect(CTL_BTN_ADD_GENE, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(PanelLimits::OnAddGene));
    
	this->SetSizer(this->sizer_top);
	this->sizer_top->Fit(this);
	this->redraw();
}


void PanelLimits::set_investigation(std::string investigation_name){
    this->lbx_limits->Clear();
    this->investigation->set_current(investigation_name);
    this->redraw();
}


void PanelLimits::disable_gene_limits(){
    this->btn_add_gene->Enable(false);
    this->btn_add_gene->Show(false);
    this->sizer_top->Fit(this);
}

void PanelLimits::redraw(){
	this->btn_remove->Enable(false);
	this->lbx_limits->Enable(false);
	if( this->lbx_limits->GetCount()>0 ){
		this->lbx_limits->Enable(true);
	}
	for( int i=0; i<int(this->lbx_limits->GetCount()); i++){
		if( this->lbx_limits->IsSelected(i) ){
			this->btn_remove->Enable(true);
		}
	}
	wxCommandEvent evt( wxEVT_COMMAND_TEXT_UPDATED, this->GetId() );
	evt.SetEventObject(this);
	GetEventHandler()->ProcessEvent(evt);
}

void PanelLimits::OnSelected(wxCommandEvent& WXUNUSED(event)){
	this->redraw();
}


std::string PanelLimits::get_limits(){
	wxArrayString items = this->lbx_limits->GetStrings();
	std::stringstream ss;
	for(int i=0; i<(int)items.size(); i++){
		ss << items[i].ToAscii();
		if( i<(int)items.size()-1)
			ss << ",";
	}
	return ss.str();
}


void PanelLimits::OnRemove(wxCommandEvent& WXUNUSED(event)){
	std::vector<wxString> keep;
	for( int i=0; i<int(this->lbx_limits->GetCount()); i++){
		if( !this->lbx_limits->IsSelected(i) ){
			keep.push_back( this->lbx_limits->GetString(i) );
		}
	}
	this->lbx_limits->Clear();
	for(int i=0; i<(int)keep.size(); i++){
		this->lbx_limits->Append(keep.at(i));
	}
    this->redraw();
}

void PanelLimits::OnAddSample(wxCommandEvent& WXUNUSED(event)){
	SampleFilterDialog dialog(this->investigation);
	if (dialog.ShowModal() == wxID_OK){
		this->lbx_limits->Append( wxString::FromAscii(dialog.limit.c_str()) );
	}
	this->redraw();
}

void PanelLimits::OnAddGene(wxCommandEvent& WXUNUSED(event)){
	GeneFilterDialog dialog(this->investigation);
	if (dialog.ShowModal() == wxID_OK){
		this->lbx_limits->Append( wxString::FromAscii(dialog.limit.c_str()) );
	}
	this->redraw();
}

void PanelLimits::Clear(){
	this->set_limits(std::string(""));
}

void PanelLimits::set_limits(std::string limits){
	this->lbx_limits->Clear();
	std::vector<std::string> lim;
	alg::split(lim, limits, alg::is_any_of(",") );
	for(int i=0; i<(int)lim.size(); i++){
		if( lim.at(i).size()>0 ){
			this->lbx_limits->Append( wxString::FromAscii( lim.at(i).c_str()) );
		}
	}
    this->redraw();
}
