#define _CRT_SECURE_NO_DEPRECATE 1
#define _SCL_SECURE_NO_DEPRECATE 1
#include <time.h>
#include <wx/dialog.h>
#include <wx/sizer.h>
#include <wx/stattext.h>
#include <wx/button.h>
#include <wx/textctrl.h>
#include <wx/radiobox.h>
#include <wx/statbox.h>
#include <wx/choice.h>
#include <wx/font.h>
#include <wx/wx.h>
#include <wx/process.h>
#include <wx/dynarray.h>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/random.hpp>
#include <boost/thread/mutex.hpp>
namespace alg = boost::algorithm;
#include <fstream>
#include <string>
#include <vector>
#include <set>
#include "carmen/src/DataStructures.h"
#include "carmen/src/Attributes.h"
#include "carmen/src/Rawdata.h"
#include "carmen/src/ClassMinerOptions.h"
#include "carmen/src/Discretize.h"
#include "carmen/src/Dataset.h"
#include "carmen/src/graph.h"
#include "carmen/src/Parser.h"
#include "carmen/src/spear.h"
#include "Investigation.h"
using namespace boost;
using namespace std;
#include "PanelLimits.h"
#include "DlgProgress.h"
#include "DlgCorrelationNetwork.h"

IMPLEMENT_CLASS( CorrelationNetworkDialog, wxDialog )

CorrelationNetworkDialog::CorrelationNetworkDialog( Investigation* investigation, std::string fn_spear ){
	this->investigation = investigation;
	this->fn_spear = fn_spear;
    SetExtraStyle(wxWS_EX_BLOCK_EVENTS|wxDIALOG_EX_CONTEXTHELP);
	wxDialog::Create( NULL, wxID_ANY, wxT("Generate Correlation Network"), wxDefaultPosition, wxDefaultSize, wxCAPTION|wxRESIZE_BORDER|wxSYSTEM_MENU );
	CreateControls();
	GetSizer()->Fit(this); 
	GetSizer()->SetSizeHints(this); 
	Centre();
}

void CorrelationNetworkDialog::CreateControls(){
	char timebuf[80];
	struct tm* newtime;
	time_t long_time;
	time( &long_time );
	newtime = localtime( &long_time );
	strftime(timebuf, 80, "%Y_%m_%d_%H_%M_%S", newtime);
	std::string gene_name_col = this->investigation->ga->get_gene_name_column();
    for(int i=0; i<(int)this->investigation->ga->identifiers.size(); i++){
		std::string identifier = this->investigation->ga->identifiers.at(i);
		std::string gene = this->investigation->ga->prop_for_identifier(identifier, gene_name_col);
		boost::algorithm::to_lower(gene);
		this->g_p.push_back( new std::pair<std::string, std::string>(gene, identifier) );
	}
    wxStaticBoxSizer* static_correlation_sizer = new wxStaticBoxSizer(wxVERTICAL, this, wxString::FromAscii("Correlation Settings") );
    wxStaticBoxSizer* static_seed_sizer = new wxStaticBoxSizer(wxVERTICAL, this, wxString::FromAscii("Seed Probes") );
    wxStaticBoxSizer* static_sample_sizer = new wxStaticBoxSizer(wxVERTICAL, this, wxString::FromAscii("Sample Limits") );
    wxStaticBoxSizer* static_eqtl_sizer = new wxStaticBoxSizer(wxVERTICAL, this, wxString::FromAscii("Expression QTL Settings") );
    wxStaticBoxSizer* static_output_sizer = new wxStaticBoxSizer(wxVERTICAL, this, wxString::FromAscii("Output Settings") );

    wxStaticText* lbl_file_name = new wxStaticText( this, wxID_ANY, wxT("File Name:") );
    this->txt_file_name = new wxTextCtrl( this, ID_CORR_NET_TXT_FILE_NAME, wxString::FromAscii(timebuf), wxDefaultPosition, wxSize(200, 22) );
    this->txt_file_name->ChangeValue(wxString::FromAscii(timebuf));

	wxArrayString ars_new_file;
	ars_new_file.Add(wxT( "Generate new correlations") );
	if( this->fn_spear.size() > 0)
		ars_new_file.Add(wxT( "Use selected correlation file") );
	this->cho_new_file = new wxChoice( this, ID_CORR_NEW_FILE, wxDefaultPosition, wxDefaultSize, ars_new_file );
	this->cho_new_file->SetSelection(0);
	if( this->fn_spear.size() == 0){
		this->cho_new_file->Show(false);
	}

	wxStaticText* lbl_min_corr = new wxStaticText( this, wxID_ANY, wxT("Minimum correlation: +/-"));
	this->txt_min_corr = new wxTextCtrl( this, ID_CORR_NET_TXT_MIN_CORR, wxT("0.7"), wxDefaultPosition, wxSize(60, 22) );
	this->txt_min_corr->SetValue( wxT("0.7") );
	wxStaticText* lbl_min_clique = new wxStaticText( this, wxID_ANY, wxT("Minimum clique size:"));
	wxArrayString ars_min_clique;
	ars_min_clique.Add(wxT("1"));
	ars_min_clique.Add(wxT("2"));
	ars_min_clique.Add(wxT("3"));
	ars_min_clique.Add(wxT("4"));
	this->cho_min_clique = new wxChoice( this, ID_CORR_NET_CHO_MIN_CLIQUE, wxDefaultPosition, wxDefaultSize, ars_min_clique );
	this->cho_min_clique->SetSelection(0);

	this->txt_seeds = new wxTextCtrl( this, ID_CORR_NET_TXT_SEEDS, wxT(""), wxDefaultPosition, wxSize(200, 70), wxTE_MULTILINE );
	this->btn_GO = new wxButton( this, ID_CORR_NET_BTN_GO, wxT("&Add by Ontology"), wxDefaultPosition, wxSize(130,22));
	this->txt_gene = new wxTextCtrl( this, ID_CORR_NET_TXT_GENE, wxT(""), wxDefaultPosition, wxSize(200, 22));
	this->btn_add_gene = new wxButton( this, ID_CORR_NET_BTN_ADD_GENE, wxT("&Add Seed"), wxDefaultPosition, wxSize(100,22));
	wxStaticText* lbl_limit_seeds = new wxStaticText( this, wxID_ANY, wxT("Limit network to seeds"));
	this->chk_limit_seeds = new wxCheckBox(this, ID_CORR_NET_CHK_LIMIT_SEEDS, wxString::FromAscii(""));

	this->limit_A = new PanelLimits(this, ID_CORR_NET_PANEL_A, this->investigation);
	this->limit_A->disable_gene_limits();

	wxStaticText* lbl_min_present = new wxStaticText( this, wxID_ANY, wxT("Minimum % Present:"));
	this->txt_min_present = new wxTextCtrl( this, ID_CORR_NET_TXT_MIN_PRESENT, wxT("90"), wxDefaultPosition, wxSize(60, 22) );
	this->txt_min_present->SetValue( wxT("90") );
	this->lbl_eqtl_file_name = new wxStaticText( this, wxID_ANY, wxT("eQTL Source File:") );
	this->lbl_eqtl_file_name_value = new wxStaticText( this, wxID_ANY, wxT("(none)") );
	this->btn_eqtl = new wxButton( this, ID_CORR_NET_BTN_EQTL, wxT("&Browse"), wxDefaultPosition, wxSize(70,22));
	this->lbl_maximum_pvalue= new wxStaticText( this, wxID_ANY, wxT("Maximum eQTL P value:") );
	this->txt_maximum_pvalue= new wxTextCtrl( this, ID_CORR_NET_TXT_MAXIMUM_PVALUE, wxString::FromAscii("1.0"), wxDefaultPosition, wxSize(50, 22) );
	this->chk_require_eqtl = new wxCheckBox(this, ID_CORR_NET_CHK_REQUIRE_EQTL, wxString::FromAscii(""));
	this->lbl_require_eqtl= new wxStaticText( this, wxID_ANY, wxT("Require that gene has eQTL"));
	this->chk_allow_uncorrelated = new wxCheckBox(this, ID_CORR_NET_CHK_ALLOW_UNCORRELATED, wxString::FromAscii(""));
	this->lbl_allow_uncorrelated = new wxStaticText( this, wxID_ANY, wxT("Allow genes with eQTL but no correlations"));

	this->lbl_maximum_pvalue->Show(false);
	this->txt_maximum_pvalue->Show(false);
	this->chk_require_eqtl->Show(false);
	this->lbl_require_eqtl->Show(false);
	this->chk_allow_uncorrelated->Show(false);
	this->lbl_allow_uncorrelated->Show(false);

	wxStaticText* lbl_attributes = new wxStaticText( this, wxID_ANY, wxT("Extra attribute:"));
	wxArrayString ars_extra;
	ars_extra.Add(wxT("NONE"));
	for(int i=1; i<(int)this->investigation->ga->attrib.size(); i++){
		ars_extra.Add(wxString::FromAscii( this->investigation->ga->attrib.at(i).c_str() ) );
	}
	this->cho_attributes = new wxChoice( this, ID_CORR_NET_CHO_ATTRIBUTES, wxDefaultPosition, wxDefaultSize, ars_extra );
	this->cho_attributes->SetSelection(0);

	wxStaticText* lbl_annotation = new wxStaticText( this, wxID_ANY, wxT("Annotation:"));
	typedef std::pair<std::string, std::string> prop_pair;
	std::vector< prop_pair* >* annotations = new std::vector<prop_pair*>();
	this->investigation->cp->get_section(std::string("Annotations"), annotations);
	wxArrayString ars_annot;
	int species_idx = 0;
	for(int i=0; i<(int)annotations->size(); i++){
	   	if( annotations->at(i)->first.compare("go") == 0 ){
	   		this->fn_go = annotations->at(i)->second;
	   	}
	   	else{
	   		ars_annot.Add(wxString::FromAscii(annotations->at(i)->first.c_str()));
	   		this->fn_annotations.push_back( annotations->at(i)->second );
	   		if( annotations->at(i)->first.compare( this->investigation->species ) == 0 )
	   			 species_idx = ars_annot.Count()-1;
	   	}
	   	delete annotations->at(i);
	}
	this->cho_annotation = new wxChoice( this, ID_CORR_NET_CHO_ANNOTATION, wxDefaultPosition, wxDefaultSize, ars_annot );
	if( ars_annot.size()>0 )
		this->cho_annotation->SetSelection(species_idx);

	wxSizer* sizer_ok = CreateButtonSizer(wxOK|wxCANCEL);

	if(investigation->is_mac){
		wxFont font(lbl_file_name->GetFont());
		font.SetPointSize(wxSMALL_FONT->GetPointSize());
		cho_new_file->SetFont(font);
		lbl_min_corr->SetFont(font);
		this->txt_min_corr->SetFont(font);
		lbl_min_clique->SetFont(font);
		this->cho_min_clique->SetFont(font);
		this->txt_seeds->SetFont(font);
		this->btn_GO->SetFont(font);
		this->txt_gene->SetFont(font);
		this->btn_add_gene->SetFont(font);
		lbl_limit_seeds->SetFont(font);
		this->chk_limit_seeds->SetFont(font);
		lbl_min_present->SetFont(font);
		this->txt_min_present->SetFont(font);
		this->lbl_eqtl_file_name->SetFont(font);
		this->lbl_eqtl_file_name_value->SetFont(font);
		this->btn_eqtl->SetFont(font);
		this->lbl_maximum_pvalue->SetFont(font);
		this->txt_maximum_pvalue->SetFont(font);
		this->chk_require_eqtl->SetFont(font);
		this->lbl_require_eqtl->SetFont(font);
		this->chk_allow_uncorrelated->SetFont(font);
		this->lbl_allow_uncorrelated->SetFont(font);
		lbl_file_name->SetFont(font);
		this->txt_file_name->SetFont(font);
		lbl_attributes->SetFont(font);
		this->cho_attributes->SetFont(font);
		lbl_annotation->SetFont(font);
		this->cho_annotation->SetFont(font);
	}
	Connect(wxID_OK, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(CorrelationNetworkDialog::OnClickOk));
    Connect(ID_CORR_NET_BTN_ADD_GENE, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(CorrelationNetworkDialog::OnClickAddGene));
    Connect(ID_CORR_NET_BTN_GO, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(CorrelationNetworkDialog::OnClickGO));
    Connect(ID_CORR_NET_BTN_EQTL,  wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(CorrelationNetworkDialog::OnClickEQTL));

    this->sizer_top = new wxBoxSizer(wxVERTICAL);
    wxFlexGridSizer* sizer_filename = new wxFlexGridSizer(1, 2, 0, 0);
    sizer_filename->Add(lbl_file_name, 0, FLAGS, BORDER_PXL);
    sizer_filename->Add(this->txt_file_name, 0, FLAGS, BORDER_PXL);
    sizer_top->Add(sizer_filename, 0, FLAGS, BORDER_PXL);

	// CORRELATION SETTINGS
	static_correlation_sizer->Add(this->cho_new_file, 0, FLAGS, BORDER_PXL);

	wxFlexGridSizer* sizer_min_corr = new wxFlexGridSizer(2, 2, 0, 0);
	sizer_min_corr->Add(lbl_min_corr, 0, FLAGS, BORDER_PXL);
	sizer_min_corr->Add(this->txt_min_corr, 0, FLAGS, 0);
	sizer_min_corr->Add(lbl_min_clique, 0, FLAGS, BORDER_PXL);
	sizer_min_corr->Add(this->cho_min_clique, 0, FLAGS, BORDER_PXL);

	static_correlation_sizer->Add(sizer_min_corr, 0, FLAGS, BORDER_PXL);
	this->sizer_top->Add( static_correlation_sizer, 0, FLAGS, BORDER_PXL);

	wxBoxSizer* sizer_seeds= new wxBoxSizer(wxHORIZONTAL);
	sizer_seeds->Add(this->txt_seeds, 0, FLAGS, BORDER_PXL);
	sizer_seeds->Add(this->btn_GO, 0, FLAGS_TOP_LEFT, BORDER_PXL);
	static_seed_sizer->Add(sizer_seeds, 0, FLAGS, BORDER_PXL);

	wxBoxSizer* sizer_add_gene= new wxBoxSizer(wxHORIZONTAL);
	sizer_add_gene->Add(this->txt_gene, 0, FLAGS, BORDER_PXL);
	sizer_add_gene->Add(this->btn_add_gene, 0, FLAGS_TOP_LEFT, BORDER_PXL);
	static_seed_sizer->Add(sizer_add_gene, 0, FLAGS, BORDER_PXL);

	wxBoxSizer* sizer_limit = new wxBoxSizer(wxHORIZONTAL);
	sizer_limit->Add(this->chk_limit_seeds, 0, FLAGS, BORDER_PXL);
	sizer_limit->Add(lbl_limit_seeds, 0, FLAGS, BORDER_PXL);
	static_seed_sizer->Add(sizer_limit, 0, FLAGS, BORDER_PXL);
	this->sizer_top->Add( static_seed_sizer, 0, FLAGS, BORDER_PXL);

	// SAMPLE SETTINGS
	static_sample_sizer->Add(this->limit_A, 0, FLAGS, BORDER_PXL);
	wxBoxSizer* sizer_min_present= new wxBoxSizer(wxHORIZONTAL);
	sizer_min_present->Add(lbl_min_present, 0, FLAGS, BORDER_PXL);
	sizer_min_present->Add(this->txt_min_present, 0, FLAGS, BORDER_PXL);
	static_sample_sizer->Add(sizer_min_present, 0, FLAGS, BORDER_PXL);
	this->sizer_top->Add( static_sample_sizer, 0, FLAGS, BORDER_PXL);

	// EQTL SETTINGS
	wxFlexGridSizer* sizer_eqtl_top = new wxFlexGridSizer(2, 2, 0, 0);
	sizer_eqtl_top->Add(lbl_eqtl_file_name, 0, FLAGS, BORDER_PXL);
	wxBoxSizer* sizer_eqtl_topright = new wxBoxSizer(wxHORIZONTAL);
	sizer_eqtl_topright->Add(this->lbl_eqtl_file_name_value, 0, FLAGS, BORDER_PXL);
	sizer_eqtl_topright->Add(this->btn_eqtl, 0, FLAGS, BORDER_PXL);
	sizer_eqtl_top->Add(sizer_eqtl_topright, 0, FLAGS, BORDER_PXL);
	sizer_eqtl_top->Add(lbl_maximum_pvalue, 0, FLAGS, BORDER_PXL);
	sizer_eqtl_top->Add(txt_maximum_pvalue, 0, FLAGS, BORDER_PXL);

	wxFlexGridSizer* sizer_eqtl_bottom = new wxFlexGridSizer(2, 2, 0, 0);
	sizer_eqtl_bottom->Add(this->chk_require_eqtl, 0, FLAGS, BORDER_PXL);
	sizer_eqtl_bottom->Add(this->lbl_require_eqtl, 0, FLAGS, BORDER_PXL);
	sizer_eqtl_bottom->Add(this->chk_allow_uncorrelated, 0, FLAGS, BORDER_PXL);
	sizer_eqtl_bottom->Add(this->lbl_allow_uncorrelated, 0, FLAGS, BORDER_PXL);

	static_eqtl_sizer->Add(sizer_eqtl_top, 0, FLAGS, BORDER_PXL);
	static_eqtl_sizer->Add(sizer_eqtl_bottom, 0, FLAGS, BORDER_PXL);
	this->sizer_top->Add( static_eqtl_sizer, 0, FLAGS, BORDER_PXL);

	wxFlexGridSizer* sizer_other = new wxFlexGridSizer(2, 2, 0, 0);
	sizer_other->Add(lbl_attributes, 0, FLAGS, BORDER_PXL);
	sizer_other->Add(this->cho_attributes, 0, FLAGS, BORDER_PXL);
	sizer_other->Add(lbl_annotation, 0, FLAGS, BORDER_PXL);
	sizer_other->Add(this->cho_annotation, 0, FLAGS, BORDER_PXL);
	static_output_sizer->Add(sizer_other, 0, FLAGS, BORDER_PXL);

	this->sizer_top->Add(static_output_sizer, 0, FLAGS, BORDER_PXL);
	this->sizer_top->Add( sizer_ok, 0, FLAGS_RIGHT, BORDER_PXL);
	SetSizer(this->sizer_top);
}

void CorrelationNetworkDialog::OnClickEQTL(wxCommandEvent& WXUNUSED(event)){
	wxString fn = wxFileSelector(wxString::FromAscii("Select file"), wxString::FromAscii(""), wxString::FromAscii(""), wxString::FromAscii("txt"), wxString::FromAscii("*.txt"), wxFD_OPEN, this);
	if( fn.Len() > 0 ){
		this->fn_eqtl = std::string(fn.ToAscii());
		if( int(fn.size()) > 25 ){
			std::stringstream ss;
			ss << "..." << fn.Right(25).ToAscii();
			this->lbl_eqtl_file_name_value->SetLabel(wxString::FromAscii(ss.str().c_str()));
		}
		else
			this->lbl_eqtl_file_name_value->SetLabel(fn);
		this->lbl_maximum_pvalue->Show(true);
		this->txt_maximum_pvalue->Show(true);
		this->chk_require_eqtl->Show(true);
		this->lbl_require_eqtl->Show(true);
		this->chk_allow_uncorrelated->Show(true);
		this->lbl_allow_uncorrelated->Show(true);
		this->sizer_top->Fit(this);
	}
}

void CorrelationNetworkDialog::OnClickGO(wxCommandEvent& WXUNUSED(event)){
	GOFilterDialog dialog(this->investigation);
	if (dialog.ShowModal() == wxID_OK){
		this->txt_seeds->AppendText( wxString::FromAscii(dialog.get_probe_list().c_str()) );
	}
}

void CorrelationNetworkDialog::OnClickAddGene(wxCommandEvent& WXUNUSED(event)){

    wxString search_w = this->txt_gene->GetValue();
    std::string searchterm(search_w.ToAscii());
	boost::algorithm::to_lower(searchterm);
    wxArrayString aProbes;
	for(int i=0; i<(int)this->g_p.size(); i++)
		if( this->g_p.at(i)->first.compare(searchterm)==0 )
            aProbes.Add( wxString::FromAscii( this->g_p.at(i)->second.c_str()) );
    if( aProbes.size()==0 ){
        wxMessageBox(wxString::FromAscii("ERROR: Gene symbol not found in this investigation."), wxString::FromAscii("ERROR: gene not found"), wxOK | wxICON_EXCLAMATION, this);
		return;
    }
    wxString probe_id;
    if( aProbes.size()==1 ){
        probe_id = aProbes.Item(0);
    }
    else{
        probe_id = wxGetSingleChoice(wxString::FromAscii("Choose a probe"), wxString::FromAscii("Choose a probe"), aProbes);
    }
    if( probe_id.size()>0 ){
        wxString items = this->txt_seeds->GetValue();
        if( items.size()==0 )
            this->txt_seeds->SetValue( probe_id );
        else{
            std::string probes_raw( items.ToAscii() );
            std::vector<std::string> probes;
            alg::split(probes, probes_raw, alg::is_any_of("\n") );
            probes.push_back( std::string( probe_id.ToAscii() ) );
            std::string probe_list;
            stringhasher::join( probe_list, &probes, "\n" );
            this->txt_seeds->SetValue( wxString::FromAscii(probe_list.c_str()) );
        }
    }
}


void CorrelationNetworkDialog::OnClickOk(wxCommandEvent& event){
	std::stringstream ss;
	if(this->investigation->is_mac)
		ss << this->investigation->dir_results << '/' << this->txt_file_name->GetValue().ToAscii();
	else
		ss << this->investigation->dir_results << '\\' << this->txt_file_name->GetValue().ToAscii();
    std::string base_location( ss.str() );
    std::string f( this->fn_spear );

    double min_abs_d, percent_required_d;
    std::string min_abs(this->txt_min_corr->GetValue().ToAscii());
    std::string percent_required( this->txt_min_present->GetValue().ToAscii() );
	try{
		min_abs_d = boost::lexical_cast<double>( min_abs );
	}
	catch(boost::bad_lexical_cast blc){
		wxMessageBox(wxString::FromAscii("Minimum correlation must be a number; stopping analysis."), _T("ERROR: bad minimum correlation"), wxOK | wxICON_EXCLAMATION, this);
		return;
	}
	try{
		percent_required_d = boost::lexical_cast<double>(percent_required);
	}
	catch(boost::bad_lexical_cast blc){
		wxMessageBox(wxString::FromAscii("Minimum percent required must be a number between 0 and 100; stopping analysis."), _T("ERROR: bad minimum percent required"), wxOK | wxICON_EXCLAMATION, this);
		return;
	}
	if( percent_required_d<0 || percent_required_d>100){
		wxMessageBox(wxString::FromAscii("Minimum percent required must be a number between 0 and 100; stopping analysis."), _T("ERROR: bad minimum percent required"), wxOK | wxICON_EXCLAMATION, this);
		return;
	}
    std::string p = "";
    wxString items = this->txt_seeds->GetValue();
	std::vector<std::string> probes_final, probes_unknown, probes_passed; 
	if( (int)items.size()==0 ){
        if( wxMessageBox(wxString::FromAscii("WARNING: No seed probes selected.\nAnalysis could take a very long time."), _T("No seed probes"), wxOK|wxCANCEL, this) ==wxCANCEL)
			return;
	}
    else{
        std::vector<std::string> cmd;
        std::string probes_raw( items.ToAscii() );
        alg::split(probes_passed, probes_raw, alg::is_any_of("\n") );
        for(int i=0; i<(int)probes_passed.size(); i++){
            if( (int)probes_passed.at(i).size()>0 ){
                if( this->investigation->ga->identifier2idx.find(probes_passed.at(i)) == this->investigation->ga->identifier2idx.end() ){
					probes_unknown.push_back(probes_passed.at(i));
                }
				else{
					probes_final.push_back(probes_passed.at(i));
				}
			}
        }
		if( probes_unknown.size()>0 ){
			if( (int)probes_unknown.size() > 10 ){
				wxMessageBox(wxString::FromAscii("ERROR: more than 10 unrecognizable probes"), _T("ERROR: unknown probe"), wxOK | wxICON_EXCLAMATION, this);
			}
			else{
	 			std::string s("ERROR: could not recoginize the following probe(s):\n");
				for(int i=0; i<(int)probes_unknown.size(); i++){
					s += probes_unknown.at(i) += " ";
				}
				wxMessageBox(wxString::FromAscii(s.c_str()), _T("ERROR: unknown probe"), wxOK | wxICON_EXCLAMATION, this);
			}
			if(probes_final.size()==0){
				wxMessageBox(wxString::FromAscii("No probes were recognized; stopping analysis."), _T("ERROR: unknown probe"), wxOK | wxICON_EXCLAMATION, this);
				return;		
			}
		}
        stringhasher::join( p, &probes_final, "," );
    }
	std::string fn_cytoscape( this->investigation->cp->get(std::string("External"), std::string("cytoscape")) );
    std::string fn_props( this->investigation->cp->get(std::string("Internal"), std::string("cytoscape_vizmap")) );
    std::vector<std::string> cmd;
	cmd.push_back( "CorrelationNetwork" );
	cmd.push_back("-g" + this->investigation->fn_ga );
	if( this->cho_new_file->GetSelection()==0 ){
		cmd.push_back("-d" + this->investigation->fn_expr );
		cmd.push_back("-f" + this->investigation->fn_sa );
		cmd.push_back("-a" + this->limit_A->get_limits() );
	}
	else{
		cmd.push_back("-r" + f );
	}
	if( this->fn_eqtl.size()>0 ){
		cmd.push_back("-q" + this->fn_eqtl);
		double max_perm_p;
		try{
			max_perm_p = boost::lexical_cast<double>(this->txt_maximum_pvalue->GetValue().ToAscii());
		}
		catch(boost::bad_lexical_cast blc){
			wxMessageBox(wxString::FromAscii("Maximum P-value must be between 0 and 1; stopping analysis."), _T("ERROR: bad maximum P-value"), wxOK | wxICON_EXCLAMATION, this);
			return;
		}
		if( max_perm_p >1 || max_perm_p < 0 ){
			wxMessageBox(wxString::FromAscii("Maximum P-value must be between 0 and 1; stopping analysis."), _T("ERROR: bad maximum P-value"), wxOK | wxICON_EXCLAMATION, this);
			return;
		}
		cmd.push_back("-k" + std::string( this->txt_maximum_pvalue->GetValue().ToAscii() ) );
		if( this->chk_require_eqtl->IsChecked() )
			cmd.push_back("-jT");
		else
			cmd.push_back("-jF");
		if( this->chk_allow_uncorrelated->IsChecked() )
			cmd.push_back("-uT");
		else
			cmd.push_back("-uF");
	}
	cmd.push_back( "-p" + p); // seed probes
	cmd.push_back( "-n" + boost::lexical_cast<std::string>(percent_required_d / 100) );
	cmd.push_back( "-s" + min_abs );
	cmd.push_back( "-eT" ); // include correlations between neighbors, required to avoid "star" network
	int idx_clique = this->cho_min_clique->GetSelection();
	if( idx_clique == 0 )
		cmd.push_back("-c1");
	else if( idx_clique == 1 )
		cmd.push_back("-c2");
	else if( idx_clique == 2 )
		cmd.push_back("-c3");
	else if( idx_clique == 3 )
		cmd.push_back("-c4");
	if( this->chk_limit_seeds->IsChecked() )
		cmd.push_back("-lT");
	else
		cmd.push_back("-lF");
	cmd.push_back("-y" + this->investigation->ga->get_gene_name_column()  );
	cmd.push_back("-m0"); // minimum variance
	cmd.push_back("-vT"); // verbose
	int idx_a = this->cho_attributes->GetSelection();
	if( idx_a > 0 )
		cmd.push_back( "-x" + this->investigation->ga->attrib.at(idx_a));
	cmd.push_back("-t" + fn_cytoscape);
	cmd.push_back("-b" + fn_props);
	cmd.push_back("-h" + fn_go);
	cmd.push_back("-i" + this->fn_annotations.at( this->cho_annotation->GetSelection() ));
	cmd.push_back("-o" + base_location );

	ProgressDialog dialog(this, cmd, this->investigation, true);
	if (dialog.ShowModal() == wxID_OK){
		if( this->investigation->is_mac)
			this->batch_filename = std::string(base_location + ".sh");
		else
			this->batch_filename = std::string(base_location + ".bat");
		this->AcceptAndClose();
	}
}
