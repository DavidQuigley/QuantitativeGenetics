#include <wx/wx.h>
#include <wx/stattext.h>
#include <wx/textctrl.h>
#include <wx/sizer.h>
#include <wx/choice.h>
#include <wx/button.h>
#include <wx/arrstr.h>
#include <wx/dialog.h>
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
#include "DlgGenesToProbes.h"

IMPLEMENT_CLASS( GenesToProbesDialog, wxDialog )

GenesToProbesDialog::GenesToProbesDialog( Investigation* investigation ){
	this->investigation = investigation;
	SetExtraStyle(wxWS_EX_BLOCK_EVENTS|wxDIALOG_EX_CONTEXTHELP);
	wxDialog::Create( NULL, wxID_ANY, wxT("Convert genes to probes"), wxDefaultPosition, wxDefaultSize, wxCAPTION|wxRESIZE_BORDER|wxSYSTEM_MENU );
	CreateControls();
	GetSizer()->Fit(this);
	GetSizer()->SetSizeHints(this);
	Centre();
}

void GenesToProbesDialog::CreateControls(){

	wxStaticText* lbl_genes = new wxStaticText( this, wxID_ANY, wxT("Genes:") );
	this->txt_genes= new wxTextCtrl( this, ID_G2P_TXT_GENES, wxT(""), wxDefaultPosition, wxSize(150, 80), wxTE_MULTILINE );
	this->btn_convert = new wxButton( this, ID_G2P_BTN_CONVERT, wxT("&Convert"), wxDefaultPosition, wxSize(130,20));
	wxStaticText* lbl_probes = new wxStaticText( this, wxID_ANY, wxT("Probes:") );
	this->txt_probes = new wxTextCtrl( this, ID_G2P_TXT_PROBES, wxT(""), wxDefaultPosition, wxSize(150, 80), wxTE_MULTILINE );

	wxSizer* sizer_ok = CreateButtonSizer(wxOK|wxCANCEL);
	if( this->investigation->is_mac ){
		wxFont font(lbl_genes->GetFont());
		font.SetPointSize(wxSMALL_FONT->GetPointSize());
		lbl_genes->SetFont(font);
		this->txt_genes->SetFont(font);
		this->btn_convert->SetFont(font);
		this->txt_probes->SetFont(font);
	}
	Connect(ID_G2P_BTN_CONVERT, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(GenesToProbesDialog::OnClickConvert));

	wxBoxSizer* sizer_top = new wxBoxSizer(wxVERTICAL);
	wxBoxSizer* sizer_upper = new wxBoxSizer(wxHORIZONTAL);
	wxBoxSizer* sizer_left = new wxBoxSizer(wxVERTICAL);
	sizer_left->Add(lbl_genes, 0, FLAGS, BORDER_PXL);
	sizer_left->Add(this->txt_genes, 0, FLAGS, BORDER_PXL);

	sizer_upper->Add(sizer_left, 0, FLAGS, BORDER_PXL);
	sizer_upper->Add(btn_convert, 0, FLAGS_CENTER, BORDER_PXL);

	wxBoxSizer* sizer_right = new wxBoxSizer(wxVERTICAL);
	sizer_right->Add(lbl_probes, 0, FLAGS, BORDER_PXL);
	sizer_right->Add(this->txt_probes, 0, FLAGS, BORDER_PXL);
	sizer_upper->Add(sizer_right, 0, FLAGS, BORDER_PXL);

	sizer_top->Add(sizer_upper, 0, FLAGS, BORDER_PXL);

	sizer_top->Add(sizer_ok, 0, FLAGS_RIGHT, BORDER_PXL);
	this->SetSizer(sizer_top);
}

void GenesToProbesDialog::OnClickConvert(wxCommandEvent& WXUNUSED(event)){
	std::vector<std::pair<std::string, std::string>*> g_p;
	std::string gene_name_column = this->investigation->ga->get_gene_name_column();
	for(int i=0; i<(int)this->investigation->ga->identifiers.size(); i++){
    	std::string identifier = this->investigation->ga->identifiers.at(i);
    	std::string gene = this->investigation->ga->prop_for_identifier(identifier, gene_name_column);
    	boost::algorithm::to_lower(gene);
    	g_p.push_back( new std::pair<std::string, std::string>(gene, identifier) );
    }
	wxString genes_w = this->txt_genes->GetValue();
    std::string genes_string( genes_w.ToAscii() );
    boost::algorithm::to_lower(genes_string);
    std::vector<std::string> genes, probes;
    bool gene_found;
    int n_not_found=0;
    std::stringstream genes_not_found;
    genes_not_found << "Gene symbols not found in this investigation:\n";
    alg::split(genes, genes_string, alg::is_any_of("\n") );
    for(int i=0; i<(int)genes.size(); i++){
    	gene_found=false;
    	if( genes.at(i).size()==0 )
    		continue;
    	for(int j=0; j<int(g_p.size()); j++){
			if( g_p.at(j)->first.compare(genes.at(i))==0 ){
				probes.push_back(g_p.at(j)->second);
				gene_found=true;
			}
    	}
    	if( !gene_found ){
    		genes_not_found << "'" << genes.at(i) << "' ";
    		n_not_found++;
    		if(n_not_found >1 && n_not_found % 5 == 0)
    			genes_not_found << "\n";
    	}
    }

    if( n_not_found>0 ){
    	wxMessageBox(wxString::FromAscii(genes_not_found.str().c_str()), _T("ERROR: gene not found"), wxOK | wxICON_EXCLAMATION, this);
    }
    std::string probe_list;
    stringhasher::join( probe_list, &probes, "\n" );
    this->txt_probes->SetValue( wxString::FromAscii(probe_list.c_str()) );
    for(int i=0; i<int(g_p.size()); i++){
    	delete g_p.at(i);
    }
}
