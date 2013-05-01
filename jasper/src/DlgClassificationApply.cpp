#define _CRT_SECURE_NO_DEPRECATE 1
#define _SCL_SECURE_NO_DEPRECATE 1
#include <stdio.h>
#include <string>
#include <algorithm>
#include <iostream>
#include <math.h>
#include <fstream>
#include <vector>
using namespace std;
#include <wx/dialog.h>
#include <wx/sizer.h>
#include <wx/stattext.h>
#include <wx/arrstr.h>
#include <wx/choice.h>
#include <wx/textctrl.h>
#include <wx/panel.h>
#include <wx/process.h>
#include <wx/grid.h>
#include <wx/wx.h>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/random.hpp>
using namespace boost;
namespace alg = boost::algorithm;
#include "carmen/src/DataStructures.h"
#include "carmen/src/Attributes.h"
#include "carmen/src/Rawdata.h"
#include "carmen/src/Parser.h"
#include "Investigation.h"
#include "PanelDiscretization.h"
#include "PanelLimits.h"
#include "DlgProgress.h"
#include "DlgClassificationApply.h"

IMPLEMENT_CLASS( ClassificationApplyDialog, wxDialog )

ClassificationApplyDialog::ClassificationApplyDialog( std::string filename, Investigation* investigation){
	this->investigation = investigation;
	SetExtraStyle(wxWS_EX_BLOCK_EVENTS|wxDIALOG_EX_CONTEXTHELP);
	wxDialog::Create( NULL, wxID_ANY, wxT("Apply Classifier"), wxDefaultPosition, wxSize(400,350), wxCAPTION|wxRESIZE_BORDER|wxSYSTEM_MENU );
	CreateControls(filename);
	GetSizer()->Fit(this); 
	GetSizer()->SetSizeHints(this); 
	Centre();
}

void ClassificationApplyDialog::CreateControls(std::string filename){
    
    this->filename = filename;
    wxStaticText* lbl_file_name = new wxStaticText( this, wxID_ANY, wxT("File Name:") );
	this->txt_file_name = new wxTextCtrl( this, ID_CLASS_A_TXT_FILE_NAME, wxT(""), wxDefaultPosition, wxSize(200, 20) );

    wxStaticText* lbl_investigations = new wxStaticText( this, wxID_ANY, wxT("Choose an Investigation"), wxDefaultPosition, wxDefaultSize, 0 );
	std::string investigations_str = this->investigation->cp->get("Filesets", "names");
	alg::split( this->investigations, investigations_str, alg::is_any_of(",") );
    wxArrayString invests;
	for(int i=0; i<(int)this->investigations.size(); i++)
		invests.Add(wxString::FromAscii(this->investigations.at(i).c_str()));
	this->cho_investigations = new wxChoice( this, ID_CLASS_A_CHO_INVESTIGATIONS, wxDefaultPosition, wxDefaultSize, invests );
	if( invests.size() > 0){
		this->cho_investigations->SetSelection(0);
	}

    std::vector<std::string> path_parts;
	alg::split( path_parts, filename, alg::is_any_of("\\") );
	std::string fn_ruleset = path_parts.at( path_parts.size()-1 );
    alg::replace_all(fn_ruleset, ".classifier", "");
	alg::replace_all(fn_ruleset, ".labels", "");

    this->txt_file_name->SetValue(wxString::FromAscii(fn_ruleset.c_str() ));
    
    wxStaticText* lbl_disc = new wxStaticText( this, wxID_ANY, wxT("Discretization:") );
	this->discretization = new PanelDiscretization(this, ID_CLASS_A_DISC, this->investigation);

    wxStaticBox* static_a = new wxStaticBox(this, wxID_ANY, wxString::FromAscii("Class A Definition"));
    wxStaticBoxSizer* static_class_a_sizer = new wxStaticBoxSizer(static_a, wxHORIZONTAL);
	wxStaticBox* static_b = new wxStaticBox(this, wxID_ANY, wxString::FromAscii("Class B Definition"));
    wxStaticBoxSizer* static_class_b_sizer = new wxStaticBoxSizer(static_b, wxHORIZONTAL);
    this->limit_A = new PanelLimits(this, ID_CLASS_A_LIMIT_A, this->investigation);
	this->limit_B = new PanelLimits(this, ID_CLASS_A_LIMIT_B, this->investigation);

    std::string limit_a_str, limit_b_str, disc_str, disc_mult, disc_u_str, disc_l_str;
    if( alg::find_first(filename, ".labels")){
        CarmenParser cp(filename);
        cp.get_value("Class_A_test", 0, limit_a_str);
        cp.get_value("Class_B_test", 0, limit_b_str);
        this->limit_A->set_limits( limit_a_str );
        this->limit_B->set_limits( limit_b_str );
        cp.get_value("Discretization", 0, disc_str);
        if( disc_str.compare("none")!=0 && disc_str.compare("None")!=0){
            cp.get_value(std::string("Discretization"), 3, disc_mult);
		    std::vector<std::string> mu;
		    alg::split(mu, disc_mult, alg::is_any_of(",") );
		    disc_l_str = mu.at(0);
	    	disc_u_str = mu.at(1);
    		this->discretization->set_discretization(disc_str, disc_l_str, disc_u_str);
        }
        std::string fn_raw_tag, fn_raw_file, fn_data;
        cp.get_value("Data_file", 0, fn_data);
        for(int i=0;i<(int)this->investigations.size();i++){
            fn_raw_tag = this->investigations.at(i) + "<>raw";
            fn_raw_file = this->investigation->cp->get("Filesets", fn_raw_tag);
            if(fn_raw_file.compare(fn_data)==0){
                this->cho_investigations->SetSelection(i);
                break;
            }
        }
        
    }

    wxSizer* sizer_ok = CreateButtonSizer(wxOK|wxCANCEL);

    Connect(wxID_OK, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(ClassificationApplyDialog::OnClickOk));
    Connect(ID_CLASS_A_CHO_INVESTIGATIONS, wxEVT_COMMAND_CHOICE_SELECTED, wxCommandEventHandler(ClassificationApplyDialog::OnChooseInvestigation));

    wxFlexGridSizer* sizer_upper = new wxFlexGridSizer(3, 2, 0, 0);
    sizer_upper->Add(lbl_file_name, 0, FLAGS, BORDER_PXL);
	sizer_upper->Add(this->txt_file_name, 0, FLAGS, BORDER_PXL);
    sizer_upper->Add(lbl_investigations, 0, FLAGS, BORDER_PXL);
	sizer_upper->Add(this->cho_investigations, 0, FLAGS, BORDER_PXL);
    sizer_upper->Add(lbl_disc, 0, FLAGS, BORDER_PXL);
	sizer_upper->Add(this->discretization, 0, FLAGS, BORDER_PXL);
    
    wxBoxSizer* sizer_top = new wxBoxSizer(wxVERTICAL);
	sizer_top->Add(sizer_upper, 0, FLAGS, BORDER_PXL);

    wxBoxSizer* sizer_limits = new wxBoxSizer(wxHORIZONTAL);
	static_class_a_sizer->Add(this->limit_A, 0, FLAGS, BORDER_PXL);
	static_class_b_sizer->Add(this->limit_B, 0, FLAGS, BORDER_PXL);
	sizer_limits->Add(static_class_a_sizer, 0, FLAGS, BORDER_PXL);
	sizer_limits->Add(static_class_b_sizer, 0, FLAGS, BORDER_PXL);

    sizer_top->Add(sizer_limits, 0, FLAGS_RIGHT, BORDER_PXL);
	sizer_top->Add(sizer_ok, 0, FLAGS_RIGHT, BORDER_PXL);
	
	this->SetSizer(sizer_top);
}


void ClassificationApplyDialog::reload_investigation(){
    this->limit_A->set_investigation( this->investigations.at( this->cho_investigations->GetCurrentSelection() ) );
    this->limit_B->set_investigation( this->investigations.at( this->cho_investigations->GetCurrentSelection() ) );
}

void ClassificationApplyDialog::OnClickOk(wxCommandEvent& WXUNUSED(event)){
    
    std::string d( this->limit_A->investigation->fn_expr );
	std::string f( this->limit_A->investigation->fn_sa );
	std::string g( this->limit_A->investigation->fn_ga );
    std::string o( this->txt_file_name->GetValue().ToAscii() );
    std::string m, l, u;
    std::string a = this->limit_A->get_limits();
	std::string b = this->limit_B->get_limits();

    if(a.size()==0 && b.size()==0){
		wxMessageBox(wxString::FromAscii("Please select a limit for Class A or Class B."), _T("Not Ready"), wxOK | wxICON_EXCLAMATION, this);
		return;
	}
    l = "";
    u = "";
    int d_idx = this->discretization->get_discretization_code();
    if( d_idx==DISC_NONE )
		m = "none";
	else if( d_idx==DISC_SD ){
		m="SD";
	}
    else if( d_idx==DISC_MAD ){
		m = "MAD";
	}
    else if( d_idx==DISC_ABS ){
		m = "abs";
	}
    else if( d_idx==DISC_PER ){
		m="per";
	}
	else if( d_idx==DISC_SD_SAMPLES){
		m = "SD_SAMPLES";
	}
	if( d_idx != DISC_NONE ){
		this->discretization->get_discretization_bounds(l, u);
	}

    if( o.size()==0 ){
		wxMessageBox(wxString::FromAscii("Please enter a name for the classifier file that will be generated"), _T("Not Ready"), wxOK | wxICON_EXCLAMATION, this);
		return;
	}
	if( ! alg::find_first( o, ".labels") )
		o = o + ".labels";
	std::string attempted_filename = this->investigation->dir_results + '\\' + o;

    std::vector<std::string> cmd;
	cmd.push_back( "classify" );
	cmd.push_back("-d" + d );
	cmd.push_back("-f" + f );
	cmd.push_back("-g" + g );
    cmd.push_back("-a" + a );
    if((int)b.size()>0)
        cmd.push_back("-b" + b );
	cmd.push_back("-c" + this->filename);
    cmd.push_back("-m" + m );
    cmd.push_back("-l" + l );
    cmd.push_back("-u" + u );
    cmd.push_back("-vT");
    cmd.push_back("-o" + attempted_filename );
	ProgressDialog dialog(this, cmd, this->investigation, true);
	if (dialog.ShowModal() == wxID_OK){
        this->new_filename = attempted_filename;
		this->AcceptAndClose();
	}
}


void ClassificationApplyDialog::OnChooseInvestigation( wxCommandEvent& event ){
    this->reload_investigation();
}
