#define _CRT_SECURE_NO_DEPRECATE 1
#define _SCL_SECURE_NO_DEPRECATE 1
#include <wx/dialog.h>
#include <wx/sizer.h>
#include <wx/stattext.h>
#include <wx/button.h>
#include <wx/textctrl.h>
#include <wx/panel.h>
#include <wx/radiobox.h>
#include <wx/statbox.h>
#include <wx/wx.h>
#include <wx/process.h>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/random.hpp>
namespace alg = boost::algorithm;
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
#include "PanelLimits.h"
#include "DlgProgress.h"
#include "DlgTmev.h"

IMPLEMENT_CLASS( TmevDialog, wxDialog )

TmevDialog::TmevDialog( Investigation* investigation, std::string filename ){
	this->investigation = investigation;
	
	SetExtraStyle(wxWS_EX_BLOCK_EVENTS|wxDIALOG_EX_CONTEXTHELP);
	wxDialog::Create( NULL, wxID_ANY, wxT("Create Dataset file"), wxDefaultPosition, wxDefaultSize, wxCAPTION|wxRESIZE_BORDER|wxSYSTEM_MENU );
	CreateControls(filename);

	GetSizer()->Fit(this); // This fits the dialog to the minimum size dictated by the sizers
	GetSizer()->SetSizeHints(this); // This ensures that the dialog cannot be sized smaller than the minimum size
	Centre();
}

void TmevDialog::CreateControls(std::string rs_filename)
{
	// CONTROLS
	//wxStaticText* lbl_file_name = new wxStaticText( this, wxID_ANY, wxT("File Name:") );
	//this->txt_file_name = new wxTextCtrl( this, ID_TMEV_TXT_FILE_NAME, wxT("expr_tmev.txt"), wxDefaultPosition, wxSize(150, 20) );
	//this->txt_file_name->ChangeValue(wxString::FromAscii("expr_tmev.txt"));

	wxStaticBoxSizer* static_sample_sizer = new wxStaticBoxSizer(wxVERTICAL, this, wxString::FromAscii("Sample Limits") );
	wxStaticBoxSizer* static_probe_sizer = new wxStaticBoxSizer(wxVERTICAL, this, wxString::FromAscii("Probe Limits") );

	this->limit_A = new PanelLimits(this, ID_TMEV_LIMIT_A, this->investigation);
    this->limit_A->disable_gene_limits();
    this->rdo_all = new wxRadioButton(this, ID_TMEV_RDO_ALL, wxString::FromAscii("Use all probes"), wxDefaultPosition, wxDefaultSize, wxRB_GROUP);
    this->rdo_file = new wxRadioButton(this, ID_TMEV_RDO_FILE, wxString::FromAscii("Probes in currently selected file"));
    this->rdo_paste = new wxRadioButton(this, ID_TMEV_RDO_PASTE, wxString::FromAscii("Paste probes into box"));
    this->rs_filename = rs_filename;
    if( rs_filename.size()==0 ){
        this->rdo_file->Enable(false);
    }
    this->rdo_all->SetValue(true);
    this->txt_probes = new wxTextCtrl( this, ID_TMEV_TXT_PROBES, wxT(""), wxDefaultPosition, wxSize(150, 80), wxTE_MULTILINE );
    this->btn_GO = new wxButton( this, ID_TMEV_BTN_GO, wxT("&Add by Ontology"), wxDefaultPosition, wxSize(150,22));
    this->txt_probes->Show(false);
    this->btn_GO->Show(false);

	wxStaticText* lbl_min_percent_required = new wxStaticText( this, wxID_ANY, wxT("Minimum % called present:"));
	this->txt_min_percent_present = new wxTextCtrl( this, ID_TMEV_TXT_MIN_PERCENT_PRESENT, wxT("0"), wxDefaultPosition, wxSize(40, 22) );
	this->txt_min_percent_present->SetValue( wxT("0") );

	wxSizer* sizer_ok = CreateButtonSizer(wxOK|wxCANCEL);

	Connect(wxID_OK, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(TmevDialog::OnClickOk));
	Connect(ID_TMEV_RDO_ALL, wxEVT_COMMAND_RADIOBUTTON_SELECTED, wxCommandEventHandler(TmevDialog::OnChangeAttribute));
	Connect(ID_TMEV_RDO_FILE, wxEVT_COMMAND_RADIOBUTTON_SELECTED, wxCommandEventHandler(TmevDialog::OnChangeAttribute));
	Connect(ID_TMEV_RDO_PASTE, wxEVT_COMMAND_RADIOBUTTON_SELECTED, wxCommandEventHandler(TmevDialog::OnChangeAttribute));
	Connect(ID_TMEV_BTN_GO, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(TmevDialog::OnClickGO));
	
    this->sizer_top = new wxBoxSizer(wxVERTICAL);
	
    //wxBoxSizer* sizer_filename = new wxBoxSizer(wxHORIZONTAL);
	//sizer_filename->Add(lbl_file_name, 0, FLAGS, BORDER_PXL);
	//sizer_filename->Add(this->txt_file_name, 0, FLAGS, BORDER_PXL);
    //this->sizer_top->Add(sizer_filename, 0, FLAGS, BORDER_PXL);

    static_sample_sizer->Add( this->limit_A, 0, FLAGS, 0);
    this->sizer_top->Add(static_sample_sizer, 0, FLAGS, BORDER_PXL);

    wxBoxSizer* sizer_probe = new wxBoxSizer(wxVERTICAL);
    sizer_probe->Add(this->rdo_all, 0, FLAGS, BORDER_PXL);
    sizer_probe->Add(this->rdo_file, 0, FLAGS, BORDER_PXL);
    sizer_probe->Add(this->rdo_paste, 0, FLAGS, BORDER_PXL);
    wxBoxSizer* sizer_probe_box = new wxBoxSizer(wxVERTICAL);
    sizer_probe_box->Add(this->txt_probes, 0, FLAGS, BORDER_PXL);
    sizer_probe_box->Add(this->btn_GO, 0, FLAGS_TOP_LEFT, BORDER_PXL);
    sizer_probe->Add(sizer_probe_box, 0, FLAGS, BORDER_PXL);
    wxBoxSizer* sizer_percent= new wxBoxSizer(wxHORIZONTAL);
    sizer_percent->Add(lbl_min_percent_required , 0, FLAGS, BORDER_PXL);
	sizer_percent->Add(this->txt_min_percent_present, 0, FLAGS, BORDER_PXL);
	sizer_probe->Add(sizer_percent, 0, FLAGS, BORDER_PXL);
    static_probe_sizer->Add( sizer_probe, 0, FLAGS, 0);
    this->sizer_top->Add(static_probe_sizer, 0, FLAGS, BORDER_PXL);

    this->sizer_top->Add(sizer_ok, 0, FLAGS_RIGHT, BORDER_PXL);
    if(this->investigation->is_mac){
    	// mac uses giant fonts by default; PC is better.
    	wxFont font(this->rdo_all->GetFont());
    	font.SetPointSize(wxSMALL_FONT->GetPointSize());
        this->rdo_all->SetFont( font );
        this->rdo_file->SetFont( font );
        this->rdo_paste->SetFont( font );
        //lbl_file_name->SetFont( font );
        //this->txt_file_name->SetFont( font );
        lbl_min_percent_required->SetFont( font );
        this->txt_min_percent_present->SetFont( font );
   	}
	this->SetSizer(this->sizer_top);
}


void TmevDialog::OnClickGO(wxCommandEvent& WXUNUSED(event)){
	GOFilterDialog dialog(this->investigation);
	if (dialog.ShowModal() == wxID_OK){
		this->txt_probes->AppendText( wxString::FromAscii(dialog.get_probe_list().c_str()) );
	}
}


void TmevDialog::OnChangeAttribute(wxCommandEvent& WXUNUSED(event)){
	if( this->rdo_paste->GetValue() ){
		this->txt_probes->Show(true);
		this->btn_GO->Show(true);
	}
	else{
		this->txt_probes->Show(false);
		this->txt_probes->Clear();
		this->btn_GO->Show(true);
	}
	this->sizer_top->Fit(this);
}

void TmevDialog::OnClickOk( wxCommandEvent& WXUNUSED(event) ){
	std::string limits = this->limit_A->get_limits();
	float percent_required = boost::lexical_cast<float>( this->txt_min_percent_present->GetValue().ToAscii() );
    std::string f = "";
    std::vector<std::string> cmd;
	std::vector<std::string> probes, probes_validated;
	if( percent_required < 1 && percent_required != 0 ){
		wxMessageBox(wxString::FromAscii("Percent required should be between 0 and 100"), _T("ERROR: incorrect percent required"), wxOK | wxICON_EXCLAMATION, this);
		return;
	}
	percent_required = percent_required / 100; // expecting a value between 0 and 1

	if(this->rdo_file->GetValue()){
        f = this->rs_filename; 
		CarmenParser cp;
		cp.Load(f);
		cp.PrepareToReadValues();
		cp.ExtractProbes(probes);
    }
	else if(this->rdo_paste->GetValue() ){
        wxString items = this->txt_probes->GetValue();
        std::string probes_raw( items.ToAscii() );
        alg::split(probes, probes_raw, alg::is_any_of("\n") );
    }
    for(int i=0; i<(int)probes.size(); i++){
    	if( (int)probes.at(i).size()>0 ){
    		if( this->investigation->ga->identifier2idx.find(probes.at(i)) == this->investigation->ga->identifier2idx.end() ){
    			std::string s("ERROR: unknown probe ");
                s += probes.at(i);
                wxMessageBox(wxString::FromAscii(s.c_str()), _T("ERROR: unknown probe"), wxOK | wxICON_EXCLAMATION, this);
    		    return;
    		}
    		else{
    			probes_validated.push_back(probes.at(i));
    		}
        }
    }
	wxString default_location = wxString::FromAscii(this->investigation->dir_results.c_str());
	wxString default_name = wxString::FromAscii("dataset.txt");
	wxString default_extension = wxString::FromAscii("txt");
	wxString filename = wxFileSelector(wxString::FromAscii("Save To"), default_location, default_name, default_extension, wxString::FromAscii("*.*"), wxFD_SAVE, this);
	if( filename.size()==0 ){
		return;
	}
	wxBeginBusyCursor();
	std::string fn_out( filename.ToAscii() );
	try{
		std::vector<std::string> samples;
		if(limits.size()==0)
			limits = std::string("IDENTIFIER!NULL");
		this->investigation->sa->find_identifiers_in_class(limits, &samples);
		this->investigation->write_dataset( fn_out, probes_validated, samples, percent_required);
		this->new_filename = fn_out;
		this->AcceptAndClose();
	}
	catch( std::string err ){
		wxEndBusyCursor();
		std::stringstream ss;
		ss << "ERROR writing dataset: " << err;
		wxMessageBox( wxString::FromAscii(ss.str().c_str()));
	}
	wxEndBusyCursor();
}
