#define _CRT_SECURE_NO_DEPRECATE 1
#include <string>
#include <vector>
#include <sstream>
#include <wx/wx.h>
#include <wx/dialog.h>
#include <wx/panel.h>
#include <wx/sizer.h>
#include <wx/stattext.h>
#include <wx/button.h>
#include <wx/textctrl.h>
#include <wx/choice.h>
#include <wx/dynarray.h>
#include <wx/arrstr.h>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/random.hpp>
namespace alg = boost::algorithm;

using namespace std;
#include "carmen/src/DataStructures.h"
#include "carmen/src/Attributes.h"
#include "carmen/src/Rawdata.h"
#include "carmen/src/Parser.h"
#include "Investigation.h"
#include "PanelDiscretization.h"
#include "PanelLimits.h"
#include "DlgDiscretization.h"

IMPLEMENT_CLASS( DiscretizationDialog, wxDialog )

DiscretizationDialog::DiscretizationDialog( wxWindow* parent, Investigation* investigation ){
	
	SetExtraStyle(wxWS_EX_BLOCK_EVENTS|wxDIALOG_EX_CONTEXTHELP);
    wxDialog::Create( parent, wxID_ANY, wxT("Choose a discretization"));
	
    this->investigation = investigation;
	wxStaticText* lbl_discretization = new wxStaticText(this, wxID_ANY, wxString::FromAscii("Discretization:"));
    this->panel_disc = new PanelDiscretization(this, DISC_DLG_DISC_PANEL, this->investigation);
    
    //wxStaticText* lbl_file_name = new wxStaticText(this, wxID_ANY, wxString::FromAscii("File Name:"));
    //std::stringstream ss;
	//ss << investigation->fn_expr << ".discretized";
    //this->txt_file_name = new wxTextCtrl( this, DISC_DLG_FILE_NAME, wxString::FromAscii(ss.str().c_str()), wxDefaultPosition, wxSize(250, 20) );
    //this->txt_file_name->ChangeValue(wxString::FromAscii(ss.str().c_str()));
    
    wxStaticText* lbl_limits = new wxStaticText(this, wxID_ANY, wxString::FromAscii("Sample Limits:"));
    this->limit_A = new PanelLimits(this, DISC_DLG_DISC_LIMIT_A, this->investigation);
    this->limit_A->disable_gene_limits();

    wxSizer* sizer_buttons = CreateButtonSizer(wxOK|wxCANCEL);

    Connect(wxID_OK, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(DiscretizationDialog::OnClickOk));
    
    wxFlexGridSizer* sizer_upper = new wxFlexGridSizer(2, 2, 0, 0);
	sizer_upper->Add(lbl_discretization, 0, FLAGS, BORDER_PXL);
    sizer_upper->Add(this->panel_disc, 0, FLAGS, BORDER_PXL);

    //sizer_upper->Add(lbl_file_name, 0, FLAGS, BORDER_PXL);
    //sizer_upper->Add(this->txt_file_name, 0, FLAGS, BORDER_PXL);
    sizer_upper->Add(lbl_limits, 0, FLAGS, BORDER_PXL);
    sizer_upper->Add(this->limit_A, 0, FLAGS, BORDER_PXL);

    wxBoxSizer* sizer_top = new wxBoxSizer(wxVERTICAL);
    sizer_top->Add(sizer_upper, 0, FLAGS, BORDER_PXL);
    sizer_top->Add(sizer_buttons, 0, FLAGS_RIGHT, BORDER_PXL);

    this->SetSizer(sizer_top);
    GetSizer()->Fit(this);
}


void DiscretizationDialog::OnClickOk( wxCommandEvent& evt ){
    this->panel_disc->get_discretization_bounds(this->disc_lower, this->disc_upper);
    this->panel_disc->get_discretization_method_name(this->disc_method);
    this->limits = this->limit_A->get_limits();

    wxString default_location = wxString::FromAscii(this->investigation->dir_results.c_str());
    wxString default_name = wxString::FromAscii("expr.discretized.txt");
    wxString default_extension = wxString::FromAscii("txt");
    wxString filename = wxFileSelector(wxString::FromAscii("Save To"), default_location, default_name, default_extension, wxString::FromAscii("*.*"), wxFD_SAVE, this);
    if( this->investigation->cp->get(std::string("Internal"), std::string("debug")).compare(std::string("true"))==0)
    	std::cout << "Writing to: " << filename << "\n";
    std::string fn( filename.ToAscii() );
    if( fn.size()==0 ){
        wxMessageBox(wxString::FromAscii("Please enter a name for the discretized file that will be generated"), _T("Not Ready"), wxOK | wxICON_EXCLAMATION, this);
		return;
    }
    this->fn_discretization = fn;
    this->AcceptAndClose();
}
