#define _SCL_SECURE_NO_DEPRECATE 1
#include <wx/dialog.h>
#include <wx/sizer.h>
#include <wx/stattext.h>
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
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
using namespace boost;
using namespace std;
#include "carmen/src/DataStructures.h"
#include "carmen/src/Parser.h"
#include "DlgInvestigationChooser.h"

IMPLEMENT_CLASS( InvestigationChooserDialog, wxDialog )

InvestigationChooserDialog::InvestigationChooserDialog( ){
	ConfigParser* cp = new ConfigParser();
	Init(cp);
	delete cp;
}

InvestigationChooserDialog::InvestigationChooserDialog( wxWindow* parent, ConfigParser* cp, wxWindowID id, const wxString& caption, const wxPoint& pos, const wxSize& size, long style )
{
	Init(cp);
	Create(parent, id, caption, pos, size, style);
}

void InvestigationChooserDialog::Init( ConfigParser* cp )
{
	std::string investigations_str = cp->get("Filesets", "names");
	alg::split( this->investigations, investigations_str, alg::is_any_of(",") );
	this->investigation = std::string("");
}

bool InvestigationChooserDialog::Create( wxWindow* parent, wxWindowID id, const wxString& caption, const wxPoint& pos, const wxSize& size, long style )
{
	SetExtraStyle(wxWS_EX_BLOCK_EVENTS|wxDIALOG_EX_CONTEXTHELP);
	if (!wxDialog::Create( parent, id, caption, pos, size, style ))
		return false;
	CreateControls();
	GetSizer()->Fit(this); // This fits the dialog to the minimum size dictated by the sizers
	GetSizer()->SetSizeHints(this); // This ensures that the dialog cannot be sized smaller than the minimum size
	Centre(); 
	return true;
}


void InvestigationChooserDialog::CreateControls()
{
	// CONTROLS
	wxStaticText* lbl_investigations = new wxStaticText( this, wxID_ANY, wxT("Open:"));
	wxStaticText* lbl_sortby = new wxStaticText( this, wxID_ANY, wxT("Sort by:"));

    wxArrayString invests, sort_orders;
	for(int i=0; i<(int)this->investigations.size(); i++){
		if( this->investigations.at(i).length()>0 )
            this->investigations_current_sorted_order.push_back( this->investigations.at(i));
			invests.Add(wxString::FromAscii(this->investigations.at(i).c_str()));
	}
    sort_orders.Add( wxString::FromAscii("Most recent"));
    sort_orders.Add( wxString::FromAscii("Alphabetical"));

	this->cho_investigations = new wxChoice( this, ID_CHO_INVESTIGATIONS, wxDefaultPosition, wxDefaultSize, invests );
    this->cho_sortby = new wxChoice( this, ID_CHO_SORTBY, wxDefaultPosition, wxDefaultSize, sort_orders );

	if( invests.size() > 0){
		this->cho_investigations->SetSelection(0);
		this->investigation = this->investigations.at(0);
	}
	wxSizer* sizer_ok = CreateButtonSizer(wxOK|wxCANCEL);
#ifdef WIN32
	bool is_mac=false;
#else
	bool is_mac=true;
#endif
	if( is_mac ){
		wxFont font(lbl_investigations->GetFont());
		font.SetPointSize(wxSMALL_FONT->GetPointSize());
		cho_investigations->SetFont(font);
	}
	Connect(ID_CHO_INVESTIGATIONS, wxEVT_COMMAND_CHOICE_SELECTED, wxCommandEventHandler(InvestigationChooserDialog::OnChangeInvestigation));
    Connect(ID_CHO_SORTBY, wxEVT_COMMAND_CHOICE_SELECTED, wxCommandEventHandler(InvestigationChooserDialog::OnChangeSortBy));

	Connect(wxEVT_CHAR, wxKeyEventHandler(InvestigationChooserDialog::OnChar));
	// POSITIONING
	wxFlexGridSizer* sizer_upper = new wxFlexGridSizer(2, 2, 0, 0);
	sizer_upper->Add(lbl_investigations, 0, FLAGS, BORDER_PXL);
	sizer_upper->Add(this->cho_investigations, 0, FLAGS, BORDER_PXL);
	sizer_upper->Add(lbl_sortby, 0, FLAGS, BORDER_PXL);
	sizer_upper->Add(this->cho_sortby, 0, FLAGS, BORDER_PXL);

	wxBoxSizer* sizer_top = new wxBoxSizer(wxVERTICAL);
	sizer_top->Add(sizer_upper, 0, FLAGS, BORDER_PXL);
	sizer_top->Add(sizer_ok, 0, FLAGS_RIGHT, BORDER_PXL);
	this->SetSizer(sizer_top);
}

void InvestigationChooserDialog::OnChar(wxKeyEvent& event){
	if (event.GetKeyCode() == WXK_ESCAPE){
		this->EndModal(wxCANCEL);
	}
	event.Skip();
};

void InvestigationChooserDialog::OnChangeInvestigation( wxCommandEvent& event ){
	this->investigation = this->investigations_current_sorted_order.at(  this->cho_investigations->GetCurrentSelection() );
}

void InvestigationChooserDialog::OnChangeSortBy( wxCommandEvent& event ){
    if( this->cho_sortby->GetCurrentSelection()==1 ){
        std::sort( this->investigations_current_sorted_order.begin(), this->investigations_current_sorted_order.end() );
    }
    else{
        this->investigations_current_sorted_order.clear();
        for(int i=0; i<(int)this->investigations.size(); i++){
            if( this->investigations.at(i).length()>0 ){
                this->investigations_current_sorted_order.push_back( this->investigations.at(i) );
            }
        }
        
    }
    this->cho_investigations->Clear();
    for(int i=0; i<int(this->investigations_current_sorted_order.size()); i++){
        this->cho_investigations->Append(wxString::FromAscii(this->investigations_current_sorted_order.at(i).c_str()));
    }
    if( this->investigations_current_sorted_order.size() > 0){
		this->cho_investigations->SetSelection(0);
		this->investigation = this->investigations_current_sorted_order.at(0);
	}
}

