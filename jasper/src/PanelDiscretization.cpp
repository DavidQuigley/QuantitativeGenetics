#include <wx/sizer.h>
#include <wx/stattext.h>
#include <wx/button.h>
#include <wx/panel.h>
#include <wx/dialog.h>
#include <wx/choice.h>
#include <wx/dynarray.h>
#include <wx/arrstr.h>
#include <wx/wx.h>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>
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
#include "PanelDiscretization.h"


PanelDiscretization::PanelDiscretization(wxWindow* parent, int id, Investigation* investigation) : wxPanel(parent, id){
	this->parent = parent;
	
	wxArrayString ars_disc;
	ars_disc.Add(wxString::FromAscii("No Discretization"));
	ars_disc.Add(wxString::FromAscii("Standard Deviation"));
	ars_disc.Add(wxString::FromAscii("Median Deviation"));
	ars_disc.Add(wxString::FromAscii("Absolute Cutoffs"));
	ars_disc.Add(wxString::FromAscii("Percentile"));
	ars_disc.Add(wxString::FromAscii("Standard Deviation by Sample"));
	
    this->method_abbreviations.push_back("none");
    this->method_abbreviations.push_back("SD");
    this->method_abbreviations.push_back("MAD");
    this->method_abbreviations.push_back("abs");
    this->method_abbreviations.push_back("per");
    this->method_abbreviations.push_back("SD_SAMPLES");
    
    this->cho_disc = new wxChoice( this, ID_CHOICE_DISC, wxDefaultPosition, wxDefaultSize, ars_disc );
	this->cho_disc->SetSelection(0);
	this->bounds.push_back( std::pair<std::string, std::string>("plus/minus 0.2", "0.2" ) );
	this->bounds.push_back( std::pair<std::string, std::string>("plus/minus 0.25", "0.25" ) );
	this->bounds.push_back( std::pair<std::string, std::string>("plus/minus 0.3", "0.3" ) );
	this->bounds.push_back( std::pair<std::string, std::string>("plus/minus 0.5", "0.5" ) );
	this->bounds.push_back( std::pair<std::string, std::string>("plus/minus 0.7", "0.7" ) );
	this->bounds.push_back( std::pair<std::string, std::string>("plus/minus 1.0", "1.0" ) );
	this->bounds.push_back( std::pair<std::string, std::string>("plus/minus 1.3", "1.3" ) );
	this->bounds.push_back( std::pair<std::string, std::string>("plus/minus 1.5", "1.5" ) );
	this->bounds.push_back( std::pair<std::string, std::string>("plus/minus 1.7", "1.7" ) );
	this->bounds.push_back( std::pair<std::string, std::string>("plus/minus 2.0", "2.0" ) );
	this->bounds.push_back( std::pair<std::string, std::string>("plus/minus 2.5", "2.5" ) );
	this->bounds.push_back( std::pair<std::string, std::string>("plus/minus 3.0", "3.0" ) );

	this->bounds_per.push_back( std::pair<std::string, std::string>("0.3", "0.7" ) );
	this->bounds_per.push_back( std::pair<std::string, std::string>("0.4", "0.6" ) );

	wxArrayString ars_disc_bounds;
	for(int i=0; i<(int)bounds.size(); i++)
		ars_disc_bounds.Add(wxString::FromAscii(this->bounds.at(i).first.c_str()));

	this->cho_disc_bounds = new wxChoice( this, ID_CHOICE_DISC_BOUNDS, wxDefaultPosition, wxSize(120,22), ars_disc_bounds );
	this->cho_disc_bounds->SetSelection(0);
	this->cho_disc_bounds->Enable(false);

	wxBoxSizer* sizer_top = new wxBoxSizer(wxHORIZONTAL);
	sizer_top->Add(this->cho_disc, 0, FLAGS, BORDER_PXL);
	sizer_top->Add(this->cho_disc_bounds, 0, FLAGS, BORDER_PXL);
	if(investigation->is_mac){
		// mac uses giant fonts by default; PC is better.
		wxFont font(this->cho_disc_bounds->GetFont());
		font.SetPointSize(wxSMALL_FONT->GetPointSize());
		this->cho_disc_bounds->SetFont( font );
		this->cho_disc->SetFont( font );
	}
	Connect(ID_CHOICE_DISC, wxEVT_COMMAND_CHOICE_SELECTED, wxCommandEventHandler(PanelDiscretization::OnChooseDisc));

	this->SetSizer(sizer_top);
	sizer_top->Fit(this);
}


void PanelDiscretization::OnChooseDisc(wxCommandEvent& evt){
	redraw();
}

int PanelDiscretization::get_discretization_code(){
	return this->cho_disc->GetCurrentSelection();	
}

void PanelDiscretization::get_discretization_method_name(std::string& method){
    method = this->method_abbreviations.at( this->cho_disc->GetCurrentSelection() );
}

void PanelDiscretization::get_discretization_bounds(std::string& lower, std::string& upper){
	int disc = this->cho_disc->GetCurrentSelection();
	if( disc == DISC_NONE ){
		lower = std::string("");
		upper = std::string("");
	}
	else if( disc == DISC_PER ){
		lower = this->bounds_per.at( this->cho_disc_bounds->GetCurrentSelection() ).first;
		upper = this->bounds_per.at( this->cho_disc_bounds->GetCurrentSelection() ).second;
	}
    else if(disc == DISC_ABS ){
        std::string neg("-");
        lower = neg + this->bounds.at( this->cho_disc_bounds->GetCurrentSelection() ).second;
        upper = this->bounds.at( this->cho_disc_bounds->GetCurrentSelection() ).second;
    }
	else{
		lower = this->bounds.at( this->cho_disc_bounds->GetCurrentSelection() ).second;
		upper = this->bounds.at( this->cho_disc_bounds->GetCurrentSelection() ).second;
	}
}


void PanelDiscretization::set_discretization(std::string& disc, std::string& lower, std::string& upper){
	
	this->cho_disc_bounds->Clear();
	int d_bounds = 0;
	this->cho_disc_bounds->Enable(false);
	int disc_code = DISC_NONE;
	alg::to_lower(disc);
	if( disc.compare("sd")==0 )
		disc_code = DISC_SD;
	else if(disc.compare("mad")==0)
		disc_code = DISC_MAD;
	else if(disc.compare("abs")==0)
		disc_code = DISC_ABS;
	else if(disc.compare("per")==0)
		disc_code = DISC_PER;
	else if(disc.compare("sd_samples")==0)
		disc_code = DISC_SD_SAMPLES;
	this->cho_disc->SetSelection(disc_code);
	float lower_f, upper_f;
	
	if( disc_code==DISC_SD || disc_code==DISC_SD_SAMPLES || disc_code==DISC_MAD || disc_code==DISC_ABS ){
		lower_f = boost::lexical_cast<float>(lower);
		upper_f = boost::lexical_cast<float>(upper);
		this->cho_disc_bounds->Enable(true);
		for(int i=0; i<(int)this->bounds.size();i++){
			this->cho_disc_bounds->Append( wxString::FromAscii(this->bounds.at(i).first.c_str() ) );
			if( boost::lexical_cast<float>( this->bounds.at(i).second ) ==lower_f )
				d_bounds = i;
		}
		this->cho_disc_bounds->SetSelection(d_bounds);
	}
	else if(disc_code==DISC_PER ){
		lower_f = boost::lexical_cast<float>(lower);
		upper_f = boost::lexical_cast<float>(upper);
		this->cho_disc_bounds->Enable(true);
		for(int i=0; i<(int)this->bounds_per.size(); i++){
			std::stringstream ss;
			ss << this->bounds_per.at(i).first << ", " << this->bounds_per.at(i).second;
			this->cho_disc_bounds->Append( wxString::FromAscii( ss.str().c_str() ) );
			if( boost::lexical_cast<float>( this->bounds_per.at(i).first) == lower_f && boost::lexical_cast<float>( this->bounds_per.at(i).second ) == upper_f )
				d_bounds = i;
		}
		this->cho_disc_bounds->SetSelection(d_bounds);
	}
}


void PanelDiscretization::redraw(){
	this->cho_disc_bounds->Enable(false);
	int d_sel = this->cho_disc->GetCurrentSelection();
	int d_bounds = this->cho_disc_bounds->GetCurrentSelection();
	this->cho_disc_bounds->Clear();

	if( d_sel==DISC_SD || d_sel==DISC_SD_SAMPLES || d_sel==DISC_MAD || d_sel==DISC_ABS ){
		this->cho_disc_bounds->Enable(true);
		for(int i=0; i<(int)this->bounds.size();i++)
			this->cho_disc_bounds->Append( wxString::FromAscii(this->bounds.at(i).first.c_str() ) );
		if(d_bounds==-1)
			d_bounds = 0;
		this->cho_disc_bounds->SetSelection(d_bounds);
	}
	else if(d_sel==DISC_PER ){
		this->cho_disc_bounds->Enable(true);
		for(int i=0; i<(int)this->bounds_per.size();i++){
			std::stringstream ss;
			ss << this->bounds_per.at(i).first << ", " << this->bounds_per.at(i).second ;
			this->cho_disc_bounds->Append( wxString::FromAscii( ss.str().c_str() ) );
		}
		if(d_bounds==-1 || d_bounds > (int)this->bounds_per.size() - 1)
			d_bounds = 0;
		this->cho_disc_bounds->SetSelection(d_bounds);
	}
	wxCommandEvent evt( wxEVT_COMMAND_CHOICE_SELECTED, this->GetId() );
	evt.SetEventObject(this);
	GetEventHandler()->ProcessEvent(evt);
	this->parent->Refresh();
	this->parent->Update();
}
