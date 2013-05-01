#define _CRT_SECURE_NO_DEPRECATE 1
#include <stdio.h>
#include <string>
#include <algorithm>
#include <iostream>
#include <math.h>
#include <fstream>
#include <vector>
using namespace std;
#include <wx/wx.h>
#include <wx/grid.h>
#include <wx/choicdlg.h>
#include <wx/button.h>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/random.hpp>
using namespace boost;
namespace alg = boost::algorithm;
#include "carmen/src/DataStructures.h"
#include "carmen/src/Attributes.h"
#include "carmen/src/Rawdata.h"
#include "carmen/src/Parser.h"
#include "carmen/src/graph.h"
#include "Investigation.h"
#include "GridClassifierResults.h"
#include "DlgClassificationViewer.h"

IMPLEMENT_CLASS( ClassificationViewerDialog, wxDialog )

ClassificationViewerDialog::ClassificationViewerDialog( std::string filename, Investigation* investigation){
	this->investigation = investigation;
	SetExtraStyle(wxWS_EX_BLOCK_EVENTS|wxDIALOG_EX_CONTEXTHELP);
	wxDialog::Create( NULL, wxID_ANY, wxT("Classification Performance"), wxDefaultPosition, wxSize(400,350), wxCAPTION|wxRESIZE_BORDER|wxSYSTEM_MENU );
	CreateControls(filename);
	GetSizer()->Fit(this); 
	GetSizer()->SetSizeHints(this); 
	Centre();
}

void ClassificationViewerDialog::CreateControls(std::string filename){
    this->grid_classification = new GridClassifierResults(this, filename, ID_CLASS_V_GRID_CORR, wxSize(300,300));
    this->btn_close = new wxButton( this, wxID_OK, wxT("Close") );
    this->btn_close->SetDefault();

    wxBoxSizer* sizer_top = new wxBoxSizer(wxVERTICAL);
    sizer_top->Add(this->grid_classification, 0, FLAGS, BORDER_PXL);
    sizer_top->Add(this->btn_close, 0, FLAGS, BORDER_PXL);
    SetSizer(sizer_top);
}
