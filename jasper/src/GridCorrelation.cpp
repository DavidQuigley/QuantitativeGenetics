#define _CRT_SECURE_NO_DEPRECATE 1
#include <stdio.h>
#include <string>
#include <algorithm>
#include <iostream>
#include <math.h>
#include <fstream>
#include <vector>
#include <set>
using namespace std;
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/random.hpp>
using namespace boost;
namespace alg = boost::algorithm;
#include <wx/wx.h>
#include <wx/grid.h>
#include <wx/choicdlg.h>
#include "carmen/src/DataStructures.h"
#include "carmen/src/Attributes.h"
#include "carmen/src/Rawdata.h"
#include "carmen/src/graph.h"
#include "carmen/src/Parser.h"
#include "Investigation.h"
#include "GridCorrelation.h"


GridCorrelation::GridCorrelation(wxWindow* parent, Investigation* investigation, std::string filename, wxWindowID id, wxSize size) : 
	wxGrid(parent, id, wxPoint(-1,-1), size){
	this->m_parent = parent;
	this->investigation = investigation;
	//this->read_file_into_graph(filename);
	this->G.read(filename);
	std::string gn("Gene Name");
	for(HASH_S_I::iterator itr=investigation->ga->identifier2idx.begin(); itr != investigation->ga->identifier2idx.end(); itr++){
		this->id2p[(*itr).second] = (*itr).first;
		this->p2g[ (*itr).first ] = investigation->ga->prop_for_identifier( (*itr).first, gn );
	}
	this->CreateGrid( 0, 0 );
	this->SetRowLabelSize(0);
    this->SetColLabelSize(15);
	this->Clear();
	this->AppendCols(3);
	this->SetColLabelValue(0, wxString::FromAscii("Name"));
	this->SetColLabelValue(1, wxString::FromAscii("ID"));
	this->SetColLabelValue(2, wxString::FromAscii("r-value"));
	this->read_file_into_graph(filename);
	this->filename = filename;
    this->EnableEditing(false);
	this->Fit();
	this->m_parent->FitInside();
}


void GridCorrelation::Clear(){
	int n_rows = this->GetNumberRows();
	if( n_rows>0 )
		this->DeleteRows(0, n_rows);
	this->Layout();
}

void GridCorrelation::read_file_into_graph(std::string filename){
	wxBeginBusyCursor();

	CarmenParser cp(filename);
	cp.PrepareToReadValues();
	std::vector<std::string> values;
	double rho;
	Attributes* ga = this->investigation->ga;
	for(HASH_S_I::iterator itr=ga->identifier2idx.begin(); itr != ga->identifier2idx.end(); itr++){
		this->id2p[(*itr).second] = (*itr).first;		
	}
	while( cp.ReadNextValue(values) ){
		this->p2g[ values.at(1) ] = values.at(0);
		this->p2g[ values.at(3) ] = values.at(2);
		rho = boost::lexical_cast<double>(values.at(4));
		this->G.add_edge( ga->identifier2idx[values.at(1)], ga->identifier2idx[values.at(3)], rho );
	}
	wxEndBusyCursor();
}


bool GridCorrelation::write_to_file(std::string file_path){
    std::ofstream f_out(file_path.c_str());
	if( !f_out.is_open() ){
        wxMessageBox(wxString::FromAscii("ERROR: unable to open file for writing."), wxString::FromAscii("Error writing file"), wxOK|wxICON_EXCLAMATION, this);
        return false;
    }
    f_out << "Correlation with probe " << this->current_probe_id.c_str() << "\n";
    f_out << "Gene\tProbe\tRho\n";
    int n_col = this->GetNumberCols();
    for(int i=0; i<this->GetNumberRows(); i++){
        for(int j=0; j<n_col; j++){
            f_out << this->GetCellValue(i, j).ToAscii();
            if( j<n_col-1 )
                f_out << "\t";
        }
        f_out << "\n";
    }
    f_out.close();
    return true;
}


bool weight_sort(const std::pair<double,int> *a, const std::pair<double,int> *b){
	return a->first > b->first;
}

void GridCorrelation::Update(std::string gene_name, double min_abs_corr){
	std::string probe_id;
	this->Clear();
	this->current_probe_id = "";
	if( this->p2g.find(gene_name) != this->p2g.end() ){
		probe_id = gene_name;
	}
	else{
		// not a probe id; look in genes
		Attributes* ga = this->investigation->ga;
		std::vector<std::string> probes;
		std::string gn("Gene Name");
		alg::to_lower(gene_name);
		for(int i=0; i<(int)ga->identifiers.size(); i++){
			std::string g = ga->prop_for_identifier(ga->identifiers.at(i), gn);
			alg::to_lower(g);
			if( g.compare(gene_name)==0 ){
				probes.push_back(ga->identifiers.at(i));
			}
		}
		if( probes.size()==0 ){
			wxMessageBox(wxString::FromAscii("Gene not found in experiment"), wxString::FromAscii("Gene not found"), wxOK|wxICON_EXCLAMATION, this);
			return;
		}
		else if((int)probes.size() == 1){
			probe_id = probes.at(0);
		}
		else{
			wxArrayString aChoices;
			for(int i=0; i<(int)probes.size(); i++){
				aChoices.Add( wxString::FromAscii(probes.at(i).c_str()) );
			}
			wxString probe_wx = wxGetSingleChoice(wxString::FromAscii("Choose the probe:"),wxString::FromAscii("Choose probe"), aChoices );
			probe_id = probe_wx.ToAscii();
		}
	}
	if( this->investigation->ga->identifier2idx.find(probe_id) != this->investigation->ga->identifier2idx.end() ){
		int root_id = this->investigation->ga->identifier2idx[probe_id];
		std::vector<int> ids;
		std::vector<double> weights;
		this->G.neighbors(root_id, ids, weights);
		
		char rho[10];
		std::vector< std::pair<double, int>* > returned;
		std::pair<double,int>* p;
		for(int i=0;i<(int)weights.size(); i++){
			if( abs(weights.at(i)) >= min_abs_corr ){
				p = new std::pair<double, int>();
				p->first = weights.at(i);
				p->second = ids.at(i);
				returned.push_back( p );
			}
		}
		std::sort( returned.begin(), returned.end(), weight_sort );
		this->AppendRows((int)returned.size());
		for(int i=0; i<(int)returned.size(); i++){
			sprintf(rho, "%2.3f", returned.at(i)->first );
			this->SetCellValue(i, 0, wxString::FromAscii( (p2g[ this->id2p[ returned.at(i)->second ] ]).c_str() ) );
			this->SetCellValue(i, 1, wxString::FromAscii( (this->id2p[ returned.at(i)->second ]).c_str() ) );
			this->SetCellValue(i, 2, wxString::FromAscii( rho )  );
			delete returned.at(i);
		}
		
		if( (int)returned.size()==0 ){
			wxMessageBox(wxString::FromAscii("No probes are correlated at that stringency."), wxString::FromAscii("No probes found"), wxOK|wxICON_INFORMATION, this);
		}
		this->current_probe_id = probe_id;
	}

	this->SetColSize(0, 150);
	this->AutoSizeColumn(1);
	this->AutoSizeColumn(2);
	this->ForceRefresh();
	this->UpdateWindowUI();

}
		
