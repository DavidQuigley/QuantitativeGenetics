#define _CRT_SECURE_NO_DEPRECATE 1
#define _SCL_SECURE_NO_DEPRECATE 1
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <wx/artprov.h>
#include <wx/treectrl.h>
#include <wx/dialog.h>
#include <wx/sizer.h>
#include <wx/imaglist.h>
#include <wx/stattext.h>
#include <wx/button.h>
#include <wx/textctrl.h>
#include <wx/wx.h>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/random.hpp>
using namespace boost;
namespace alg = boost::algorithm;
#include "carmen/src/DataStructures.h"
#include "carmen/src/Attributes.h"
#include "carmen/src/Rawdata.h"
#include "carmen/src/ClassMinerOptions.h"
#include "carmen/src/Discretize.h"
#include "carmen/src/Dataset.h"
#include "carmen/src/chi2.h"
#include "carmen/src/Rule.h"
#include "carmen/src/Parser.h"
#include "Investigation.h"
#include "DlgMergeRules.h"


IMPLEMENT_CLASS( MergeRulesDialog, wxDialog )

MergeRulesDialog::MergeRulesDialog( Investigation* investigation, std::string filename ){
	this->cancel_load = false;
	this->saved_new_file = false;
	this->investigation = investigation;
	SetExtraStyle(wxWS_EX_BLOCK_EVENTS|wxDIALOG_EX_CONTEXTHELP);
	wxDialog::Create( NULL, wxID_ANY, wxT("Merge Redundant Rules"), wxDefaultPosition, wxDefaultSize, wxCAPTION|wxRESIZE_BORDER|wxSYSTEM_MENU );
	if( !CreateControls(filename) )
		this->cancel_load = true;
	else{	
		GetSizer()->Fit(this); 
		GetSizer()->SetSizeHints(this);
		Centre();
	}
}


bool MergeRulesDialog::CreateControls(std::string filename){
	std::cout << "1\n";
	this->filename = filename;
	wxStaticText* lbl_threshold = new wxStaticText( this, wxID_ANY, wxT("Merge rules that are:") );
	this->txt_threshold = new wxTextCtrl( this, ID_TXT_THRESHOLD, wxT("100"), wxDefaultPosition, wxSize(50, 20) );
	wxButton* btn_refresh = new wxButton( this, ID_BTN_REFRESH, wxT("Merge") );
	wxStaticText* lbl_percent = new wxStaticText( this, wxID_ANY, wxT("% similar") );
	this->lbl_results = new wxStaticText( this, ID_LBL_RESULTS, wxT("Results:") );
	std::cout << "2\n";
	btn_refresh->SetDefault();
	std::cout << "3\n";
	this->ruleset = new RuleSet(filename);
	std::cout << "4\n";
	this->r_by_s = new Matrix<double>();
	ClassifierDataset* data = new ClassifierDataset();
	ClassMinerOptions * cmo = new ClassMinerOptions();
	cmo->file_name_dis = ruleset->file_name_dis;	
	cmo->file_name_sa = ruleset->file_name_sa;
	cmo->file_name_ga = ruleset->file_name_ga;
	cmo->discretization = ruleset->disc;
	cmo->disc_lower = (float)ruleset->disc_lower;
	cmo->disc_upper = (float)ruleset->disc_upper;
	cmo->gene_limit = ruleset->gene_limit;
	cmo->class_a = ruleset->class_a;
	cmo->class_b = ruleset->class_b;
	std::cout << "5\n";
	try{
		data->load(this->investigation->sa, this->investigation->ga, cmo);
	}
	catch( std::string msg ){
		std::stringstream ss;
		ss << "ERROR loading data:" << msg;
		wxMessageBox(wxString::FromAscii(ss.str().c_str()), _T("Error loading data"), wxOK | wxICON_EXCLAMATION, this);
		delete this->ruleset;
		delete this->r_by_s;
		delete data;
		delete cmo;
		
		return false;
	}
	std::cout << "Loaded\n";
	if( this->ruleset->rules.size() > 1000 ){
		if( wxMessageBox( wxString::FromAscii("Warning: large rulesets may take a long time to calculate.  Continue?"), wxString::FromAscii("Warning: lengthy calculation"), wxYES_NO | wxICON_QUESTION) == wxNO ){
			delete this->ruleset;
			delete this->r_by_s;
			delete data;
			delete cmo;
			return false;
		}
	}

	this->ruleset->build_rules_by_samples(this->r_by_s, data);  
	this->tree = new wxTreeCtrl(this, ID_TXT_MERGE_RULE_TREE, wxPoint(-1, -1), wxSize(500,500), wxTR_HAS_BUTTONS);
	this->imagelist.Create(16,16);
	this->img_file_normal = this->imagelist.Add(wxArtProvider::GetBitmap(wxART_NORMAL_FILE, wxART_OTHER, wxSize(16,16))); //wxART_NORMAL_FILE
	this->img_folder = this->imagelist.Add(wxArtProvider::GetBitmap(wxART_FOLDER, wxART_OTHER, wxSize(16,16))); //wxART_FOLDER
	this->img_file_open = this->imagelist.Add(wxArtProvider::GetBitmap(wxART_FILE_OPEN, wxART_OTHER, wxSize(16,16))); //wxART_FILE_OPEN
	this->tree->SetImageList( &this->imagelist );
	this->root = this->tree->AddRoot(_T("Merged Rules"));
	this->tree->SetItemImage(this->root, this->img_folder, wxTreeItemIcon_Normal);
	this->tree->SetItemImage(this->root, this->img_file_open, wxTreeItemIcon_Expanded);
	this->merge_rules(1);
	this->redraw();
	this->tree->Expand(this->root);
	wxButton* btn_ok = new wxButton( this, wxID_OK, wxT("&Close"));
	wxButton* btn_save_merged = new wxButton( this, ID_BTN_SAVE_MERGED, wxT("Save Merged") );
	wxBoxSizer* sizer_threshold = new wxBoxSizer(wxHORIZONTAL);
	sizer_threshold->Add(lbl_threshold, 0, FLAGS, BORDER_PXL);
	sizer_threshold->Add(this->txt_threshold, 0, FLAGS, BORDER_PXL);
	sizer_threshold->Add(lbl_percent, 0, FLAGS, BORDER_PXL);
	sizer_threshold->Add(btn_refresh, 0, FLAGS, BORDER_PXL);

	wxBoxSizer* sizer_top = new wxBoxSizer(wxVERTICAL);
	sizer_top->Add(sizer_threshold, 0, FLAGS, BORDER_PXL);
	sizer_top->Add(this->lbl_results, 0, FLAGS, BORDER_PXL);
	sizer_top->Add(this->tree, 0, FLAGS, BORDER_PXL);
	wxBoxSizer* sizer_bottom = new wxBoxSizer(wxHORIZONTAL);
	sizer_bottom->Add(btn_save_merged, 0, FLAGS, BORDER_PXL);
	sizer_bottom->Add(btn_ok, 0, FLAGS_RIGHT, BORDER_PXL);
	sizer_top->Add(sizer_bottom, 0, FLAGS, BORDER_PXL);
	Connect(wxID_OK, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(MergeRulesDialog::OnClickOk));
	Connect(ID_BTN_REFRESH, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(MergeRulesDialog::OnClickRefresh));
	Connect(ID_BTN_SAVE_MERGED, wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(MergeRulesDialog::OnClickSaveMerged));
	this->SetSizer(sizer_top);

	delete data;
	delete cmo;

	return true;
}


void MergeRulesDialog::redraw(){
	
	this->tree->DeleteChildren(this->root);
	int n_post_merge = 0;
	std::string label;
	int idx_merged;
	for(int i=0; i<(int)this->rules.size(); i++){	
		if( this->rules.at(i)->size()> 0 ){
			label = this->ruleset->rules.at(i)->get_friendly_label();
			wxTreeItemId r = this->tree->AppendItem( this->root, wxString::FromAscii(label.c_str()), this->img_folder, this->img_folder);
			for( int j=1; j<(int)this->rules.at(i)->size(); j++ ){
				idx_merged = this->rules.at(i)->at(j);
				label = this->ruleset->rules.at(idx_merged)->get_friendly_label();
				this->tree->AppendItem( r, wxString::FromAscii(label.c_str()), this->img_folder, this->img_folder);
			}
			n_post_merge += 1;
		}
	}
	std::stringstream ss;
	ss << "Originally ";
	ss << this->rules.size();
	ss << " rules, merged to ";
	ss << n_post_merge;
	ss << ".";
	this->lbl_results->SetLabel( wxString::FromAscii( ss.str().c_str() ) );
	this->tree->Expand(this->root);
}


bool comp_rules(const int* a, const int* b){
	return a[2] < b[2];
}



void MergeRulesDialog::merge_rules( double threshold ){
	this->threshold = threshold;
	wxBeginBusyCursor();
	int n_rules = r_by_s->rows();
	int n_samples = r_by_s->cols();
	enum {I=0, J, DIST};
	std::vector<int*> rule_sorter;
	rules.clear();
	std::vector<int> valid_rows, valid_cols;
	std::vector<int>* rule_idx_source;
	int i,j,r,n_diff, ctr=0;
	for(i=0; i<n_rules; i++){
		rule_idx_source = new std::vector<int>();
		rule_idx_source->push_back(i);
		this->rules.push_back( rule_idx_source );
		valid_rows.push_back(i);
		valid_cols.push_back(i);
	}
	Matrix<int> dist(n_rules, n_rules);
	int* rule;
	double one_minus_threshold = 1.0 - threshold;
	for( i=0; i<n_rules; i++){
		for( j=1; j<n_rules; j++){
			if( i<j ){
                n_diff = 0;
				for( r=0;r<n_samples; r++){
                    if( r_by_s->arr[i][r] != r_by_s->arr[j][r] )
                        n_diff += 1;
				}
				dist.arr[i][j]=n_diff;
				rule = new int[3];
				rule[I] = i; rule[J] = j; rule[DIST] = n_diff;
				rule_sorter.push_back(rule);
				ctr += 1;
			}
		}
	}
	std::sort(rule_sorter.begin(), rule_sorter.end(), comp_rules);	
	int n_sorted_rules = (int)rule_sorter.size();
	bool merged = true;
	int min_dist=0, row=0, col=0, min_row=0, min_col=0, min_row_idx=0, min_col_idx=0, idx_source=0, idx_dest=0;
	Rule* rule_min_row;
	Rule* rule_min_col;
	
	int idx_lowest=0;
	while(merged){
		merged = false;
		for(i=idx_lowest; i<n_sorted_rules; i++){
			row = rule_sorter.at(i)[I];
			col = rule_sorter.at(i)[J];
			if( valid_rows.at(row)>-1 && valid_cols.at(col)>-1){
				min_row = row;
				min_col = col;
				min_dist = dist.arr[row][col];
				min_row_idx  = row;
				min_col_idx = col;
				idx_lowest = i;
				break;
			}
		}
		
		if( (double)min_dist / (double)n_samples <= one_minus_threshold ){
			rule_min_row = ruleset->rules.at(min_row);
			rule_min_col = ruleset->rules.at(min_col);
			if( rule_min_row->chi2 < rule_min_col->chi2){
				idx_source = min_col;
				idx_dest = min_row;
				valid_rows.at(min_col_idx) = -1;
				valid_cols.at(min_col_idx) = -1;
			}
			else{
				idx_source = min_row;
				idx_dest = min_col;
				valid_rows.at(min_row_idx) = -1;
				valid_cols.at(min_row_idx) = -1;
			}
			rule_idx_source = rules.at(idx_source);
			while( rule_idx_source->size()>0 ){
				rules.at(idx_dest)->push_back( rules.at(idx_source)->at( rules.at(idx_source)->size()-1 ) );
				rules.at(idx_source)->pop_back();
			}

			merged = true;
			min_dist = n_samples+1;
		}
	}

	for(i=0; i<n_sorted_rules; i++)
		delete [] rule_sorter.at(i);
	wxEndBusyCursor();
}


void MergeRulesDialog::OnClickRefresh( wxCommandEvent& WXUNUSED(event) ){
	std::string thresh( this->txt_threshold->GetValue().ToAscii() );
	double threshold;
	try{
		threshold = boost::lexical_cast<double>(thresh);
	}
	catch( ... ){
		wxMessageBox(wxString::FromAscii("Threshold must be a value between 0 and 100"), _T("Not Ready"), wxOK | wxICON_EXCLAMATION, this);
		return;
	}
	if( threshold > 100 || threshold < 0 ){
		wxMessageBox(wxString::FromAscii("Threshold must be a value between 0 and 100"), _T("Not Ready"), wxOK | wxICON_EXCLAMATION, this);
		return;
	}
	threshold = threshold / 100.0;
	this->merge_rules(threshold);
	this->redraw();
}


void MergeRulesDialog::OnClickSaveMerged(wxCommandEvent& evt){
	// save foo.ruleset as foo.merged.ruleset
	// should note in header that ruleset is the product of a merge
	//this->ruleset->write(modified_filename);
	std::stringstream file_name;
	std::vector<std::string> file_name_parts;
	alg::split( file_name_parts, this->filename, alg::is_any_of(".") );
	for(int i=0; i<(int)file_name_parts.size()-1; i++){
		file_name << file_name_parts.at(i);
		if(i<(int)file_name_parts.size()-2)
			file_name << ".";
	}
	file_name << ".merged.ruleset";
	std::string fn = file_name.str();
	wxString fileloc = wxGetTextFromUser(wxString::FromAscii("Merged Ruleset file location:"), wxString::FromAscii("Choose file location"), wxString::FromAscii(fn.c_str() ), this);
	if( fileloc.size()==0 )
		return;
	std::ofstream f_out( fileloc.ToAscii() );
	if( !f_out.is_open() ){
		std::stringstream ss;
		ss << "Unable to open file for writing: " << fileloc.ToAscii();
		wxMessageBox(wxString::FromAscii(ss.str().c_str() ), wxString::FromAscii("Error writing rules"));
		return;
	}
	int n_top_level = 0;
	for(int i=0; i<(int)this->rules.size(); i++){	
		if( this->rules.at(i)->size()> 0 ){
			n_top_level += 1;
		}
	}
	this->ruleset->write_header(f_out, n_top_level, this->threshold*100);
	for(int i=0; i<(int)this->rules.size(); i++){	
		if( this->rules.at(i)->size()> 0 ){
			this->ruleset->write_rule(f_out, i);
		}
	}
	f_out.close();
    this->new_filename = fileloc.ToAscii();
	this->saved_new_file = true;
}


void MergeRulesDialog::OnClickOk( wxCommandEvent& WXUNUSED(event) ){
	
	delete this->ruleset;
	delete this->r_by_s;
	for(int i=0; i<(int)this->rules.size(); i++)
		delete this->rules[i];
	this->AcceptAndClose();
}
