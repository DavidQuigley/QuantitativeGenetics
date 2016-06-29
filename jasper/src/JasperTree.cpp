#define _CRT_SECURE_NO_DEPRECATE 1
#define _SCL_SECURE_NO_DEPRECATE 1
#include <string>
#include <vector>
#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "wx/treectrl.h"
#include "wx/artprov.h"
#include "wx/imaglist.h"
#include "wx/msgdlg.h"
#include "wx/log.h"
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>
#include <boost/random.hpp>
using namespace boost;
namespace alg = boost::algorithm;
namespace fs = boost::filesystem;
namespace alg = boost::algorithm;
using namespace std;
#include "carmen/src/DataStructures.h"
#include "carmen/src/Attributes.h"
#include "carmen/src/Rawdata.h"
#include "carmen/src/Parser.h"
#include "Investigation.h"
#include "JasperTree.h"


JasperTree::JasperTree(wxWindow* parent, wxWindowID id, Investigation* investigation, wxSize size):
	wxTreeCtrl(parent, id, wxPoint(-1, -1), size, wxTR_HAS_BUTTONS){
	this->parent = parent;
	this->cp = investigation->cp;
	this->investigation = investigation;
	this->imagelist.Create(16,16);
	this->img_file_normal = this->imagelist.Add(wxArtProvider::GetBitmap(wxART_NORMAL_FILE, wxART_OTHER, wxSize(16,16))); //wxART_NORMAL_FILE
	this->img_folder = this->imagelist.Add(wxArtProvider::GetBitmap(wxART_FOLDER, wxART_OTHER, wxSize(16,16))); //wxART_FOLDER
    this->img_file_open = this->imagelist.Add(wxArtProvider::GetBitmap(wxART_FILE_OPEN, wxART_OTHER, wxSize(16,16))); //wxART_FILE_OPEN
    if(this->investigation->is_mac){
    	// mac uses giant fonts by default; PC is better.
    	wxFont font(this->GetFont());
    	font.SetPointSize(wxSMALL_FONT->GetPointSize());
    	this->SetFont(font);
    }
    this->SetImageList( &this->imagelist );
	this->refresh_from_directory();
}

void JasperTree::get_selected_label(std::string& label){
	wxTreeItemId current = this->GetSelection();
	wxString wx_lbl;
	if( current.IsOk() ){
		wx_lbl = this->GetItemText(current);
	}
	else
		wx_lbl = wxEmptyString;
	label = std::string(wx_lbl.ToAscii());
}


void JasperTree::open_tree_to_file(std::string filename){
    wxTreeItemId cur_label, cur_file_type, cur_file;
    TreeData* d;
    this->CollapseAll();
    this->Expand(this->root);
    wxTreeItemIdValue token1, token2, token3;
	std::string f("/");
	std::string b("\\");
	boost::algorithm::replace_all(filename,b,f);
    //std::cout << "Seeking " << filename << "\n";
    cur_label = this->GetFirstChild( this->root, token1);
    while( cur_label.IsOk() ){
    	//std::cout << ((TreeData*)this->GetItemData(cur_label))->original_label << "\n";
        cur_file_type = this->GetFirstChild( cur_label, token2 );
        while( cur_file_type.IsOk() ){
        	//std::cout << "\t" << ((TreeData*)this->GetItemData(cur_file_type))->original_label << "\n";
            cur_file = this->GetFirstChild( cur_file_type, token3 );
            while( cur_file.IsOk() ){
            	d = (TreeData*)this->GetItemData(cur_file);
                //std::cout << "\t\tChecking " << d->str << "\n";
				std::string check(d->str);
				boost::algorithm::replace_all(check, b,f);
				if( check.compare(filename) == 0 ){
                	this->Expand(cur_label);
                    this->Expand(cur_file_type);
                    this->SelectItem(cur_file);
                    return;
                }
                cur_file = this->GetNextChild( cur_file_type, token3 );
            }
            cur_file_type = this->GetNextChild( cur_label, token2);
        }
        cur_label = this->GetNextChild( this->root, token1 );
    }
}


void JasperTree::get_selected_original_label(std::string& label){
	wxTreeItemId current = this->GetSelection();
	if( current.IsOk() ){
		TreeData* d = (TreeData*) this->GetItemData(current);
		label = d->original_label;
	}
	else
		label = "";
}

int JasperTree::get_selected_leaftype(){
	wxTreeItemId current = this->GetSelection();
	if( current.IsOk() ){
		TreeData* d = (TreeData*) this->GetItemData(current);
		return d->leaf_type;
	}
	return LEAF_ERROR;
}

void JasperTree::get_selected_filename(std::string& s){
	wxTreeItemId current = this->GetSelection();
	if( current.IsOk() ){
		TreeData* d = (TreeData*) this->GetItemData(current);
		if( d->leaf_type == LEAF_FILE)
			s = d->str;
		else
			s = "";
	}
	else
		s = "";
}

void JasperTree::get_limit_and_label(std::string filepath, std::string extension, std::string& limit, std::string& label){
	CarmenParser parser(filepath);
	std::stringstream lbl, file_name;
	std::vector<std::string> path_parts, file_name_parts;

	std::string df;
	if( filepath.find("/") != std::string::npos )
		alg::split( path_parts, filepath, alg::is_any_of("/") );
	else
		alg::split( path_parts, filepath, alg::is_any_of("\\") );
	alg::split( file_name_parts, path_parts.at(path_parts.size()-1), alg::is_any_of(".") );
	
	for(int i=0; i<(int)file_name_parts.size()-1; i++){
		file_name << file_name_parts.at(i);
		if(i<(int)file_name_parts.size()-2)
			file_name << ".";
	}

	std::string A("");
	std::string B("");
	std::string bang("!");
	std::string eq("=");
	std::string isnot(" is not ");
	std::string is(" is ");
	std::string vs(" vs. ");
	std::string comma(",");
	std::string slash("/");
	if( extension.compare("classifier")==0){
		parser.get_value("Class_A_train", 0, A);
		parser.get_value("Class_B_train", 0, B);
	}
    else if( extension.compare("labels")==0){
        std::string A_test, B_test;
        parser.get_value("Class_A_train", 0, A);
		parser.get_value("Class_B_train", 0, B);
        parser.get_value("Class_A_test", 0, A_test);
		parser.get_value("Class_B_test", 0, B_test);
        A = A + " > " + A_test;
        B = B + " > " + B_test;
    }
	else{
		parser.get_value("Class_A", 0, A);
		parser.get_value("Class_B", 0, B);
	}
	alg::replace_all(A, bang, isnot);
	alg::replace_all(A, eq, is);
	alg::replace_all(B, bang, isnot);
	alg::replace_all(B, eq, is);
	if(A.size()==0){
		limit = std::string(B);
	}
	else if(B.size()==0){
		limit = std::string(A);
	}
	else{
		std::stringstream ss;
		ss << A << vs << B;
		limit = ss.str();
	}
	
	if( extension.compare("ruleset")==0){    
		std::string sci_str, disc, mine_type, n_rules, max_depth;
		parser.get_value("Mine_Type", 0, mine_type);
		parser.get_value("Generated", 0, n_rules);
		parser.get_value("Max_Depth", 0, max_depth);
		if( parser.get_value("S_C_I", 0, sci_str) ){
			alg::replace_all(sci_str, comma, slash);
			parser.get_value("Discretization", 0, disc);
			if(disc.compare("None") != 0){
				std::string mult;
				parser.get_value("Discretization", 3, mult);
				if( disc.compare("per")==0)
					disc = disc + " +/- " + mult.substr(0, mult.find(comma)) + "/" + mult.substr( mult.find(comma)+1 );
				else
					disc = disc + " +/- " + mult.substr( mult.find(comma)+1 );
			}
			if( mine_type.compare("Core Rules")==0)
				lbl << sci_str << "  " << disc << "  (" << n_rules << " core rules, depth " << max_depth << ")";
			else
				lbl << sci_str << "  " << disc << "  (" << n_rules << " rules, depth " << max_depth << ")";
		}
		else{
			lbl << " " << mine_type << "  (" << n_rules << " rules, depth " << max_depth << ")";
		}
	}
	else if(extension.compare("classifier")==0){
        std::string method;
        std::string acc_all, acc_a, acc_b;
        parser.get_value("Method", 0, method);
		parser.get_value("Acc", 0, acc_all);
        parser.get_value("Acc", 1, acc_a);
        parser.get_value("Acc", 2, acc_b);
        
        lbl << method << " Acc: ";
        lbl << acc_all;
        lbl << " (";
        lbl << acc_a;
        lbl << " | " ;
        lbl << acc_b;
        lbl << ")";
		
	}
	else if(extension.compare("difference")==0){
		std::string n_perms, max_p_value;
		parser.get_value("N_Permutations", 0, n_perms);
		parser.get_value("Max_p_value", 0, max_p_value);
		lbl << n_perms << " perms, max p-value " << max_p_value;
	}
    else if( extension.compare("labels")==0 ){
		std::string method;
		std::string acc_all, acc_a, acc_b;
		parser.get_value("Method", 0, method);
		parser.get_value("Method", 0, method);
		parser.get_value("Acc", 3, acc_all);
        parser.get_value("Acc", 4, acc_a);
        parser.get_value("Acc", 5, acc_b);
        lbl << method << " Acc Tst: ";
        lbl << acc_all;
        lbl << " (";
        lbl << acc_a;
        lbl << " | " ;
        lbl << acc_b;
        lbl << ")";
    }
	else if(extension.compare("spear")==0){
		std::string Min_var, Delta_rho, Abs_rho, probe_id;
		parser.get_value("Min_var", 0, Min_var);
		parser.get_value("Delta_rho", 0, Delta_rho);
		parser.get_value("Abs_rho", 0, Abs_rho);
		parser.get_value("Focus_Probe", 0, probe_id);
		if(probe_id.size()>0 ){
			lbl << "Probe: " << probe_id << " ";
		}
		lbl << "min Var/abs(r)/change(r): " << Min_var << "/" << Abs_rho << "/" << Delta_rho;
	}
	lbl << ". File: " << file_name.str();
	label = std::string(lbl.str().c_str());
}


wxTreeItemId JasperTree::add_label(wxTreeItemId item_q, std::string filename, std::string extension, std::string label){
	TreeData* td = new TreeData(filename, label, LEAF_FILE);
	if( extension.compare("ruleset")==0 )
		return this->AppendItem( find_child_named(item_q, "Rulesets"), wxString::FromAscii(label.c_str()), this->img_folder, this->img_folder, td );
	else if(extension.compare("classifier")==0)
		return this->AppendItem( find_child_named(item_q, "Classifiers"), wxString::FromAscii(label.c_str()), this->img_folder, this->img_folder, new TreeData(filename, label, LEAF_FILE)  );
	else if( extension.compare("difference")==0 )
		return this->AppendItem( find_child_named(item_q, "Differences"), wxString::FromAscii(label.c_str()), this->img_folder, this->img_folder, new TreeData(filename, label, LEAF_FILE)  );
	else if(extension.compare("spear")==0 )
		return this->AppendItem( find_child_named(item_q, "Correlations"), wxString::FromAscii(label.c_str()), this->img_folder, this->img_folder, new TreeData(filename, label, LEAF_FILE)  );
	else if(extension.compare("labels")==0 )
		return this->AppendItem( find_child_named(item_q, "Classified Datasets"), wxString::FromAscii(label.c_str()), this->img_folder, this->img_folder, new TreeData(filename, label, LEAF_FILE)  );
	else
		delete td;
	wxTreeItemId dummy;
	return dummy;
}

wxTreeItemId JasperTree::create_query(std::string label_to_show, std::string original_label){
	
	wxTreeItemId new_q = append_child_with_icons( this->root, label_to_show, new TreeData("", original_label, LEAF_QUERY));
	//append_child_with_icons(new_q, "Rulesets", new TreeData("", "Rulesets", LEAF_SUBHEAD));
	//append_child_with_icons(new_q, "Classifiers", new TreeData("", "Classifiers", LEAF_SUBHEAD));
	append_child_with_icons(new_q, "Differences", new TreeData("", "Differences", LEAF_SUBHEAD));
	append_child_with_icons(new_q, "Correlations", new TreeData("", "Correlations", LEAF_SUBHEAD));
	//append_child_with_icons(new_q, "Classified Datasets", new TreeData("", "Classified Datasets", LEAF_SUBHEAD));
    return new_q;
}

wxTreeItemId JasperTree::append_child_with_icons(wxTreeItemId parent, std::string label, TreeData* treedata){
	wxTreeItemId added = this->AppendItem(parent, wxString::FromAscii(label.c_str()), -1, -1, treedata );
	this->SetItemImage(added, this->img_folder, wxTreeItemIcon_Normal);
	this->SetItemImage(added, this->img_file_open, wxTreeItemIcon_Expanded);
	return added;
}


wxTreeItemId JasperTree::find_child_named(wxTreeItemId item_top, std::string query){
	wxTreeItemIdValue cookie;
	wxTreeItemId item = this->GetFirstChild( item_top, cookie );
	while( item.IsOk() ){
		wxString sData = this->GetItemText(item);
		if( sData.Cmp( wxString::FromAscii(query.c_str()) ) == 0 ){
			return item;
		}
		item = this->GetNextChild( item_top, cookie);
	}
	wxTreeItemId dummy;
	return dummy;
}


TreeData::TreeData(std::string str, std::string original_label, int leaf_type){
	this->str = str;
	this->leaf_type = leaf_type;
	this->original_label = original_label;
}

void JasperTree::refresh_from_directory(){
	// read the investigation directory ("dir" parameter) and create leaves in the
	// tree for the following file extensions:
	// .ruleset,  .difference,  .classifier,  .spear, .labels
	std::string investigation_name = this->cp->get("History", "most_recent_fileset");
	if( this->cp->get(std::string("Internal"), std::string("debug") ).compare(std::string("true")) != 0 )
		std::cout << "Reading investigation: " << investigation_name << "\n";
	if( investigation_name.size()==0 ){
		this->DeleteAllItems(); // special case: no investigations, so no leaves
		return;
	}
	try{
		// The next call will load the sample attributes and gene attributes data
		this->investigation->set_current( investigation_name );
	}
	catch(std::string errmsg){
		wxMessageBox(wxString::FromAscii(errmsg.c_str()), wxString::FromAscii("Error in properties file"), wxOK | wxICON_EXCLAMATION, this);
		this->investigation->investigation_name = std::string(""); // this will disable most features of Carmen
		if( wxIsBusy() )
			wxEndBusyCursor();
		return;
	}
	if( this->root.IsOk() ){
		this->DeleteChildren(this->root);
		this->SetItemText(this->root, wxString::FromAscii(investigation_name.c_str()));
	}
	else
		this->root = this->AddRoot(wxString::FromAscii(investigation_name.c_str()));
	this->SetItemImage(this->root, this->img_folder, wxTreeItemIcon_Normal);
	this->SetItemImage(this->root, this->img_file_open, wxTreeItemIcon_Expanded);
	this->SetItemData(this->root, new TreeData(investigation_name, investigation_name, LEAF_ROOT) );

	this->SetItemBold(this->root);
	if( investigation_name.size()==0 )
		return;
	std::string file_name, extension, label, query;
	
	std::string posix(investigation->dir_results.c_str());
	for(int i=0; i<(int)posix.size(); i++){
		if( posix.at(i)== '\\'){
			posix.at(i) = '/';
		}
	}
	if( this->cp->get(std::string("Internal"), std::string("debug") ).compare(std::string("true")) != 0 )
		std::cout <<  "From directory " << investigation->dir_results << "\n";
	std::vector<std::string> file_name_parts;
	boost::filesystem::path results(posix.c_str());
	
	if( fs::exists( results ) ){
		fs::directory_iterator end;
		for( fs::directory_iterator iter(results) ; iter != end ; ++iter ){
			if ( !fs::is_directory( iter->status()) ){
				file_name = iter->path().filename().string();
				alg::split( file_name_parts, file_name, alg::is_any_of(".") );
				extension = file_name_parts.at( file_name_parts.size()-1 );				
				if( extension.compare("ruleset")==0 || 
					extension.compare("classifier")==0 ||
					extension.compare("difference")==0 ||
					extension.compare("spear")==0 || 
                    extension.compare("labels")==0 ){
					try{
						get_limit_and_label(iter->path().string(), extension, query, label);
						std::string nickname(query);
						std::string query_LC( query );
						boost::algorithm::to_lower(query_LC);
						if( investigation->nicknames.find(query_LC) != investigation->nicknames.end() )
							nickname = investigation->nicknames[query_LC];
						
						wxTreeItemId id_query = find_child_named(this->root, nickname); 
						if( !id_query.IsOk() ){
							id_query = create_query(nickname, query);
						}
						add_label(id_query, iter->path().string(), extension, label);
					}
					catch( std::string errormessage){
						std::stringstream ss;
						ss << "Error parsing file at\n" << iter->path().filename() << "\n" << errormessage;
						wxMessageBox(wxString::FromAscii(errormessage.c_str()), wxString::FromAscii(ss.str().c_str()) , wxOK | wxICON_EXCLAMATION, this);
					}
					catch( ... ){
						std::stringstream ss;
						ss << "Unhandled error reading file " << file_name;
						wxMessageBox(wxString::FromAscii(ss.str().c_str()), wxString::FromAscii(ss.str().c_str()) , wxOK | wxICON_EXCLAMATION, this);
					}
				}		
			}
		}	
	}
	this->Expand(this->root);
}
