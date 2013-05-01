#define _CRT_SECURE_NO_DEPRECATE 1
#define _SCL_SECURE_NO_DEPRECATE 1
#include <string>
#include <vector>
#include <fstream>
#include <set>
using namespace std;
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/random.hpp>
namespace alg = boost::algorithm;
#include <wx/wx.h>
#include "carmen/src/DataStructures.h"
#include "carmen/src/Attributes.h"
#include "carmen/src/Rawdata.h"
#include "carmen/src/Parser.h"
#include "Investigation.h"

Investigation::Investigation(){
	this->dir_results = std::string("");
	this->fn_expr = std::string("");
	this->fn_ga = std::string("");
	this->fn_sa = std::string("");
	this->species = std::string("");
	this->data_type = std::string("");
	this->gene_name_column = std::string("Gene Name");
	this->cp = NULL;
    this->sa = NULL;
    this->ga = NULL;
}


void Investigation::clone_from(Investigation *source){
    this->dir_results = source->dir_results;
    this->fn_expr = source->fn_expr;
    this->fn_ga = source->fn_ga;
    this->fn_sa = source->fn_sa;
    this->data_type = source->data_type;
    this->species = source->species;
    this->gene_name_column = source->gene_name_column;
    if( this->sa != NULL )
        delete this->sa;
    if(this->ga != NULL )
        delete this->ga;
    this->sa = new Attributes("NA");
	this->ga = new Attributes("NA");
    this->fn_properties = source->fn_properties;
    if( (int)this->fn_properties.size() > 0 )
        read_from_file(this->fn_properties);
    this->investigation_name = source->investigation_name;
    if( (int)this->investigation_name.size()>0)
        set_current(this->investigation_name);
    this->is_mac = source->is_mac;
}


Investigation::Investigation(std::string fn_properties){
	this->dir_results = std::string("");
	this->fn_expr = std::string("");
	this->fn_ga = std::string("");
	this->fn_sa = std::string("");
	this->species = std::string("");
	this->data_type = std::string("");
	this->gene_name_column = std::string("Gene Name");
	this->cp = NULL;
	this->sa = new Attributes("NA");
	this->ga = new Attributes("NA");
#ifdef WIN32
	this->is_mac = false;
#else
	this->is_mac = true;
#endif
	this->read_from_file(fn_properties);
}


void Investigation::read_from_file(std::string fn_properties){
	if( this->cp == NULL )
		this->cp = new ConfigParser(fn_properties);
	else{
		delete this->cp;
		this->cp = new ConfigParser(fn_properties);
	}
#ifdef WIN32
	this->is_mac = false;
#else
	this->is_mac = true;
#endif
	this->fn_properties = fn_properties;
}

void Investigation::set_current(std::string investigation_name){
	
	if( investigation_name.size()==0)
		throw string("No investigations found in properties file.");
	std::string names = this->cp->get("Filesets", "names");
	alg::to_lower(names);
	std::string i_name_lower( investigation_name );
	alg::to_lower(i_name_lower);
	if( !boost::algorithm::find_first(names, i_name_lower) )
		throw string("ERROR: Requested investigation not found in Filesets->names");

	this->investigation_name = i_name_lower;

	this->dir_results = cp->get("Filesets", i_name_lower + "<>dir");
	this->fn_expr = cp->get("Filesets", i_name_lower + "<>raw");
	this->fn_sa = cp->get("Filesets", i_name_lower + "<>sa");
	this->fn_ga = cp->get("Filesets", i_name_lower + "<>ga");
	this->species = cp->get("Filesets", i_name_lower + "<>species");
	this->data_type = cp->get("Filesets", i_name_lower + "<>data_type");
	this->gene_name_column = cp->get("Filesets", i_name_lower + "<>gene_name_column");
	if( this->gene_name_column.size()==0 )
		this->gene_name_column = std::string("Gene Name");
	if(this->dir_results.size()==0){ throw string("Missing dir element in properties file for investigation " + this->investigation_name); }
	if(this->fn_expr.size()==0){ throw string("Missing raw element in properties file for investigation " + this->investigation_name); }
	if(this->fn_sa.size()==0){ throw string("Missing sa element in properties file for investigation " + this->investigation_name); }
	if(this->fn_ga.size()==0){ throw string("Missing ga element in properties file for investigation " + this->investigation_name); }
	cp->set("History", "most_recent_fileset", investigation_name);
	delete this->sa;
	delete this->ga;
	this->sa = new Attributes("NA");
	this->ga = new Attributes("NA");
	wxBeginBusyCursor();
    try{ this->sa->load(this->fn_sa); }
    catch( ... ){
        wxEndBusyCursor();
		std::stringstream ss;
		ss << "ERROR: Error reading sample attributes file, attempting to load\n" << this->fn_sa; 
		throw ss.str();
    }
    try{
    	this->ga->load(this->fn_ga);
    }
    catch( std::string err){
    	std::stringstream ss;
    	ss << "ERROR: Error reading gene attributes file, attempting to load\n" << this->fn_ga << "\n" << err;
    	throw ss.str();
    }
    try{
    	if( this->gene_name_column.size()>0 )
    		this->ga->set_gene_name_column( this->gene_name_column);
    }
    catch( std::string err ){
    	// attempt to rescue wayward symbol column before blowing up
    	if( this->ga->attribute2idx.find( std::string("Gene Name") ) != this->ga->attribute2idx.end() ){
    		this->ga->set_gene_name_column( std::string("Gene Name") );
    		this->gene_name_column = std::string("Gene Name");
    	}
    	else if( this->ga->attribute2idx.find( std::string("symbol") ) != this->ga->attribute2idx.end() ){
    		this->ga->set_gene_name_column( std::string("symbol") );
    		this->gene_name_column = std::string("symbol");
    	}
    	else{
    		throw std::string("Specified gene name column " + this->gene_name_column + " was not found in gene attributes file");
    	}
    }
    wxEndBusyCursor();
	int idx = 0;
	std::vector<prop_pair*> V;
	std::string key, value, name_cur;
	std::string colon_word("<colon>");
	std::string colon(":");
    std::string angle("<");

	cp->get_section("Nicknames", &V);
	this->nicknames.clear();
	for(int i=0; i<(int)V.size(); i++){
		key = V.at(i)->first;
		value = V.at(i)->second;
		idx = key.find(angle);
		name_cur = key.substr(1, idx-1); // strip leading quote, end before first angle
		if( name_cur.compare(i_name_lower)==0 ){
			// for key, val strip trailing quote and strip leading, trailing quotes
			key = key.substr(idx+2, key.size()-idx-3);
			value = value.substr(1, value.size()-2);
			boost::algorithm::replace_all(key, colon_word, colon);
			this->nicknames[ key ] = value;
		}
		delete V.at(i);
	}
    // push most recently loaded investigation to front of list
    vector<string> names_cur;
    stringstream names_new;
    names_new << i_name_lower;
    boost::algorithm::split(names_cur, names, boost::algorithm::is_any_of(",") );
    for(int i=0; i<(int)names_cur.size(); i++){
        if( names_cur.at(i).compare(i_name_lower) != 0)
            names_new << "," << names_cur.at(i);
    }
    cp->set("Filesets", "names", names_new.str() );
    cp->write();
}


void Investigation::probes_with_percent_req(Rawdata& expr, boost::unordered_map<int,int>& probe_idx_present, float fraction_required, std::vector<int>& idx_expr){
	//for each probe in Rawdata expr, identify index of all probes with sufficient present values
	if( fraction_required > 1)
		throw std::string("ERROR: fraction_required must be betwen 0 and 1");
	if( fraction_required==0 || !expr.data->has_missing_values ){
		for(int i=0; i<(int)expr.data->rows(); i++){
			probe_idx_present[i]=1;
		}
		return;
	}
	int n_req = (int)((double)idx_expr.size() * fraction_required );
	int n_present;
	for(int i=0; i<(int)expr.data->rows(); i++){
		n_present = 0;
		if(expr.data->row_has_missing_value(i)){
			for(int a=0;a<(int)idx_expr.size(); a++){
                if(!expr.data->is_missing(i,idx_expr.at(a))){
                    n_present++;
                }
            }
			if( n_present >= n_req )
				probe_idx_present[i]=1;
		}
		else
			probe_idx_present[i]=1;
	}
}


void Investigation::write_dataset( std::string fn_out, std::vector<std::string> probes, std::vector<std::string> samples, float percent_req ){

	Rawdata expr;
	expr.load(this->fn_expr);
	std::vector<int> idx_expr;
	int n_ids = (int)samples.size();
	for(int i=0; i<n_ids; i++){
		idx_expr.push_back( expr.sample_name2idx[samples.at(i)] );
	}
	std::vector<int> valid_probe_idx;
	boost::unordered_map<int,int> probe_idx_present;	
	probes_with_percent_req(expr, probe_idx_present, percent_req, idx_expr);
	if( probes.size() > 0 ){
		for( int i=0; i<(int)probes.size(); i++){
			try{
				if( probe_idx_present.find(this->ga->identifier2idx[ probes.at(i) ]) != probe_idx_present.end() )
					valid_probe_idx.push_back( this->ga->identifier2idx[ probes.at(i) ] );
			}
			catch( ... ){
				std::stringstream err;
				err << "Error limiting dataset by probes; probe " << probes.at(i) << " not found in gene attributes.";
				throw std::string(err.str());
			}
		}
	}
	else{
		for(int i=0;i<(int)this->ga->identifiers.size();i++){
			if( probe_idx_present.find(i) != probe_idx_present.end() )
				valid_probe_idx.push_back(i);
		}
	}
	
	std::ofstream f_out(fn_out.c_str());
	if( !f_out.is_open() ) 
		throw std::string("Unable to open for writing: " + fn_out);
	f_out << "PROBE_ID\tName";
	for(int i=0; i<n_ids; i++)
		f_out << "\t" << samples.at(i);
	f_out << "\n";
	for(int i=0; i<(int)this->sa->attrib.size(); i++){
		f_out << "\t" << this->sa->attrib.at(i);
		for( int j=0; j<n_ids; j++)
			f_out << "\t" + this->sa->prop_for_identifier(samples.at(j),this->sa->attrib.at(i));
		f_out << "\n";
	}
	if( this->ga->identifiers.size() != expr.identifiers.size())
		throw std::string( "Gene attribute file and expression data file do not have\nthe same number of probes.  Check source data." );
	std::string probe_id;
	std::string gene_name_column = this->ga->get_gene_name_column();
	int row=0;
	int col=0;
	for(int i=0; i<(int)valid_probe_idx.size(); i++ ){
		probe_id = expr.identifiers.at(valid_probe_idx.at(i));
		f_out << probe_id << "\t" << ga->prop_for_identifier(probe_id, gene_name_column);
		row = valid_probe_idx.at(i);
		if( expr.data->row_has_missing_value(row ) ){
			for( int j=0; j<n_ids; j++ ){
				col = idx_expr.at(j);
				if(expr.data->is_missing(row, col) ){
					f_out << "\tNA";
				}
				else{
					f_out << "\t" << expr.data->arr[valid_probe_idx.at(i)][col];
				}
			}
		}
		else{
			for( int j=0; j<n_ids; j++){
				f_out << "\t" << expr.data->arr[valid_probe_idx.at(i)][idx_expr.at(j)];
			}
		}
		f_out << "\n";
	}
	f_out.close();
}

Investigation::~Investigation(){
	delete this->cp;
	this->cp = NULL;
}
