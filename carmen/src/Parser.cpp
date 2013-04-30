#define _SCL_SECURE_NO_DEPRECATE 1
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/random.hpp>
using namespace boost;
namespace alg = boost::algorithm;
using namespace std;
#include "DataStructures.h"
#include "Parser.h"


GOAnnotation::GOAnnotation(){
	this->idx = -1;
	this->go_id = std::string("");
	this->branch = GOAnnotation::GO_BP;
	this->description = std::string("");
}


GOAnnotation::GOAnnotation( int idx, std::string go_id, GO_branch branch, std::string description ){
	this->idx = idx;
	this->go_id = go_id;
	this->branch = branch;
	this->description = description;
}


GOAnnotationParser::GOAnnotationParser(){
	// allow for instantiation without load
}

void GOAnnotationParser::load(std::string fn){

	this->go_id2idx.clear();
	for(int i=0; i<(int)this->GO.size(); i++){
		delete this->GO.at(i);
	}
	this->GO.clear();

	ifstream f_lines(fn.c_str());
	if( !f_lines.is_open() ){
		throw string("Unable to open file for reading: " + fn );
	}
	bool has_CRLF = false;
	std::string line;
	getline(f_lines, line); // header
	if( line.at( line.size()-1 ) == '\r' )
		has_CRLF=true; // Windows-originated file
	getline(f_lines, line);
	std::vector<std::string> holder;
	int idx = 0;
	GOAnnotation::GO_branch branch = GOAnnotation::GO_BP;
	while( !f_lines.eof() ){
		if( (int)line.size()>0 ){
			if( has_CRLF )
				line = line.substr(0,line.size()-1);
			holder.clear();
			alg::split(holder, line, alg::is_any_of("\t"));
			if(holder.size() != 3){
				throw std::string("Parse error (incorrect number of elements) reading " + line);
			}
			if(holder.at(1).compare("BP")==0)
				branch = GOAnnotation::GO_BP;
			else if(holder.at(1).compare("CC")==0)
				branch = GOAnnotation::GO_CC;
			else if(holder.at(1).compare("MF")==0 )
				branch = GOAnnotation::GO_MF;
			else
				throw std::string("Parse error reading branch in GO annotation: " + line);
			this->GO.push_back( new GOAnnotation( idx, holder.at(0), branch, holder.at(2)) );
			go_id2idx[holder.at(0)] = idx;
			idx++;
		}
		getline(f_lines, line);
	}
}

GOAnnotationParser::GOAnnotationParser(std::string fn){
	this->load(fn);
}

int GOAnnotationParser::get_idx(std::string go_id){
	if( this->go_id2idx.find( go_id ) == this->go_id2idx.end() )
		return -1;
	else
		return this->go_id2idx[ go_id ];
}

int GOAnnotationParser::size(){
	return this->go_id2idx.size();
}

void GOAnnotationParser::get_idx_matching_description(std::string searchterm, std::vector<int> &results){
	results.clear();
	GOAnnotation* g;
	boost::algorithm::to_lower(searchterm);
	HASH_S_I desc2idx;
	std::vector<std::string> sorter;
	for(int i=0; i<(int)this->GO.size(); i++){
		g = this->GO.at(i);
		std::string desc_lower = g->description;
		boost::algorithm::to_lower(desc_lower);
		if( !boost::find_first( desc_lower, searchterm).empty()  ){
			desc2idx[ desc_lower ] = g->idx;
			sorter.push_back(desc_lower);
		}
	}
	std::sort(sorter.begin(), sorter.end() );
	for(int i=0; i<(int)sorter.size(); i++){
		results.push_back( desc2idx[sorter.at(i)] );
	}
}

GOAnnotation* GOAnnotationParser::get_annotation(int idx ){
	if( idx < 0 || idx >= (int)this->GO.size() )
		throw std::string("GO idx request is out of range");
	return this->GO.at(idx);
}

GOAnnotationParser::~GOAnnotationParser(){
	for(int i=0;i<(int)this->GO.size(); i++)
		delete this->GO.at(i);
}


GeneAnnotation::GeneAnnotation(std::string symbol, std::string chr, std::string ucsc_chr, std::string full_name, std::string loc_start, std::vector<int>& GO ){
	this->symbol = symbol;
	this->chr = chr;
	this->ucsc_chr = ucsc_chr;
	this->full_name = full_name;
	this->loc_start = loc_start;
	for(int i=0; i<(int)GO.size(); i++)
		this->GO_annotations.push_back(GO.at(i));
}


GeneAnnotationParser::GeneAnnotationParser(){
}

GeneAnnotationParser::GeneAnnotationParser(std::string fn, GOAnnotationParser& GOP){
	load(fn, GOP);
}

void GeneAnnotationParser::load(std::string fn, GOAnnotationParser& GOP){
	// TODO: Clear hash of symbol2annotation, deleting .second
	ifstream f_lines(fn.c_str());
	if( !f_lines.is_open() ){
		throw string("Unable to open file for reading: " + fn );
	}
	bool has_CRLF = false;
	std::string line;
	getline(f_lines, line); // header
	if( line.at( line.size()-1 ) == '\r' )
		has_CRLF=true; // Windows-originated file
	getline(f_lines, line);
	std::vector<std::string> holder, go_holder;
	int idx;
	std::string symbol;
	while( !f_lines.eof() ){
		if( has_CRLF )
			line = line.substr(0,line.size()-1);
		holder.clear();
		alg::split(holder, line, alg::is_any_of("\t"));
		//symbol	chr	ucsc_chr	loc_start	fullname	GO
		std::vector<int> GO;
		symbol = holder.at(0);
		boost::algorithm::to_lower(symbol);
		if(holder.size() != 6)
			throw std::string("Parse error reading gene annotation: " + line);
		if( holder.at(5).compare("NA") != 0){
			alg::split(go_holder, holder.at(5), alg::is_any_of(","));
			for( int i=0; i<(int)go_holder.size(); i++){
				idx = GOP.get_idx(go_holder.at(i) );
				if( idx >= 0 )
					GO.push_back( idx );
			}
		}
		GeneAnnotation* g = new GeneAnnotation(symbol, holder.at(1), holder.at(2), holder.at(4), holder.at(3), GO );
		this->symbol2annotation[ symbol ] = g;
		getline(f_lines, line);
	}
}

GeneAnnotation* GeneAnnotationParser::get_annotation(std::string symbol){
	boost::algorithm::to_lower(symbol);
	if( this->symbol2annotation.find(symbol) == this->symbol2annotation.end() )
		throw std::string("Symbol not found");
	else{
		return this->symbol2annotation[symbol];
	}
}


void GeneAnnotationParser::get_symbols_with_GO_idx(int target_go_idx, std::vector<std::string>& symbols){
	symbols.clear();
	GeneAnnotation* g;
	for (boost::unordered_map <std::string,GeneAnnotation*>::iterator i=this->symbol2annotation.begin(); i!=symbol2annotation.end(); ++i){
		g = (GeneAnnotation*)i->second;
		for(int j=0; j<(int)g->GO_annotations.size(); j++){
			if( g->GO_annotations.at(j)==target_go_idx ){
				symbols.push_back(g->symbol);
			}
		}
	}
}


KeyValuesParser::KeyValuesParser(std::string fn){
    
	ifstream f_lines(fn.c_str());
	if( !f_lines.is_open() )
		throw string("Unable to open file for reading:" + fn );
    
    this->key2values.clear();
	boost::char_separator<char> sep("\t");	
    string line;
    
    this->has_CRLF=false;
    getline(f_lines, line);
	if(line.size()>0 ){
		if( line.at( line.size()-1 ) == '\r' ){
			this->has_CRLF=true;    
        }
    }
    int ctr;
    std::string header;
    typedef boost::tokenizer< boost::char_separator<char> > tokenizer;
    while(line.size()>0){
        if( this->has_CRLF ){
            line = line.substr(0,line.size()-1);
        }
        ctr=0;
        tokenizer tok(line, sep);
        for(tokenizer::iterator beg=tok.begin(); beg!=tok.end(); ++beg){
            if(ctr==0){
                header = *beg;
                if( this->key2values.find(header) != this->key2values.end() ){
                    std::stringstream ss;
                    ss << "Duplicate key value: " << header;
                    throw ss.str();
                }
                this->key2values[ header ] = new std::vector<std::string>();
                this->keys_in_order.push_back( header );
                ctr++;
            }
            else {
                this->key2values[header]->push_back(*beg);
            }
        }
        getline(f_lines, line);
    }
	f_lines.close();
}

void KeyValuesParser::keys(std::vector<std::string>& k){
    // prefer to preserve original order of keys.
    k.clear();
    for(int i=0; i<int(keys_in_order.size()); i++)
        k.push_back(this->keys_in_order.at(i) );
}



bool KeyValuesParser::has_key(std::string key){
    return !( this->key2values.find(key)==this->key2values.end() );
}

void KeyValuesParser::get(std::string key, std::vector<std::string>& values){
    values.clear();
    if( this->key2values.find(key) != this->key2values.end() ){
        std::vector<std::string>* sp = this->key2values[key];
        for( int i=0; i<int(sp->size()); i++)
            values.push_back(sp->at(i));
    }
}


ConfigParser::ConfigParser(){
	this->fn_loaded = std::string("");
}

ConfigParser::ConfigParser( std::string fn ){
	this->load(fn);
}

void ConfigParser::clear(){
	std::vector< prop_pair* >* pairs;
	for(boost::unordered_map<std::string, std::vector< prop_pair* >* >::iterator itr=this->contents.begin(); itr != this->contents.end(); itr++){
		pairs = (*itr).second;
		for(int i=0; i<(int)pairs->size(); i++)
			delete pairs->at(i);
	}
	this->sections.clear();
	this->contents.clear();
}

void ConfigParser::load(std::string fn){
	// if already loaded with other information, clear this cleanly.
	this->clear();
	ifstream f_lines(fn.c_str());
	if( !f_lines.is_open() ){	
		throw string("Unable to open file for reading: " + fn );
	}
	this->fn_loaded = fn;
	string line;
	string eq(" = ");
	string current_header;
	getline(f_lines, line);
	int idx=0;
	while( !f_lines.eof() ){
		if( (int)line.size()>0 ){
			if(  line.at(0)=='['){
				current_header = line.substr(1, line.size() - 2);
				this->sections.push_back(current_header);
				this->contents[ current_header ] = new std::vector<prop_pair*>();
			}	
			else{
				idx=(int)line.find(eq);
				prop_pair* pp = new prop_pair();
				pp->first = line.substr(0, idx);
				pp->second = line.substr(idx+3,line.size());
				this->contents[current_header]->push_back( pp );
			}
		}
		getline(f_lines, line);
	}
	f_lines.close();
}


void ConfigParser::add_section(std::string section){
	if( this->contents.find(section) == this->contents.end() ){
		this->contents[ section ] = new std::vector<prop_pair*>();
		this->sections.push_back(section);
	}
}


void ConfigParser::create(std::string file_name){
	std::ofstream f_out(file_name.c_str());
	if( !f_out.is_open() )
		throw std::string("Unable to open for writing: " + file_name);
	f_out.close();
	this->fn_loaded = file_name;
}


void ConfigParser::set(std::string section, std::string A, std::string B){
	
	if( this->contents.find(section) == this->contents.end() ){	
		this->contents[ section ] = new std::vector<prop_pair*>();
		this->sections.push_back(section);
	}
	for(int i=0; i<(int)this->contents[section]->size(); i++){
		if( this->contents[section]->at(i)->first.compare(A)==0){
			this->contents[section]->at(i)->second = B;
			return;
		}
	}
	prop_pair* pp = new prop_pair();
	pp->first = A;
	pp->second = B;
	this->contents[section]->push_back( pp );
}


std::string ConfigParser::get(std::string section, std::string A){
	if( this->contents.find(section) != this->contents.end() ){
		for(int i=0; i<(int)this->contents[section]->size(); i++){
			if( this->contents[section]->at(i)->first.compare(A)==0){
				return this->contents[section]->at(i)->second;
			}
		}
	}
	return std::string();
}

void ConfigParser::remove(std::string section_name, std::string A){
	// if we find the item to be removed, 
	//		delete the pair holding the item
	//		store the remaining pairs in remember
	//		pop pairs from removed to end
	//		push stored items from remember 
	if( this->contents.find(section_name) != this->contents.end() ){
		std::vector< prop_pair* >* section = this->contents[section_name];
		std::vector< prop_pair* > remember;
		bool store = false;
		for(int i=0; i<(int)section->size(); i++){
			if(store){
				remember.push_back(section->at(i));
			}
			if( section->at(i)->first.compare(A)==0){
				store = true;
				delete section->at(i);
			}
		}
		if(store){
			for(int i=0; i<(int)remember.size() + 1; i++)
				section->pop_back();
			for(int i=0; i<(int)remember.size(); i++)
				section->push_back( remember.at(i) );
		}
	}
}

void ConfigParser::get_section(std::string section, std::vector<prop_pair*>* V){
	V->clear();
	if( this->contents.find(section) != this->contents.end() ){
		for(int i=0; i<(int)this->contents[section]->size(); i++){
			prop_pair* pp = new prop_pair();
			pp->first = this->contents[section]->at(i)->first;
			pp->second = this->contents[section]->at(i)->second;
			V->push_back(pp);
		}
	}
}

void ConfigParser::write(){
	write(this->fn_loaded.c_str() );
}

void ConfigParser::write(std::string fn){
	std::ofstream f_out(fn.c_str());
	if( !f_out.is_open() ){
		throw string( "Unable to open file for writing:" + fn );
	}
	std::vector< prop_pair* >* pairs;
	for(int i=0; i<(int)this->sections.size(); i++){
		f_out << "[" << this->sections.at(i) << "]\n";
		pairs = this->contents[this->sections.at(i)];
		for(int i=0; i<(int)pairs->size(); i++)
			f_out << pairs->at(i)->first << " = " << pairs->at(i)->second << "\n"; 
		f_out << "\n";
	}
	f_out.close();
}

ConfigParser::~ConfigParser(){
	this->clear();
}

// CarmenParser: Read CARMEN-style files (ruleset, difference, spear, classifier)

CarmenParser::CarmenParser(){
}


CarmenParser::CarmenParser(std::string filename){
	this->file_has_values = false;
	this->has_CRLF = false;
	Load(filename);
}

int CarmenParser::CountValueLines(){
	// Counts the number of lines that do not start with #
	// Used by methods that require me to know this before I can 
	// call them (e.g. grid CreateRows()
	int n_lines=0;
	std::string line;
	ifstream f_lines(this->fn_loaded.c_str());
	if( !f_lines.is_open() )
		throw string("Unable to open file for reading:" + this->fn_loaded );
	getline(f_lines, line);
	while( !f_lines.eof() ){
		if( (int)line.size()>0 ){
			if(  line.at(0)!='#'){
				n_lines++;
			}
		}
		getline(f_lines, line);
	}
	f_lines.close();
	return n_lines;
}

void CarmenParser::Load(std::string fn){
	// Loads the header into header2values

	std::vector<std::string> file_name_parts;
	boost::algorithm::split(file_name_parts, fn, boost::algorithm::is_any_of(".") );
	this->extension = file_name_parts.at( file_name_parts.size()-1 );

	ifstream f_lines(fn.c_str());
	if( !f_lines.is_open() )
		throw string("Unable to open file for reading:" + fn );
	this->fn_loaded = fn;
	this->header2values.clear();
	std::string line;
	std::string space(" ");
	std::string current_header;
	if(f_lines.eof() ){
		f_lines.close();
		return;
	}
	getline(f_lines, line);
	if(line.size()>0 )
		if( line.at( line.size()-1 ) == '\r' )
			this->has_CRLF=true; //Windows-originated file
	int idx=0;//, token_idx=0;
	typedef boost::tokenizer< boost::char_separator<char> > tokenizer;
	boost::char_separator<char> sep(" ");
	std::vector<std::string> holder;
	while( !f_lines.eof() ){
		if( (int)line.size()>0 ){
			if( line.at(0)=='#'){
				if( this->has_CRLF )
					line = line.substr(0,line.size()-1);
				idx=(int)line.find( space,2 );
				current_header = line.substr(2, idx - 2);
				std::vector<std::string>* values = new std::vector<std::string>();
				holder.clear();
				alg::split(holder, line, alg::is_any_of(" "));
				for(int i=2; i<(int)holder.size(); i++)
					values->push_back(holder.at(i));
				this->header2values[ current_header ] = values;
				if(int(line.size())>=3){
					current_header = line.substr( 2, line.size() ); // cut off leading #
					alg::split(this->column_names, current_header, alg::is_any_of("\t") );
				}
			}	
			else{
				// current value of line will be used for column names
				break;
			}
		}
		getline(f_lines, line);
	}
	f_lines.close();	
}


void CarmenParser::write_compressed_spear(std::string fn_out, HASH_S_I& p2id){
	std::string dummy;
	if( !this->get_value("Abs_rho", 0, dummy) )
		throw std::string("This function can only be called on a spear file.");
	PrepareToReadValues();
	std::vector<std::string> values;
	std::ofstream f_out(fn_out.c_str());
	if( !f_out.is_open() ){
		throw std::string( "Unable to open requested file for writing");
	}
	while( ReadNextValue(values) ){
		f_out << p2id[ values.at(1) ] << " " << p2id[values.at(3) ] << " " << values.at(4) << "\n";
	}
	f_out.close();
}


void CarmenParser::get_column_names( std::vector<std::string> & cn ){
	// return a line from the file header, if present.
	cn.clear();
	for(int i=0; i<int(this->column_names.size()); i++){
		cn.push_back( this->column_names.at(i) );
	}
}

bool CarmenParser::get_value_vector( std::string header_keyword, std::vector<std::string>* header_values ){
	// return a line from the file header, if present.
	if(this->header2values.find(header_keyword)==this->header2values.end()){
		header_values->clear(); // header_keyword not present
		return false;
	}
	header_values = this->header2values[header_keyword];
	return true;
}


bool CarmenParser::get_value(std::string header_keyword, int idx, std::string& s){
	// Return a string from the file header, if present. 
	// Header values are space delimited, and indexing starts after the
	// keyword (e.g. use index 0 to read the 5 in # Max_Depth 5 )
	if(this->header2values.find(header_keyword)==this->header2values.end()){
		s = ""; // header_keyword not present
		return false;
	}
	std::vector<std::string>* v = this->header2values[header_keyword];
	if(idx<0 || idx > (int)v->size()-1){
		s = ""; // index error
		return false;
	}
	s = this->header2values[header_keyword]->at(idx);
	return true;
}

void CarmenParser::ExtractEQTL( std::vector<QTL*>& qtls, double max_eqtl_pval){
	//1.50573e-13	Ets1	1422027_a_at	E09.034.366_10	0.00000	6.54406	8.72455	0.00000
	std::vector<std::string> values;
	for(int i=0; i<int(qtls.size()); i++){
		delete qtls.at(i);
	}
	qtls.clear();
	std::string probe_id, probe_name, locus_id;
	double pval;
	while( ReadNextValue(values)){
		probe_name = values.at(1);
		probe_id = values.at(2);
		locus_id = values.at(3);
		try{ pval= boost::lexical_cast<double>( values.at(4)); }
		catch( boost::bad_lexical_cast &){
			throw std::string( "Invalid pval reading QTL file");
		}
		// intentionally pushing locus_id twice, don't currently have fancy name for loci
		if( pval <= max_eqtl_pval ){
			qtls.push_back( new QTL(locus_id, locus_id, probe_id, probe_name, pval) );
		}
	}
}


void CarmenParser::ExtractProbes(std::vector<string>& probes){
	std::string delta, ids_raw, mine_type, identifier, value;
	std::vector<std::string> values;
	std::vector<std::string> ids, ids_with_is;
	HASH_S_I hsh_probes;
	if( this->extension.compare("ruleset")==0){ // ruleset
		while( ReadNextValue(values)){
			ids_raw = values.at(8);
			alg::split(ids_with_is, ids_raw, alg::is_any_of("|") );
			for(int i=0; i<(int)ids_with_is.size(); i++){
				stringhasher::get_idval(ids_with_is.at(i), identifier, value);
				hsh_probes[identifier] = 1;
			}
		}
	}
	else if( this->extension.compare("spear")==0 ){ // spear
		while( ReadNextValue(values)){
			hsh_probes[values.at(1)] = 1;
			hsh_probes[values.at(3)] = 1;
		}
	}
	else{ // just read first column
		while( ReadNextValue(values)){
			hsh_probes[values.at(0)] = 1;
		}
	}
	for( HASH_S_I::iterator iter = hsh_probes.begin(); iter != hsh_probes.end(); ++iter){
		probes.push_back( (*iter).first );
	}
}


void CarmenParser::ExtractGeneCounts(HASH_S_I& hsh_genes){
	std::string delta, ids_raw, mine_type, identifier, value;
	std::vector<std::string> values;
	std::vector<std::string> ids, ids_with_is;
	if( this->extension.compare("ruleset")==0){ // ruleset
		while( ReadNextValue(values)){
			ids_raw = values.at(7);
			alg::split(ids_with_is, ids_raw, alg::is_any_of("|") );
			for(int i=0; i<(int)ids_with_is.size(); i++){
				stringhasher::get_idval(ids_with_is.at(i), identifier, value);
				if( hsh_genes.find(identifier) == hsh_genes.end() )
					hsh_genes[identifier] = 1;
				else
					hsh_genes[identifier] = hsh_genes[identifier]+1;
			}
		}
	}
	else if( this->extension.compare("spear")==0 ){ // spear
		while( ReadNextValue(values)){
			if( hsh_genes.find(values.at(0)) == hsh_genes.end() )
				hsh_genes[values.at(0)] = 1;
			else
				hsh_genes[values.at(0)] = hsh_genes[values.at(0)]+1;

			if( hsh_genes.find(values.at(2)) == hsh_genes.end() )
				hsh_genes[values.at(2)] = 1;
			else
				hsh_genes[values.at(2)] = hsh_genes[values.at(2)]+1;
		}
	}
	else{ // just read first column
		while( ReadNextValue(values)){
			if( hsh_genes.find(values.at(0)) == hsh_genes.end() )
				hsh_genes[values.at(0)] = 1;
			else
				hsh_genes[values.at(0)] = hsh_genes[values.at(0)]+1;
		}
	}
}


void CarmenParser::ExtractGenes(std::vector<string>& genes){
	HASH_S_I hsh_genes;
	this->ExtractGeneCounts(hsh_genes);
	for( HASH_S_I::iterator iter = hsh_genes.begin(); iter != hsh_genes.end(); ++iter){
		genes.push_back( (*iter).first );
	}
}

void CarmenParser::PrepareToReadValues(){
	// figure out how many lines of # header there are.  Close the file pointer, open it again,
	// and read up to one before this so that the next call to ReadNextValue will pick
	// up the first line of values.

	if( this->fs.is_open() )
		this->fs.close();

	this->fs.open(this->fn_loaded.c_str());
	if( !this->fs.is_open() )
		throw string("Unable to open file for reading:" + this->fn_loaded );
	int first_real = 0;
	std::string line;
	getline(this->fs, line);
	while(!this->fs.eof() && line.length()>0 && line.at(0)=='#'){
		getline(this->fs, line);
		++first_real;
	}

	if(this->fs.eof() || first_real==0 || line.size()==0 )
		this->file_has_values = false;
	else
		this->file_has_values = true;
	this->fs.close();
	if(this->file_has_values){
		this->fs.open(this->fn_loaded.c_str() );
		for(int i=0; i<first_real; i++)
			getline(this->fs, line);
	}
}

bool CarmenParser::ReadNextValue( std::vector<std::string> & values){
	// Assumes that this->fs points to a valid file that has been 
	// opened using PrepareToReadValues
	if( this->file_has_values ){
		if( !this->fs.is_open() )
			throw string("ERROR: File not open.  Did you forget PrepareToReadValues()?");
		boost::char_separator<char> sep("\t");	
		string line;
		if( this->fs.eof() ){
			this->fs.close();
			return false;
		}
		getline(this->fs, line);
		if(line.size()==0){
			this->fs.close();
			return false;
		}
		//std::cout << "<" << line << ">" << "\n";
		if( this->has_CRLF ){
			line = line.substr(0,line.size()-1);
			//std::cout << "Has crlf\n";
			//std::cout << "<" << line << ">" << "\n";
		}
		tokenizer tok(line, sep);
		values.clear();
		for(tokenizer::iterator beg=tok.begin(); beg!=tok.end(); ++beg)
			values.push_back(*beg);
		return true;
	}
	else
		return false;
	
}

