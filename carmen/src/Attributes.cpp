#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <algorithm>
using namespace std;
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/random.hpp>
using namespace boost;
namespace alg = boost::algorithm;
#include "DataStructures.h"
#include "Attributes.h"

Attributes::Attributes(string no_result){
	this->no_result = no_result;
    this->file_name = std::string("");
    this->col_chr = std::string("");
    this->col_locus = std::string("");
    this->col_gene_names = std::string("Gene Name");
}

std::string Attributes::get_gene_name_column(){
	return this->col_gene_names;
}


std::string Attributes::get_chr_column(){
    return this->col_chr;
}

void Attributes::set_gene_name_column(std::string col_name){
	if( this->attribute2idx.find(col_name) == this->attribute2idx.end() )
		throw std::string("Requested gene name column not found in attributes");
	this->col_gene_names = std::string(col_name);
}

void Attributes::load(std::string file_name){
	// It is possible-- though perverse-- for a file to consist of a mix of unix-style (\n)
	// and windows-style ("\n\r") lines. That's the reason for checking the line ending on each line.
	//
	this->file_name = file_name;
    std::ifstream f_a(file_name.c_str());
	if( !f_a.is_open() ){ 
		throw std::string( "Unable to open file for reading: " + file_name );
	}
	std::string line;
		
	typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
	boost::char_separator<char> sep("\t");
	getline(f_a, line);

	if( line.at( line.size()-1 ) == '\r'  )
		line = line.substr(0,line.size()-1);

	tokenizer tok(line, sep);
	int col_idx=0;
	for(tokenizer::iterator iter=tok.begin(); iter!=tok.end(); iter++){
		std::string attribute = (*iter);
		if( col_idx==0 && attribute.compare("IDENTIFIER") != 0 )
			throw std::string("Incorrect file format in attribute file: first element of header must be IDENTIFIER; it is " + attribute);
		
		this->attrib.push_back( attribute );
		this->attribute2idx[ (*iter) ] = col_idx;
		++col_idx;
	}
	int n_attrib = (int)this->attrib.size();
	int row_idx=0;
	while( !f_a.eof() ){
		getline(f_a, line);
		if(line.length()==0)
			break;
		if( line.at( line.size()-1 ) == '\r'  )
			line = line.substr(0,line.size()-1);
		tokenizer tok(line, sep);	
		col_idx=0;
		string* val = new string[ n_attrib ];
		for(tokenizer::iterator iter=tok.begin(); iter!=tok.end(); ++iter){
			if( col_idx == n_attrib ){
				std::stringstream ss;
				ss << "incorrect number of columns at line " << row_idx+2;
				throw std::string( ss.str() );
			}
			val[col_idx] = *iter;
			if( col_idx==0 ){
				if( this->identifier2idx.find( (*iter) ) != this->identifier2idx.end() )
					throw std::string("Duplicate identifier in attributes: " + (*iter));
				this->identifiers.push_back((*iter));
				this->identifier2idx[(*iter)]=row_idx;
			}	
			++col_idx;
		}
		this->values.push_back( val );
		++row_idx;
	}
}

Attributes::~Attributes(){
	for(int i=0; i<(int)this->values.size(); i++)
		delete [] values.at(i);
}


void Attributes::restrict_identifiers( std::vector<std::string> new_identifiers ){
    // after call, identifiers should match the order of new_identifiers
    // If an identifier in target is not present, that's fatal
    // UPDATEs identifiers, values, identifier2idx
    
    std::vector<int> idx_of_orig_to_keep, idx_orig_to_remove;
    HASH_S_I hsh_target;
    for(int i=0; i<int(new_identifiers.size());i++){
        if( this->identifier2idx.find( new_identifiers.at(i) ) == this->identifier2idx.end() ){
            std::stringstream ss;
            ss << "called restrict_identifiers with identifier not found in attributes: " << new_identifiers.at(i);
            throw ss.str();
        }
        hsh_target[new_identifiers.at(i)]=1; // track seen so we know which to delete
        idx_of_orig_to_keep.push_back( this->identifier2idx[ new_identifiers.at(i) ] ); // store original index of identifier
    }
    // keep pointers to values for each new identifier
    std::vector<std::string*> new_values;
    for(int i=0; i<int(idx_of_orig_to_keep.size()); i++){
        new_values.push_back( this->values.at( idx_of_orig_to_keep.at(i)) );
    }
    // delete values that are not present in target
    for(int i=0; i<int(this->identifiers.size()); i++){
        if( hsh_target.find(this->identifiers.at(i)) == hsh_target.end() ){
            delete [] this->values.at(i);
        }
    }

    this->identifiers.clear();
    this->values.clear(); // pointers to the values that we're keeping are in new_values
    this->identifier2idx.clear();
    for(int i=0; i<int(new_identifiers.size()); i++){
        this->identifiers.push_back(new_identifiers.at(i) );
        this->values.push_back(new_values.at(i));
        this->identifier2idx[ this->identifiers.at(i) ]=i;     // rebuild identifier2idx
    }

}

bool Attributes::get_chromosome_and_locus_by_identifier( std::string identifier, std::string& chromosome, int& locus){
    // If not data are missing, return false; 
    if( this->col_chr.size()==0 || this->col_locus.size()==0 )
        throw std::string("Genomic information not loaded for this attribute object.");
    bool locus_and_chr_present = true;
    chromosome = this->prop_for_identifier(identifier, this->col_chr);
    if( chromosome.compare("NA")==0 )
        locus_and_chr_present=false;
    std::string locus_str = this->prop_for_identifier(identifier, this->col_locus);
    try{
        locus = boost::lexical_cast<int>(locus_str);
    }
    catch( boost::bad_lexical_cast & ){
        locus=-1;
        locus_and_chr_present = false;
    }
    return locus_and_chr_present;
}

void Attributes::set_chromosome_and_locus( std::string col_chr, std::string col_locus ){
    // If user tells us which column contains the chromosome and locus data
    // we can store the column names and vectors called chromosome and locus with the 
    // relevant strings.
 	if( this->attribute2idx.find(col_chr) == this->attribute2idx.end() )
		throw std::string("Requested chromosome column not found in attributes");
	if( this->attribute2idx.find(col_locus) == this->attribute2idx.end() )
		throw std::string("Requested locus column not found in attributes");
    this->col_chr = col_chr;
    this->col_locus = col_locus;
    std::string chr, loc;
    int locus;
    std::string NA("NA");
    for(int i=0; i<int(this->identifiers.size()); i++){
        chr = this->prop_for_identifier(this->identifiers.at(i) , col_chr);
        loc = this->prop_for_identifier(this->identifiers.at(i) , col_locus);
        try{
            locus = boost::lexical_cast<int>(loc);
        }
        catch( boost::bad_lexical_cast & ){
            continue;
        }
        this->chromosome.push_back(chr);
        this->locus.push_back(locus);
        this->idx_of_locus.push_back(i);
    }
}

void Attributes::get_idx_by_genomic_range( std::string chr, int loc_start, int loc_end, std::vector<int> & matching_loci){
    matching_loci.clear();
    for(int i=0; i<int(this->chromosome.size()); i++){
        if( chr.compare(this->chromosome.at(i))==0 ){
            if( loc_start <= this->locus.at(i) && loc_end >= this->locus.at(i) ){
                matching_loci.push_back(this->idx_of_locus.at(i) );
            }
        }
    }
}


void Attributes::indices_with_property(string attribute, string value, vector<int>* row_idx){
	this->indices_with_property(attribute, value, row_idx, true);
}

void Attributes::indices_with_property(string attribute, string value, vector<int>* row_idx, bool case_sensitive){
	row_idx->clear();
	int idx_attributes;
	HASH_S_I::iterator f = this->attribute2idx.find(attribute);
	if( f == this->attribute2idx.end() )
		return;
	else
		idx_attributes = (*f).second;
	
	string* row = NULL;
	std::string value_lc = value;
	boost::algorithm::to_lower(value_lc);
	std::string comparitor_lc;
	for(int i=0; i<(int)this->values.size(); i++){
		row = this->values.at(i);
		if( case_sensitive ){
			if( row[ idx_attributes ].compare( value ) == 0){
				row_idx->push_back(i);
			}
		}
		else{
			comparitor_lc = row[ idx_attributes ];
			boost::algorithm::to_lower(comparitor_lc);
			if( comparitor_lc.compare( value_lc ) == 0){
				row_idx->push_back(i);
			}
		}
	}
}


void Attributes::indices_without_property(string attribute, string value, vector<int>* row_idx){
	row_idx->clear();
	int idx_attributes;
	HASH_S_I::iterator f = this->attribute2idx.find(attribute);
	if( f == this->attribute2idx.end() )
		return;
	else
		idx_attributes = (*f).second;
	
	string* row = NULL;
	for(int i=0; i<(int)this->values.size(); i++){
		row = this->values.at(i);
		if( row[ idx_attributes ].compare( value ) != 0)
			row_idx->push_back(i);
	}
}


std::string Attributes::prop_for_identifier(string identifier, string attribute){
	// given an identifier and an attribute, return the property.  
	// If there is no such identifier, no such attribute, or the value is missing,
	// return the "no_result" string.
	int idx_attributes, idx_identifiers;
	HASH_S_I::iterator f = this->identifier2idx.find(identifier);
	if( f == this->identifier2idx.end() )
		return this->no_result;
	else
		idx_identifiers = (*f).second;
	
	f = this->attribute2idx.find(attribute);
	if( f == this->attribute2idx.end() )
		return this->no_result;
	else
		idx_attributes = (*f).second;
	
	string* row = this->values.at(idx_identifiers);
	return row[idx_attributes];
}


void Attributes::intersection( vector<int>* source, vector<int>* diff){
	// limit source to those items also in diff
	vector<int> keep = vector<int>();
	for(vector<int>::iterator iter = diff->begin(); iter != diff->end(); ++iter){
		int d = (*iter);
		for(vector<int>::iterator inner = source->begin(); inner != source->end(); ++inner){
			if( d==(*inner)){
				keep.push_back(d);
				break;
			}
		}
	}
	source->clear();
	for(vector<int>::iterator iter = keep.begin(); iter != keep.end(); ++iter)
		source->push_back(*iter);
}


void Attributes::write(std::string file_name_out){
	std::ofstream f_out(file_name_out.c_str());
	if( !f_out.is_open() ){ 
		std::stringstream ss;
		ss << "Unable to open file for writing:" << file_name_out << "\n";
		throw std::string( ss.str() );
	}
	int n_attrib = (int)this->attribute2idx.size();
	int n_id = (int)this->identifiers.size();
	int i;
	for( i=0; i < n_attrib; i++){
		f_out << this->attrib.at(i);
		if( i < n_attrib-1 )
			f_out << "\t";
		else
			f_out << "\n";
	}
	
	for( i=0; i < n_id; i++){
		f_out << this->identifiers.at(i) << "\t";
		for(int a = 0; a<n_attrib; a++){
			f_out << this->prop_for_identifier( this->identifiers.at(i), this->attrib.at(a) );
			if( a < n_attrib-1 )
				f_out << "\t";
			else
				f_out << "\n";
		}
	}
}


void Attributes::find_identifiers_in_class(std::string limit, std::vector<std::string>* identifiers_in_class){
	// Only finds identifiers based on sample attributes.
	// If user has specified gene: limits, these cannot yet be evaluated here because 
	// attributes don't know about discretized values. If (for example) class_a=gene:12312_at=1, 
	// then this will return all samples.
	// During a dataset load, these samples are later on pruned out.
	if( limit.length() == 0 )
		return;
	vector<int>* idx = new vector<int>();
	vector<string> limits;
	stringhasher::split(limit, &limits, std::string(",").c_str());
	string attribute, value;
	idx->clear();
	for(int i=0; i<(int)this->identifiers.size();i++){
		idx->push_back(i);  
	}
	vector<int>* row_idx = new vector<int>();
	bool is_equal=true;
	for(vector<string>::iterator iter=limits.begin(); iter != limits.end(); ++iter){
		stringhasher::parse_limit(*iter, attribute, value, is_equal);
		if(attribute.length()<5 || attribute.substr(0,5).compare("gene:")!=0 ){
			// ignore attributes that start with "gene:"...
			if(this->attribute2idx.find(attribute) == this->attribute2idx.end() )
				throw std::string("Sample attributes file does not contain an attribute " + attribute );
			if(is_equal)
				this->indices_with_property(attribute, value, row_idx);
			else
				this->indices_without_property(attribute, value, row_idx);
			this->intersection(idx, row_idx);
		}
	}
	for(int i=0; i<(int)idx->size(); i++){
		identifiers_in_class->push_back( this->identifiers.at( idx->at(i) ) );
	}
	delete row_idx;
	delete idx;
}
