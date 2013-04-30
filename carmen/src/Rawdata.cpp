#include <vector>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
using namespace std;
#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>
#include <boost/random.hpp>
#include "DataStructures.h"
#include "Rawdata.h"

Rawdata::Rawdata(){
    this->file_name = "";
	this->MISSING_DATA_FLOAT = -99999;
	this->data = NULL;
	verbose = false;
}

Rawdata::~Rawdata(){
	if(this->data != NULL)
		delete this->data;
}

void Rawdata::load( std::string file_name ){
	// Note that this function scans using %f because the standard missing data element is
	// -99999, and allows conversion to int later. 
	// This function stores the results in rows, n_cols, and has_missing_data; it must be called by
	// one of the derived classes's load() functions.
	const int MAX_IDENTIFIER_LENGTH = 200;
	int n_cols = -1;
	bool has_CRLF=false;
	this->file_name = file_name;

	typedef boost::tokenizer<boost::char_separator<char> > tokenizer;	
	if( file_name.size()==0 )
		throw std::string("Empty file name string passed to Rawdata::load");
	ifstream f_lines(file_name.c_str());
	if( !f_lines.is_open() )
		throw std::string("Unable to open file for reading: " + file_name );

	boost::char_separator<char> sep("\t");
	// read sample names from first line into this->sample_names, also set n_cols.
	std::string line;
	getline(f_lines, line);
	if( line.at( line.size()-1 ) == '\r' )
		has_CRLF=true; //Windows-originated file
	if(has_CRLF)
		line = line.substr(0,line.size()-1);
	tokenizer tok(line, sep);
	for(tokenizer::iterator beg=tok.begin(); beg!=tok.end(); ++beg){
		if(n_cols==-1){
			if( (*beg).compare("IDENTIFIER") != 0 )
				throw std::string("Incorrect file format in raw data: first element of header must be IDENTIFIER; it is " + (*beg));
		}
		else{
			this->sample_names.push_back( *beg );
			this->sample_name2idx[(*beg)] = n_cols;
		}
		++n_cols;
	}
	f_lines.close();
	FILE * pFile;
	bool has_missing_values=false;

	char x[MAX_IDENTIFIER_LENGTH];
	pFile = fopen(file_name.c_str(), "r");
	int out_f, out_line;	
	for(int col=0;col<n_cols+1;col++){ 
		out_f = fscanf(pFile, "%s", x);
	} // ignore first line, already read it.
	
	std::vector<std::string>::iterator iter;
	float f;
	int c;
	std::vector<float*> rows;
	HASH_S_I::iterator it;
	char identifier[MAX_IDENTIFIER_LENGTH];
	out_line = fscanf(pFile, "%s", identifier);
	int line_no = 1;
	std::vector<int> missing_r, missing_c;
	float* cols;
	while(out_line>0){
		string identifier_str = std::string(identifier);
		if(this->identifier2idx.find(identifier_str)!=this->identifier2idx.end()){
			std::stringstream ss;
			ss << "Duplicate identifier (" << identifier_str << ") on line " << line_no;
			throw std::string(ss.str());
		}
		line_no++;
		this->identifiers.push_back(identifier_str);
		this->identifier2idx[identifier_str] = (int)rows.size();

		cols = new float[n_cols];
		for(c=0; c<n_cols; c++){
			out_f = fscanf(pFile, "%f", &f);
			if( out_f==0 ){
				// couldn't parse as floating point, interpret as missing data
				f=MISSING_DATA_FLOAT;
#ifndef UBUNTU
				out_f = fscanf(pFile, "%s", x);
#endif
				has_missing_values=true;
			}
			else if( f==MISSING_DATA_FLOAT ){
				has_missing_values=true;
			}
			cols[c]=f;
		}
		rows.push_back(cols);	// push the values onto rows
		if( (int)rows.size() % 50000 == 0 ){
			if( this->verbose ){
				std::cout << "MESSAGE: Reading line " << (int)rows.size()+1 << "\n";
				std::cout.flush();
			}
		}
		out_line = fscanf(pFile, "%s", identifier);
	}
	fclose(pFile);
	int n_rows = (int)rows.size();
	int r=0;
	this->data = new Matrix<float>(n_rows, n_cols, has_missing_values);
	for(r=0; r<n_rows; r++){
		cols = rows.at(r);
		for(c=0; c<n_cols; c++){
			f = cols[c];
			if( f==MISSING_DATA_FLOAT )
				this->data->set_missing_value(r,c);
			else
				this->data->arr[r][c] = f;
		}
		delete[] rows.at(r);
	}
}

int Rawdata::index_of_identifier(std::string identifier){
	HASH_S_I::iterator it = this->identifier2idx.find(identifier);
	if( it==this->identifier2idx.end() )
		return -1;
	else
		return( (*it).second );
}

int Rawdata::index_of_sample(std::string sample_name){
	HASH_S_I::iterator it = this->sample_name2idx.find(sample_name);
	if( it==this->sample_name2idx.end() )
		return -1;
	else
		return( (*it).second );
}

