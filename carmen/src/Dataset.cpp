#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <stdio.h>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/random.hpp>
#include "boost/regex.hpp"  
#include "boost/algorithm/string/regex.hpp"
using namespace boost;

#include "DataStructures.h"
#include "ClassMinerOptions.h"
#include "Attributes.h"
#include "Rawdata.h"
#include "Discretize.h"
#include "Dataset.h"


Dataset::Dataset(){
	this->sample_names = new std::vector<std::string>();
	this->identifiers = new std::vector<std::string>();
	this->a_idx = new std::vector<int>();
	this->b_idx = new std::vector<int>();
    this->discrete = NULL;
	limits_a = "";
    limits_b = "";
	MISSING_DATA_FLOAT=-99999;
	this->raw_data = NULL;
}

Dataset::~Dataset(){
	delete this->sample_names;
	delete this->identifiers;
	delete this->a_idx;
	delete this->b_idx;
	if( this->raw_data != NULL )
		delete this->raw_data;
	if( this->discrete != NULL)
		delete this->discrete;
}

ClassifierDataset::ClassifierDataset(){
	this->features = NULL;
	this->translate = new std::vector<int*>();
}

ClassifierDataset::~ClassifierDataset(){
	for(int i=0; i<(int)this->translate->size(); i++)
		delete this->translate->at(i);
	delete this->translate;
	if( this->features != NULL)
		delete this->features;
}

void Dataset::difference( std::vector<int>* source, std::vector<int>* diff){
	// limit source to those items NOT in diff.  Slower than a hash but who cares.
	std::vector<int> keep;
	for(std::vector<int>::iterator inner = source->begin(); inner != source->end(); ++inner){
		int d = (*inner);
		bool found_in_diff = false;
		for(std::vector<int>::iterator iter = diff->begin(); iter != diff->end(); ++iter){
			if( d==(*iter)){
				found_in_diff=true;
				break;
			}
		}
		if(!found_in_diff)
			keep.push_back(d);
	}
	source->clear();
	for(std::vector<int>::iterator iter = keep.begin(); iter != keep.end(); ++iter)
		source->push_back(*iter);
}


void Dataset::intersection( std::vector<int>* source, std::vector<int>* diff){
	// limit source to those items also in diff.  Slower than a hash but who cares.
	std::vector<int> keep;
	for(std::vector<int>::iterator iter = diff->begin(); iter != diff->end(); ++iter){
		int d = (*iter);
		for(std::vector<int>::iterator inner = source->begin(); inner != source->end(); ++inner){
			if( d==(*inner)){
				keep.push_back(d);
				break;
			}
		}
	}
	source->clear();
	for(std::vector<int>::iterator iter = keep.begin(); iter != keep.end(); ++iter)
		source->push_back(*iter);
}

bool Dataset::limits_contain_gene_restriction(std::string class_limit){
	// check whether class limit contains a limit on genes
	// identified by an attribute token that starts out "gene:"
	std::vector<std::string>* limits = new std::vector<std::string>();
	stringhasher::split_comma(class_limit, limits);
	std::string attribute, value;
	bool is_equal;
	for(int i=0; i<(int)limits->size(); i++){
		stringhasher::parse_limit(limits->at(i), attribute, value, is_equal);
		if(attribute.length()>6 && attribute.substr(0,5).compare("gene:")==0 ){
			delete limits;
			return true;
		}
	}
	delete limits;
	return false;
}

void Dataset::restrict_identifiers_with_gene_disc(Rawdata* D, std::string class_limit, std::vector<std::string>* sample_identifiers, Attributes* sa){
	// Called to restrict the sample identifiers when the user has defined a class by the discretization of 
	// a gene.  This must occur after we have discretized the genes, which is the reason it isn't called
	// at the top of dataset->load().  
	if( class_limit.length()==0 )
		return;
	// convert sample_identifiers to sample_idx
	std::string attribute, value, identifier;
	std::vector<int> * sample_idx = new std::vector<int>();
	std::vector<int> * sample_idx_with_dis_value = new std::vector<int>();
	
	for(int i=0; i<(int)sample_identifiers->size(); i++){		
		sample_idx->push_back( D->index_of_sample( sample_identifiers->at(i) ) );
	}
	if( sample_idx->size()==0 ){
		// if gene:foo=1 is the ONLY restriction, then we won't have a pre-existing set of items in sample_identifiers.
		// in this case, populate them just as would have been done in Attributes::find_identifiers_in_class()
		for(int i=0; i<(int)D->identifiers.size();i++){
			sample_idx->push_back(i);
		}
	}
	std::vector<std::string>* limits = new std::vector<std::string>();
	stringhasher::split_comma(class_limit, limits);

	bool is_equal;
	std::vector<std::string> * samples_with_dis_value = new std::vector<std::string>();
	for(int i=0; i<(int)limits->size(); i++){
		sample_idx_with_dis_value->clear();
		stringhasher::parse_limit(limits->at(i), attribute, value, is_equal);
		if(attribute.length()>6 && attribute.substr(0,5).compare("gene:")==0 ){
			identifier = attribute.substr(5, attribute.length() - 5);
			int row = D->identifier2idx[identifier];
			int v = boost::lexical_cast<int>(value);
			for(int col=0; col<this->discrete->dis->cols(); col++){
				if( this->discrete->dis->arr[row][col]== v )
					sample_idx_with_dis_value->push_back(col);
			}
			intersection(sample_idx, sample_idx_with_dis_value);
		}
	}
	sample_identifiers->clear();
	for(int i=0; i<(int)sample_idx->size(); i++){
		sample_identifiers->push_back( D->sample_names.at( sample_idx->at(i) ) );
	}
	delete sample_idx;
	delete sample_idx_with_dis_value;
	delete limits;
	delete samples_with_dis_value;
}


void Dataset::find_permitted_genes(HASH_S_I* permitted_genes, ClassMinerOptions* cmo, Attributes* ga){
	std::vector<std::string>* limits = new std::vector<std::string>();
	std::vector<std::string>* identifiers_in_limits = new std::vector<std::string>();
	bool is_equal;
	std::string attribute, value, identifier;

	HASH_S_I id_in_lim;
	stringhasher::split_comma(cmo->class_a, limits);
	// check each class limit for the presence of a "gene:identifier=value"
	// if found, store in hash id_in_lim, vector identifiers_in_limits.
	for(int i=0; i<(int)limits->size(); i++){
		stringhasher::parse_limit(limits->at(i), attribute, value, is_equal);
		if(attribute.length()>5 || attribute.substr(0,5).compare("gene:")==0 ){
			identifier = attribute.substr(5, attribute.length() - 5);
			id_in_lim[identifier]=1;
			identifiers_in_limits->push_back( identifier );
		}
	}

	stringhasher::split_comma(cmo->class_b, limits);
	for(int i=0; i<(int)limits->size(); i++){
		stringhasher::parse_limit(limits->at(i), attribute, value, is_equal);
		if(attribute.length()>5 || attribute.substr(0,5).compare("gene:")==0 ){
			identifier = attribute.substr(5, attribute.length() - 5);
			id_in_lim[identifier]=1;
			identifiers_in_limits->push_back( identifier );
		}
	}

	if( cmo->gene_limit.length()>0){
		// genes limited by user with -r switch.  Include only those genes that match this limit.
		std::vector<std::string>* gene_identifiers = new std::vector<std::string>();
		ga->find_identifiers_in_class(cmo->gene_limit, gene_identifiers);
		for( std::vector<std::string>::iterator giter = gene_identifiers->begin(); giter != gene_identifiers->end(); giter++)
			(*permitted_genes)[ (*giter) ] = 1;
		delete gene_identifiers;
	}
	if(identifiers_in_limits->size()>0){
		// if user is making class definition conditional on gene discretization, we'll need to pull
		// those identifiers out of the discretization.  Otherwise, we'll just find them, which is meaningless.
		if( permitted_genes->size()>0 ){
			// if user is already restricting genes anyway, there will be values in permitted_genes.
			// remove our identifiers.
			for(int i=0; i<(int)identifiers_in_limits->size(); i++){
				identifier = identifiers_in_limits->at(i);
				permitted_genes->erase(identifier);
			}
		}
		else{
			// if user is not otherwise restricting genes, allow all genes except those found in id_in_lim.
			for(int i=0; i<(int)ga->identifiers.size(); i++){
				identifier = ga->identifiers.at(i);
				if( id_in_lim.find(identifier) == id_in_lim.end() )
					(*permitted_genes)[ identifier ] = 1;
			}
		}
	}
	delete limits;
	delete identifiers_in_limits;
}

void Dataset::check_sample_integrity(Rawdata* D, std::vector<std::string>* a, std::vector<std::string>* b){
	// ensure that each sample identified in a and b is present in raw data.
	std::string identifier;
	for(int i=0; i<(int)a->size(); i++){
		identifier = a->at(i);
		if( D->sample_name2idx.find(identifier) == D->sample_name2idx.end() ){
			throw( std::string( "Identifier in Class A not found in expression file: '" + identifier + "'") );
		}
	}
	for(int i=0; i<(int)b->size(); i++){
		identifier = b->at(i);
		if( D->sample_name2idx.find(identifier) == D->sample_name2idx.end() ){
			throw( std::string("Identifier in Class B not found in expression file: '" + identifier + "'") );
		}
	}
}


void Dataset::load(Attributes* sa, Attributes* ga, std::string fn_data, std::string limit_a, std::string limit_b){
	ClassMinerOptions* cmo = new ClassMinerOptions();
	cmo->file_name_dis = fn_data;
	cmo->class_a = limit_a;
	cmo->class_b = limit_b;
	cmo->discretization = "NONE";
	load(sa, ga, cmo);
	delete cmo;
}

void Dataset::confirm_disjoint(){
	boost::unordered_map<int, int> ids_in_a, ids_in_b;
	for(int i=0; i<(int)this->a_idx->size(); i++)
		ids_in_a[this->a_idx->at(i)] = 1;
	for(int i=0; i<(int)this->b_idx->size(); i++)
		ids_in_b[this->b_idx->at(i)] = 1;
	for( boost::unordered_map<int, int>::iterator itr = ids_in_a.begin(); itr != ids_in_a.end(); itr++ )
		if( ids_in_b.find( (*itr).first ) != ids_in_b.end() )
			throw std::string( "Element in class A also present in class B" );
	for( boost::unordered_map<int, int>::iterator itr = ids_in_b.begin(); itr != ids_in_b.end(); itr++ )
		if( ids_in_a.find( (*itr).first ) != ids_in_a.end() )
			throw std::string( "Element in class A also present in class B" );
}


void Dataset::check_gene_integrity(Rawdata* R, Attributes* ga ){
	if( ga->identifiers.size() != R->identifiers.size() ){
        std::stringstream ss;
        ss << "Gene attributes and data files do not have same number of identifiers (" << ga->identifiers.size() << " vs " <<  R->identifiers.size() << ")";
		throw std::string( ss.str() );
	}
	for(int i=0; i<(int)ga->identifiers.size(); i++){
		if( ga->identifiers.at(i).compare(R->identifiers.at(i)) != 0 ){
			std::stringstream ss;
			ss << "Gene attributes and data differ in gene order at identifier " << i+1 << " data: " << R->identifiers.at(i) << " gene " << ga->identifiers.at(i);
			throw std::string(ss.str());
		}
	}
}


void Dataset::load(Attributes* sa, Attributes* ga, ClassMinerOptions* cmo){
	this->sa = sa;
	this->ga = ga;
	// Sample limits take precedence over gene limits.  If user requests ER=0, gene:foo=1
    // as a limit, we first restrict to samples with er=0.  We then discretize these 
    // samples and further restrict to samples where gene foo has value 1.

    this->limits_a = cmo->class_a;
    this->limits_b = cmo->class_b;
	std::vector<std::string>* sample_identifiers_a = new std::vector<std::string>();
	std::vector<std::string>* sample_identifiers_b = new std::vector<std::string>();
    // Identify identifers that meet any sample limits
	sa->find_identifiers_in_class( cmo->class_a, sample_identifiers_a);
	sa->find_identifiers_in_class( cmo->class_b, sample_identifiers_b);
	
	// load raw data.  Remember that for speed, NA values should be converted to 
	// the magic value -99999; otherwise a load failure occurs.  This is only 
	// true for raw values.  We use fgets instead of parsing because it's much faster.
	this->raw_data = new Rawdata();
	this->raw_data->load(cmo->file_name_dis);
	// ensure that each sample identified in a and b is present in raw data.
	// ensure that raw data and gene attributes have the same genes, in the same order.
	this->check_sample_integrity(this->raw_data, sample_identifiers_a, sample_identifiers_b);
	this->check_gene_integrity(this->raw_data, ga );

	// Make sure we don't get confused by sample identifiers being in a different 
	// order from expression identifiers.  This led to subtle bugs in the past.  Also,
	// ensure that any overlapping class members are only included once
	boost::unordered_map<int, int> prelim_hash;
	std::vector<int> preliminary_idx;
	for(int i=0; i<(int)sample_identifiers_a->size(); i++){
		prelim_hash[ this->raw_data->index_of_sample( sample_identifiers_a->at(i) ) ] = 1;
	}
	for(int i=0; i<(int)sample_identifiers_b->size(); i++){
		prelim_hash[ this->raw_data->index_of_sample( sample_identifiers_b->at(i) ) ] = 1;
	}
	for( boost::unordered_map<int, int>::iterator itr = prelim_hash.begin(); itr != prelim_hash.end(); itr++ )
		preliminary_idx.push_back( (*itr).first );
	if( cmo->discretization.compare("none") != 0 )
		this->discrete = new Discretize(this->raw_data, preliminary_idx, cmo->discretization, cmo->disc_lower, cmo->disc_upper);	

	// User may have restricted gene list; find genes we should use.
	// If no items have been added to permitted_genes, assume no restrictions
	// and mark this by setting permitted_genes to NULL
	HASH_S_I* permitted_genes = new HASH_S_I;
	this->find_permitted_genes( permitted_genes, cmo, ga );
	if(permitted_genes->size()==0){
		delete permitted_genes;
		permitted_genes = NULL;
	}
	//int n_cols = (int)this->raw_data->sample_names.size();
	for( std::vector<std::string>::iterator iter=this->raw_data->sample_names.begin(); iter != this->raw_data->sample_names.end(); iter++)
		this->sample_names->push_back( (*iter) );
	
	// further restrict classes if the user has defined classes based on the discretization of a gene
	// (e.g., gene:foo=1).  If user has not done this, these calls have no effect.
	if(this->limits_contain_gene_restriction(cmo->class_a))
		this->restrict_identifiers_with_gene_disc( this->raw_data, cmo->class_a, sample_identifiers_a, sa);
	if(this->limits_contain_gene_restriction(cmo->class_b))
		this->restrict_identifiers_with_gene_disc( this->raw_data, cmo->class_b, sample_identifiers_b, sa);
	// Load permitted genes from raw_data->identifiers into dataset->identifiers
	for(int r=0; r<(int)this->raw_data->identifiers.size(); r++){
		if( permitted_genes==NULL || permitted_genes->find( this->raw_data->identifiers.at(r)) != permitted_genes->end() ){
			this->identifiers->push_back( this->raw_data->identifiers.at(r) );
		}
	}
	for(int i=0; i<(int)sample_identifiers_a->size(); i++)
		this->a_idx->push_back( this->raw_data->index_of_sample( sample_identifiers_a->at(i) ) );
	for(int i=0; i<(int)sample_identifiers_b->size(); i++)
		this->b_idx->push_back( this->raw_data->index_of_sample( sample_identifiers_b->at(i) ) );

	// if we have not set class B, mark that we don't want results (b_is_target=false) 
	// and include all non-A items in class B (so we can correctly calculate confidence) 
	this->a_is_target=true;
	this->b_is_target=true;
	if( cmo->class_a.length()==0 ){
		this->a_is_target=false;
		for(int i=0; i<(int)this->raw_data->sample_names.size(); i++)
			this->a_idx->push_back(i);
		difference(this->a_idx, this->b_idx);
	}
	else if(cmo->class_b.length()==0){
		this->b_is_target=false;
		for(int i=0; i<(int)this->raw_data->sample_names.size(); i++)
			this->b_idx->push_back(i);
		difference(this->b_idx, this->a_idx);
	}
	this->confirm_disjoint();
	// if any limit has * in it, convert this to !
	boost::algorithm::replace_all(cmo->class_a, std::string("*"), std::string("!"));
	boost::algorithm::replace_all(cmo->class_b, std::string("*"), std::string("!"));
	delete sample_identifiers_a;
	delete sample_identifiers_b;
	if( permitted_genes != NULL )
		delete permitted_genes;
}

void Dataset::write_discretized(std::string fn, std::string separator, std::string discretization, float disc_lower, float disc_upper){
	if( this->raw_data == NULL ){
		throw std::string("Cannot write discretized if not yet loaded");
	}
	Discretize* discrete = new Discretize(this->raw_data, *this->a_idx, discretization, disc_lower, disc_upper);
	discrete->dis->write_with_names(fn, separator, this->sa->identifiers, this->ga->identifiers );
	delete discrete;
}


void ClassifierDataset::load(Attributes* sa, Attributes* ga, ClassMinerOptions* cmo){

	// Sample limits take precedence over gene limits.  If user requests ER=0, gene:foo=1
    // as a limit, we first restrict to samples with er=0.  We then discretize these 
    // samples and further restrict to samples where gene foo has value 1.
	this->sa = sa;
	this->ga = ga;
    this->limits_a = cmo->class_a;
    this->limits_b = cmo->class_b;
	std::vector<std::string>* sample_identifiers_a = new std::vector<std::string>();
	std::vector<std::string>* sample_identifiers_b = new std::vector<std::string>();
	
    // Identify identifers that meet any sample limits
    sa->find_identifiers_in_class( cmo->class_a, sample_identifiers_a);
	sa->find_identifiers_in_class( cmo->class_b, sample_identifiers_b);

	// load raw data.  Remember that for speed, NA values should be converted to 
	// the magic value -99999; otherwise a load failure occurs.  This is only 
	// true for raw values.  We use fgets instead of parsing because it's much faster.
	this->raw_data = new Rawdata();
	this->raw_data->load(cmo->file_name_dis);

	// ensure that each sample identified in a and b is present in raw data.
	this->check_sample_integrity(this->raw_data, sample_identifiers_a, sample_identifiers_b);
	// Make sure we don't get confused by sample identifiers being in a different 
	// order from expression identifiers.  This led to subtle bugs in the past.  Also,
	// ensure that any overlapping class members are only included once
	boost::unordered_map<int, int> prelim_hash;
	std::vector<int> preliminary_idx;
	for(int i=0; i<(int)sample_identifiers_a->size(); i++){
		prelim_hash[ this->raw_data->index_of_sample( sample_identifiers_a->at(i) ) ] = 1;
	}
	for(int i=0; i<(int)sample_identifiers_b->size(); i++){
		prelim_hash[ this->raw_data->index_of_sample( sample_identifiers_b->at(i) ) ] = 1;
	}
	for( boost::unordered_map<int, int>::iterator itr = prelim_hash.begin(); itr != prelim_hash.end(); itr++ )
		preliminary_idx.push_back( (*itr).first );
	this->discrete = new Discretize(this->raw_data, preliminary_idx, cmo->discretization, cmo->disc_lower, cmo->disc_upper);	
	// User may have restricted gene list; find genes we should use.
	// If no items have been added to permitted_genes, assume no restrictions
	// and mark this by setting permitted_genes to NULL
	HASH_S_I* permitted_genes = new HASH_S_I;
	this->find_permitted_genes( permitted_genes, cmo, ga );
	if(permitted_genes->size()==0){
		delete permitted_genes;
		permitted_genes = NULL;
	}
	int n_cols = (int)this->raw_data->sample_names.size();
	int c;
	for( std::vector<std::string>::iterator iter=this->raw_data->sample_names.begin(); iter != this->raw_data->sample_names.end(); iter++)
		this->sample_names->push_back( (*iter) );
	std::vector<int*> rows;
	boost::unordered_map<int, int> unique;
	int orig_row=0;
	for(int r=0; r<(int)this->raw_data->identifiers.size(); r++){
		unique.clear();  // unique tracks the unique discretized values.
		if( permitted_genes==NULL || permitted_genes->find( this->raw_data->identifiers.at(r)) != permitted_genes->end() ){
			// if not restricting genes, or this is an allowed gene,
			// for each unique row-value combination, push an int[2] onto this->translate
			// where elem[0] is the original row and elem[1] is the value in that row
			this->identifiers->push_back( this->raw_data->identifiers.at(r) );
			int * cols = new int[n_cols];
			if( this->raw_data->data->row_has_missing_value(r) ){
				for(c=0; c<n_cols; c++){
					if( this->raw_data->data->is_missing(r,c) )
						cols[c]=0;
					else{
						cols[c]=this->discrete->dis->arr[r][c];
						if( cols[c] != 0 )
							unique[ cols[c] ]=1;
					}
				}
			}
			else{
				for(c=0; c<n_cols; c++){
					cols[c]=this->discrete->dis->arr[r][c];
					if( cols[c] != 0 )
						unique[ cols[c] ]=1;
				}
			}
			
			for( boost::unordered_map<int, int>::iterator iter = unique.begin(); iter != unique.end(); iter++ ){
				int* trans = new int[2];
				trans[0] = orig_row;
				trans[1] = (*iter).first;	// value in that row
				this->translate->push_back( trans );
			}
			++orig_row;
			rows.push_back(cols);	// push the discretized values onto rows
		}
	}
	// now that we know the number of rows in this->translate, define this->features
	// features has a row for each discretized value (which is probably much larger 
	// than the number of samples) and a 1 at i,j if that feature is present.
	int n_rows = (int)this->translate->size();
	this->features = new Matrix<int>(n_rows, n_cols, this->raw_data->data->has_missing_values);
	int filter_value;
	int* col = NULL;
	for(int dis_row=0; dis_row<n_rows; dis_row++){
		orig_row = this->translate->at(dis_row)[0];
		filter_value = this->translate->at(dis_row)[1];
		col = rows.at(orig_row);
		for(int c=0; c<n_cols; c++){
			if( col[c] == filter_value ){
				this->features->arr[dis_row][c] = 1;
			}
		}
	}

	// further restrict classes if the user has defined classes based on the discretization of a gene
	// (e.g., gene:foo=1).  If user has not done this, these calls have no effect.
	if(this->limits_contain_gene_restriction(cmo->class_a))
		this->restrict_identifiers_with_gene_disc( this->raw_data, cmo->class_a, sample_identifiers_a, sa);
	if(this->limits_contain_gene_restriction(cmo->class_b))
		this->restrict_identifiers_with_gene_disc( this->raw_data, cmo->class_b, sample_identifiers_b, sa);
	int idx;
	for(int i=0; i<(int)sample_identifiers_a->size(); i++){
		idx = this->raw_data->index_of_sample( sample_identifiers_a->at(i) );
		this->a_idx->push_back( idx );
	}
	for(int i=0; i<(int)sample_identifiers_b->size(); i++){
		idx = this->raw_data->index_of_sample( sample_identifiers_b->at(i) );
		this->b_idx->push_back( idx );
	}
	// if we have not set class B, mark that we don't want results (b_is_target=false) 
	// and include all non-A items in class B (so we can correctly calculate confidence) 
	this->a_is_target=true;
	this->b_is_target=true;
	if( cmo->class_a.length()==0 ){
		this->a_is_target=false;
		for(int i=0; i<(int)this->raw_data->sample_names.size(); i++)
			this->a_idx->push_back(i);
		difference(this->a_idx, this->b_idx);
	}
	else if(cmo->class_b.length()==0){
		this->b_is_target=false;
		for(int i=0; i<(int)this->raw_data->sample_names.size(); i++)
			this->b_idx->push_back(i);
		difference(this->b_idx, this->a_idx);
	}
	this->confirm_disjoint();

	for(int i=0; i<(int)rows.size(); i++)
		delete rows.at(i);
	delete sample_identifiers_a;
	delete sample_identifiers_b;
	if( permitted_genes != NULL )
		delete permitted_genes;
}


void ClassifierDataset::load_original(Attributes* sa, Attributes* ga, ClassMinerOptions* cmo){

	// Find samples to use in each class
	// If there are any gene: restrictions in cmo->class_a or cmo->class_b 
	// (e.g., "AURKA is up") we will adjust classes after performing the 
	// discretization.  We shouldn't try to classify the sample using both 
	// sample properties (e.g. ER is Yes) and gene: restrictions.
	this->sa = sa;
	this->ga = ga;
    this->limits_a = cmo->class_a;
    this->limits_b = cmo->class_b;
	std::vector<std::string>* sample_identifiers_a = new std::vector<std::string>();
	std::vector<std::string>* sample_identifiers_b = new std::vector<std::string>();
	sa->find_identifiers_in_class( cmo->class_a, sample_identifiers_a);
	sa->find_identifiers_in_class( cmo->class_b, sample_identifiers_b);
	
	// load raw data.  Remember that for speed, NA values should be converted to 
	// the magic value -99999; otherwise a load failure occurs.  This is only 
	// true for raw values.  We use fgets instead of parsing because it's much faster.
	this->raw_data = new Rawdata();
	this->raw_data->load(cmo->file_name_dis);

	// ensure that each sample identified in a and b is present in raw data.
	this->check_sample_integrity(this->raw_data, sample_identifiers_a, sample_identifiers_b);
	
	// Make sure we don't get confused by sample identifiers being in a different 
	// order from expression identifiers.  This led to subtle bugs in the past.  Also,
	// ensure that any overlapping class members are only included once
	boost::unordered_map<int, int> prelim_hash;
	std::vector<int> preliminary_idx;
	for(int i=0; i<(int)sample_identifiers_a->size(); i++){
		prelim_hash[ this->raw_data->index_of_sample( sample_identifiers_a->at(i) ) ] = 1;
	}
	for(int i=0; i<(int)sample_identifiers_b->size(); i++){
		prelim_hash[ this->raw_data->index_of_sample( sample_identifiers_b->at(i) ) ] = 1;
	}
	for( boost::unordered_map<int, int>::iterator itr = prelim_hash.begin(); itr != prelim_hash.end(); itr++ )
		preliminary_idx.push_back( (*itr).first );
	this->discrete = new Discretize(this->raw_data, preliminary_idx, cmo->discretization, cmo->disc_lower, cmo->disc_upper);	

	// User may have restricted gene list; find genes we should use.
	// If no items have been added to permitted_genes, assume no restrictions
	// and mark this by setting permitted_genes to NULL
	HASH_S_I* permitted_genes = new HASH_S_I;
	this->find_permitted_genes( permitted_genes, cmo, ga );
	if(permitted_genes->size()==0){
		delete permitted_genes;
		permitted_genes = NULL;
	}
	int n_cols = (int)this->raw_data->sample_names.size();
	int c;
	for( std::vector<std::string>::iterator iter=this->raw_data->sample_names.begin(); iter != this->raw_data->sample_names.end(); iter++)
		this->sample_names->push_back( (*iter) );
	std::vector<int*> rows;
	boost::unordered_map<int, int> unique;
	int orig_row=0;
	for(int r=0; r<(int)this->raw_data->identifiers.size(); r++){
		unique.clear();  // unique tracks the unique discretized values.
		if( permitted_genes==NULL || permitted_genes->find( this->raw_data->identifiers.at(r)) != permitted_genes->end() ){
			// if not restricting genes, or this is an allowed gene,
			// for each unique row-value combination, push an int[2] onto this->translate
			// where elem[0] is the original row and elem[1] is the value in that row
			this->identifiers->push_back( this->raw_data->identifiers.at(r) );
			int * cols = new int[n_cols];
			if( this->raw_data->data->row_has_missing_value(r) ){
				for(c=0; c<n_cols; c++){
					if( this->raw_data->data->is_missing(r,c) )
						cols[c]=0;
					else{
						cols[c]=this->discrete->dis->arr[r][c];
						if( cols[c] != 0 )
							unique[ cols[c] ]=1;
					}
				}
			}
			else{
				for(c=0; c<n_cols; c++){
					cols[c]=this->discrete->dis->arr[r][c];
					if( cols[c] != 0 )
						unique[ cols[c] ]=1;
				}
			}
			
			for( boost::unordered_map<int, int>::iterator iter = unique.begin(); iter != unique.end(); iter++ ){
				int* trans = new int[2];
				trans[0] = orig_row;
				trans[1] = (*iter).first;	// value in that row
				this->translate->push_back( trans );
			}
			++orig_row;
			rows.push_back(cols);	// push the discretized values onto rows
		}
	}

	// now that we know the number of rows in this->translate, define this->features
	// features has a row for each discretized value (which is probably much larger 
	// than the number of samples) and a 1 at i,j if that feature is present.
	int n_rows = (int)this->translate->size();
	this->features = new Matrix<int>(n_rows, n_cols, this->raw_data->data->has_missing_values);
	int filter_value;
	int* col = NULL;
	for(int dis_row=0; dis_row<n_rows; dis_row++){
		orig_row = this->translate->at(dis_row)[0];
		filter_value = this->translate->at(dis_row)[1];
		col = rows.at(orig_row);
		for(int c=0; c<n_cols; c++){
			if( col[c] == filter_value ){
				this->features->arr[dis_row][c] = 1;
			}
		}
	}
	
	// further restrict classes if the user has defined classes based on the discretization of a gene
	// (e.g., gene:foo=1).  If user has not done this, these calls have no effect.
	this->restrict_identifiers_with_gene_disc( this->raw_data, cmo->class_a, sample_identifiers_a, sa);
	this->restrict_identifiers_with_gene_disc( this->raw_data, cmo->class_b, sample_identifiers_b, sa);
	int idx;
	for(int i=0; i<(int)sample_identifiers_a->size(); i++){
		idx = this->raw_data->index_of_sample( sample_identifiers_a->at(i) );
		this->a_idx->push_back( idx );
	}

	for(int i=0; i<(int)sample_identifiers_b->size(); i++){
		idx = this->raw_data->index_of_sample( sample_identifiers_b->at(i) );
		this->b_idx->push_back( idx );
	}

	// if we have not set class B, mark that we don't want results (b_is_target=false) 
	// and include all non-A items in class B (so we can correctly calculate confidence) 
	this->a_is_target=true;
	this->b_is_target=true;
	if( cmo->class_a.length()==0 ){
		this->a_is_target=false;
		for(int i=0; i<(int)this->raw_data->sample_names.size(); i++)
			this->a_idx->push_back(i);
		difference(this->a_idx, this->b_idx);
	}
	else if(cmo->class_b.length()==0){
		this->b_is_target=false;
		for(int i=0; i<(int)this->raw_data->sample_names.size(); i++)
			this->b_idx->push_back(i);
		difference(this->b_idx, this->a_idx);
	}
	// last check: make sure that a_idx and b_idx are disjoint
	confirm_disjoint();

	for(int i=0; i<(int)rows.size(); i++)
		delete rows.at(i);
	delete sample_identifiers_a;
	delete sample_identifiers_b;
	if( permitted_genes != NULL )
		delete permitted_genes;
}


void ClassifierDataset::samples_with_gene_value(std::string identifier, int value, std::vector<int>* sample_idx){
	// populate sample_idx with the index of each sample that has a discretized value of "value"
	sample_idx->clear();
	int row=-1;
	int r=0;
	for(std::vector<std::string>::iterator iter=this->raw_data->identifiers.begin(); iter!=this->raw_data->identifiers.end(); ++iter){
		if( (*iter).compare(identifier)==0){
			row=r;
			break;
		}
		++r;
	}
	if( row==-1 ){
		throw std::string("No gene in dataset with name ") + identifier;
	}
	r=0;
	int dis_row=-1;
	for(std::vector<int*>::iterator iter=this->translate->begin(); iter!=this->translate->end(); ++iter){
		if((*iter)[0]==row && (*iter)[1]==value){
			dis_row = r;
			break;
		}
		++r;
	}
	if(dis_row>-1){
		for(int col=0; col < (int)this->raw_data->sample_names.size(); ++col){
			if( this->features->arr[dis_row][col] == 1 )
				sample_idx->push_back(col);
		}
	}
}

extract_sorter::extract_sorter(double val, std::string val_str, int index){
	this->val = val;
	this->val_str = val_str;
	this->idx = index;
}

bool extract_sorter_by_double_value(const extract_sorter* a, const extract_sorter* b){
	return a->val < b->val;
}


void Dataset::grouped_sort(std::vector<int>& idx_sorted, std::vector<std::string> samples, int idx_probe, std::string group_by, bool sort_within_groups){
	// This function populates idx_sorted by pushing extract_sorter pointers onto
	// the vector sorter and then sorting subsets of these vectors.
    std::vector<std::string> g_values;
	std::vector<double> g_values_as_double;
	bool all_group_labels_are_numeric=true;
	std::string v;

	// if all of the group labels that are not missing can be cast to doubles,
	// treat the string values as numeric and sort using numeric order rather than
	// alphabetical order.
	for( int j=0; j<(int)samples.size(); j++ ){
		v = sa->prop_for_identifier(samples.at(j), group_by );
		if( v.compare( sa->no_result ) == 0 ){
			g_values_as_double.push_back(-99999);
			g_values.push_back( std::string("NA") );
		}
		else{
			g_values.push_back( v );
			try{ g_values_as_double.push_back( boost::lexical_cast<double>( v ) ); }
		    catch( boost::bad_lexical_cast &){
		    	all_group_labels_are_numeric = false;
		    	g_values_as_double.push_back(0);
		    }
		}
	}
    std::vector<extract_sorter*> sorter;
	// sort results by group
	if(all_group_labels_are_numeric){
		for( int j=0; j<(int)g_values.size(); j++ ){
			int sample_idx = this->raw_data->sample_name2idx[samples.at(j)];
			extract_sorter* e = new extract_sorter( g_values_as_double.at(j), std::string(g_values.at(j).c_str()), sample_idx );
			sorter.push_back( e );
		}
		std::sort(sorter.begin(), sorter.end(), extract_sorter_by_double_value );
	}
	else{
		// XXX this nutty code works around a very subtle bug that developed when I sent
		// certain extract_sorter payloads (e.g. all M or F) to a string version of the
		// extract_sorter comparison function. Not all string payloads failed; I couldn't divine why.
		// Most likely I am at fault, but it's some pointer-related nightmare. Rather than spend
		// another day on this I coded a workaround that creates a new array 'to_sort' where each
		// item is a string made of value_to_sort<<>>index_of_that_value. This is a way to carry along
		// the index order. I am making the assumption that <<>> will never appear in a sample attribute
		// that I want to sort.
		std::vector<std::string> to_sort, getsplit;
		for(int i=0; i<(int)g_values.size(); i++){
			std::stringstream ss;
			ss << g_values.at(i) << "MATCHTHISPLEASE" << this->raw_data->sample_name2idx[samples.at(i)];
			to_sort.push_back(ss.str());
		}
		std::sort(to_sort.begin(), to_sort.end() );
        for(int i=0; i<(int)g_values.size(); i++){
			boost::algorithm::split_regex( getsplit, to_sort.at(i), boost::regex( "MATCHTHISPLEASE" ) ) ;
            //stringhasher::split(to_sort.at(i), &getsplit, prong);
            int idx_from_split = boost::lexical_cast<int>( getsplit[1] );
			extract_sorter* e = new extract_sorter( 0, getsplit[0], idx_from_split );
			sorter.push_back( e );
		}
	}

	if( sort_within_groups ){
		if((int)sorter.size()>1){
			// identify groups, sort again within each group using values
			std::string cur = sorter.at(0)->val_str;
			sorter.at(0)->val = this->raw_data->data->arr[idx_probe][sorter.at(0)->idx];
			int group_start=0;
			std::vector<int> subgroup;
			int i=1;
			while(i<(int)sorter.size()){
				if( cur.compare( sorter.at(i)->val_str )==0 )
					i++;
				else{
					for(int n=group_start; n<i; n++){
						sorter.at(n)->val = this->raw_data->data->arr[idx_probe][sorter.at(n)->idx];
					}
					extract_sorter** ptr_begin = &sorter[group_start];
					extract_sorter** ptr_end = (&sorter[i-1])+1;
					std::sort(ptr_begin, ptr_end, extract_sorter_by_double_value );
					i++;
					if( i<(int)sorter.size())
						cur = sorter.at(i)->val_str;
					group_start = i;
				}
			}

			if( group_start != (int)sorter.size() ){
				extract_sorter** ptr_begin;
				if( group_start==0 ){
					for(int n=0; n<(int)sorter.size(); n++)
						sorter.at(n)->val = this->raw_data->data->arr[idx_probe][sorter.at(n)->idx];
					ptr_begin = &sorter[0];
				}
				else{
					for(int n=group_start-1; n<(int)sorter.size(); n++)
						sorter.at(n)->val = this->raw_data->data->arr[idx_probe][sorter.at(n)->idx];
					ptr_begin = &sorter[group_start-1];
				}
				// for next line, want one after sorter.size()
				extract_sorter** ptr_end = (&sorter[ (int)sorter.size()-1 ]) + 1 ;
				std::sort(ptr_begin, ptr_end, extract_sorter_by_double_value );
			}
		}
	}
	for(int i=0; i<(int)sorter.size(); i++){
		idx_sorted.push_back( sorter.at(i)->idx );
		delete sorter.at(i);
	}
}


void Dataset::extract(std::vector<std::vector<double>*>& values,
		std::vector<std::string>& labels,
		std::vector<std::string>& genes,
		std::vector<std::string> probes,
		std::string sample_limit,
		std::string group_by,
		bool sort_by_value,
		std::string label_attribute){
	// populate values and labels with data from Dataset
	// restrict samples by sample_limit, no restriction if string of length 0
	// if group_by, sort samples into groups by value of group_by. Default to no groups.
	//   If group_by is passed, labels are group values
	// if sort_by_value, sort within groups by first probe.

	for(int i=0; i<(int)values.size(); i++)
		if( values.at(i) != NULL )
			delete values.at(i);
	values.clear();
	labels.clear();
	genes.clear();

	std::vector<std::string> samples;
	sa->find_identifiers_in_class( sample_limit, &samples);

	for( int i=0; i<(int)probes.size(); i++ ){
		if( ga->identifier2idx.find(probes.at(i)) == ga->identifier2idx.end() )
			throw std::string("probe not found:" + probes.at(i));
		genes.push_back( this->ga->prop_for_identifier(probes.at(i), ga->get_gene_name_column() ) );
	}
	if( group_by.size() > 0){
		if( sa->attribute2idx.find(group_by) == sa->attribute2idx.end() )
			throw std::string("attribute not found: " + group_by );
	}

	std::vector<int> idx;
	if( group_by.size()==0  ){
		if( sort_by_value ){
			// sort idx by values of first probe.
			// labels are still sample IDENTIFIERs
			std::vector<extract_sorter*> sorter;
			int idx_first_probe = this->raw_data->identifier2idx[probes.at(0)];
			for( int j=0; j<(int)samples.size(); j++ ){
				int sample_idx = this->raw_data->index_of_sample(samples.at(j));
				extract_sorter* e = new extract_sorter( this->raw_data->data->arr[idx_first_probe][sample_idx], std::string(""), sample_idx );
				sorter.push_back( e );
			}
			std::sort(sorter.begin(), sorter.end(), extract_sorter_by_double_value );
			idx.clear();
			for( int j=0; j<(int)sorter.size(); j++ ){
				idx.push_back(sorter.at(j)->idx);
				delete sorter.at(j);
			}
		}
		else{
			// no sorting, present in original order
			for(int i=0; i<(int)samples.size(); i++){
				idx.push_back( this->raw_data->index_of_sample( samples.at(i) ) );
			}
		}
	}
	else{
		// sort idx by groups, also potentially by first probe value within groups
		// labels are group values
		int idx_probe = this->raw_data->identifier2idx[probes.at(0)];
		std::vector<int> idx_ordered;
		grouped_sort(idx_ordered, samples, idx_probe, group_by, sort_by_value );
		for( int j=0; j<(int)idx_ordered.size(); j++ ){
			idx.push_back(idx_ordered.at(j));
		}
	}

	// common code to push values onto values stack
	for(int i=0; i<(int)probes.size(); i++ ){
		int idx_probe = this->raw_data->identifier2idx[probes.at(i)];
		std::vector<double>* V = new std::vector<double>();
		if( this->raw_data->data->row_has_missing_value(idx_probe) ){
			for( int j=0; j<(int)idx.size(); j++ ){
				if( this->raw_data->data->is_missing(idx_probe, idx.at(j) ))
					V->push_back(-99999);
				else
					V->push_back( this->raw_data->data->arr[idx_probe][idx.at(j)] );
			}
		}
		else{
			for( int j=0; j<(int)idx.size(); j++ ){
				V->push_back( this->raw_data->data->arr[idx_probe][idx.at(j)] );
			}
		}
		values.push_back(V);
	}
	for(int i=0; i<(int)idx.size(); i++){
		std::string identifier = this->raw_data->sample_names[ idx.at(i) ];
		labels.push_back( this->sa->prop_for_identifier(identifier, label_attribute ) );
	}
}


