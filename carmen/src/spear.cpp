#include <iostream>
#include <iomanip>
#include <vector>
#include <set>
#include <boost/unordered_map.hpp>
#include <boost/random.hpp>
#include <algorithm>
#include "math.h"
using namespace std;
using namespace boost;
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/thread.hpp>
#include <boost/bind.hpp>
#include <boost/thread/mutex.hpp>
#include "DataStructures.h"
#include "ClassMinerOptions.h"
#include "Rawdata.h"
#include "Attributes.h"
#include "Discretize.h"
#include "Dataset.h"
#include "graph.h"
#include "Parser.h"
#include "spear.h"

/*
 * David Quigley
 * UCSF Hellen Diller Family Comprehensive Cancer Center, 2007-2014
 *
 * Spear is used for various operations that include a measure of correlation between
 * rows of elements. When Spearman rank correlation is calculated ties are ranked using
 * the mean of the ranks that are spanned, so (1,2,2) is ranked (1, 2.5, 2.5).
 * By passing in only one class, you can find genes that are
 * correlated in a single group.  By passing in two non-intersecting classes,
 * you can find genes whose correlation changes between two groups.
 *
 * OVERVIEW OF HOW THE CODE IS CALLED
 *
 * The user instantiates a Spearman object and specifies a number of parameters by
 * calls to set_ functions which modify local variables. This includes the location of
 * data files and restrictions on the kind of correlation results to return. Two
 * basic kinds of operations can be performed: calculations of which pairs of probes
 * meet certain correlation criteria (e.g. all probes correlated with a given single
 * probe where r >= +/- 0.8) and calculations that define a significance level for
 * correlation (e.g. the GWER when using all probes in a data set).
 *
 * When calculating correlation, results are stored in Spear objects, which are simple
 * data-storage structs.
 *
 * Two output formats are commonly used: a text file with the extension .spear or
 * a set of files in the format used by Cytoscape (cytoscape.org).
 */

// CACHEHASH IS CODE THAT WAS ULTIMATELY NOT USED
//
// CacheHash is a data structure designed to store the results of recalculating
// spearman correlation ranks due to missing data. Each time we compare two probes
// we must remove the samples where either probe has missing data and then re-rank
// the remaining values. This is very expensive, so we store the rank for a given
// probe (indexed by the integer idx) with a given set of missing values (indexed by vals)
// The hope is that we will frequently see the same missing data for a given sample,
// and therefore save time.
CacheHash::CacheHash( Matrix<float>* M, std::vector<int> valid_cols ){
	// for each row, identify which columns (if any) are missing and
	// store their columns in a set. We will then use this a lookup
	if( !M->has_missing_values )
		return;
	int NC = int(valid_cols.size());
	for(int i=0; i<NC; i++)
		if( valid_cols.at(i) >= M->cols() )
			throw std::string("Out of bounds valid_cols value in CacheHash constructor");
	for(int r=0; r<M->rows(); r++){
		std::set<int>* miss = new std::set<int>();
		if( M->row_has_missing_value(r) ){
			for(int i=0; i<NC; i++){
				if( M->is_missing(r, valid_cols.at(i) ) )
					miss->insert(valid_cols.at(i));
			}
		}
		this->idx2missing[ r ] = miss;
	}
}

std::set<int>* CacheHash::missing_by_idx( int idx ){
	if( this->idx2missing.find(idx) != this->idx2missing.end() ){
		return this->idx2missing[idx];
	}
	else{
		return NULL;
	}
}

CacheHash::~CacheHash(){
	std::vector<int> to_clear;
	for( boost::unordered_map<int, HASH_SIZE_VEC*>::iterator itr=this->hsh.begin(); itr!=hsh.end(); itr++){
		to_clear.push_back( (*itr).first );
	}
	for(int i=0; i<int(to_clear.size()); i++)
		this->clear(to_clear.at(i));
	for( boost::unordered_map< int, std::set<int>* >::iterator i=this->idx2missing.begin(); i!=this->idx2missing.end(); i++){
		delete (*i).second;
	}
}

bool CacheHash::has(int idx, std::set<int>& vals){
	if( this->hsh.find(idx)==this->hsh.end() )
		return false;
	else{
		if( this->hsh[idx]->find(boost::hash_value(vals)) == this->hsh[idx]->end() )
			return false;
		else
			return true;
	}
}

std::vector<int>* CacheHash::retreive(int idx, std::set<int>& vals){
	if( this->has( idx, vals ) ){
		HASH_SIZE_VEC* P = this->hsh[idx];
		return (*P)[boost::hash_value(vals)];
	}
	else{
		return NULL;
	}
}


void CacheHash::store( int idx, std::set<int>& vals, std::vector<int>* ranks ){
	HASH_SIZE_VEC* v;
	if( this->hsh.find(idx)==this->hsh.end() ){
		v = new HASH_SIZE_VEC;
		(*v)[ boost::hash_value(vals) ] = ranks;
		hsh[idx] = v;
	}
	else{
		v = this->hsh[idx];
		(*v)[ boost::hash_value(vals) ] = ranks;
	}
}

void CacheHash::clear(int idx){
	HASH_SIZE_VEC* v;
	if( this->hsh.find(idx) != this->hsh.end() ){
		v = this->hsh[idx];
		for(HASH_SIZE_VEC::iterator i=v->begin(); i != v->end(); i++){
			delete (*i).second;
		}
		delete v;
		hsh.erase(idx);
	}
}

RewiringResult::RewiringResult(int row_1, double z_sum, double z_increase, double z_decrease, double p_perm,
                               int n_z_lt5, int n_z_5_6, int n_z_6_7, int n_z_7_8, int n_z_8_9, int n_z_gt9 ) : row1(row_1),  zsum(z_sum), zinc(z_increase), zdec(z_decrease), pperm(p_perm), z_lt5(n_z_lt5), z_5_6(n_z_5_6), z_6_7(n_z_6_7), z_7_8(n_z_7_8), z_8_9(n_z_8_9), z_gt9(n_z_gt9){
    
}

int RewiringResult::row_1(){
    return this->row1;
}

double RewiringResult::z_sum(){
    return this->zsum;
}

double RewiringResult::z_increase(){
    return this->zinc;
}

double RewiringResult::z_decrease(){
    return this->zdec;
}

int RewiringResult::n_z_lt5(){
    return this->z_lt5;
}
int RewiringResult::n_z_5_6(){
    return this->z_5_6;
}
int RewiringResult::n_z_6_7(){
    return this->z_6_7;
}
int RewiringResult::n_z_7_8(){
    return this->z_7_8;
}
int RewiringResult::n_z_8_9(){
    return this->z_8_9;
}
int RewiringResult::n_z_gt9(){
    return this->z_gt9;
}
double RewiringResult::p_perm(){
    return this->pperm;
}

void RewiringResult::increment_p_perm(double value){
    this->pperm += value;
}

/*
 * Spear objects are a storage struct with no business intelligence of their own.
 * They are accumulated over the course of a call to run().
 * Note that to save memory, correlation values are stored as integers. This matters if we
 * are storing very large numbers of Spear objects. Since real-world use indicates that our
 * correlation results are only meaningful to two or three decimal places, this should have
 * no impact on accuracy
 */
Spear::Spear(double rho_a, double rho_b, int row1, int row2, double z_score) : rhoa((int)(rho_a * 100000)), rhob((int)(rho_b * 100000)), row1(row1), row2(row2), zscore((int)(z_score* 100000)), perm_pvalue(-1) {
}

Spear::Spear(double rho_a, double rho_b, int row1, int row2, double z_score, double perm_pvalue) : rhoa((int)(rho_a * 100000)), rhob((int)(rho_b * 100000)), row1(row1), row2(row2), zscore((int)(z_score* 100000)), perm_pvalue(perm_pvalue){
}

Spear::Spear(Spear* s): rhoa((int)(s->rho_a() * 100000)), rhob((int)(s->rho_b() * 100000)), row1(s->row_1()), row2(s->row_2()), zscore((int)(zscore* 100000)), perm_pvalue(s->perm_p_value()){
}


const double Spear::rho_a(){
	return ((double)this->rhoa) / 100000; // rhoa is stored internally as an int
}

const double Spear::rho_b(){
	return ((double)this->rhob) / 100000; // rhob is stored internally as an int
}

const double Spear::z_score(){
	return ((double)this->zscore) / 100000; // z_score is stored internally as an int
}

const int Spear::row_1(){
	return this->row1;
}

const int Spear::row_2(){
	return this->row2;
}

const double Spear::perm_p_value(){
	return this->perm_pvalue;
}

const double Spear::delta(){
	double a = ((double)this->rhoa) / 100000;
	double b = ((double)this->rhob) / 100000;
	return a - b;
}

Sorter::Sorter(double val, double index){
	this->v[0] = val;
	this->v[1] = index;
	this->v[2] = 0;
}

bool sorter_by_value(const Sorter* a, const Sorter* b){
	return a->v[0] < b->v[0];
}

bool sorter_orig(const Sorter* a, const Sorter* b){
	return a->v[1] < b->v[1];
}

bool spears_sort(const Spear* a, const Spear* b){
	double da = ( ((double)a->rhoa) / 100000 ) - (((double)a->rhob) / 100000 );
	double db = ( ((double)b->rhoa) / 100000 ) - (((double)b->rhob) / 100000 );
	return abs(da) > abs(db);
}

bool spears_sort_asc(const Spear* a, const Spear* b){
	double da = abs( ((double)a->rhoa) / 100000 ) - (((double)a->rhob) / 100000 );
	double db = abs( ((double)b->rhoa) / 100000 ) - (((double)b->rhob) / 100000 );
	return abs(da) < abs(db);
}

DifferentalCorrelationResult::DifferentalCorrelationResult( std::string title, int N, double meanA, double meanB, double mean_diff_obs, double pval){
	this->title=title;
	this->N = N;
	this->meanA = meanA;
	this->meanB = meanB;
	this->meanDiff = mean_diff_obs;
	this->pval = pval;
}

ProbeSet::ProbeSet(std::string title, std::vector<std::string>& passed_probes, std::vector<std::string>& names){
	this->title = title;
	for(int i=0; i<int(passed_probes.size()); i++)
		this->probelist.push_back( passed_probes.at(i) );
	for(int i=0; i<int(names.size()); i++)
		this->probe_names.push_back( names.at(i) );
}


ProbeSets::ProbeSets(){
	this->itr_counter=0;
}


void ProbeSets::add_from_genelist( std::string set_title, Attributes* ga, std::vector<std::string> genelist, int min_probe_set_size ){
	// a genelist may contain EITHER gene names or probe identifiers.
	std::string gene_column= ga->get_gene_name_column();
	std::vector<int> v;
	std::vector<std::string> probelist;
	for(int i=0; i<int(genelist.size()); i++ ){
		if( ga->identifier2idx.find(genelist.at(i)) == ga->identifier2idx.end() ){
			ga->indices_with_property(gene_column , genelist.at(i), &v, false);
			for(int k=0; k<int(v.size()); k++ ){
				probelist.push_back( ga->identifiers.at( v.at(k) ));
			}
		}
		else{
			probelist.push_back( ga->identifiers.at( ga->identifier2idx[genelist.at(i)] ));
		}
	}
	if( int(probelist.size()) >= min_probe_set_size){
		this->probesets.push_back( new ProbeSet( set_title, probelist, genelist) );
	}
}


void ProbeSets::add_from_probelist(std::string set_title, Attributes* ga, std::vector<std::string> probelist, int min_probe_set_size ){

	for(int i=0; i<int(probelist.size()); i++ ){
		if( ga->identifier2idx.find(probelist.at(i)) == ga->identifier2idx.end() ){
			std::stringstream ss;
			ss << "Probe " << probelist.at(i) << " not found in gene attributes";
			throw std::string( ss.str() );
		}
	}
	if( int(probelist.size()) >= min_probe_set_size){
		this->probesets.push_back( new ProbeSet( set_title, probelist, probelist) );
	}

}

void ProbeSets::clear(){
	this->reset_iterator();
	for(int i=0; i<int(this->probesets.size()); i++){
		delete this->probesets.at(i);
	}
	this->probesets.clear();
}


int ProbeSets::size(){
	return int( this->probesets.size() );
}

void ProbeSets::reset_iterator(){
	this->itr_counter=0;
}


ProbeSet* ProbeSets::next_probeset(){
	if( this->itr_counter >= int(this->probesets.size() ) ){
		return NULL;
	}
	else{
		itr_counter += 1;
		return this->probesets.at(this->itr_counter - 1);
	}
}

ProbeSets::~ProbeSets(){
	this->clear();
}


Spearman::Spearman(){
	// default values
	this->fn_expr = std::string("");;
	this->fn_sa = std::string("");;
	this->fn_ga = std::string("");;
	this->class_a = std::string("");;
	this->class_b = std::string("");;
	this->corr_type = std::string("spearman");
	this->min_var = 0;
	this->corr_abs_a = 0;
	this->corr_abs_b = 0;
	this->corr_diff = 0;
	this->max_eqtl_pval = 1;
	this->percent_required = 0;
	this->fn_out = std::string("");
	this->fn_eqtl = std::string("");
	this->n_perms = 0;
    this->current_permutation=0;
    this->n_threads=1;
	this->limit_network_to_seeds = false;
	this->min_clique_size = 0;
	this->include_seed_neighbor_correlations = false;
	this->verbose = false;
	this->data = NULL;
	this->ga = NULL;
	this->sa = NULL;
    this->permutations_idx_a=NULL;
    this->permutations_idx_b=NULL;
    this->z_score_sum_obs=0;
    this->global_p_perm=0;
	this->do_distribution=false;
	this->batch_extension = std::string("bat");
	this->fn_cytoscape = std::string("");
	this->fn_cytoscape_props = std::string("");
}


Spearman::Spearman(std::string fn_in, Attributes* ga){
	// This call will load an existing .spear file rather than instantiating a blank one
	// Note that since Spear objects store probes by their index in an Attributes object,
	// a gene attributes file must be passed in order to have a meaningful index for the probes.
	// Did not originally calculate Z score, so we need to check whether it is present or not.
	this->data = NULL;
	this->sa = NULL;
	this->ga = ga;

	CarmenParser cp(fn_in);
	cp.PrepareToReadValues();
	std::vector<std::string> values;
	double rho_a=0, rho_b=0, zscore=0;
	int row1=0, row2=0, nA=0, nB=0;
	std::vector<std::string> header;
	cp.get_column_names( header );
	int zscore_index = -1;
	for(int i=0; i<int(header.size()); i++){
	    if( header.at(i).compare("z_score")==0 ){
	    	zscore_index=i;
	    }
	}
	if( zscore_index==-1){
		std::string nA_str, nB_str;
		cp.get_value("Class_A", 1, nA_str);
		cp.get_value("Class_B", 1, nB_str);
	    try{ nA = boost::lexical_cast<int>( nA_str ); }
	    catch( boost::bad_lexical_cast &){ throw( std::string( "ERROR: unable to parse number of items in class A" ) ); }
	    try{ nB = boost::lexical_cast<int>( nB_str ); }
	    catch( boost::bad_lexical_cast &){ throw( std::string( "ERROR: unable to parse number of items in class B" ) ); }
	}
	while( cp.ReadNextValue( values) ){
		row1 = boost::lexical_cast<int>( ga->identifier2idx[ values[1] ] );
		row2 = boost::lexical_cast<int>( ga->identifier2idx[ values[3] ] );
		rho_a = boost::lexical_cast<double>( values[4] ); 
		rho_b = boost::lexical_cast<double>( values[5] );
		if( zscore_index==-1 ){
			zscore = this->fisher_zscore(rho_a, rho_b, nA, nB);
		}
		this->spears.push_back(new Spear(rho_a, rho_b, row1, row2, zscore));
	}

	this->n_perms=0;
    this->current_permutation=0;
    this->n_threads=1;
	this->limit_network_to_seeds = false;
	this->min_clique_size=0;
    this->z_score_sum_obs=0;
    this->global_p_perm=0;
	this->do_distribution=false;
	this->include_seed_neighbor_correlations = false;
	this->batch_extension = std::string("bat");
	this->fn_cytoscape = std::string("");
	this->fn_cytoscape_props = std::string("");
	this->verbose = false;
}

Spearman::~Spearman(){
	if(this->data != NULL)
		delete this->data;
	if(this->ga != NULL )
		delete this->ga;
	if(this->sa != NULL)
		delete this->sa;
	for(int i=0; i<(int)spears.size(); i++)
		delete this->spears.at(i);
    for(int i=0; i<(int)rewiring_results.size(); i++)
        delete this->rewiring_results.at(i);
    if( this->permutations_idx_a != NULL )
        delete this->permutations_idx_a;
    if( this->permutations_idx_b != NULL )
        delete this->permutations_idx_b;
	this->data = NULL;
	this->ga = NULL;
	this->sa = NULL;
}


bool Spearman::get_success(){
	return this->success;
}


int Spearman::get_n_probes_to_process(){
	return this->idx.size();
}

void Spearman::set_allow_uncorrelated_loci_with_eQTL(bool allow){
	this->allow_uncorrelated_loci_with_eQTL = allow;
}


void Spearman::set_GWER_method(std::string method){
	this->GWER_method = method;
}


void Spearman::set_do_distribution(bool do_distribution){
	this->do_distribution = do_distribution;
}

void Spearman::set_input_files(std::string fn_expr, std::string fn_sa, std::string fn_ga){
	this->fn_expr = fn_expr;
	this->fn_sa = fn_sa;
	this->fn_ga = fn_ga;
}

void Spearman::set_gene_name_column(std::string gene_name_column){
	this->gene_name_column = gene_name_column;
}

void Spearman::set_limit_a(std::string class_a){
	this->class_a = class_a;
}

void Spearman::set_limit_b(std::string class_b){
	this->class_b = class_b;
}

void Spearman::set_max_eqtl_pval(double pval){
	this->max_eqtl_pval = pval;
}

void Spearman::set_min_clique_size(int min){
	this->min_clique_size = min;
}

void Spearman::set_method(std::string corr_method){
	if(corr_type.compare("spearman")!=0 && corr_type.compare("pearson")!=0)
		throw std::string("corr_type must be one of {spearman, pearson}");
	this->corr_type = corr_method;
}

void Spearman::set_min_zscore(double min_zscore){
	this->min_zscore = min_zscore;
}

void Spearman::set_seeds(std::vector<std::string>& seeds){
    this->seeds.clear();
	for(int i=0; i<(int)seeds.size(); i++)
		this->seeds.push_back( seeds.at(i) );
}

void Spearman::set_limit_network_to_seeds( bool limit_to_seeds ){
	this->limit_network_to_seeds = limit_to_seeds;
	if(limit_to_seeds)
		this->include_seed_neighbor_correlations = false;
}

void Spearman::set_min_var(double min_var){
	if( min_var < 0)
		throw std::string( "min_var cannot be less than 0");
	this->min_var = min_var;
}

void Spearman::set_percent_required(double percent_required){
	// expects a value betwen 0 and 1
	if( percent_required > 1)
		throw std::string( "percent_present cannot be greater than 1");
	if( percent_required < 0)
		throw std::string( "percent_present cannot be less than 0");
	this->percent_required = percent_required;
}

void Spearman::set_corr_diff(double corr_diff){
	// correlation is a value between -1 and 1, so max( abs(X)-abs(Y) ) <= 2
	if( corr_diff > 2)
		throw std::string( "corr_diff cannot be greater than 2");
	if( corr_diff < 0)
		throw std::string( "corr_diff cannot be less than 0");
	this->corr_diff = corr_diff;
}

void Spearman::set_corr_abs_a(double corr_abs){
	if( corr_abs > 1)
		throw std::string( "corr_abs cannot be greater than 1");
	if( corr_abs < 0)
		throw std::string( "corr_abs cannot be less than 0");
	this->corr_abs_a = corr_abs;
}

void Spearman::set_corr_abs_b(double corr_abs){
	if( corr_abs > 1)
		throw std::string( "corr_abs cannot be greater than 1");
	if( corr_abs < 0)
		throw std::string( "corr_abs cannot be less than 0");
	this->corr_abs_b = corr_abs;
}

void Spearman::set_n_permutations(int n_perms){
	if( n_perms<0 )
		throw std::string("n_perms must be greater than 0");
	this->n_perms = n_perms;
}

void Spearman::set_include_seed_neighbor_correlations(bool include_seed_neighbor_correlations){
	this->include_seed_neighbor_correlations = include_seed_neighbor_correlations;
}

void Spearman::set_batch_extension(std::string ext){
	this->batch_extension = ext;
}

void Spearman::set_fn_cytoscape(std::string fn){
	this->fn_cytoscape = fn;
}

void Spearman::set_fn_cytoscape_props(std::string fn){
	this->fn_cytoscape_props = fn;
}

void Spearman::set_fn_eqtl( std::string fn){
	this->fn_eqtl = fn;
}

void Spearman::set_fn_out( std::string fn_out){
	this->fn_out = fn_out;
}

void Spearman::set_n_threads( int n_threads ){
    this->n_threads=n_threads;
}

void Spearman::set_column_extra_attributes(std::string col){
	this->col_extra_attributes = col;
}

void Spearman::set_verbose( bool is_verbose ){
	this->verbose = is_verbose;
}

void Spearman::set_require_eqtl(bool req){
	this->require_eqtl = req;
}

double Spearman::find_var(std::vector<double>* C){
	// Calculate the sample variance of a vector of doubles C.
	// WARNING: Does not check for missing data
	int N = (int)C->size();
	if(N<=1)
		return 0.0; // Don't want to deal with NaN
	else{
		double sum=0.0, mean=0.0, sum_of_diff_squared=0.0, diff=0.0;
		for(int i=0; i<N; i++)
			sum += C->at(i);
		mean = sum/N;
		for(int i=0; i<N; i++){
			diff = C->at(i) - mean;
			sum_of_diff_squared += diff * diff;
		}
		return ( sum_of_diff_squared / (N-1) );
	}
}


void Spearman::limit_ids_by_var(){
	// calculate variance for each probe currently specified in idx
	// if variance >= this->min_var, keep the probe
	// If class A and class B are both specified, variance across both classes is used.
	// Missing data are excluded and do not affect variance
	if(this->min_var==0){
		return;
	}
	std::vector<int> keep;
	std::vector<double>* C = new std::vector<double>();
	int n, c, val;
	for(int r=0; r<(int)this->idx.size(); r++){
		val = this->idx.at(r);
		C->clear();
		if( this->data->raw_data->data->has_missing_values ){
			for(n=0; n<(int)this->data->a_idx->size(); n++){
				c = this->data->a_idx->at(n);
				if( !this->data->raw_data->data->is_missing(val, c) ){
					C->push_back( this->data->raw_data->data->arr[val][c] );
				}
			}
			for(n=0; n<(int)this->data->b_idx->size(); n++){
				c = this->data->b_idx->at(n);
				if( !this->data->raw_data->data->is_missing(val, c) ){
					C->push_back( this->data->raw_data->data->arr[val][c] );
				}
			}
		}
		else{
			for(n=0; n<(int)this->data->a_idx->size(); n++){
				C->push_back( this->data->raw_data->data->arr[val][this->data->a_idx->at(n)] );
			}
			for(n=0; n<(int)this->data->b_idx->size(); n++){
				C->push_back( this->data->raw_data->data->arr[val][this->data->b_idx->at(n)] );
			}
		}
		if( int(C->size()) > 1 ){
			if( this->find_var(C) >= this->min_var )
				keep.push_back(val);
		}
	}
	this->idx.clear();
	for(int i=0; i<(int)keep.size(); i++)
		idx.push_back( keep.at(i) );
	delete C;
}

void Spearman::limit_ids_by_NA(){
    // idx is a vector of the index of probes we will consider
	// restrict idx by elminating any rows where class (a,b)
    // has fewer than (n_req_a, n_req_b) samples present
	// We require at least two present values per probeset no matter what is specified.

    std::vector<int> keep;
    int n_in_a, n_in_b, a, b, val;
    if( this->data->raw_data->data->has_missing_values == false )
        return;
	//if( this->percent_required==0 )
	//	return;
	int n_req_a = (int)((double)this->data->a_idx->size() * this->percent_required );
    int n_req_b = (int)((double)this->data->b_idx->size() * this->percent_required );
    if( n_req_a < 2 )
    	n_req_a = 2;
    if( n_req_b < 2 )
    	n_req_b = 2;
    for(int i=0; i<(int)this->idx.size(); i++){
        n_in_a = 0;
        n_in_b = 0;
        val = this->idx.at(i);
        if(this->data->raw_data->data->row_has_missing_value(val)){
            for(a=0;a<(int)this->data->a_idx->size(); a++){
                if(!this->data->raw_data->data->is_missing(val,data->a_idx->at(a))){
                    n_in_a++;
                }
            }
            if(this->data->b_is_target){
                for(b=0;b<(int)this->data->b_idx->size();b++){
                    if(!this->data->raw_data->data->is_missing(val,data->b_idx->at(b))){
                        n_in_b++;
                    }
                }
                if( n_in_a>=n_req_a && n_in_b>=n_req_b )
                    keep.push_back(i);
            }
            else{
                if( n_in_a>=n_req_a ){
                    keep.push_back(val);
                }
            }
        }
        else{
            keep.push_back(val);
        }
    }
    idx.clear();
	for(int i=0; i<(int)keep.size(); i++)
		idx.push_back( keep.at(i) );
}


void Spearman::limit_ids_by_seeds(){
	// user may set a vector probe IDs in this->seeds
	// restrict analysis to these probes.
	if(this->seeds.size()==0)
		return; // sanity check
	std::vector<int> keep;
	HASH_S_I hsh_seeds;
	for(int i=0; i<(int)this->seeds.size();i++){
		hsh_seeds[this->seeds.at(i)] = 1;
	}
	for( int i=0; i<(int)this->idx.size(); i++){
		if( hsh_seeds.find(this->data->identifiers->at(i)) != hsh_seeds.end() ){
			keep.push_back(i);
		}
	}
    idx.clear();
	for(int i=0; i<(int)keep.size(); i++)
		idx.push_back( keep.at(i) );
}


double Spearman::mean(std::vector<double>& M){
	double sum=0;
	int N = (int)M.size();
	if( N==0 )
		return 0;
	for(int i=0; i<N; i++){
		sum += M.at(i);
	}
	return sum / N;
}


void Spearman::find_ranks( Matrix<float>* raw_data, std::vector<int>& sample_idx, Matrix<double>* ranks){
	// Each sample has n_genes measurements that are a floating-point value.
    // We will convert these floating point measurements into ranks for use by the Spearman
    // rank correlation.  It is efficient to pre-calculate all of the ranks in the raw_data because
    // when we compare one gene to thousands of others we don't want to needlessly calculate this
    // ranking over and over.
	//
	// Sorter is a three-element int array: val, true index, rank.
    // The first two are given; Rank is calculated.
	// For each identifier:
	//		generate ranks (N samples)
	//		store ranks
    // INPUT:
    //   Matrix raw_data  {n_genes, n_samples} with doubles, sample_idx vector indicating which samples to use.
    // OUTPUT:
    //   The matrix ranks {n_genes, n_samples} contains the rank of each gene within a sample
	int n_genes = (int)raw_data->rows();
	int N = (int)sample_idx.size();
	int i,j;
	ranks->resize(n_genes, N);
	std::vector< Sorter* > sorters_a;
	for(i=0; i<N; i++)
		sorters_a.push_back( new Sorter(0, i) );
	for(int idx_gene=0; idx_gene<n_genes; idx_gene++){
		// write values and sort in value-order.
		for(j=0; j<N; j++){
			sorters_a.at(j)->v[0] = raw_data->arr[idx_gene][ sample_idx.at(j) ];
			sorters_a.at(j)->v[1] = j; // true index
		}
		std::sort(sorters_a.begin(), sorters_a.end(), sorter_by_value );
		// break ties with mean of tied ranks, 0 0 0 1 ->
		double mean_val;
		int rank_start, rank_end;
		std::vector<double> ranks_for_mean;
		for(i=0; i<N; i++){
			if(i==N-1 || sorters_a.at(i)->v[0] != sorters_a.at(i+1)->v[0] ){
				sorters_a.at(i)->v[2] = i+1; // last item (i-1 was not a tie) or no tie between i and i+1
			}
			else{
				ranks_for_mean.clear();
				ranks_for_mean.push_back(i+1);
				ranks_for_mean.push_back(i+2);
				rank_start = i;
				rank_end = i+2;
				while(rank_end < N && sorters_a.at(rank_start)->v[0] == sorters_a.at(rank_end)->v[0] ){
					ranks_for_mean.push_back(rank_end+1);
					rank_end++;
				}
				mean_val = mean(ranks_for_mean);
				for( i=i; i<rank_end;i++)
					sorters_a.at(i)->v[2]=mean_val;
				i--; // because i is incremented in for loop
			}
		}
		// write ranks
		for(j=0; j<N; j++){
			ranks->arr[idx_gene][ int(sorters_a.at(j)->v[1]) ] = sorters_a.at(j)->v[2];
		}
	}
	for(int i=0; i<int(sorters_a.size());i++)
		delete sorters_a.at(i);

	if( raw_data->has_missing_values ){
		for(int r=0; r<ranks->rows(); r++){
			if( raw_data->row_has_missing_value(r) ){
				for(int c=0; c<ranks->cols(); c++){
					if( raw_data->is_missing(r, sample_idx.at(c) ) ){
						ranks->set_missing_value(r, c );
					}
				}
			}
		}
	}
}

/*
 * This python code correctly (to my understanding) re-ranks two rows
 * and is much more efficient than re-calculating ranks from scratch,
 * under the assumption that missing values are a small percentage of
 * the total number of samples.
def rerank(A_orig, B_orig, is_miss_A, is_miss_B):
    N = len(A_orig)
    min_miss_A, min_miss_B = 99999, 99999
    n_miss_A=0
    A = [x for x in A_orig]
    B = [x for x in B_orig]
    m_A, m_B, is_missing_either = [], [], []
    for i in range(0,len(is_miss_A) ):
        if is_miss_A[i] or is_miss_B[i]:
            m_A.append( A[i] )
            m_B.append( B[i] )
            if A[i] < min_miss_A:
                min_miss_A = A[i]
            if B[i] < min_miss_B:
                min_miss_B = B[i]
            is_missing_either.append(True)
        else:
            is_missing_either.append(False)
    for idx in range(0,N):
        if A[idx] >= min_miss_A:
            me = A[idx]
            n_miss_lt_me=0
            for m in m_A:
                if m < me:
                    n_miss_lt_me += 1
            A[idx] -= n_miss_lt_me
        if B[idx] >= min_miss_B:
            me = B[idx]
            n_miss_lt_me=0
            for m in m_B:
                if m < me:
                    n_miss_lt_me += 1
            B[idx] -= n_miss_lt_me
    for idx, missing in enumerate(is_missing_either):
        if missing:
            A[idx] = 0
            B[idx] = 0
    return A, B
 */

double Spearman::find_spearman_from_ranks( Matrix<double>* ranks, int row1, int row2, std::vector<int>& cols1, std::vector<int>& cols2 ){
	// Use ranks (n_genes x n_samples_in_class) to calculate Spearman rank correlation (rho)
	// of any two genes. Normally ranks1 and ranks2 are the SAME pointer. If we are performing a permutation
	// analysis ranks1 is observed (the true values loaded from files) and ranks2 is a shuffled permutation.
	// The key statistic to calculate is 1.0 - (6 * SUM(D)^2 / ( (N*N*N) - N ) );
	//
	// If there are missing data, we must generate updated ranks (A and B) every time we do a comparison.
	// When doing a non-trivial number of comparisons it is prohibitively expensive to do this
	// by re-calculating from raw data.
	//
	// This version uses a single matrix of ranks and allows the user to specify both which rows and which columns
	// to use in each comparison.

	int c1, c2;
	if( cols1.size() != cols2.size() )
		throw std::string("Processing error in find_spearman_from_ranks: column lists do not have same number of columns");

	if( ranks->row_has_missing_value(row1) || ranks->row_has_missing_value(row2) ){
		// mark rows that are present in both
		int N = int(cols1.size());
		std::vector<int> m_A, m_B; // missing ranks in A, B
		boost::unordered_map<int, int> is_missing_either;
		for(int i=0; i<int(cols1.size()); i++){
			c1 = cols1.at(i);
			c2 = cols2.at(i);
			if( ranks->is_missing(row1, c1) || ranks->is_missing(row2, c2) ){
				m_A.push_back( ranks->arr[row1][c1] );
				m_B.push_back( ranks->arr[row2][c2] );
				is_missing_either[i] = 1;
			}
		}
		std::sort( m_A.begin(), m_A.end() );
		std::sort( m_B.begin(), m_B.end() );
		double sigma_d2=0, d=0;
		double a,b,v1, v2;
		int offset=0, N_present=0;
		for(int i=0; i<N; i++){
			c1 = cols1.at(i);
			c2 = cols2.at(i);
			if( is_missing_either.find(i) == is_missing_either.end() ){
				N_present += 1;
				a = ranks->arr[row1][c1];
				offset=0;
				for(int j=0; j<int(m_A.size()); j++){
					if( a < m_A.at(j) )
		                break;
				    offset += 1;
				}
		        v1 = a - offset;
		        b = ranks->arr[row2][c2];
		        offset=0;
		        for(int j=0; j<int(m_B.size()); j++){
		        	if( b < m_B.at(j) )
		        		break;
		        	offset += 1;
		        }
		        v2 = b - offset;
		        d = v1-v2;
				sigma_d2 += (d * d);
			}
		}
		return 1.0 - (6 * sigma_d2 / ( (N_present*N_present*N_present) - N_present ) );
	}
	else{
		// no missing values.
		double sigma_d2=0, d=0;
		int N = int(cols1.size());
		for( int i=0; i<N; i++ ){
			c1 = cols1.at(i);
			c2 = cols2.at(i);
			d = ranks->arr[row1][c1] - ranks->arr[row2][c2];
			sigma_d2 += (d * d);
		}
		return 1.0 - (6 * sigma_d2 / ( (N*N*N) - N ) );
	}
}

double Spearman::find_spearman_from_ranks( Matrix<double>* ranks1, int row1, Matrix<double>* ranks2, int row2){
	// Use ranks (n_genes x n_samples_in_class) to calculate Spearman rank correlation (rho)
	// of any two genes. Normally ranks1 and ranks2 are the SAME pointer. If we are performing a permutation
	// analysis ranks1 is observed (the true values loaded from files) and ranks2 is a shuffled permutation.
	// The key statistic to calculate is 1.0 - (6 * SUM(D)^2 / ( (N*N*N) - N ) );
	//
	// If there are missing data, we must generate updated ranks (A and B) every time we do a comparison.
	// When doing a non-trivial number of comparisons it is prohibitively expensive to do this
	// by re-calculating from raw data.
    //
    // This version uses two separate matrices of ranks, which must have the same dimensions
	//
	if( ranks1->rows() != ranks2->rows() )
		throw std::string("Processing error in find_spearman: ranks1 and rank2 do not have same number of rows");
	if( ranks1->cols() != ranks2->cols() )
		throw std::string("Processing error in find_spearman: ranks1 and rank2 do not have same number of columns");
	if( ranks1->row_has_missing_value(row1) || ranks2->row_has_missing_value(row2) ){
		// mark rows that are present in both
		int N = ranks1->cols();
		std::vector<int> m_A, m_B; // missing ranks in A, B
		boost::unordered_map<int, int> is_missing_either;
		for(int c=0; c<ranks1->cols(); c++){
			if( ranks1->is_missing(row1, c) || ranks2->is_missing(row2, c) ){
				m_A.push_back( ranks1->arr[row1][c] );
				m_B.push_back( ranks2->arr[row2][c] );
				is_missing_either[c] = 1;
			}
		}
		std::sort(m_A.begin(), m_A.end() );
		std::sort(m_B.begin(), m_B.end() );
		double sigma_d2=0, d=0;
		double a,b,v1, v2;
		int offset=0, N_present=0;
		for(int i=0; i<N; i++){
			if( is_missing_either.find(i) == is_missing_either.end() ){
				N_present += 1;
				a = ranks1->arr[row1][i];
				offset=0;
				for(int j=0; j<int(m_A.size()); j++){
					if( a < m_A.at(j) )
		                break;
				    offset += 1;
				}
		        v1 = a - offset;
		        b = ranks2->arr[row2][i];
		        offset=0;
		        for(int j=0; j<int(m_B.size()); j++){
		        	if( b < m_B.at(j) )
		        		break;
		        	offset += 1;
		        }
		        v2 = b - offset;
		        d = v1-v2;
				sigma_d2 += (d * d);
			}
		}
		return 1.0 - (6 * sigma_d2 / ( (N_present*N_present*N_present) - N_present ) );
	}
	else{
		// no missing values.
		double sigma_d2=0, d=0;
		int N = ranks1->cols();
		for(int i=0; i<N; i++){
			d = ranks1->arr[row1][i] - ranks2->arr[row2][i];
			sigma_d2 += (d * d);
		}
		return 1.0 - (6 * sigma_d2 / ( (N*N*N) - N ) );
	}
}


void Spearman::calculate_ranks( Matrix<float>* raw_data, int row, std::vector<int>* idx, std::vector<int>& ranks){
	// given a row of raw_data, calculate the ranks of columns indexed by idx
	// store results in ranks
	// !!! assumes that all values are present (i.e. not marked missing)
	std::vector< Sorter* > sorters;
	int N = (int)idx->size();
	int i=0;
	for(i=0; i<N; i++)
		sorters.push_back( new Sorter(0, i) );

	for(i=0; i<N; i++){
		sorters.at(i)->v[0] = raw_data->arr[row][ idx->at(i) ];
		sorters.at(i)->v[1] = i; // true index
	}
	std::sort(sorters.begin(), sorters.end(), sorter_by_value );
	// break ties with mean of tied ranks, 0 0 0 1 ->
	double mean_val;
	int rank_start, rank_end;
	std::vector<double> ranks_for_mean;
	for(i=0; i<N; i++){
		if(i==N-1 || sorters.at(i)->v[0] != sorters.at(i+1)->v[0] ){
			sorters.at(i)->v[2] = i+1; // last item (i-1 was not a tie) or no tie between i and i+1
		}
		else{
			ranks_for_mean.clear();
			ranks_for_mean.push_back(i+1);
			ranks_for_mean.push_back(i+2);
			rank_start = i;
			rank_end = i+2;
			while(rank_end < N && sorters.at(rank_start)->v[0] == sorters.at(rank_end)->v[0] ){
				ranks_for_mean.push_back(rank_end+1);
				rank_end++;
			}
			mean_val = mean(ranks_for_mean);
			for( i=i; i<rank_end;i++)
				sorters.at(i)->v[2]=mean_val;
			i--; // because i is incremented in for loop
		}
	}
	// sort back to original order
	ranks.clear();
	ranks.reserve(N);
	for(i=0; i<N; i++)
		ranks.push_back(0);
	for(i=0; i<N; i++){
		ranks.at( sorters.at(i)->v[1] ) = sorters.at(i)->v[2] ;
	}
	for(i=0; i<N; i++)
		delete sorters.at(i);
}





bool Spearman::get_mean_difference(double& meanA, double& meanB){
	// calculates the mean correlation coefficients in class A compared to class B
	double sumA=0, sumB = 0;
	if( this->spears.size() == 0 ){
		meanA=0;
		meanB=0;
		return false;
	}
	for(int i=0; i<int(this->spears.size()); i++){
		sumA += this->spears.at(i)->rho_a();
		sumB += this->spears.at(i)->rho_b();
	}
	meanA = sumA / double(this->spears.size());
	meanB = sumB / double(this->spears.size());
	return true;
}


void Spearman::permute_group_labels_persistent( std::vector<int>* a_idx, std::vector<int>* b_idx, std::vector<int>& a_idx_perm, std::vector<int>& b_idx_perm, boost::mt19937& rng ){
	// Given true class labels passed in: a_idx and b_idx,
	// Mixes sample labels across two groups with results stored in a_idx_perm and b_idx_perm.
	std::vector<int> perm;
	for(int i=0; i<int( a_idx->size() ); i++)
		perm.push_back( a_idx->at(i) );
	if( b_idx != NULL ){
		for(int i=0; i<int( b_idx->size() ); i++)
			perm.push_back( b_idx->at(i) );
	}
    //http://stackoverflow.com/questions/4778797/setting-seed-boostrandom
    
    boost::random_number_generator<boost::mt19937> generator( rng );
    
	std::random_shuffle(perm.begin(), perm.end(), generator );
	a_idx_perm.clear();
	b_idx_perm.clear();
	for(int i=0; i<int( a_idx->size() ); i++)
		a_idx_perm.push_back(perm.at(i));
	if(b_idx != NULL){
		for(int i=int( a_idx->size() ); i<int( perm.size() ); i++ )
			b_idx_perm.push_back(perm.at(i));
	}
}



void Spearman::permute_group_labels( std::vector<int>* a_idx, std::vector<int>* b_idx, std::vector<int>& a_idx_perm, std::vector<int>& b_idx_perm ){
	// Given true class labels passed in: a_idx and b_idx,
	// Mixes sample labels across two groups with results stored in a_idx_perm and b_idx_perm.
	std::vector<int> perm;
	for(int i=0; i<int( a_idx->size() ); i++)
		perm.push_back( a_idx->at(i) );
	if( b_idx != NULL ){
		for(int i=0; i<int( b_idx->size() ); i++)
			perm.push_back( b_idx->at(i) );
	}
    //http://stackoverflow.com/questions/4778797/setting-seed-boostrandom
    
    boost::mt19937 rng;
    rng.seed(std::time(0));
    boost::random_number_generator<boost::mt19937> generator( rng );

	std::random_shuffle(perm.begin(), perm.end(), generator );
	a_idx_perm.clear();
	b_idx_perm.clear();
	for(int i=0; i<int( a_idx->size() ); i++)
		a_idx_perm.push_back(perm.at(i));
	if(b_idx != NULL){
		for(int i=int( a_idx->size() ); i<int( perm.size() ); i++ )
			b_idx_perm.push_back(perm.at(i));
	}
}


void Spearman::load_data_for_differential_correlation(){
	// called by the various DC functions
	if( this->fn_expr.length()==0){
		throw std::string( "Must set fn_expr to a valid file path");
	}
	if( this->fn_sa.length()==0){
		throw std::string( "Must set fn_sa to a valid file path");
	}
	if( this->fn_ga.length()==0){
		throw std::string( "Must set fn_ga to a valid file path");
	}
	if( this->class_a.length()==0 || this->class_b.length()==0){
		throw std::string( "Must set both class A and class B" );
	}
	if( this->data != NULL )
		delete this->data;
	this->data = new Dataset();

	if(this->verbose){
		std::cout << "MESSAGE: Loading data set...\n";
		std::cout.flush();
	}
	if( this->sa != NULL )
		delete this->sa;
	if( this->ga != NULL )
		delete this->ga;

	this->sa = new Attributes("NA");
	this->ga = new Attributes("NA");
	try{
		this->sa->load(this->fn_sa);
	}
	catch(std::string e){
		std::cout << "ERROR: " << e << "\n";
		return;
	}
	try{
		this->ga->load(this->fn_ga);
		if( this->gene_name_column.size() > 0 )
			this->ga->set_gene_name_column(this->gene_name_column);
	}
	catch(std::string e){
		std::cout << "ERROR: " << e << "\n";
		return;
	}
	if(this->verbose){
		std::cout << "MESSAGE: Class A defined as: " << this->class_a << "\n";
		   std::cout << "MESSAGE: Class B defined as: " << this->class_b <<"\n";
	   }
	try{
		data->load(sa, ga, this->fn_expr, this->class_a, this->class_b);
	}
	catch(std::string e){
		std::cout << "ERROR: " << e << "\n";
		return;
	}
	if(this->verbose){
		std::cout << "MESSAGE: Loaded expression file "<< this->fn_expr << "\n";
		std::cout << "MESSAGE: Loaded sample file "<< this->fn_sa << "\n";
		std::cout << "MESSAGE: Loaded gene file "<< this->fn_ga << "\n";
		std::cout << "MESSAGE: Loaded data set.  Class A has " << data->a_idx->size() << " members, Class B has " << data->b_idx->size() << " members\n";
		std::cout << "MESSAGE: Calculating " << n_perms << " permutations\n";
		std::cout.flush();
	}
}


double Spearman::find_spearman( Matrix<float>* raw_data1, int row1, std::vector<int>* sample_idx1, Matrix<float>* raw_data2, int row2, std::vector<int>* sample_idx2){
	// Calculate ranks for row1, row2, and then calculate Spearman rank correlation (rho)
	// of any two genes
	// Normally ranks1 and ranks2 are the SAME pointer.
    // If we are shuffling, one is observed (true) and one is shuffled
    // ranks1 and ranks2 must have the same dimensions!
	if( raw_data1->rows() != raw_data2->rows() )
		throw std::string("Processing error in find_spearman: raw_data1 and raw_data2 do not have same number of rows");
	if( raw_data1->cols() != raw_data2->cols() )
		throw std::string("Processing error in find_spearman: raw_data1 and raw_data2 do not have same number of columns");
	if( sample_idx1->size() != sample_idx2->size() )
		throw std::string("Preprocessing error in find_spearman: sample_idx1 and sample_idx2 do not have same number of elements");
	std::vector<int> ranks1, ranks2;
	int N = (int)sample_idx1->size();

	if( raw_data1->row_has_missing_value(row1) || raw_data2->row_has_missing_value(row2) ){
		std::vector<int>* idx1 = new vector<int>();
		std::vector<int>* idx2 = new vector<int>();
		int i1,i2;
		// reset sample_idx to those cases where both row1 and row2 are present
		for(int i=0; i<N; i++){
			i1 = sample_idx1->at(i);
			i2 = sample_idx2->at(i);
			if( !raw_data1->is_missing( row1, i1 ) && !raw_data2->is_missing( row2, i2 ) ){
				idx1->push_back(i1);
				idx2->push_back(i2);
			}
		}
		N = (int)idx1->size(); // reset N to the number present in both
		calculate_ranks( raw_data1, row1, idx1, ranks1 );
		calculate_ranks( raw_data2, row2, idx2, ranks2 );
		delete idx1;
		delete idx2;
	}
	else{
		calculate_ranks( raw_data1, row1, sample_idx1, ranks1 );
		calculate_ranks( raw_data2, row2, sample_idx2, ranks2 );
	}
	double sigma_d2=0, d=0;
	if(N>0){
		double N_dbl = (double)N; // avoid overflow problems incurred with very large N (>15000)
        for(int i=0; i<N; i++){
			d = ranks1.at(i) - ranks2.at(i);
			sigma_d2 += (d * d);
		}
        double Ns = (N_dbl*N_dbl*N_dbl)-N_dbl;
		return 1.0 - (6.0 * sigma_d2 / Ns );
	}
	else
		return 0; // XXX really should return a flag to indicate NA
}


double Spearman::find_correlation_r(Matrix<float>* raw_data1, int row1, std::vector<int>* idx_1, Matrix<float>* raw_data2, int row2, std::vector<int>* idx_2){
	// Use raw data (n_probes x n_samples_in_class) to calculate Pearson correlation (r)
	// Normally raw_data1 and raw_data2 are the SAME pointer.
    // If we are shuffling, one is observed (true) and one is shuffled
    // ranks1 and ranks2 must have the same dimensions!
	if( raw_data1->rows() != raw_data2->rows() )
		throw std::string("Processing error in find_correlation_r: raw_data1 and raw_data2 do not have same number of rows");
	if( raw_data1->cols() != raw_data2->cols() )
		throw std::string("Processing error in find_correlation_r: raw_data1 and raw_data2 do not have same number of columns");
	if( idx_1->size() != idx_2->size() )
		throw std::string("Processing error in find_correlation_r: idx_1 not the same size as idx_2");
	// need sum of squares for X, Y, X*Y
    float SSx = 0, SSy = 0, SSxy = 0;
	int N = (int)idx_1->size();
	float x,y, sum_x=0, sum_y=0;
	int both_present=0;
	if( raw_data1->row_has_missing_value(row1) || raw_data2->row_has_missing_value(row2) ){
		for(int i=0; i<N; i++){
			if( !raw_data1->is_missing(row1,idx_1->at(i)) && !raw_data2->is_missing(row2,idx_2->at(i) ) ){
				x = raw_data1->arr[row1][ idx_1->at(i) ];
				y = raw_data2->arr[row2][ idx_2->at(i) ];
				sum_x += x;  sum_y += y;  SSx += x*x;  SSy += y*y;  SSxy += x*y;
				both_present += 1;
			}
		}
	}
	else{
		both_present=N;
		for(int i=0; i<N; i++){
			x = raw_data1->arr[row1][ idx_1->at(i) ];
			y = raw_data2->arr[row2][ idx_2->at(i) ];
			sum_x += x;  sum_y += y;  SSx += x*x;  SSy += y*y;  SSxy += x*y;
		}
	}
	float obs_cov = SSxy - (sum_x * sum_y / (float)both_present);
	SSx -= (sum_x*sum_x / (float)both_present );
	SSy -= (sum_y*sum_y / (float)both_present );
	float max_possible_pos_cov = SSx * SSy;
	if( max_possible_pos_cov == 0 )
		return 0;
	else
		return obs_cov / sqrt(max_possible_pos_cov);
}


double Spearman::fisher_zscore(double rhoA, double rhoB, int nA, int nB){
	if( rhoA == 1 )
		rhoA = 0.99999;
	if( rhoB == 1 )
		rhoB = 0.99999;
	if( rhoA == -1 )
		rhoA = -0.99999;
	if( rhoB == -1 )
		rhoB = -0.99999;
	double transA = 0.5 * log( (1+rhoA)/(1-rhoA) );
	double transB = 0.5 * log( (1+rhoB)/(1-rhoB) );
	return (transA - transB) / ( sqrt( (1/double(nA-3) ) + (1/double(nB-3)) ) );
}


//***************************************************
// Thread worker functions
//***************************************************


int Spearman::request_permutation_number(){
    // If we're done, return -1
    // We're returning the index into idx, not the actual value of idx.
    // This is so individual threads can report where we are in the vector;
    // otherwise we wouldn't be able to see how far we've gotten
    
    // MUST be locked to continue
    int next_permutation=-1;
    boost::mutex::scoped_lock lock( this->thread_iter_mutex);
    if( this->current_permutation < this->n_perms ){
        next_permutation=this->current_permutation;
        this->current_permutation++;
    }
    return next_permutation;
}

void Spearman::process_rewiring_in_thread(int thread_id ){
    
    // permutation indexes are pre-calculated, use the row that corresponds to the current permutation.
    {
        boost::mutex::scoped_lock lock( this->io_mutex );
        std::cout << "MESSAGE: Initiated permutation thread " << thread_id+1 << ".\n";
    }
    // seed_idx is either set by user, or all valid probesets
    std::vector<int>* seed_idx;
    if( this->seeds.size()>0 ){
        seed_idx = new std::vector<int>();
        for(int i=0; i<(int)this->seeds.size(); i++){
            seed_idx->push_back( this->data->raw_data->identifier2idx[seeds.at(i)] );
        }
    }
	else{
		seed_idx = &(this->idx);
	}
    int N = (int)seed_idx->size();
    
    std::vector<int> perm_a_idx, perm_b_idx;
    Matrix<double>* ranks_a_perm = new Matrix<double>();
	Matrix<double>* ranks_b_perm = new Matrix<double>();
    double rho_a, rho_b;
    bool is_spearman = true;
	if(this->corr_type.compare("spearman") != 0){
		is_spearman = false;
    }
    
    int row1, row2;
    double z_sum_perm, z_val, total_z_sum_perm; // z_sum_perm is one probe, total_z_sum_perm is for all probes
    
    int permutation_number = request_permutation_number();
    while( permutation_number != -1 ){
        // Fill permutation_a_idx, permutation_b_idx with correct scrambled values
        perm_a_idx.clear();
        perm_b_idx.clear();
        for(int i=0; i<this->permutations_idx_a->cols(); i++){
            perm_a_idx.push_back( this->permutations_idx_a->arr[permutation_number][i] );
        }
        for(int i=0; i<this->permutations_idx_b->cols(); i++){
            perm_b_idx.push_back( this->permutations_idx_b->arr[permutation_number][i] );
        }
        if( is_spearman ){
            this->find_ranks(this->data->raw_data->data, perm_a_idx, ranks_a_perm); // recalculate ranking. Expensive.
            this->find_ranks(this->data->raw_data->data, perm_b_idx, ranks_b_perm);
        }
        
        total_z_sum_perm=0.0;
        for(int i=0; i<N; i++){
            row1 = seed_idx->at(i);
            z_sum_perm=0.0;
            for(int j=0; j<(int)idx.size(); j++){ // intersect(seed_idx, idx) may be (1) a subset of idx or (2) equal to idx.
                row2 = idx.at(j);
                if(row1==row2)
                    continue;
                if( is_spearman ){
                    rho_a = (float)this->find_spearman_from_ranks( ranks_a_perm, row1, ranks_a_perm, row2);
                    rho_b = (float)this->find_spearman_from_ranks( ranks_b_perm, row1, ranks_b_perm, row2);
                }
                else{
                    rho_a = (float)this->find_correlation_r( this->data->raw_data->data, row1, &perm_a_idx, this->data->raw_data->data, row2, &perm_a_idx); // note permutation idx
                    rho_b = (float)this->find_correlation_r( this->data->raw_data->data, row1, &perm_b_idx, this->data->raw_data->data, row2, &perm_b_idx);
                }
                z_val = this->fisher_zscore(rho_a, rho_b, this->permutations_idx_a->cols(), this->permutations_idx_b->cols() );
                if(z_val<0){
                    z_val = -1.0 * z_val; // had trouble with abs()
                }
                z_sum_perm += z_val; // correlation magnitude increased
            }
            z_sum_perm = z_sum_perm / double( idx.size() ); // normalize by number of probes
            total_z_sum_perm += z_sum_perm;
            if( z_sum_perm > rewiring_results.at(i)->z_sum() ){
                boost::mutex::scoped_lock lock( this->results_mutex );
                rewiring_results.at(i)->increment_p_perm( 1.0 / double(n_perms) );
            }
            if(this->verbose & (i>1) & (i % 5000 == 0) ){
                boost::mutex::scoped_lock lock( this->io_mutex );
                char timebuf[80];
                struct tm* newtime;
                time_t long_time;
                time( &long_time );
                newtime = localtime( &long_time );
                strftime(timebuf, 80, "%H_%M_%S", newtime);
                std::string timeout(timebuf);
                std::cout << "MESSAGE: " << timeout << " thread " << thread_id+1 << " completed row " << i+1 << " of " << N << " in permutation " << permutation_number+1 << "\n";
                std::cout.flush();
            }
        }
        
        if( total_z_sum_perm > this->z_score_sum_obs ){
            boost::mutex::scoped_lock lock( this->results_mutex );
            this->global_p_perm  += 1.0 / double( this->n_perms );
        }
        {
            boost::mutex::scoped_lock lock( this->io_mutex );
            std::cout << "MESSAGE: Thread " << thread_id+1 << " completed permutation " << permutation_number+1 << ".\n";
            std::cout.flush();
        }
        permutation_number = request_permutation_number();
    }
}


void Spearman::calculate_rewiring_coefficient(){
    // positive z_sum indicates that on the whole correlation is being created in B compared to A.
    // negative z_sum indicates that on the whole correlation is being removed in B compared to A.
    // the null hypothesis is that z_sum is 0.
    
    // Initialize. Store values in this->rewiring_results.
	this->success=false;
	for(int i=0; i<(int)rewiring_results.size(); i++){
		delete rewiring_results.at(i);
	}
	rewiring_results.clear();
    load_data_for_run();
	prune_data_by_seeds();
    if( this->verbose ){
        std::cout << "MESSAGE: Data load complete\n";
        std::cout.flush();
    }
    
    int i, j=0, row1, row2;
	double rho_a, rho_b;
	Matrix<double>* ranks_a = new Matrix<double>(); // stored ranks for spearman rank correlation
	Matrix<double>* ranks_b = new Matrix<double>();
    std::vector<int> perm_a_idx, perm_b_idx;
    int n_A = int(this->data->a_idx->size());
	int n_B = int(this->data->b_idx->size());
    double z_sum = 0, z_increase=0, z_decrease=0, z_val=0, p_perm=0;
    int n_z_lt5=0, n_z_5_6=0, n_z_6_7=0, n_z_7_8=0, n_z_8_9=0, n_z_gt9=0;
    
    std::vector<int>* seed_idx;
    
    // seed_idx is either set by user, or all valid probesets
    if( this->seeds.size()>0 ){
        seed_idx = new std::vector<int>();
        for(int i=0; i<(int)this->seeds.size(); i++){
            seed_idx->push_back( this->data->raw_data->identifier2idx[seeds.at(i)] );
        }
    }
	else{
		seed_idx = &(this->idx);
	}
    int N = (int)seed_idx->size();
    bool is_spearman = true;
	if(this->corr_type.compare("spearman") != 0){
		is_spearman = false;
    }
    else{
		this->find_ranks(this->data->raw_data->data, *(this->data->a_idx), ranks_a);
        this->find_ranks(this->data->raw_data->data, *(this->data->b_idx), ranks_b);
    }
    std::vector<thread*> threads;
      
    for(i=0; i<N; i++){
        row1 = seed_idx->at(i);
        z_sum=0;
        z_increase=0;
        z_decrease=0;
        n_z_lt5=0; n_z_5_6=0; n_z_6_7=0; n_z_7_8=0; n_z_8_9=0; n_z_gt9=0;
        for(j=0; j<(int)idx.size(); j++){ // intersect(seed_idx, idx) may be (1) a subset of idx or (2) equal to idx.
            row2 = idx.at(j);
            if(row1==row2)
                continue;
            if( is_spearman ){
                rho_a = (float)this->find_spearman_from_ranks( ranks_a, row1, ranks_a, row2);
                rho_b = (float)this->find_spearman_from_ranks( ranks_b, row1, ranks_b, row2);
            }
            else{
                rho_a = (float)this->find_correlation_r( this->data->raw_data->data, row1, this->data->a_idx, this->data->raw_data->data, row2, this->data->a_idx);
            	rho_b = (float)this->find_correlation_r( this->data->raw_data->data, row1, this->data->b_idx, this->data->raw_data->data, row2, this->data->b_idx);
            }
            z_val = this->fisher_zscore(rho_a, rho_b, n_A, n_B);
            if(z_val<0.0){
                z_val = -1.0 * z_val; // had trouble with abs()
            }
            if( abs(rho_a)<abs(rho_b) ){
                z_increase += z_val;
            }
            else{
                z_decrease += z_val;
            }
            z_sum += z_val; // total change
            if(z_val<=5){      n_z_lt5++; }
            else if(z_val<=6){ n_z_5_6++; }
            else if(z_val<=7){ n_z_6_7++; }
            else if(z_val<=8){ n_z_7_8++; }
            else if(z_val<=9){ n_z_8_9++; }
            else{              n_z_gt9++; }
        }
        z_sum = z_sum / double( idx.size() ); // normalize by number of probes
        z_increase = z_increase / double( idx.size() ); // normalize by number of probes
        z_decrease = z_decrease / double( idx.size() ); // normalize by number of probes
        this->rewiring_results.push_back(new RewiringResult(row1, z_sum, z_increase, z_decrease, p_perm, n_z_lt5, n_z_5_6, n_z_6_7, n_z_7_8, n_z_8_9, n_z_gt9 ) );
        z_score_sum_obs += z_sum;
        if(this->verbose & (i % 100 == 0) ){
			std::cout << "MESSAGE: Completed row " << i+1 << " of " << N << "...\n";
			std::cout.flush();
		}
    }
    
    // to assess statistical significance of a disruption, permute category assignments (a vs b)
    // calculating ranking is expensive, so I only want to do it n_perms times.
    // I have already calculated the observed Z scores and stored them in this->rewiring_results
    // for each permutation,
    //    randomly shuffle the assignment of samples into class_a and class_b
    //    find_ranks in these new groups
    //    calculate Z score for each gene pair; if it exceeds z_observed, increment p value
    
    boost::mt19937 rng;
    rng.seed(std::time(0)); // seed RNG for permutations
    this->permutations_idx_a = new Matrix<int>(n_perms, n_A, false, 0);
    this->permutations_idx_b = new Matrix<int>(n_perms, n_B, false, 0);
    for(int perm=0; perm<n_perms; perm++){
        this->permute_group_labels_persistent( this->data->a_idx, this->data->b_idx, perm_a_idx, perm_b_idx, rng );
        for(int i=0; i<int(perm_a_idx.size()); i++){
            this->permutations_idx_a->arr[perm][i] = perm_a_idx.at(i);
        }
        for(int i=0; i<int(perm_b_idx.size()); i++){
            this->permutations_idx_b->arr[perm][i] = perm_b_idx.at(i);
        }
    }
    
    if( this->n_perms>0 ) {
        this->current_permutation=0;
        for( int thread_id=0; thread_id<this->n_threads; thread_id++){
            threads.push_back(new thread( boost::bind( &Spearman::process_rewiring_in_thread, this, thread_id)) );
        }
        for( int thread_id=0; thread_id<this->n_threads; thread_id++){
            threads.at(thread_id)->join();
        }
        for( int thread_id=0; thread_id<this->n_threads; thread_id++){
            delete threads.at(thread_id);
        }
    }
    delete ranks_a;
    delete ranks_b;
    this->success=true;
}

void Spearman::find_spears(){
	int n_A = int(this->data->a_idx->size());
	int n_B = int(this->data->b_idx->size());
	int i, j=0, row1, row2;
	double rho_a, rho_b, z_score;
	boost::unordered_map<int,int> hsh_neighbors;
	bool is_spearman = true;
	if(this->corr_type.compare("spearman") != 0)
		is_spearman = false;
	Matrix<double>* ranks_a = new Matrix<double>();
	Matrix<double>* ranks_b = new Matrix<double>();

	// If user did not specify seed probes, consider {all probes} vs {all probes}.
	//   In this case, seed_idx is set to idx, the set of all valid probe index values.
	// If user specified seeds, consider {seeds} vs. {all probes}
	//   In this case, seed_idx is set to the index values of the seeds.
	//   If the user has specified set_limit_network_to_seeds(true), idx has already
	//     been restricted to the seed pool.
	//   This will generate a "star shaped" network of correlations.
	//   If the user has specified set_include_seed_neighbor_correlations(true), we will
	//     go on to calculate correlations between seed neighbor probes.
	//
	// MISSING DATA
	// If row1 and row2 are all present, use pre-calculated ranks
	// If there is a single missing value, we must re-calculate ranks.
	// Ideally, this would be "a single missing value given the constraints of limits placed
	// on group 1 and group 2", but that would require we make two clones of the matrix
	// subsetted by column; this would make permutation very expensive. Could implement later if
	// speed is an issue.
	
    
    // This graph is used to avoid adding redundant comparisons (e.g. {12, 32} and {32, 12} )
    // required in the case where we (1) pass seeds and (2) do not restrict the comparison just to
    // those seeds, because the seed list is not going to be in the same order as the full list of genes.
    
    std::vector<int>* seed_idx;
    
    if( this->seeds.size()>0 ){
        seed_idx = new std::vector<int>();
        for(int i=0; i<(int)this->seeds.size(); i++){
            seed_idx->push_back( this->data->raw_data->identifier2idx[seeds.at(i)] );
        }
    }
	else{
		seed_idx = &(this->idx); // set seed_idx to all valid probesets
	}
    int N = (int)seed_idx->size();

    if( is_spearman ){
		std::vector<int> valid_cols;
		this->find_ranks(this->data->raw_data->data, *(this->data->a_idx), ranks_a);
		if(n_B>0){
			find_ranks(this->data->raw_data->data, *(this->data->b_idx), ranks_b);
		}
	}
    // Using seen to keep track of which pairs we've already tested. Changed to matrix to keep memory footprint constant.
    // If you set the size of seen to int(idx.size()) the code blows up if you restrict to seeds; seen is set to be
    // could be 3x3 for 3 seeds, while row1 and row2 are indexed into the gene attributes and full data tables and will be
    // in the thousands
    Matrix<bool> seen( int(this->data->raw_data->identifiers.size()),
                       int(this->data->raw_data->identifiers.size()),
                      false, false);
	for(i=0; i<N; i++){
        row1 = seed_idx->at(i);
        hsh_neighbors[row1] = 1;
        for(j=0; j<(int)idx.size(); j++){ // intersect(seed_idx, idx) may be (1) a subset of idx or (2) equal to idx.
            row2 = idx.at(j);
            if( row1==row2 || seen.arr[row1][row2] )
                continue;
            else{
                seen.arr[row2][row1] = true;
                seen.arr[row1][row2] = true;
            }
            if( is_spearman )
                rho_a = this->find_spearman_from_ranks( ranks_a, row1, ranks_a, row2);
            else
                rho_a = this->find_correlation_r( this->data->raw_data->data, row1, this->data->a_idx, this->data->raw_data->data, row2, this->data->a_idx);
            if( n_B == 0 ){
				z_score = this->fisher_zscore(rho_a, 0, n_A, n_A);
				if( abs(rho_a) >= this->corr_abs_a && abs(z_score) >= this->min_zscore ){
					this->spears.push_back(new Spear(rho_a, 0, row1, row2, z_score ) );
					hsh_neighbors[row2] = 1;
				}
            }
            else{
				if( is_spearman )
					rho_b = this->find_spearman_from_ranks( ranks_b, row1, ranks_b, row2);
				else
					rho_b = this->find_correlation_r( this->data->raw_data->data, row1, this->data->b_idx, this->data->raw_data->data, row2, this->data->b_idx);
				z_score = this->fisher_zscore(rho_a, rho_b, n_A, n_B);
				if( abs(rho_a) >= this->corr_abs_a && abs(rho_b) >= this->corr_abs_b ){
					if( abs( rho_a - rho_b ) >= this->corr_diff && abs(z_score) >= this->min_zscore ){
						this->spears.push_back(new Spear(rho_a, rho_b, row1, row2, z_score ));
						hsh_neighbors[row2]=1;                        
					}
				}
            }
		}
		if(this->verbose && i % 100==0){
			std::cout << "MESSAGE: Completed row " << i+1 << " of " << N << ", " << this->spears.size() << " stored...\n";
			std::cout.flush();
		}
	}

	if( this->seeds.size()>0 && include_seed_neighbor_correlations ){
        std::cout << "MESSAGE: Calculating neighbor correlations\n"; 
		std::vector<int> neighbors;
		for(boost::unordered_map<int, int>::iterator iter = hsh_neighbors.begin(); iter != hsh_neighbors.end(); iter++){
			neighbors.push_back( iter->first );
		}
		for(int i=0; i<(int)neighbors.size(); i++){
			row1 = neighbors.at(i);
			for(int j=i+1; j<(int)neighbors.size(); j++){
				row2 = neighbors.at(j);
				if( row1==row2 || seen.arr[row1][row2] )
                    continue;
                else{
                    seen.arr[row2][row1] = true;
                    seen.arr[row1][row2] = true;
                }
				if( is_spearman )
					rho_a = this->find_spearman_from_ranks( ranks_a, row1, ranks_a, row2);
				else
					rho_a = this->find_correlation_r( this->data->raw_data->data, row1, this->data->a_idx, this->data->raw_data->data, row2, this->data->a_idx);
				if( n_B == 0 ){
					z_score = this->fisher_zscore(rho_a, 0, n_A, n_A);
					if( abs(rho_a) >= this->corr_abs_a && abs(z_score) >= this->min_zscore ){
						spears.push_back(new Spear(rho_a, 0, row1, row2, z_score));
					}
				}
				else{
					if( is_spearman )
						rho_b = this->find_spearman_from_ranks( ranks_b, row1, ranks_b, row2);
					else
						rho_b = this->find_correlation_r( this->data->raw_data->data, row1, this->data->b_idx, this->data->raw_data->data, row2, this->data->b_idx);
					z_score = this->fisher_zscore(rho_a, rho_b, n_A, n_B);
					if( abs(rho_a) >= this->corr_abs_a || abs(rho_b) >= this->corr_abs_b ){
						if( abs( rho_a - rho_b ) >= this->corr_diff && abs(z_score) >= this->min_zscore ){
							this->spears.push_back(new Spear(rho_a, rho_b, row1, row2, z_score));
						}
					}
				}
			}
		}
        std::cout << "MESSAGE: Neighbor correlations complete\n";         
	}
	delete ranks_a;
	delete ranks_b;
	
    if( this->seeds.size()>0 )
        delete seed_idx; // if seeds.size()==0, this is a pointer to this->idx and shouldn't be deleted.
    
}

void Spearman::find_DC_distribution(){
	this->rho_distribution.clear();
	for(int i=0; i<10000; i++){
		this->rho_distribution.push_back(0);
	}
	int i, j=0, row1, row2;
	double rho_a, rho_b;
	bool is_spearman = true;
	if(this->corr_type.compare("spearman") != 0)
		is_spearman = false;
	Matrix<double>* ranks_a = new Matrix<double>();
	Matrix<double>* ranks_b = new Matrix<double>();
	std::vector<int>* seed_idx;
	if( this->seeds.size()>0 ){
		seed_idx = new std::vector<int>();
		for(int i=0; i<(int)this->seeds.size(); i++){
			seed_idx->push_back( this->data->raw_data->identifier2idx[seeds.at(i)] );
		}
	}
	else{
		seed_idx = &(this->idx); // set seed_idx to all valid probesets
	}
	int N = (int)seed_idx->size();
	if( is_spearman ){
		std::vector<int> valid_cols;
		this->find_ranks(this->data->raw_data->data, *(this->data->a_idx), ranks_a);
		this->find_ranks(this->data->raw_data->data, *(this->data->b_idx), ranks_b);
	}
	int dist_idx;
	for(i=0; i<N; i++){
		row1 = seed_idx->at(i);
		for(j=i+1; j<(int)idx.size(); j++){
			row2 = idx.at(j);
			if(row1==row2)
				continue; // required because seeds is a subset of idx, so row1 could be equal to row2
			if( is_spearman ){
				rho_a = this->find_spearman_from_ranks( ranks_a, row1, ranks_a, row2);
				rho_b = this->find_spearman_from_ranks( ranks_b, row1, ranks_b, row2);
			}
			else{
				rho_a = this->find_correlation_r( this->data->raw_data->data, row1, this->data->a_idx, this->data->raw_data->data, row2, this->data->a_idx);
				rho_b = this->find_correlation_r( this->data->raw_data->data, row1, this->data->b_idx, this->data->raw_data->data, row2, this->data->b_idx);
			}
			if( abs(rho_a) >= this->corr_abs_a || abs(rho_b) >= this->corr_abs_b ){
				if( abs( rho_a - rho_b ) >= this->corr_diff ){
					dist_idx = floor( ((rho_a-rho_b)+2) * 2500);
					if(dist_idx==10000)
						dist_idx = 9999;
					this->rho_distribution.at( dist_idx ) = this->rho_distribution.at( dist_idx )+1;
				}
			}
		}
		if(this->verbose && i % 100==0){
			std::cout << "MESSAGE: Completed row " << i+1 << " of " << N << "...\n";
			std::cout.flush();
		}
	}
	delete ranks_a;
	delete ranks_b;
}


void Spearman::find_DC_experimentwise(){
	int n_A = int(this->data->a_idx->size());
	int n_B = int(this->data->b_idx->size());
	int i, j=0, row1, row2;
	double rho_a, rho_b, z_score;
	Matrix<double>* ranks_obs_a = new Matrix<double>();
	Matrix<double>* ranks_obs_b = new Matrix<double>();
	Matrix<double>* ranks_perm_a=NULL, *ranks_perm_b=NULL;
	boost::mt19937 state;
	RNG generator(state);
	bool is_spearman = false;
	std::vector<int> valid_columns_a, valid_columns_b;
	std::vector< std::pair<int, int>* > idx_pairs;

	if(this->corr_type.compare("spearman") == 0){
		is_spearman = true;
		this->find_ranks(this->data->raw_data->data, *(this->data->a_idx), ranks_obs_a);
		ranks_perm_a = new Matrix<double>(ranks_obs_a->rows(), ranks_obs_a->cols(), ranks_obs_a->has_missing_values);
		ranks_obs_a->clone(ranks_perm_a);
		this->find_ranks( this->data->raw_data->data, *(this->data->b_idx), ranks_obs_b );
		ranks_perm_b = new Matrix<double>( ranks_obs_b->rows(), ranks_obs_b->cols(), ranks_obs_b->has_missing_values );
		ranks_obs_b->clone(ranks_perm_b);
	}

	// identify indices of probe pairs where abs(rho_A) >= min_rho_A OR abs(rho_B) >= min_rho_B
	// limit analyis to these probes. If min_rho_A and min_rho_B not set, this defaults to full analysis
	int N = (int)this->idx.size();

	for(i=0;i<N; i++){
		if( this->verbose ){
		    std::cout << "MESSAGE: Identifying significant indices in row " << i+1 << " of " << N << "...\n";
		    std::cout.flush();
		}
		for(j=i+1; j<N; j++){
			row1 = idx.at(i);
			row2 = idx.at(j);
			if(row1==row2)
				continue; // required because seeds is a subset of idx, so row1 could be equal to row2
			if( is_spearman ){
				rho_a = this->find_spearman_from_ranks( ranks_obs_a, row1, ranks_obs_a, row2);
				rho_b = this->find_spearman_from_ranks( ranks_obs_b, row1, ranks_obs_b, row2);
			}
			else{
				rho_a = this->find_correlation_r( this->data->raw_data->data, row1, this->data->a_idx, this->data->raw_data->data, row2, this->data->a_idx);
				rho_b = this->find_correlation_r( this->data->raw_data->data, row1, this->data->b_idx, this->data->raw_data->data, row2, this->data->b_idx);
			}
			z_score = this->fisher_zscore(rho_a, rho_b, n_A, n_B);
			if( (abs(rho_a)>=this->corr_abs_a && abs(rho_b)<this->corr_abs_b ) || (abs(rho_b) >= this->corr_abs_b && abs(rho_a) < this->corr_abs_a  )){
				if( abs(z_score) >= this->min_zscore ){
					this->spears.push_back( new Spear(rho_a, rho_b, row1, row2, z_score, n_perms) );
				}
			}
		}
	}
	N = (int)this->spears.size();
	if( this->verbose ){
		std::cout << "MESSAGE: Performing DC experimentwise with " << N << " pairs of probes\n";
		std::cout.flush();
	}

	std::vector<int> valid_columns_perm, perm_a_idx, perm_b_idx;
	for( int perm=0; perm<n_perms; perm++){
		if( is_spearman ){
			ranks_perm_a->shuffle_by_rows(valid_columns_a);
			ranks_perm_b->shuffle_by_rows(valid_columns_b);
		}
		else{
			this->permute_group_labels( this->data->a_idx, this->data->b_idx, perm_a_idx, perm_b_idx );
		}
		double max_value=0;

		for(i=0; i<N; i++){
			row1 = this->spears.at(i)->row1;
			row2 = this->spears.at(i)->row2;
			if( is_spearman ){
				// ranks_perm_a and ranks_perm_b have been shuffled
				rho_a = this->find_spearman_from_ranks( ranks_perm_a, row1, ranks_perm_a, row2);
				rho_b = this->find_spearman_from_ranks( ranks_perm_b, row1, ranks_perm_b, row2);
			}
			else{
				// only need to use shuffled index vectors (perm_a_idx, perm_b_idx)
				rho_a = this->find_correlation_r( this->data->raw_data->data, row1, this->data->a_idx, this->data->raw_data->data, row2, &perm_a_idx);
				rho_b = this->find_correlation_r( this->data->raw_data->data, row1, this->data->b_idx, this->data->raw_data->data, row2, &perm_b_idx);
			}
			if( abs(rho_a-rho_b) < abs(this->spears.at(i)->delta()) ){
				this->spears.at(i)->perm_pvalue--;
			}
			if(verbose && i% 1000==0){
				std::cout << "MESSAGE: Permutation " << perm+1 << " of " << n_perms << " completed row " << i+1 << " of " << N << " max " << max_value << "\n";
				std::cout.flush();
			}
		}
		if(verbose){
			std::cout << "MESSAGE: Completed permutation " << perm+1 << " of " << n_perms << "...\n";
			std::cout.flush();
		}
	}
	for(i=0; i<N; i++){
		this->spears.at(i)->perm_pvalue = this->spears.at(i)->perm_pvalue / double(n_perms);
	}
	delete ranks_obs_a;
	delete ranks_obs_b;
	delete ranks_perm_a;
	delete ranks_perm_b;
}
/*
void Spearman::find_DC_GWER_by_probe(){
	// calculate genome wide error rate for differential correlation by probe
	int i, j=0, row1, row2;
	double rho_a, rho_b, abs_a, abs_b;
	Matrix<double>* ranks_obs_a = new Matrix<double>();
	Matrix<double>* ranks_obs_b = new Matrix<double>();
	Matrix<double>* ranks_perm_a=NULL, *ranks_perm_b=NULL;
	boost::mt19937 state;
	RNG generator(state);
	bool is_spearman = false;
	std::vector<int> valid_columns_a, valid_columns_b;
	std::vector< std::pair<int, int>* > idx_pairs;

	if(this->corr_type.compare("spearman") == 0){
		is_spearman = true;
		this->find_ranks(this->data->raw_data->data, *(this->data->a_idx), ranks_obs_a);
		ranks_perm_a = new Matrix<double>(ranks_obs_a->rows(), ranks_obs_a->cols(), ranks_obs_a->has_missing_values);
		ranks_obs_a->clone(ranks_perm_a);
		this->find_ranks( this->data->raw_data->data, *(this->data->b_idx), ranks_obs_b );
		ranks_perm_b = new Matrix<double>( ranks_obs_b->rows(), ranks_obs_b->cols(), ranks_obs_b->has_missing_values );
		ranks_obs_b->clone(ranks_perm_b);
	}

	// identify indices of probe pairs where abs(rho_A) >= min_rho_A OR abs(rho_B) >= min_rho_B
	// limit GWER analyis to these probes. If min_rho_A and min_rho_B not set, this defaults to full analysis
	int N = (int)this->idx.size();
	std::vector<int> valid_columns_perm, perm_a_idx, perm_b_idx;

	for(i=0;i<N; i++){
		for(j=i+1; j<N; j++){
			row1 = idx.at(i);
			row2 = idx.at(j);
			if( is_spearman ){
				ranks_perm_a->shuffle_by_rows(valid_columns_a);
				ranks_perm_b->shuffle_by_rows(valid_columns_b);
			}
			else{
				this->permute_group_labels( this->data->a_idx, this->data->b_idx, perm_a_idx, perm_b_idx );
			}
			double max_value=0;

			for(i=0; i<N; i++){
				row1 = idx_pairs.at(i)->first;
				row2 = idx_pairs.at(i)->second;
				if( is_spearman ){
					// ranks_perm_a and ranks_perm_b have been shuffled
					rho_a = this->find_spearman_from_ranks( ranks_perm_a, row1, ranks_perm_a, row2);
					rho_b = this->find_spearman_from_ranks( ranks_perm_b, row1, ranks_perm_b, row2);
				}
				else{
					// only need to use shuffled index vectors (perm_a_idx, perm_b_idx)
					rho_a = this->find_correlation_r( this->data->raw_data->data, row1, this->data->a_idx, this->data->raw_data->data, row2, &perm_a_idx);
					rho_b = this->find_correlation_r( this->data->raw_data->data, row1, this->data->b_idx, this->data->raw_data->data, row2, &perm_b_idx);
				}
				abs_a = abs(rho_a);
				abs_b = abs(rho_b);
				if( abs(rho_a-rho_b) >= max_value ){
					max_value = abs(rho_a-rho_b);
					if( this->verbose )
						std::cout << "New max: " << max_value << "\n";
				}
				if(verbose && i% 10==0){
					std::cout << "MESSAGE: Permutation " << perm+1 << " of " << n_perms << " completed row " << i+1 << " of " << N << " max " << max_value << "\n";
					std::cout.flush();
				}
			}
			this->spears.push_back(new Spear(max_value, 0, 0, 0));
			if(verbose){
				std::cout << "MESSAGE: Completed permutation " << perm+1 << " of " << n_perms << "...\n";
				std::cout.flush();
			}

		}
	}
	for(int i=0; i<int(idx_pairs.size()); i++)
		delete idx_pairs.at(i);
	delete ranks_obs_a;
	delete ranks_obs_b;
	delete ranks_perm_a;
	delete ranks_perm_b;
}
*/

void Spearman::find_DC_GWER(){
	// calculate experiment-wide genome wide error rate for differential correlation
	// Works. Very stringent. Can be constrained by corr_abs_a and corr_abs_b
	int n_A = int(this->data->a_idx->size());
	int n_B = int(this->data->b_idx->size());
	int i, j=0, row1, row2;
	double rho_a, rho_b, abs_a, abs_b;
	Matrix<double>* ranks_obs_a = new Matrix<double>();
	Matrix<double>* ranks_obs_b = new Matrix<double>();
	Matrix<double>* ranks_perm_a=NULL, *ranks_perm_b=NULL;
	boost::mt19937 state;
	RNG generator(state);
	bool is_spearman = false;
	std::vector<int> valid_columns_a, valid_columns_b;

	std::vector< std::pair<int, int>* > idx_pairs;

	if(this->corr_type.compare("spearman") == 0){
		is_spearman = true;
		this->find_ranks(this->data->raw_data->data, *(this->data->a_idx), ranks_obs_a);
		ranks_perm_a = new Matrix<double>(ranks_obs_a->rows(), ranks_obs_a->cols(), ranks_obs_a->has_missing_values);
		ranks_obs_a->clone(ranks_perm_a);
		this->find_ranks( this->data->raw_data->data, *(this->data->b_idx), ranks_obs_b );
		ranks_perm_b = new Matrix<double>( ranks_obs_b->rows(), ranks_obs_b->cols(), ranks_obs_b->has_missing_values );
		ranks_obs_b->clone(ranks_perm_b);
		for(int i=0; i<int(this->data->a_idx->size()); i++){
			valid_columns_a.push_back( i );
		}
		for(int i=0; i<int(this->data->b_idx->size()); i++){
			valid_columns_b.push_back( i );
		}
	}

	// identify indices of probe pairs where abs(rho_A) >= min_rho_A OR abs(rho_B) >= min_rho_B
	// limit GWER analyis to these probes. If min_rho_A and min_rho_B not set, this defaults to full analysis
	int N = (int)this->idx.size();
	std::pair<int, int>* p;
	for(i=0;i<N; i++){
		if( this->verbose ){
		    std::cout << "MESSAGE: Identifying significant indices in row " << i+1 << " of " << N << "...\n";
		    std::cout.flush();
		}
		for(j=i+1; j<N; j++){
			row1 = idx.at(i);
			row2 = idx.at(j);
			if(row1==row2)
				continue; // required because seeds is a subset of idx, so row1 could be equal to row2
			if( is_spearman ){
				rho_a = this->find_spearman_from_ranks( ranks_obs_a, row1, ranks_obs_a, row2);
				rho_b = this->find_spearman_from_ranks( ranks_obs_b, row1, ranks_obs_b, row2);
			}
			else{
				rho_a = this->find_correlation_r( this->data->raw_data->data, row1, this->data->a_idx, this->data->raw_data->data, row2, this->data->a_idx);
				rho_b = this->find_correlation_r( this->data->raw_data->data, row1, this->data->b_idx, this->data->raw_data->data, row2, this->data->b_idx);
			}
			if( (abs(rho_a)>=this->corr_abs_a && abs(rho_b)<this->corr_abs_b ) || (abs(rho_b) >= this->corr_abs_b && abs(rho_a) < this->corr_abs_a  )){
				if( abs( this->fisher_zscore(rho_a, rho_b, n_A, n_B) ) >= this->min_zscore ){
					p = new std::pair<int, int>();
					p->first = row1;
					p->second = row2;
					idx_pairs.push_back(p);
				}
			}
		}
	}
	N = (int)idx_pairs.size();
	if( this->verbose ){
		std::cout << "MESSAGE: Performing DC GWER with " << N << " pairs of probes\n";
		std::cout.flush();
	}

	std::vector<int> valid_columns_perm, perm_a_idx, perm_b_idx;
    for( int perm=0; perm<n_perms; perm++){
   		if( is_spearman ){
   			ranks_perm_a->shuffle_by_rows(valid_columns_a);
   			ranks_perm_b->shuffle_by_rows(valid_columns_b);
   		}
   		else{
   			this->permute_group_labels( this->data->a_idx, this->data->b_idx, perm_a_idx, perm_b_idx );
   		}
        double max_value=0;

		for(i=0; i<N; i++){
			row1 = idx_pairs.at(i)->first;
			row2 = idx_pairs.at(i)->second;
			if( is_spearman ){
				// ranks_perm_a and ranks_perm_b have been shuffled
				rho_a = this->find_spearman_from_ranks( ranks_perm_a, row1, ranks_perm_a, row2);
				rho_b = this->find_spearman_from_ranks( ranks_perm_b, row1, ranks_perm_b, row2);
			}
			else{
				// only need to use shuffled index vectors (perm_a_idx, perm_b_idx)
				rho_a = this->find_correlation_r( this->data->raw_data->data, row1, this->data->a_idx, this->data->raw_data->data, row2, &perm_a_idx);
				rho_b = this->find_correlation_r( this->data->raw_data->data, row1, this->data->b_idx, this->data->raw_data->data, row2, &perm_b_idx);
			}
			if( abs(rho_a-rho_b) >= max_value ){
				max_value = abs(rho_a-rho_b);
				if( this->verbose )
					std::cout << "New max: " << max_value << "\n";
			}
			if(verbose && i% 10==0){
				std::cout << "MESSAGE: Permutation " << perm+1 << " of " << n_perms << " completed row " << i+1 << " of " << N << " max " << max_value << "\n";
				std::cout.flush();
			}
		}
        this->spears.push_back(new Spear(max_value, 0, 0, 0, 0));
        if(verbose){
		    std::cout << "MESSAGE: Completed permutation " << perm+1 << " of " << n_perms << "...\n";
		    std::cout.flush();
		}
    }
    for(int i=0; i<int(idx_pairs.size()); i++)
    	delete idx_pairs.at(i);
	delete ranks_obs_a;
	delete ranks_obs_b;
	delete ranks_perm_a;
	delete ranks_perm_b;
}



void Spearman::find_GWER(){
	// If there data->b_idx has no values, we are finding the GWER for individual probes.
	// If there are values in data->b_idx, we are finding the GWER for differential correlation.

	int N = (int)this->idx.size();
	int i, j=0, row1, row2;
	double rho_a, rho_b, abs_a, abs_b;
	Matrix<double>* ranks_obs_a = new Matrix<double>();
	Matrix<double>* ranks_obs_b = new Matrix<double>();
	Matrix<double>* ranks_perm_a=NULL, *ranks_perm_b=NULL;

	bool has_B = this->data->b_idx->size()>0;
	boost::mt19937 state;
	RNG generator(state);
	bool is_spearman = false;
	std::vector<int> valid_columns_a, valid_columns_b;

	if(this->corr_type.compare("spearman") == 0){
		is_spearman = true;
		this->find_ranks(this->data->raw_data->data, *(this->data->a_idx), ranks_obs_a);
		ranks_perm_a = new Matrix<double>(ranks_obs_a->rows(), ranks_obs_a->cols(), ranks_obs_a->has_missing_values);
		ranks_obs_a->clone(ranks_perm_a);
		for(int i=0; i<int(this->data->a_idx->size()); i++){
			valid_columns_a.push_back( i );
		}
		if( has_B ){
			this->find_ranks( this->data->raw_data->data, *(this->data->b_idx), ranks_obs_b );
			ranks_perm_b = new Matrix<double>( ranks_obs_b->rows(), ranks_obs_b->cols(), ranks_obs_b->has_missing_values );
			ranks_obs_b->clone(ranks_perm_b);
			for(int i=0; i<int(this->data->b_idx->size()); i++){
				valid_columns_b.push_back( i );
			}
		}
	}

	std::vector<int> valid_columns_perm, perm_a_idx, perm_b_idx;
    for( int perm=0; perm<n_perms; perm++){
   		if( is_spearman ){
   			ranks_perm_a->shuffle_by_rows(valid_columns_a);
   			if( has_B ){
   				ranks_perm_b->shuffle_by_rows(valid_columns_b);
   			}
   		}
   		else{
   			this->permute_group_labels( this->data->a_idx, this->data->b_idx, perm_a_idx, perm_b_idx );
   		}
        double max_value=0;

		for(i=0;i<N; i++){
			for(j=i+1; j<N; j++){
				row1 = idx.at(i);
				row2 = idx.at(j);
				if( has_B ){
					if( is_spearman ){
						// ranks_perm_a and ranks_perm_b have been shuffled
						rho_a = this->find_spearman_from_ranks( ranks_perm_a, row1, ranks_perm_a, row2);
						rho_b = this->find_spearman_from_ranks( ranks_perm_b, row1, ranks_perm_b, row2);
					}
					else{
						// only need to use shuffled index vectors (perm_a_idx, perm_b_idx)
						rho_a = this->find_correlation_r( this->data->raw_data->data, row1, this->data->a_idx, this->data->raw_data->data, row2, &perm_a_idx);
						rho_b = this->find_correlation_r( this->data->raw_data->data, row1, this->data->b_idx, this->data->raw_data->data, row2, &perm_b_idx);
					}
					abs_a = abs(rho_a);
					abs_b = abs(rho_b);

					if( abs(rho_a-rho_b) >= max_value ){
						max_value = abs(rho_a-rho_b);
						if( this->verbose )
							std::cout << "New max: " << max_value << "\n";
					}
				}
				else{
					// permutation for correlation, use shuffled data (or ranks derived from shuffled data)
					if( is_spearman )
						rho_a = this->find_spearman_from_ranks( ranks_perm_a, row1, ranks_perm_a, row2);
					else
						rho_a = this->find_correlation_r( this->data->raw_data->data, row1, this->data->a_idx, this->data->raw_data->data, row2, &perm_a_idx);
					abs_a = abs(rho_a);
					if( abs_a >= max_value ){
						max_value = abs_a;
						if( this->verbose )
							std::cout << "New max: " << max_value << "\n";
					}
				}
			}
			if(verbose && i% 10==0){
				std::cout << "MESSAGE: Permutation " << perm+1 << " of " << n_perms << " completed row " << i+1 << " of " << N << " max " << max_value << "\n";
				std::cout.flush();
			}
		}
        this->spears.push_back(new Spear(max_value, 0, 0, 0, 0));
        if(verbose){
		    std::cout << "MESSAGE: Completed permutation " << perm+1 << " of " << n_perms << "...\n";
		    std::cout.flush();
		}
    }
	delete ranks_obs_a;
	delete ranks_obs_b;
	delete ranks_perm_a;
	delete ranks_perm_b;
}

double Spearman::correlation(Matrix<float>* A, Matrix<float>* B, std::string method){
	// expects Matrix A and B to be identical dimensions, one row
	// used for unit tests; not called by find_spears()
	bool is_spearman;
	double result;
	if(method.compare("spearman")==0 )
		is_spearman=true;
	else if( method.compare("pearson")==0 )
		is_spearman=false;
	else
		throw std::string("Method must be one of {spearman,pearson}");
	if( A->rows() != 1 )
		throw std::string("Matrix A should have exactly one row");
	if( B->rows() != 1 )
		throw std::string("Matrix B should have exactly one row");
	if( A->cols() != B->cols() )
		throw std::string("Matrix A and B must have the same number of columns");
	std::vector<int>* sample_idx = new std::vector<int>();
	for(int i=0;i<A->cols(); i++)
		sample_idx->push_back(i);
	if( is_spearman ){
		result = (this->find_spearman( A, 0, sample_idx, B, 0, sample_idx) );
	}
	else{
		result = this->find_correlation_r( A, 0, sample_idx, B, 0, sample_idx);
	}
	delete sample_idx;
	return result;
}


void Spearman::print_spears(){
	std::string header;
	this->generate_result_header(header);
	std::cout << header;
	Spear* spear;
	std::string id1, id2;
	std::cout << "# name_1\tid_1\tname_2\tid_2\trho_a\trho_b\tdiff\tperm_pval\tz_score\n";
	std::string gene_name_col = this->ga->get_gene_name_column();
	for(int i=0; i<(int)this->spears.size(); i++){
		spear = this->spears.at(i);
		id1 = this->data->raw_data->identifiers.at( spear->row_1() );
		id2 = this->data->raw_data->identifiers.at( spear->row_2() );
		std::cout << this->ga->prop_for_identifier( id1, gene_name_col ) << "\t" << id1 << "\t";
		std::cout << this->ga->prop_for_identifier( id2, gene_name_col ) << "\t" << id2 << "\t";
		std::cout << spear->rho_a() << "\t" << spear->rho_b() << "\t" << spear->delta() << "\t" << spear->perm_p_value() << "\t" << spear->z_score() << "\n";
	}
}

void Spearman::write_distribution(){
	if(this->verbose){
		std::cout << "MESSAGE: Writing to " << this->fn_out << "\n";
		std::cout.flush();
	}
	std::ofstream f_out(this->fn_out.c_str());
	if( !f_out.is_open() )
		throw std::string( "Unable to open file for writing:" + this->fn_out );
	std::string header;
	this->generate_result_header(header);
	f_out << header;
	f_out << "rho\tcount\n";
	for(int i=0; i<int(this->rho_distribution.size()); i++){
		f_out << (float(i)/2500.0) - 2.0 << '\t' << this->rho_distribution.at(i) << '\n';
	}
	f_out.close();
}

void Spearman::print_distribution(){
	std::string header;
	this->generate_result_header(header);
	std::cout << header;
	std::cout << "rho\tcount\n";
	for(int i=0; i<int(this->rho_distribution.size()); i++){
		std::cout << (float(i)/2500.0) - 2.0 << '\t' << this->rho_distribution.at(i) << '\n';
	}
}


void Spearman::write_rewiring_coefficient(){
    if( this->verbose ){
        std::cout << "MESSAGE: Writing to " << this->fn_out << "\n";
        std::cout.flush();
    }
    std::ofstream f_out(this->fn_out.c_str());
	if( !f_out.is_open() )
		throw std::string( "Unable to open file for writing:" + this->fn_out );
	std::string header;
	this->generate_result_header(header);
	f_out << header;
    std::string id1;
	RewiringResult* result;
    std::string gene_name_col = this->ga->get_gene_name_column();
    f_out << "# Analysis rewiring_coefficient\n";
    f_out << "#probe\tsymbol\tz_sum\tz_increase\tz_decrease\tp.perm\tglobal.p.perm\tz.lt5\tz.5_6\tz.6_7\tz.7_8\tz.8_9\tz.gt9\n";
    for(int i=0; i<(int)this->rewiring_results.size(); i++){
        result = this->rewiring_results.at(i);
        id1 = this->ga->identifiers.at( result->row_1() );
        f_out << id1 << "\t" << this->ga->prop_for_identifier( id1, gene_name_col ) << "\t";
        f_out << std::fixed << std::setprecision(3) << result->z_sum() << "\t" << result->z_increase() << "\t" << result->z_decrease() << "\t";
        f_out << std::fixed << std::setprecision(5) << result->p_perm() << "\t";
        f_out << std::fixed << std::setprecision(5) << this->global_p_perm << "\t";
        f_out << result->n_z_lt5() << "\t" << result->n_z_5_6() << "\t" << result->n_z_6_7() << "\t" << result->n_z_7_8() << "\t" << result->n_z_8_9() << "\t" << result->n_z_gt9() << "\n";
    }
    f_out.close();
}

void Spearman::print_rewiring_coefficient(){
	std::string header;
	this->generate_result_header(header);
    std::cout << header;
    std::string id1;
	RewiringResult* result;
    std::string gene_name_col = this->ga->get_gene_name_column();
    std::cout << "# Analysis rewiring_coefficient\n";
    std::cout << "#probe\tsymbol\tz_sum\tz_increase\tz_decrease\tp_perm\tglobal.p.perm\tz.lt5\tz.5_6\tz.6_7\tz.7_8\tz.8_9\tz.gt9\n";
    for(int i=0; i<(int)this->rewiring_results.size(); i++){
        result = this->rewiring_results.at(i);
        id1 = this->ga->identifiers.at( result->row_1() );
        std::cout << id1 << "\t" << this->ga->prop_for_identifier( id1, gene_name_col ) << "\t";
        std::cout << std::fixed << std::setprecision(3) << result->z_sum() << "\t" << result->z_increase() << "\t" << result->z_decrease() << "\t";
        std::cout << std::fixed << std::setprecision(5) << result->p_perm() << "\t";
        std::cout << std::fixed << std::setprecision(5) << this->global_p_perm << "\t";
        std::cout << result->n_z_lt5() << "\t" << result->n_z_5_6() << "\t" << result->n_z_6_7() << "\t" << result->n_z_7_8() << "\t" << result->n_z_8_9() << "\t" << result->n_z_gt9() << "\n";
    }
}

void Spearman::write_spears(){

	if(this->verbose){
		std::cout << "MESSAGE: Writing to " << this->fn_out << "\n";
		std::cout.flush();
	}
	std::ofstream f_out(this->fn_out.c_str());
	if( !f_out.is_open() )
		throw std::string( "Unable to open file for writing:" + this->fn_out );
	std::string header;
	this->generate_result_header(header);
	f_out << header;
	Spear* spear;
	std::string id1, id2;
	if( this->n_perms > 0 && this->GWER_method.compare("probe")==0 ){
		f_out << "# id_1\tid_2\tnumber_times_perm_gt_obs\n";
		for(int i=0; i<(int)this->spears.size(); i++){
			spear = this->spears.at(i);
			id1 = this->data->ga->identifiers.at( spear->row_1() );
			id2 = this->data->ga->identifiers.at( spear->row_2() );
			f_out << id1 << "\t" << id2 << "\t" << spear->rho_a() << "\n";
		}
	}
	else{
		f_out << "# name_1\tid_1\tname_2\tid_2\trho_a\trho_b\tdiff\tperm_pval\tz_score\n";
		std::string gene_name_col = this->ga->get_gene_name_column();
		for(int i=0; i<(int)this->spears.size(); i++){
			spear = this->spears.at(i);
			id1 = this->ga->identifiers.at( spear->row_1() );
			id2 = this->ga->identifiers.at( spear->row_2() );
			f_out << this->ga->prop_for_identifier( id1, gene_name_col ) << "\t" << id1 << "\t";
			f_out << this->ga->prop_for_identifier( id2, gene_name_col ) << "\t" << id2 << "\t";
			f_out << spear->rho_a() << "\t" << spear->rho_b() << "\t" << spear->delta() << "\t" << spear->perm_p_value() << "\t" << spear->z_score() << "\n";
		}
	}
	f_out.close();
}

void Spearman::load_data_for_run(){
	// sanity checks for parameters. Load attributes and data.
	// If user requests, reduce idx by variance filter and percent present filter
	if( this->fn_expr.length()==0){
		throw std::string( "Must set fn_expr to a valid file path");
	}
	if( this->fn_sa.length()==0){
		throw std::string( "Must set fn_sa to a valid file path");
	}
	if( this->fn_ga.length()==0){
		throw std::string( "Must set fn_ga to a valid file path");
	}
	if( this->class_a.length()==0 && this->class_b.length()==0){
		if(this->verbose){
			std::cout << "MESSAGE: No sample restriction specified; using all samples.\n";
			this->class_a = "IDENTIFIER!-99999";
			std::cout.flush();
		}
	}
	if( this->limit_network_to_seeds && this->seeds.size()<2 ){
		throw std::string("Limiting network to seeds, but less than two seeds were specified.");
	}
	if( this->limit_network_to_seeds && this->include_seed_neighbor_correlations ){
		throw std::string("Cannot specify that I limit network to seeds AND include seed neighbors.");
	}
	if( this->data != NULL )
		delete this->data;
	this->data = new Dataset();

	if(this->verbose){
		std::cout << "MESSAGE: Loading data set...\n";
		std::cout.flush();
	}
	//if( this->sa != NULL ){
		delete this->sa;
        this->sa = NULL;
    //}
	//if( this->ga != NULL ){
		delete this->ga;
        this->ga = NULL;
    //}

	this->sa = new Attributes("NA");
	this->ga = new Attributes("NA");
	this->sa->load(this->fn_sa);
	this->ga->load(this->fn_ga);
	if( this->gene_name_column.size() > 0 )
		this->ga->set_gene_name_column(this->gene_name_column);

	if(this->verbose){
		std::cout << "MESSAGE: Class A defined as: " << this->class_a << "\n";
		std::cout << "MESSAGE: Class B defined as: " << this->class_b <<"\n";
	}
	data->load(sa, ga, this->fn_expr, this->class_a, this->class_b);
	// If class_B == "", Dataset by default puts all samples not in class_A into class_B.
	// This is correct for data mining, but in this context we want to ignore those samples.
	if( this->class_b.length()==0 ){
		if( this->corr_diff != 0 ){
			throw std::string( "Should not set minimum correlation difference without setting class_b_limit" );
		}
		data->b_is_target = false;
		data->b_idx->clear();
	}
	if(this->verbose){
		std::cout << "MESSAGE: Loaded expression file "<< this->fn_expr << "\n";
		std::cout << "MESSAGE: Loaded sample file "<< this->fn_sa << "\n";
		std::cout << "MESSAGE: Loaded gene file "<< this->fn_ga << "\n";
		std::cout << "MESSAGE: Loaded data set.  Class A has " << data->a_idx->size() << " members, Class B has " << data->b_idx->size() << " members\n";
		std::cout.flush();
	}
	this->idx.clear();
	for(int i=0; i<int(this->data->raw_data->identifiers.size()); i++){
		this->idx.push_back(i);
	}
	int n_full = int(this->idx.size());
	this->limit_ids_by_NA();
	int n_after_NA = int(idx.size());
	this->limit_ids_by_var();
	int n_after_var= int(idx.size());

	if(this->verbose){
		std::cout << "MESSAGE: Full data set has " << n_full << " identifiers.\n";
		std::cout << "MESSAGE: Requiring " << this->percent_required*100 << " percent of samples present for each group.\n";
		if( this->percent_required == 0 ){
			std::cout << "MESSAGE: Minimum of 2 samples present in each group still enforced.\n";
		}
		std::cout << "MESSAGE: Minimum change in correlation between classes: " << this->corr_diff << "\n";
		std::cout << "MESSAGE: Minimum absolute value of correlation class A: " << this->corr_abs_a << "\n";
		std::cout << "MESSAGE: Minimum absolute value of correlation class B: " << this->corr_abs_b << "\n";
		std::cout << "MESSAGE: Minimum variance: " << this->min_var << "\n";
        std::cout << "MESSAGE: Minimum Z score: " << this->min_zscore << "\n";
		std::cout << "MESSAGE: Removed " << n_full - n_after_NA << " identifiers for insufficient number of data points.\n";
		std::cout << "MESSAGE: Removed " << n_after_NA - n_after_var << " additional identifiers for insufficient variance.\n";
		if( this->include_seed_neighbor_correlations )
			std::cout << "MESSAGE: Calculating seed neighbor correlations\n";
		else
			std::cout << "MESSAGE: Not calculating neighbor correlations\n";
        std::cout << "MESSAGE: Calculating with " << this->n_threads << " threads.\n";
        std::cout << "MESSAGE: Performing " << this->n_perms << " permutations to assess statistical significance.\n";
		std::cout.flush();
	}
}


void Spearman::prune_data_by_seeds(){
	// called by run() and run_DC()
	// reduce the set of probe indexes in this->idx according to rules specified by
	// seeds. Throws error if unidentifiable seeds is passed, or if all seeds are legal but
	// are not in idx (because, for example, they're not present in the minimal number of samples).
	// If limit_network_to_seeds is true, will further reduce the composition of idx to those
	// seeds which pass muster
	//
	// Should be called AFTER load_data_for_run()
    
	if( seeds.size() > 0 ){
        if(this->verbose){
            std::cout << "MESSAGE: Called with "  << seeds.size() << " seed probes.\n";
			std::cout.flush();
        }
		// remove any empty probes, and check that probes specified by user exist in raw data file
		std::vector<std::string> seedcopy;
		for(int i=0; i<(int)this->seeds.size();i++)
			seedcopy.push_back(this->seeds.at(i));
		this->seeds.clear();
		int n_seeds_removed = 0;
		boost::unordered_map<int,int> idx_lookup;
		for(int i=0; i<int(idx.size()); i++){
			idx_lookup[idx.at(i)] = 1;
		}
		for(int i=0; i<(int)seedcopy.size();i++){
			if( seedcopy.at(i).size()>0 ){
				if( data->raw_data->identifier2idx.find( seedcopy.at(i) ) == data->raw_data->identifier2idx.end() ){
					std::stringstream err; err << "Seed probe " << seedcopy.at(i) << " not found in raw data";
					throw( err.str() );
				}
				if( idx_lookup.find( data->raw_data->identifier2idx[seedcopy.at(i)] ) == idx_lookup.end() ){
					// seed is a valid probe id, but was excluded from analysis by other filters
					n_seeds_removed += 1;
				}
				else{
					this->seeds.push_back(seedcopy.at(i));
				}
			}
		}
		if( this->verbose && n_seeds_removed > 0 ){
			std::cout << "MESSAGE: Removed "  << n_seeds_removed << " seed probes because of filters.\n";
			std::cout.flush();
		}
		if( this->seeds.size()==0){
			// it is legal but very confusing to specify probes and then remove them because of other filters
			throw(std::string("Removed ALL seed probes because of filters."));
		}
		if( this->limit_network_to_seeds){
			if(this->verbose){
				std::cout << "MESSAGE: Restricting network to "  << this->seeds.size() << " seed probes.\n";
				std::cout.flush();
			}
			this->limit_ids_by_seeds();
			if(this->verbose){
				std::cout << "MESSAGE: After restricting to seeds, " << this->idx.size() << " probes will be used for calculation\n";
				std::cout.flush();
			}
		}
	}
}


Rowpair_iterator::Rowpair_iterator(int N, int initial_idx, int maximum_idx){
	// if initial_idx >0, iterate through valid pairs of i,j to the
	// initial_idx value. Used to skip to tranche X of Y.
	this->i=0;
	this->j=i+1;
	this->N = N;
	this->ctr=0;
	this->maximum_idx = maximum_idx;
	if( initial_idx >0 ){
		for(i=0; i<N; i++){
			for(j=i+1; j<N; j++){
				ctr++;
				if( ctr==initial_idx )
					break;
				if( ctr > maximum_idx ){
					i=N; j=N;
					break;
				}
			}
		}
	}
}

Rowpair_iterator::Rowpair_iterator(std::vector< std::pair<int, int>* > & inbound_pairs){
	this->ctr=0;
	// Instead of iterating through a loop, this subclass is passed a list of integer pairs
	// corresponding to the index of probe pairs in the expression dataset
	for(int i=0; i< int(inbound_pairs.size()); i++){
		std::pair<int, int>* p = new std::pair<int, int>();
		p->first = inbound_pairs.at(i)->first;
		p->second = inbound_pairs.at(i)->second;
		this->pairs.push_back( p );
	}
	this->ctr = 0;
}

bool Rowpair_iterator::next_rowpair(int& row1, int& row2 ){
	if( this->pairs.size()==0 ){
		// This unpacks a simple nested loop:
		// 	for(i=0; i<N; i++){
		//     for(j=i+1; j<N; j++){
		//		  ctr++;

		this->ctr++;
		if( this->ctr > this->maximum_idx )
			return false;
		if( this->j==this->N ){
			this->i++;
			this->j=i+1;
		}
		if( this->i==this->N )
			return false;
		row1 = this->i;
		row2 = this->j;
		return true;
	}
	else{
		if( this->ctr == int(this->pairs.size() ) )
			return false;
		this->i = this->pairs[ctr]->first;
		this->j = this->pairs[ctr]->second;
		ctr++;
		return true;
	}
}

int Rowpair_iterator::get_ctr(){
	return this->ctr;
}

Rowpair_iterator::~Rowpair_iterator(){
	for(int i=0; i< int(this->pairs.size()); i++){
		delete this->pairs.at(i);
	}
}


void Spearman::run_DC(Rowpair_iterator* RPI, std::vector<int>* pvalue_distribution, double early_stop_pvalue){
	// ranks_obs_a and ranks_obs_b have the same number of columns as size(a_idx), size(b_idx)
	// we pre-calculate n_perms shuffled datasets and store those in vectors ranks_perms_a, ranks_perms_b
	//
	// We will generate an element-wise p-value for probe pair. Instead of storing the astronomical
	// number of spears, most of which are useless, we instead store the DISTRIBUTION of p-values.
	// Since p-values are calculated by permutation, each value has an exact value and the number of
	// possible distinct p-values is equal to the number of permutations performed.
	//
	// Only those pairs which have a p-value lower than this->max_eqtl_pval will be stored as spear
	// objects.

	// initialize p-value distribution
	if(pvalue_distribution==NULL)
		throw std::string("Must pass in vector for stored p-values");
	pvalue_distribution->clear();
	for(int i=0; i<=this->n_perms; i++)
		pvalue_distribution->push_back(0);
	for(int i=0; i<(int)spears.size(); i++){
		delete spears.at(i);
	}
	spears.clear();

	if( early_stop_pvalue<0 || early_stop_pvalue>1)
		throw std::string("early stop p value not between 0 and 1");
	else{
		if(this->verbose){
			std::cout << "MESSAGE: Stopping p-value permutations for pairs that exceed P = " << early_stop_pvalue << "\n";
			std::cout.flush();
		}
	}

	int N_A = int( this->data->a_idx->size() );
	int N_B = int( this->data->b_idx->size() );
	int row1=-1, row2=-1, pvalue_int, perm;
	double rho_a, rho_b, diff_obs, pvalue, perm_rho_a, perm_rho_b, diff_perm, z_score;
	Matrix<double>* ranks_obs_a = new Matrix<double>();
	Matrix<double>* ranks_obs_b = new Matrix<double>();
	Matrix<double>* ranks_perm_combined = new Matrix<double>();

	boost::mt19937 state;
	RNG generator(state);
	bool is_spearman = false;
	std::vector< Matrix<double>* > ranks_perms_a;
	std::vector< Matrix<double>* > ranks_perms_b;
	std::vector<int> valid_columns, perm_idx_a, perm_idx_b;

	if(this->corr_type.compare("spearman") == 0){
		is_spearman = true;
		// Find true ranks, stored in ranks_obs_a, ranks_obs_b
		this->find_ranks(this->data->raw_data->data, *(this->data->a_idx), ranks_obs_a);
		this->find_ranks( this->data->raw_data->data, *(this->data->b_idx), ranks_obs_b );
		for(int i=0; i<int(this->data->a_idx->size()); i++){ valid_columns.push_back( this->data->a_idx->at(i) ); }
		for(int i=0; i<int(this->data->b_idx->size()); i++){ valid_columns.push_back( this->data->b_idx->at(i) ); }
		if( this->verbose ){
			std::cout << "MESSAGE: Calculating permutation correlation ranks...\n";
			std::cout.flush();
		}
		for(int i=0; i<this->n_perms; i++){
			// shuffle true indices, distribute proportionately into perm_idx_a, perm_idx_b
			std::random_shuffle(valid_columns.begin(), valid_columns.end(), generator );
			perm_idx_a.clear();
			perm_idx_b.clear();
			for(int idx=0; idx<N_A; idx++){ perm_idx_a.push_back( valid_columns.at(idx) ); }
			for(int idx=0; idx<N_B; idx++){	perm_idx_b.push_back( valid_columns.at(N_A + idx) ); }
			// calculate permutation ranks based on shuffled data locations, stored in ranks_perms_a, ranks_perms_b
			Matrix<double>* M_A = new Matrix<double>();
			Matrix<double>* M_B = new Matrix<double>();
			this->find_ranks(this->data->raw_data->data, perm_idx_a, M_A );
			this->find_ranks( this->data->raw_data->data, perm_idx_b, M_B );
			ranks_perms_a.push_back( M_A );
			ranks_perms_b.push_back( M_B );
		}
	}

	int t_begin=0;

	int max_p_counts = floor(early_stop_pvalue * this->n_perms);
	t_begin = time(NULL);

	while( RPI->next_rowpair( row1, row2 ) ){
		//if( this->verbose && (RPI.get_ctr() % 10000==0 ) ){
		//	double one_operation= double((time(NULL)-t_begin)) / double(RPI.get_ctr()-idx_begin);
		//	double timeleft = one_operation * (idx_end-RPI.get_ctr()) / 3600;
		//	std::cout << "MESSAGE: Pair " << RPI.get_ctr() - idx_begin + 1 << " of " << idx_end-idx_begin << " (i,j)=" << row1 << "," << row2 << ", ETA " << std::setprecision(2) << timeleft << " hours\n";
		//	std::cout.flush();
		//}
		if( is_spearman ){
			rho_a = this->find_spearman_from_ranks( ranks_obs_a, row1, ranks_obs_a, row2 );
			rho_b = this->find_spearman_from_ranks( ranks_obs_b, row1, ranks_obs_b, row2 );
		}
		else{
			rho_a = this->find_correlation_r( this->data->raw_data->data, row1, this->data->a_idx, this->data->raw_data->data, row2, this->data->a_idx);
			rho_b = this->find_correlation_r( this->data->raw_data->data, row1, this->data->b_idx, this->data->raw_data->data, row2, this->data->b_idx);
		}
		diff_obs = rho_a - rho_b;
		if( abs(diff_obs) >= this->corr_diff && ( abs(rho_a) >= this->corr_abs_a || abs(rho_b) >= this->corr_abs_b ) ){
			pvalue_int = this->n_perms; // Start at 100 or 1000 and count down if diff_obs>=diff_perm
			for(perm=0; perm< this->n_perms; perm++){
				if( is_spearman ){
					perm_rho_a = this->find_spearman_from_ranks( ranks_perms_a.at(perm), row1, ranks_perms_a.at(perm), row2 );
					perm_rho_b = this->find_spearman_from_ranks( ranks_perms_b.at(perm), row1, ranks_perms_b.at(perm), row2 );
				}
				else{
					this->permute_group_labels( this->data->a_idx, this->data->b_idx, perm_idx_a, perm_idx_b );
					perm_rho_a = this->find_correlation_r( this->data->raw_data->data, row1, this->data->a_idx, this->data->raw_data->data, row2, &perm_idx_a);
					perm_rho_b = this->find_correlation_r( this->data->raw_data->data, row1, this->data->b_idx, this->data->raw_data->data, row2, &perm_idx_b);
				}
				diff_perm = perm_rho_a - perm_rho_b;
				if( abs(diff_obs) >= abs(diff_perm) ){
					pvalue_int--;
				}
				if( this->n_perms - pvalue_int > max_p_counts ){
					pvalue_int = -1; // mark early stop by setting pvalue_int to -1
					break;
				}
			}
			if( pvalue_int >= 0 ){
				// skip pairs where we hit the early stop (pvalue_int is -1 in those cases)
				pvalue_distribution->at(pvalue_int) += 1;
				pvalue = double(pvalue_int)/this->n_perms;
				if( pvalue <= this->max_eqtl_pval ){
					z_score = this->fisher_zscore(rho_a, rho_b, N_A, N_B);
					this->spears.push_back(new Spear(rho_a, rho_b, row1, row2, z_score, pvalue) );
				}
			}
		}
	}
	delete ranks_obs_a;
	delete ranks_obs_b;
	delete ranks_perm_combined;
	for(int i=0; i<this->n_perms; i++){
		delete ranks_perms_a.at(i);
		delete ranks_perms_b.at(i);
	}
}

void Spearman::run(){
    std::cout.flush();
    
	this->success=false;
	for(int i=0; i<(int)spears.size(); i++){
		delete spears.at(i);
	}
	spears.clear();
	// Following two initialization functions are shared with run_DC()
	load_data_for_run();
	prune_data_by_seeds();
    if( this->verbose ){
        std::cout << "MESSAGE: Data load complete\n";
        std::cout.flush();
    }
    if( n_perms>0 ){
        // perform analysis on shuffled data to identify strongest rho.  Used for permutation analysis
        if( this->verbose ){
        	std::cout << "MESSAGE: Performing " << this->n_perms << " permutations to calculate GWER.\n";
        	std::cout.flush();
        }
        if( this->class_b.size()==0 ){
        	find_GWER();
        }
        else{
        	find_DC_GWER();
        	//find_DC_experimentwise();
        }
    }
    else{
    	if( this->do_distribution )
    		this->find_DC_distribution();
    	else
    		this->find_spears();
	}
    std::sort(spears.begin(), spears.end(), spears_sort);
	if(this->verbose){
		std::cout << "MESSAGE: Found " << spears.size() << " relevant identifier pairs.\n";
		std::cout.flush();
	}
	this->success=true;
}



void Spearman::convert_spears_to_graph(Graph* G, Attributes* ga, boost::unordered_map<int, std::vector<std::string>* >& g2p, double min_abs){
	// Used by write_to_cytoscape (called from prepare_graphs)
	// Resulting graph (returned as G) is indexed by gene names, not by probes.
	// This merges multiple probes with the same gene name.
	bool has_gene_name = true;
	if( ga==NULL )
		throw std::string("Passed null ga into convert_spears_to_graph");
	if( G == NULL )
		throw std::string("Passed null G (graph object) into convert_spears_to_graph");
	if( this->verbose ){
		std::cout << "MESSAGE: Gene Symbols come from " << this->ga->get_gene_name_column() << "\n";
		std::cout.flush();
	}
	std::string gene_name_col = this->ga->get_gene_name_column();
	if( ga->attribute2idx.find(gene_name_col) == ga->attribute2idx.end() )
		has_gene_name=false;
	std::string g1, g2;
	int idx = 0;
	int g1_idx, g2_idx;
	HASH_S_I str2idx, seen;
	std::string probe1, probe2;
	std::vector<std::string>* v;
    
	for(int i=0; i<(int)this->spears.size(); i++){
        if( abs( this->spears.at(i)->rho_a() ) < min_abs )
			continue;
		probe1 = ga->identifiers.at( this->spears.at(i)->row_1() );
		probe2 = ga->identifiers.at( this->spears.at(i)->row_2() );
		if( has_gene_name ){
			g1 = ga->prop_for_identifier( probe1, gene_name_col );
			g2 = ga->prop_for_identifier( probe2, gene_name_col );
		}
		else{
			g1 = probe1;
			g2 = probe2;
		}
		if( g1.compare(g2)==0 )
			continue;
		if( str2idx.find(g1) == str2idx.end() ){
			str2idx[g1] = idx;
			g1_idx = idx;
			idx++;
		}
		else
			g1_idx = str2idx[g1];

		if( str2idx.find(g2) == str2idx.end() ){
			str2idx[g2] = idx;
			g2_idx = idx;
			idx++;
		}
		else
			g2_idx = str2idx[g2];
		// reverse-lookup gene -> vector of probes
		if( seen.find(probe1)==seen.end() ){
			seen[probe1] = 1;
			if( g2p.find( g1_idx ) == g2p.end() ){
				v = new std::vector<std::string>();
				v->push_back( probe1 );
				g2p[g1_idx] = v;
			}
			else{
				g2p[g1_idx]->push_back(probe1);
			}
		}
		if( seen.find(probe2)==seen.end() ){
			seen[probe2] = 1;
			if( g2p.find( g2_idx ) == g2p.end() ){
				v = new std::vector<std::string>();
				v->push_back( probe2 );
				g2p[g2_idx] = v;
			}
			else{
				g2p[g2_idx]->push_back(probe2);
			}
		}
		G->add_edge(g1_idx, g2_idx, this->spears.at(i)->rho_a());
		G->set_node_label(g1_idx, g1);
		G->set_node_label(g2_idx, g2);
	}
}

void Spearman::prepare_graphs( Graph* G, Graph* G_QTL, HASH_I_VECTOR_STR& g2p, double min_abs, Attributes* ga, std::string fn_eQTL, double max_perm_pval, bool require_eQTL){
	if( ga == NULL )
		throw std::string("Passed null ga object into prepare_graphs");

	// Read eQTL file fn_eQTL and create graph G_QTL of eQTL results
	if( fn_eQTL.size()>0){
		CarmenParser cp(fn_eQTL);
		std::vector< QTL* > qtls;
		cp.PrepareToReadValues();
		cp.ExtractEQTL( qtls, this->max_eqtl_pval );
		G_QTL->add_qtls(qtls);
		if(this->verbose){
			std::cout << "MESSAGE: loaded eQTL file, keeping " << qtls.size() << " eQTL at max pval " << this->max_eqtl_pval << "\n";
			std::cout.flush();
		}
		for(int i=0; i<int(qtls.size()); i++)
			delete qtls.at(i);
	}

	if( this->verbose && this->min_clique_size > 2 ){
		std::cout << "MESSAGE: Reducing network to cliques of size " << this->min_clique_size << " or greater.\n";
		std::cout.flush();
	}
    
	convert_spears_to_graph(G, ga, g2p, min_abs);

	if( this->min_clique_size == 3 ){
		Graph g_triangles;
		g_triangles.set_verbose(this->verbose);
		G->reduce_to_triangles(g_triangles);
		G->copy(&g_triangles);
	}
	else if( this->min_clique_size == 4 ){
		Graph g_quad;
		g_quad.set_verbose(this->verbose);
		G->reduce_to_quads(g_quad);
		G->copy(&g_quad);
	}
	std::vector<int> nodes_eqtl, nodes_gene, nodes_keep_gene, nodes_keep_eqtl, neighbors;
	G_QTL->nodes(nodes_eqtl);
	if( this->require_eqtl && !(this->allow_uncorrelated_loci_with_eQTL) ){
		// mark genes with a locus in G_QTL
		// also mark those loci
        if(this->verbose){
            std::cout << "MESSAGE: Requiring eQTL but not correlated partner node for all nodes\n";
            std::cout.flush();
        }
		int idx_to_keep;
		for(int i=0; i<int(nodes_eqtl.size()); i++){
			idx_to_keep=G->get_node_index_by_label( G_QTL->get_node_label_by_index(nodes_eqtl.at(i)) );
			if(idx_to_keep>=0){
				nodes_keep_gene.push_back(idx_to_keep);
				nodes_keep_eqtl.push_back(nodes_eqtl.at(i));
				G_QTL->neighbors(nodes_eqtl.at(i), neighbors);
				for(int j=0; j<int(neighbors.size()); j++){
					if( G_QTL->node_type(neighbors.at(j)) == NODE_LOCUS){
						nodes_keep_eqtl.push_back(neighbors.at(j));
					}
				}
			}
		}
		G->subgraph(nodes_keep_gene);
		G_QTL->subgraph(nodes_keep_eqtl);
	}
	else if( nodes_eqtl.size() > 0 ){
		if( !this->allow_uncorrelated_loci_with_eQTL ){
			// only report genes with a locus that also meet our correlation criteria
			// mark loci associated with a probe present in G
            if(this->verbose){
                std::cout << "MESSAGE: Requiring eQTL and correlation for all nodes\n";
                std::cout.flush();
            }
            G->nodes(nodes_gene);
			int idx_to_keep;
            std::cout << nodes_gene.size() << "\n";
			for(int i=0; i<int(nodes_gene.size()); i++){
				idx_to_keep = G_QTL->get_node_index_by_label( G->get_node_label_by_index( nodes_gene.at(i) ) );
				if(idx_to_keep >= -1){
					nodes_keep_eqtl.push_back(idx_to_keep);
					G_QTL->neighbors(idx_to_keep, neighbors);
					for( int j=0; j<int(neighbors.size()); j++){
						if( G_QTL->node_type(neighbors.at(j))==NODE_LOCUS ){
							//nodes_keep_eqtl.push_back(neighbors.at(i)); // changed this to j 2/29/2012
							nodes_keep_eqtl.push_back(neighbors.at(j));                            
						}
					}
				}
			}
			G_QTL->subgraph(nodes_keep_eqtl);
		}
		if( this->seeds.size()>0 && this->limit_network_to_seeds ){
           HASH_S_I seed_hsh;
			for(int i=0; i<int(seeds.size()); i++)
				seed_hsh[seeds.at(i)]=1;
			for(int i=0; i<int(nodes_eqtl.size()); i++){
				if( G_QTL->node_type(neighbors.at(i)) == NODE_LOCUS)
					nodes_keep_eqtl.push_back(i);
				else if( seed_hsh.find( G_QTL->get_node_label_by_index(i) ) != seed_hsh.end() )
					nodes_keep_eqtl.push_back(i);
			}
			G_QTL->subgraph(nodes_keep_eqtl);
		}
	}
}

void Spearman::write_to_cytoscape(double min_abs, std::string fn_base, Attributes* ga){
	this->write_to_cytoscape(min_abs, fn_base, ga, NULL, NULL, std::string(""), 1, false);
}

void Spearman::write_to_cytoscape(double min_abs, std::string fn_base, Attributes* ga, GOAnnotationParser* go, GeneAnnotationParser* gene_parser){
	this->write_to_cytoscape(min_abs, fn_base, ga, go, gene_parser, std::string(""), 1, false);
}

void Spearman::generate_result_header(std::string& header){
	std::stringstream h;
	char timebuf[80];
	struct tm* newtime;
	time_t long_time;
	time( &long_time );
	newtime = localtime( &long_time );
	strftime(timebuf, 80, "%Y_%m_%d_%H_%M_%S", newtime);
	h << "# Source CARMEN 1.0\n";
	h << "# Date " << timebuf << "\n";
	int a_idx_size=0;
	int b_idx_size = 0;
	if( this->data != NULL){
		a_idx_size = this->data->a_idx->size();
		b_idx_size = this->data->b_idx->size();
	}
	h << "# Class_A " << this->class_a << " " << a_idx_size << " members.\n";
	h << "# Class_B " << this->class_b << " " << b_idx_size << " members.\n";
	h << "# Focus_Probe";
	if( this->seeds.size()==0 )
		h << " \n";
	else{
		for(int i=0;i<(int)this->seeds.size();i++)
			h << " " << this->seeds.at(i);
		h  << "\n";
	}
	h << "# Min_var " << this->min_var << "\n";
	h << "# Delta_rho " << this->corr_diff << "\n";
	h << "# Abs_rho " << this->corr_abs_a << "\n";
	h << "# Abs_rho_Class_B " << this->corr_abs_b << "\n";
	h << "# Data_File " << this->fn_expr << "\n";
	h << "# Sample_File " << this->fn_sa << "\n";
	h << "# Genes_File " << this->fn_ga << "\n";
	if( this->fn_eqtl.size()>0 ){
		h << "# EQTL_File " << this->fn_eqtl << "\n";
		if( this->require_eqtl )
			h << "# Require_EQTL True\n";
		else
			h << "# Require_EQTL False\n";
		h << "# Max_eqtl_pval " << this->max_eqtl_pval << "\n";
		if( this->allow_uncorrelated_loci_with_eQTL )
			h << "# Allow_uncorrelated_loci_with_eqtl True\n";
		else
			h << "# Allow_uncorrelated_loci_with_eqtl False\n";
	}
	h << "# Max_pval_eqtl " << max_eqtl_pval << "\n";
	h << "# Rho_difference " << "\n";
	h << "# Min_clique_size " << this->min_clique_size << "\n";
	h << "# Min_percent_present " << this->percent_required << "\n";
	h << "# Gene_name_column " << this->gene_name_column << "\n";
	if( this->include_seed_neighbor_correlations )
		h << "# Include_seed_neighbor_correlations True\n";
	else
		h << "# Include_seed_neighbor_correlations False\n";
	if( this->limit_network_to_seeds )
		h << "# Limit_network_to_seeds True\n";
	else
		h << "# Limit_network_to_seeds False\n";
	h << "# Min_clique_size " << this->min_clique_size << "\n";
	h << "# N_perms " << this->n_perms << "\n";
	h << "# Percent_required " << this->percent_required << "\n";
	header = h.str();
}



void Spearman::write_GO_terms_for_cytoscape( GeneAnnotationParser* gene_parser, GOAnnotationParser* go, std::string fn_bp, std::string fn_mf, std::string fn_cc, Graph* G ){
    // Only called once, but need to make write_to_cytoscape more reasonable in size
    
    std::ofstream f_bp(fn_bp.c_str());
    std::ofstream f_mf(fn_mf.c_str());
    std::ofstream f_cc(fn_cc.c_str());
    if( !f_bp.is_open() )
        throw std::string( "Unable to open file for writing: " + fn_bp);
    if( !f_mf.is_open() )
        throw std::string( "Unable to open file for writing: " + fn_mf);
    if( !f_cc.is_open() )
        throw std::string( "Unable to open file for writing: " + fn_cc);
    f_bp << "Biological_process (java.lang.String)\n";
    f_mf << "Molecular_function (java.lang.String)\n";
    f_cc << "Cellular_component (java.lang.String)\n";
    std::vector<int> node_idx;
    G->nodes(node_idx);
    GeneAnnotation* gene_annot;
    GOAnnotation* go_annot;
    std::string gene;
    std::vector<std::string> bp_list,mf_list,cc_list;
    for(int i=0; i<int(node_idx.size());i++){
        try{
            gene = G->get_node_label_by_index( node_idx.at(i) );
            bp_list.clear(); mf_list.clear(); cc_list.clear();
            gene_annot = gene_parser->get_annotation( gene );
            for(int g=0; g<int(gene_annot->GO_annotations.size() ); g++){
                go_annot = go->get_annotation( gene_annot->GO_annotations.at(g) );
                if(      go_annot->branch == GOAnnotation::GO_BP ){ bp_list.push_back(go_annot->description);}
                else if( go_annot->branch == GOAnnotation::GO_MF ){ mf_list.push_back(go_annot->description);}
                else if( go_annot->branch == GOAnnotation::GO_CC ){ cc_list.push_back(go_annot->description);}
            }
            f_bp << gene << " = ("; f_mf << gene << " = ("; f_cc << gene << " = (";
            for(int j=0; j<int(bp_list.size()); j++){ f_bp << bp_list.at(j); if( j<( int(bp_list.size()) - 1)) { f_bp << "::";} }
            for(int j=0; j<int(mf_list.size()); j++){ f_mf << mf_list.at(j); if( j<( int(mf_list.size()) - 1)) { f_mf << "::";} }
            for(int j=0; j<int(cc_list.size()); j++){ f_cc << cc_list.at(j); if( j<( int(cc_list.size()) - 1)) { f_cc << "::";} }
            f_bp << ")\n"; f_mf << ")\n"; f_cc << ")\n";
        }
        catch(std::string err){
            // symbol not found
            continue;
        }
    }
    f_bp.close();
    f_mf.close();
    f_cc.close();    
}

void Spearman::write_extra_attribute_for_cytoscape( Attributes* ga, std::string fn_extra, Graph* G, HASH_I_VECTOR_STR& g2p ){
    // Only called once, but need to make write_to_cytoscape more reasonable in size

    double d;
    if( ga->attribute2idx.find(this->col_extra_attributes)==ga->attribute2idx.end() )
        throw new std::string( "Gene attributes does not have column: " + this->col_extra_attributes);
    std::ofstream f_extra(fn_extra.c_str());
    if( !f_extra.is_open() )
        throw new std::string( "Unable to open file for writing: " + fn_extra );
    f_extra << "NodeAttribute (java.lang.Double)\n";
    f_extra << "Cytoscape_bug_workaround = 1.0\n"; // Cytoscape gets confused about data type without this if first element has no decimal
    std::vector<int> nodes;
    G->nodes(nodes);
    std::string val_extra_attribute;
    for(int i=0; i<(int)nodes.size(); i++){
        std::vector<std::string>* v;
        if( g2p.find( nodes.at(i) ) != g2p.end() ){
            v =g2p[nodes.at(i)];
            for(int j=0; j<(int)v->size(); j++){
                val_extra_attribute = ga->prop_for_identifier(v->at(j), this->col_extra_attributes);
                if( val_extra_attribute.compare("NA") == 0 )
                    f_extra << G->get_node_label_by_index(nodes.at(i)) << " = NA\n";
                else{
                    try{
                        d = boost::lexical_cast<double>(val_extra_attribute);
                        f_extra.precision(2);
                        f_extra << G->get_node_label_by_index(nodes.at(i)) << " = " << d << "\n";                        
                    }
                    catch(  boost::bad_lexical_cast &){
                        std::cout << "ERROR: extra attribute column value could not be converted to floating-point number: " << val_extra_attribute << "\n";
                    }
                }
            }
        }
    }
    f_extra.close();
}

void Spearman::get_gene_name(Graph* G, int idx, std::string& g1){
    
    std::vector<std::string>* splitter = new std::vector<std::string>();
	char space[1];
	space[0] = ' ';  
    
    g1 = G->get_node_label_by_index(idx);
    if( g1.find(" ") != g1.npos ){
        stringhasher::split(g1, splitter, space);
        g1 = splitter->at(0);
    }
    delete splitter;
}

void Spearman::get_go_attributes_into_vectors( GOAnnotationParser* go, GeneAnnotationParser* gene_parser, std::string gene, std::vector<std::string> &bp, std::vector<std::string> &mf, std::vector<std::string> &cc ){
    // populates vectors bp, mf, cc with appropriate GO annotations
    bp.clear();
    mf.clear();
    cc.clear();
    if( gene_parser == NULL)
        return;
    GeneAnnotation* gene_annot;
    GOAnnotation* go_annot;
    try{
        gene_annot = gene_parser->get_annotation( gene );
    
        for(int g=0; g<int(gene_annot->GO_annotations.size() ); g++){
            go_annot = go->get_annotation( gene_annot->GO_annotations.at(g) );
            if(      go_annot->branch == GOAnnotation::GO_BP ){
                bp.push_back(go_annot->description);
            }
            else if( go_annot->branch == GOAnnotation::GO_MF ){
                mf.push_back(go_annot->description);
            }
            else if( go_annot->branch == GOAnnotation::GO_CC ){
                cc.push_back(go_annot->description);
            }
        }
    }
    catch( std::string str){
        // pass
    }
}

void Spearman::write_to_cytoscape_v3(double min_abs, std::string fn_base, Attributes* ga, GOAnnotationParser* go, GeneAnnotationParser* gene_parser, std::string fn_eQTL, double max_perm_pval, bool require_eQTL){
    // Cytoscape version 3 can be automated by passing in a script that reads table and
    // network files. This is distinct from version 2, which takes a network file and a
    // set of node/edge description files.
    
    Graph* G = new Graph();
	Graph* G_QTL = new Graph();
	HASH_I_VECTOR_STR g2p;
	if( this->verbose ){
		std::cout << "MESSAGE: Writing Cytoscape files to base " << fn_base << "\n";
		std::cout.flush();
	}
	prepare_graphs( G, G_QTL, g2p, min_abs, ga, fn_eQTL, max_perm_pval, require_eQTL);
	this->fn_out = fn_base + "_summary.spear";
	this->write_spears();
    
    std::string fn_net(fn_base + "_network.txt");
	std::string fn_noa(fn_base + ".noa");
    std::ofstream f_net(fn_net.c_str());
    std::ofstream f_noa(fn_noa.c_str());
    if( !f_net.is_open() )
		throw std::string( "Unable to open file for writing: " + fn_net);
    if( !f_noa.is_open() )
		throw std::string( "Unable to open file for writing: " + fn_noa);
    std::vector<std::vector<int>*>* e = new std::vector<std::vector<int>*>();
	G->edges(e);
	std::string g1, g2;
	double rho, pval;
	int idx1, idx2;
    std::string Dir;
    //gene1	gene2	interaction	rho pval
    //a	b	gg	0.5 0
    int n_cor_edges_written=0, n_eqtl_edges_written=0, n_genes_written=0, n_loci_written=0;
	n_cor_edges_written = int(e->size());
    
    f_net << "gene1\tinteraction\tgene2\tRhoA\tDir\tpval\n";
    for(int i=0; i<(int)e->size(); i++){
		idx1 = e->at(i)->at(0);
		idx2 = e->at(i)->at(1);
        get_gene_name(G, idx1, g1);
        get_gene_name(G, idx2, g2);
        rho = G->edge_weight( idx1, idx2 );
        if( rho >= 0 )
            Dir = "direct";
        else
            Dir = "inverse";
        if( rho==1 ){ // work around Cytoscape bug
            f_net << g1 << "\tgg\t" << g2 << "\t1.0\t" << Dir << "\tNA\n";
        }
        else if( rho == -1 ){
            f_net << g1 << "\tgg\t" << g2 << "\t-1.0\t" << Dir << "\tNA\n";
        }
        else{
            f_net << g1 << "\tgg\t" << g2 << "\t" << rho << "\t" << Dir << "\tNA\n";
        }
		n_cor_edges_written += 1;
	}
    G_QTL->edges(e);
	n_eqtl_edges_written = int(e->size());
	for(int i=0; i<(int)e->size(); i++){
		idx1 = e->at(i)->at(0);
		idx2 = e->at(i)->at(1);
        get_gene_name( G_QTL, idx1, g1 );
        get_gene_name( G_QTL, idx2, g2 );
        pval = G_QTL->edge_weight( idx1, idx2 );
        if( pval==1 ) // work around Cytoscape bug
            f_net << g1 << "\tgs\t" << g2 << "\tNA\tNA\t1.0\n";
        else if(pval==0 )
            f_net << g1 << "\tgs\t" << g2 << "\tNA\tNA\t0.0\n";
        else
            f_net << g1 << "\tgs\t" << g2 << "\tNA\tNA\t" << pval << "\n";
        delete e->at(i);
	}
	delete e;
    
    if( this->col_extra_attributes.size()>0 ){
        if( ga->attribute2idx.find(this->col_extra_attributes)==ga->attribute2idx.end() )
            throw new std::string( "Gene attributes does not have column: " + this->col_extra_attributes);
    }
    std::string val_extra_attribute;
    // write out node types {gene,locus} and possible extra attribute
    f_noa << "source\tNodeType\tNodeAttribute\tBiological_process\tMolecular_function\tCellular_component\n";
	std::vector<int> nodes;
	G->nodes(nodes);
	n_genes_written = int(nodes.size());
    HASH_S_I gene_seen;
    std::vector<std::string> bp_list,mf_list,cc_list;
    std::stringstream ss_bp, ss_mf, ss_cc;

	for(int i=0; i<int(nodes.size()); i++){
		get_gene_name( G, nodes.at(i), g1 );
        val_extra_attribute = ga->prop_for_identifier( ga->identifiers.at( nodes.at(i) ), this->col_extra_attributes);
        
        this->get_go_attributes_into_vectors( go, gene_parser, g1, bp_list, mf_list, cc_list );
        if( int( bp_list.size() )==0 ){ bp_list.push_back( std::string("NA") ); }
        if( int( mf_list.size() )==0 ){ mf_list.push_back( std::string("NA") ); }
        if( int( cc_list.size() )==0 ){ cc_list.push_back( std::string("NA") ); }
        ss_bp.str( "" ); ss_mf.str( "" ); ss_cc.str("");
        ss_bp.clear(); ss_mf.clear(); ss_cc.clear();
        for(int j=0; j<int(bp_list.size()); j++){ ss_bp << bp_list.at(j); if( j<( int(bp_list.size()) - 1)) { ss_bp << "::";} }
        for(int j=0; j<int(mf_list.size()); j++){ ss_mf << mf_list.at(j); if( j<( int(mf_list.size()) - 1)) { ss_mf << "::";} }
        for(int j=0; j<int(cc_list.size()); j++){ ss_cc << cc_list.at(j); if( j<( int(cc_list.size()) - 1)) { ss_cc << "::";} }
        f_noa << g1 << "\tgene\t" << val_extra_attribute << "\t" << ss_bp.str() << "\t" << ss_mf.str() << "\t" << ss_cc.str() << "\n";
        gene_seen[g1]=1;
	}
	G_QTL->nodes(nodes);
	for(int i=0; i<int(nodes.size()); i++){
		if( G_QTL->node_type( nodes.at(i) )==NODE_LOCUS ){
            get_gene_name( G_QTL, nodes.at(i), g1 );
            f_noa << g1 << "\tlocus\tNA\tNA\tNA\tNA\n";
			n_loci_written += 1;
		}
        else{
            get_gene_name( G_QTL, nodes.at(i), g1 );
            if( gene_seen.find(g1) == gene_seen.end() ){
                val_extra_attribute = ga->prop_for_identifier( ga->identifiers.at( nodes.at(i) ), this->col_extra_attributes);
                this->get_go_attributes_into_vectors( go, gene_parser, g1, bp_list, mf_list, cc_list );
                if( int(bp_list.size() )==0 ){ bp_list.push_back(std::string("NA")); }
                if( int(mf_list.size() )==0 ){ mf_list.push_back(std::string("NA")); }
                if( int(cc_list.size() )==0 ){ cc_list.push_back(std::string("NA")); }
                ss_bp.str( "" ); ss_mf.str( "" ); ss_cc.str("");
                ss_bp.clear(); ss_mf.clear(); ss_cc.clear();
                for(int j=0; j<int(bp_list.size()); j++){ ss_bp << bp_list.at(j); if( j<( int(bp_list.size()) - 1)) { ss_bp << "::";} }
                for(int j=0; j<int(mf_list.size()); j++){ ss_mf << mf_list.at(j); if( j<( int(mf_list.size()) - 1)) { ss_mf << "::";} }
                for(int j=0; j<int(cc_list.size()); j++){ ss_cc << cc_list.at(j); if( j<( int(cc_list.size()) - 1)) { ss_cc << "::";} }
                f_noa << g1 << "\tgene\t" << val_extra_attribute << "\t" << ss_bp.str() << "\t" << ss_mf.str() << "\t" << ss_cc.str() << "\n";
                gene_seen[g1]=1;
            }
        }
	}
	f_noa.close();

    
    // write command script
    std::string fn_cmd(fn_base + "_cytoscape_commands.txt");
	std::ofstream f_cmd(fn_cmd.c_str());
	if( !f_cmd.is_open() )
		throw new std::string( "Unable to open file for writing: " + fn_cmd);
    f_cmd << "network import file file=\"" << fn_net << "\"";
    f_cmd << " firstRowAsColumnNames=true startLoadRow=2 indexColumnSourceInteraction=1";
    f_cmd << " indexColumnTargetInteraction=3 indexColumnTypeInteraction=2\n";
    f_cmd << "table import file file=\"" << fn_noa << "\"";
    f_cmd << " firstRowAsColumnNames=true keyColumnIndex=1 startLoadRow=1\n";
    f_cmd << "vizmap load file file=\"" << fn_cytoscape_props << "\"\n";
    f_cmd << "vizmap apply styles=\"Correlation\"\n";
    f_cmd.close();
    
    delete G_QTL;
	delete G;
	for( boost::unordered_map<int, std::vector<std::string>* >::iterator i=g2p.begin(); i != g2p.end(); i++){
		delete i->second;
	}
	if( this->verbose ){
		std::cout << "MESSAGE: Wrote " << n_cor_edges_written << " gene-gene edges\n";
		std::cout << "MESSAGE: Wrote " << n_eqtl_edges_written << " gene-locus edges\n";
		std::cout << "MESSAGE: Wrote " << n_genes_written << " genes\n";
		std::cout << "MESSAGE: Wrote " << n_loci_written << " loci\n";
	}
    
}

void Spearman::write_to_cytoscape(double min_abs, std::string fn_base, Attributes* ga, GOAnnotationParser* go, GeneAnnotationParser* gene_parser, std::string fn_eQTL, double max_perm_pval, bool require_eQTL){

	Graph* G = new Graph();
	Graph* G_QTL = new Graph();
	HASH_I_VECTOR_STR g2p;
	if( this->verbose ){
		std::cout << "MESSAGE: Writing Cytoscape files to base " << fn_base << "\n";
		std::cout.flush();
	}
	prepare_graphs( G, G_QTL, g2p, min_abs, ga, fn_eQTL, max_perm_pval, require_eQTL);
	this->fn_out = fn_base + "_summary.spear";
	this->write_spears();
	std::string fn_sif(fn_base + ".sif");
	std::string fn_rho(fn_base + ".eda");
	std::string fn_pval(fn_base + "_pval.eda");
	std::string fn_type(fn_base + ".noa");
	std::string fn_dir(fn_base + "_dir.eda");
	std::ofstream f_sif(fn_sif.c_str());
	std::ofstream f_rho(fn_rho.c_str());
	std::ofstream f_dir(fn_dir.c_str());
	std::ofstream f_pval(fn_pval.c_str());
	std::ofstream f_type(fn_type.c_str());
	if( !f_sif.is_open() )
		throw std::string( "Unable to open file for writing: " + fn_sif);
	if( !f_rho.is_open() )
		throw std::string( "Unable to open file for writing: " + fn_rho);
	if( !f_pval.is_open() )
		throw std::string( "Unable to open file for writing: " + fn_pval);
	if( !f_type.is_open() )
		throw std::string( "Unable to open file for writing: " + fn_type);
	if( !f_dir.is_open() )
		throw std::string( "Unable to open file for writing: " + fn_dir);
	if( this->fn_cytoscape.size()==0 )
		throw std::string( "Must set fn_cytoscape");
	f_rho << "RhoA (java.lang.Double)\n";
	f_dir << "Dir (java.lang.String)\n";
	f_type << "NodeType (java.lang.String)\n";
	f_pval << "pvalue (java.lang.Double)\n";
	std::vector<std::vector<int>*>* e = new std::vector<std::vector<int>*>();
	G->edges(e);
	std::string g1, g2;
	double rho, pval;
	int idx1, idx2;
  
	std::string fn_bp(fn_base + "_bp.nda");
	std::string fn_mf(fn_base + "_mf.nda");
	std::string fn_cc(fn_base + "_cc.nda");

	if(gene_parser != NULL){
        write_GO_terms_for_cytoscape( gene_parser, go, fn_bp, fn_mf, fn_cc, G );
	}
	int n_cor_edges_written=0, n_eqtl_edges_written=0, n_genes_written=0, n_loci_written=0;
	n_cor_edges_written = int(e->size());

	for(int i=0; i<(int)e->size(); i++){
		idx1 = e->at(i)->at(0);
		idx2 = e->at(i)->at(1);
        get_gene_name(G, idx1, g1);
        get_gene_name(G, idx2, g2);
        rho = G->edge_weight( idx1, idx2 );
		f_sif << g1 << " gg " << g2 << "\n";
		f_rho << g1 << " (gg) " << g2 << " = ";
        if( rho == 1 )  // work around Cytoscape bug
			f_rho << "1.0\n";
		else if( rho == -1 )
			f_rho << "-1.0\n";
		else
			f_rho << rho << "\n";
		if( rho > 0 )
			f_dir << g1 << " (gg) " << g2 << " = direct\n";
		else
			f_dir << g1 << " (gg) " << g2 << " = inverse\n";
		n_cor_edges_written += 1;
	}
	G_QTL->edges(e);
	n_eqtl_edges_written = int(e->size());
	for(int i=0; i<(int)e->size(); i++){
		idx1 = e->at(i)->at(0);
		idx2 = e->at(i)->at(1);
        get_gene_name( G_QTL, idx1, g1 );
        get_gene_name( G_QTL, idx2, g2 );        
        pval = G_QTL->edge_weight( idx1, idx2 );
		f_sif << g1 << " gs " << g2 << "\n";
		f_pval << g1 << " (gs) " << g2 << " = ";
		if( pval == 1 )  // work around Cytoscape bug
			f_pval << "1.0\n";
		else if( pval == 0 )
            f_pval << "0.0\n";
        else
			f_pval << pval << "\n";
		delete e->at(i);
	}
	delete e;
	f_sif.close();
	f_rho.close();
	f_dir.close();
	// write out node types {gene,locus}
	std::vector<int> nodes;
	G->nodes(nodes);
	n_genes_written = int(nodes.size());
    HASH_S_I gene_seen;
	for(int i=0; i<int(nodes.size()); i++){
		get_gene_name( G, nodes.at(i), g1 );
		f_type << g1 << " = gene\n";
        gene_seen[g1]=1;
	}
	G_QTL->nodes(nodes);
	for(int i=0; i<int(nodes.size()); i++){
		if( G_QTL->node_type( nodes.at(i) )==NODE_LOCUS ){
            get_gene_name( G_QTL, nodes.at(i), g1 );
            f_type << g1 << " = locus\n";
			n_loci_written += 1;
		}
        else{
            get_gene_name( G_QTL, nodes.at(i), g1 );
            if( gene_seen.find(g1) == gene_seen.end() ){
                f_type << g1 << " = gene\n";
                gene_seen[g1]=1;
            }
        }
	}

	std::string fn_extra(fn_base + "_extra.noa"); // may not be used
	if( this->col_extra_attributes.size()>0 ){
		write_extra_attribute_for_cytoscape( ga, fn_extra, G, g2p );
	}

	std::string fn_bat(fn_base + "." + this->batch_extension);
	std::ofstream f_bat(fn_bat.c_str());
	if( !f_bat.is_open() )
		throw new std::string( "Unable to open file for writing: " + fn_bat);

	std::string result_header;
	generate_result_header(result_header);
	f_bat << result_header;

	if(this->limit_network_to_seeds)
		f_bat << "# Limit_to_seeds True\n";
	else
		f_bat << "# Limit_to_seeds False\n";
	f_bat << "# Seeds";
	for(int i=0; i<(int)this->seeds.size(); i++)
		f_bat << " " << seeds.at(i);
	f_bat << "\n";
	if(this->fn_cytoscape.at(0)=='\"')
		f_bat << this->fn_cytoscape;
	else
		f_bat << "\"" << this->fn_cytoscape << "\"";
	f_bat << " -N ";
	if(fn_sif.at(0)=='\"')
		f_bat << fn_sif;
	else
		f_bat << "\"" << fn_sif << "\"";
	f_bat << " -e ";
	if(fn_rho.at(0)=='\"')
		f_bat << fn_rho;
	else
		f_bat << "\"" << fn_rho << "\"";
	f_bat << " -e ";
	if(fn_dir.at(0)=='\"')
		f_bat << fn_dir;
	else
		f_bat << "\"" << fn_dir << "\"";
	if( this->col_extra_attributes.size()>0 ){
		f_bat << " -n ";
		if(fn_extra.at(0)=='\"')
			f_bat << fn_extra;
		else
			f_bat << "\"" << fn_extra << "\"";
	}
	if( this->fn_cytoscape_props.size()>0 ){
		f_bat << " -V ";
		if(this->fn_cytoscape_props.at(0)=='\"')
			f_bat << this->fn_cytoscape_props;
		else
			f_bat << "\"" << fn_cytoscape_props << "\"";
	}
	if( go != NULL ){
		f_bat << " -n ";
		if(fn_bp.at(0)=='\"'){ f_bat << fn_bp; } else { f_bat << "\"" << fn_bp << "\""; }
		f_bat << " -n ";
		if(fn_mf.at(0)=='\"'){ f_bat << fn_mf; } else { f_bat << "\"" << fn_mf << "\""; }
		f_bat << " -n ";
		if(fn_cc.at(0)=='\"'){ f_bat << fn_cc; } else { f_bat << "\"" << fn_cc << "\""; }
	}
	f_bat << " -n ";
	if(fn_type.at(0)=='\"')
		f_bat << fn_type;
	else
		f_bat << "\"" << fn_type << "\"";
	f_bat << " -e ";
	if(fn_pval.at(0)=='\"')
		f_bat << fn_pval;
	else
		f_bat << "\"" << fn_pval << "\"";

	f_bat << "\n";
	f_bat.close();
	delete G_QTL;
	delete G;
	for( boost::unordered_map<int, std::vector<std::string>* >::iterator i=g2p.begin(); i != g2p.end(); i++){
		delete i->second;
	}
	if( this->verbose ){
		std::cout << "MESSAGE: Wrote " << n_cor_edges_written << " gene-gene edges\n";
		std::cout << "MESSAGE: Wrote " << n_eqtl_edges_written << " gene-locus edges\n";
		std::cout << "MESSAGE: Wrote " << n_genes_written << " genes\n";
		std::cout << "MESSAGE: Wrote " << n_loci_written << " loci\n";
	}
}

void Spearman::results(std::vector< Spear* >& spears){
	spears.clear();
	for(int i=0; i<(int)this->spears.size(); i++){
		Spear* s = new Spear(this->spears.at(i));
		spears.push_back(s);
	}
}

int Spearman::n_spears(){
	return (int)this->spears.size();
}



/*
* THE FOLLOWING FUNCTIONS ARE EXPERIMENTAL CODE THAT IS NOT CALLED ANYWHERE OR MAINTAINED
*
*
void Spearman::find_DC_experimentwise_GWER(){
	// This fundamentally works differently from find_GWER, which is designed to test for correlation GWER
	// Here we will calculate element-wise GWER for each pair and store it.
	// The intended use is to calculate a FDR on the distribution of observed p-values
	// PROBLEM: the correlation structure is complex and not independent.
	int N = (int)this->idx.size();

	Matrix<double>* ranks_obs_a = new Matrix<double>();
	Matrix<double>* ranks_obs_b = new Matrix<double>();
	Matrix<double>* ranks_perm_a=NULL, *ranks_perm_b=NULL;
	Matrix<float>* data_perm = new Matrix<float>();

	boost::mt19937 state;
	RNG generator(state);
	bool is_spearman = false;
	std::vector<int> valid_columns_a, valid_columns_b;

	if(this->corr_type.compare("spearman") == 0){
		is_spearman = true;
		this->find_ranks(this->data->raw_data->data, *(this->data->a_idx), ranks_obs_a);
		ranks_perm_a = new Matrix<double>(ranks_obs_a->rows(), ranks_obs_a->cols(), ranks_obs_a->has_missing_values);
		for(int i=0; i<int(this->data->a_idx->size()); i++)
			valid_columns_a.push_back( i );
		ranks_obs_a->clone(ranks_perm_a);
		this->find_ranks(this->data->raw_data->data, *(this->data->b_idx), ranks_obs_b);
		ranks_perm_b = new Matrix<double>(ranks_obs_b->rows(), ranks_obs_b->cols(), ranks_obs_b->has_missing_values);
		for(int i=0; i<int(this->data->b_idx->size()); i++)
			valid_columns_b.push_back( i );
		ranks_obs_b->clone(ranks_perm_b);
	}

	int i, j=0, row1, row2;
	double rho_a, rho_b;
	std::vector<int> valid_columns_perm, perm_a_idx, perm_b_idx;
	std::vector<int> pvals;
	std::vector<double> dc_obs;
	int total_N = N * (N-1) / 2;
	pvals.reserve(total_N);
	dc_obs.reserve(total_N);

	int ctr=0;
	// Find observed differential correlation values, initialize pvals at n_perms
	if( this->verbose ){
		std::cout << "MESSAGE: Calculating observed differential correlation...\n";
		std::cout.flush();
	}

	for(i=0;i<N; i++){
		for(j=i+1; j<N; j++){
			row1 = idx.at(i);
			row2 = idx.at(j);
			if( is_spearman ){
				rho_a = this->find_spearman_from_ranks( ranks_obs_a, row1, ranks_obs_a, row2);
				rho_b = this->find_spearman_from_ranks( ranks_obs_b, row1, ranks_obs_b, row2);
			}
			else{
				rho_a = this->find_correlation_r( this->data->raw_data->data, row1, this->data->a_idx, this->data->raw_data->data, row2, this->data->a_idx);
				rho_b = this->find_correlation_r( this->data->raw_data->data, row1, this->data->b_idx, this->data->raw_data->data, row2, this->data->b_idx);
			}
			dc_obs.push_back( abs(rho_a-rho_b) );
			pvals.push_back( n_perms );
			ctr++;
		}
	}
	if( this->verbose ){
		std::cout << "MESSAGE: Calculated observed differential correlation for " << ctr+1 << " probe pairs\n";
		std::cout.flush();
	}
	// For each permutation, if perm_dc <= obs_dc than decrement pval
	for( int perm=0; perm<n_perms; perm++){
		if( is_spearman ){
            ranks_perm_a->shuffle_by_rows(valid_columns_a);
            ranks_perm_b->shuffle_by_rows(valid_columns_b);
        }
		else{
			permute_group_labels( this->data->a_idx, this->data->b_idx, perm_a_idx, perm_b_idx);
		}
		ctr=0;
		for(i=0;i<N; i++){
			for(j=i+1; j<N; j++){
				row1 = idx.at(i);
				row2 = idx.at(j);
				if( is_spearman ){
					// ranks_perm_a and ranks_perm_b have been shuffled
					rho_a = this->find_spearman_from_ranks( ranks_perm_a, row1, ranks_perm_a, row2);
					rho_b = this->find_spearman_from_ranks( ranks_perm_b, row1, ranks_perm_b, row2);
				}
				else{
					// only need to use shuffled index vectors (perm_a_idx, perm_b_idx)
					rho_a = this->find_correlation_r( this->data->raw_data->data, row1, &perm_a_idx, this->data->raw_data->data, row2, &perm_a_idx);
					rho_b = this->find_correlation_r( this->data->raw_data->data, row1, &perm_b_idx, this->data->raw_data->data, row2, &perm_b_idx);
				}
				if( abs(rho_a-rho_b) <= dc_obs.at(ctr) ){
					pvals.at(ctr) = pvals.at(ctr) - 1;
				}
				ctr++;
			}
			if(verbose && i% 100==0){
				std::cout << "MESSAGE: Completed row " << i+1 << " of " << N << "...\n";
				std::cout.flush();
			}
		}
		if(verbose){
			std::cout << "MESSAGE: Completed permutation " << perm+1 << " of " << n_perms << "...\n";
			std::cout.flush();
		}
	}
	ctr=0;
	for(i=0;i<N; i++){
		for(j=i+1; j<N; j++){
			row1 = idx.at(i);
			row2 = idx.at(j);
			this->spears.push_back(new Spear(pvals.at(ctr), 0, row1, row2));
			ctr++;
		}
	}
	delete ranks_obs_a;
	delete ranks_obs_b;
	delete ranks_perm_a;
	delete ranks_perm_b;
	delete data_perm;
}

double Spearman::calculate_differential_correlations_probeset_pval(double mean_diff_obs){
	if( this->n_perms==0)
		return 1;
	int N = int( this->idx.size() );
	int i, j=0, row1, row2;
	double rho_a, rho_b;
	Matrix<double>* ranks_obs_a = new Matrix<double>();
	Matrix<double>* ranks_obs_b = new Matrix<double>();

	bool is_spearman = false;

	if(this->corr_type.compare("spearman") == 0){
		is_spearman = true;
		this->find_ranks(this->data->raw_data->data, *(this->data->a_idx), ranks_obs_a);
		this->find_ranks(this->data->raw_data->data, *(this->data->b_idx), ranks_obs_b);
		std::vector<int> valid_cols;
		for(int i=0; i<int(this->data->a_idx->size()) ; i++)
			valid_cols.push_back(this->data->a_idx->at(i) );
		for(int i=0; i<int(this->data->b_idx->size()); i++)
			valid_cols.push_back(this->data->b_idx->at(i) );
	}
	Matrix<double>* ranks_perm_a = new Matrix<double>(ranks_obs_a->rows(), ranks_obs_a->cols(), ranks_obs_a->has_missing_values);
	Matrix<double>* ranks_perm_b = new Matrix<double>(ranks_obs_b->rows(), ranks_obs_b->cols(), ranks_obs_b->has_missing_values);

	std::vector<int> a_idx_perm, b_idx_perm;
	double sum_diff_perm;
	double n_comparisons=0.0;
	double p = n_perms;

    for( int perm=0; perm<n_perms; perm++){
    	sum_diff_perm=0;
    	n_comparisons=0;
    	if( is_spearman ){
    		//permute_ranks_between_groups( ranks_obs_a, ranks_obs_b, a_idx_perm, b_idx_perm, ranks_perm_a, ranks_perm_b );
    		shuffle_ranks_between_labels(ranks_obs_a, ranks_obs_b, ranks_perm_a, ranks_perm_b );
    		for(i=0;i<N; i++){
    		    for(j=i+1; j<N; j++){
    			    row1 = idx.at(i);
    			    row2 = idx.at(j);
    			    n_comparisons += 1;
    			    rho_a = this->find_spearman_from_ranks( ranks_perm_a, row1, ranks_obs_a, row2);
    			    rho_b = this->find_spearman_from_ranks( ranks_perm_b, row1, ranks_obs_b, row2);
    			    sum_diff_perm += abs( rho_a - rho_b );
    		    }
    	    }
    	}
    	else{
    		permute_group_labels( this->data->a_idx, this->data->b_idx, a_idx_perm, b_idx_perm ); // swap labels between groups, not just their order within groups
    		for(i=0;i<N; i++){
    			for(j=i+1; j<N; j++){
    				row1 = idx.at(i);
    				row2 = idx.at(j);
    				n_comparisons += 1;
				    rho_a = this->find_correlation_r( this->data->raw_data->data, row1, &a_idx_perm, this->data->raw_data->data, row2, &a_idx_perm);
				    rho_b = this->find_correlation_r( this->data->raw_data->data, row1, &b_idx_perm, this->data->raw_data->data, row2, &b_idx_perm);
				    sum_diff_perm += abs( rho_a - rho_b );
    			}
		    }
	    }
	    if( abs(mean_diff_obs) >= (sum_diff_perm / n_comparisons) )
	    	p = p-1;
    }
    p = p / n_perms;
	delete ranks_obs_a;
	delete ranks_obs_b;
	delete ranks_perm_a;
	delete ranks_perm_b;
	return p;
}

void Spearman::calculate_differential_correlation_by_probesets(ProbeSets& ps, std::vector<DifferentalCorrelationResult*>& DCRs){
	// Calculate change in correlation by calculating correlation between
	for(int i=0; i<int(DCRs.size()); i++)
		delete DCRs.at(i);
	DCRs.clear();
	load_data_for_run();
	this->set_corr_abs(0.0);
	this->set_corr_diff(0.0);
	double meanA, meanB, perm_pval;

	ps.reset_iterator();
	ProbeSet* PS;
	PS = ps.next_probeset();
	int ctr=0;
	while( PS != NULL ) {
		ctr += 1;
		this->seeds.clear();
		this->idx.clear();
		for(int i=0; i<int(PS->probelist.size()); i++){
			this->seeds.push_back(PS->probelist.at(i));
			this->idx.push_back( int(this->data->raw_data->identifier2idx[ PS->probelist.at(i) ] ) );
		}
		this->find_spears();
		if( this->get_mean_difference(meanA, meanB) ){
			if( this->verbose ){
				std::cout << "MESSAGE: Calculating permutations for set " << ctr << " of " << ps.size() << ", " << this->spears.size() << " pairs\n";
				std::cout.flush();
			}
			perm_pval = calculate_differential_correlations_probeset_pval(meanA-meanB);
			DCRs.push_back( new DifferentalCorrelationResult( PS->title, int(PS->probelist.size()), meanA, meanB, (meanA-meanB), perm_pval) );
		}
		PS = ps.next_probeset();
	}
}

void Spearman::shuffle_ranks_between_labels(Matrix<double>* ranks_obs_a, Matrix<double>* ranks_obs_b, Matrix<double>* ranks_perm_a, Matrix<double>* ranks_perm_b ){
	std::vector< std::pair<int, int>* > sample_idx;
	const int OBS_A = 1;
	const int OBS_B = 2;
	ranks_perm_a->set_all_values_present();
	ranks_perm_b->set_all_values_present();

	for(int i=0; i<ranks_obs_a->cols(); i++)
		sample_idx.push_back( new std::pair<int,int>( OBS_A, i ) );
	for(int i=0; i<ranks_obs_b->cols(); i++)
		sample_idx.push_back( new std::pair<int,int>( OBS_B, i ) );

	boost::mt19937 state;
	RNG generator(state);

	std::random_shuffle(sample_idx.begin(), sample_idx.end(), generator);
	int current_idx, current_col=0;
	Matrix<double>* current_source;
	for(int i=0; i<ranks_obs_a->cols(); i++){
		current_idx=i;
		current_col = sample_idx.at(current_idx)->second;
		if( sample_idx.at(current_idx)->first==OBS_A )
			current_source = ranks_obs_a;
		else
			current_source = ranks_obs_b;
		for(int r=0; r<ranks_obs_a->rows(); r++){
			ranks_perm_a->arr[r][i] = current_source->arr[r][current_col];
			if( current_source->is_missing(r, current_col) )
				ranks_perm_a->set_missing_value(r, i);
		}
	}
	for(int i=0; i<ranks_obs_b->cols(); i++){
		current_idx = ranks_obs_b->cols() + i;
		current_col = sample_idx.at(current_idx)->second;
		if( sample_idx.at(current_idx)->first==OBS_A )
			current_source = ranks_obs_a;
		else
			current_source = ranks_obs_b;
		for(int r=0; r<ranks_obs_b->rows(); r++){
			ranks_perm_b->arr[r][i] = current_source->arr[r][current_col];
			if( current_source->is_missing(r, current_col) )
				ranks_perm_b->set_missing_value(r, i);
		}
	}
	for(int i=0; i<int(sample_idx.size() ); i++)
		delete sample_idx.at(i);
}


*/
