#include <iostream>
#include <sstream>
#include <vector>
#include <list>
#include <string>
#include <algorithm>
#include <time.h>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/random.hpp>
#include <boost/thread.hpp>
#include <boost/bind.hpp>
#include <boost/thread/mutex.hpp>
using namespace boost;
using namespace std;
#include "DataStructures.h"
#include "ClassMinerOptions.h"
#include "Rawdata.h"
#include "Attributes.h"
#include "Discretize.h"
#include "Dataset.h"
#include "ParseOptions.h"
#include "Parser.h"
#define MATHLIB_STANDALONE 1
#include "Rmath.h"
#include "QTL_Calculator.h"


rqtl::rqtl(){
    this->idx_genotype=-1;
    this->idx_probe = -1;
    this->pvalue = -1;
    this->min_permuted_pval = -1;
    this->mean_a=-99999;
    this->mean_b=-99999;
    this->mean_c=-99999;
    this->t_stat=0;
}

rqtl::rqtl(int idx_genotype, int idx_probe, double pvalue, double min_permuted_pval, double t_stat){
    this->idx_genotype = idx_genotype;
    this->idx_probe = idx_probe;
    this->pvalue = pvalue;
    this->min_permuted_pval = min_permuted_pval;
    this->mean_a=-99999;
    this->mean_b=-99999;
    this->mean_c=-99999;
    this->t_stat=t_stat;
}

bool sort_regressed_qtl(const rqtl* a, const rqtl* b){
	return a->pvalue < b->pvalue;
}

qtl_by_chrom::qtl_by_chrom(int idx_probe, int n_chrom){
    this->idx_probe = idx_probe;
    for(int i=0; i<n_chrom; i++){
        this->T_obs.push_back(0);
        this->snp_indexes.push_back(-1);
        this->is_robust_adjustment.push_back(false);
        this->pvalues_raw.push_back(0);
        this->pvalues_perm.push_back(0);
    }
}


double qtl_by_chrom::get_max_T_obs(){
    // Return the largest T statistic stored in this qtl across all chromosomes
    if( this->T_obs.size()==0 )
        return 0;
    double max_tstat= this->T_obs.at(0);
    for(int i=0; i<int(this->T_obs.size()); i++){
        if( this->T_obs.at(i) > max_tstat )
            max_tstat = this->T_obs.at(i);
    }
    return max_tstat;
}

void qtl_by_chrom::update_perm_pvalue(std::vector<double>& max_perm_tstats){
    // given a the highest t statistic from a permutation, increment all permutation pvalues for chromosomes
    // that have a lower value
    for(int i=0; i<int(this->T_obs.size()); i++){
        if( this->T_obs.at(i) <= max_perm_tstats.at(i) ){
            this->pvalues_perm.at(i) += 1;
        }
    }
}


void qtl_by_chrom::calculate_final_pvalues(int n_samples, int n_perms){
    for(int i=0; i<int(this->T_obs.size()); i++){
        this->pvalues_raw.at(i) = pt( this->T_obs.at(i), (double)n_samples-2, 0, 0)*2;
        if(n_perms>0){
            this->pvalues_perm.at(i) = pvalues_perm.at(i) / n_perms;
        }
        else{
            this->pvalues_perm.at(i) = 1;            
        }
    }
}


QTL_Calculator::QTL_Calculator(Dataset* data_expr, Attributes* sa_expr, Attributes* ga_expr, Dataset* data_snps, Attributes* sa_snps, Attributes* ga_snps, std::string shared_g_p_col){
    this->data_expr = data_expr;
    this->sa_expr = sa_expr;
    this->ga_expr = ga_expr;
    this->data_snps = data_snps;
    this->sa_snps = sa_snps;
    this->ga_snps = ga_snps;
    this->expr_snps.clear();
    this->is_robust=false;
    this->n_threads=1;
    this->n_perms=0;
    this->n_perms_max=0;
    this->index_for_next_thread=0;
    this->verbose=true;
    this->analysis_method = METHOD_REGRESSION;
    match_expr_and_snps( shared_g_p_col );
    // by default use all probesets
    for(int i=0; i<int(data_expr->raw_data->identifiers.size()); i++){
        this->idx.push_back(i);
    }
}

int QTL_Calculator::get_cis_interval(){
    // 0 for genome-wide; non-zero for base interval to examine on genome
    return this->cis_interval;
}


void QTL_Calculator::get_matched_snp_sample_idx( std::vector<int>& idx_to_fill ){
    // fill idx with the indices of the SNP samples we've matched to expression data
    idx_to_fill.clear();
    for(int i=0; i<int(this->expr_snps.size()); i++)
        idx_to_fill.push_back(this->expr_snps.at(i)->second );
}


void QTL_Calculator::get_matched_gene_sample_idx( std::vector<int>& idx_to_fill ){
    // fill idx with the indices of the expression samples we've matched to SNP data
    idx_to_fill.clear();
    for(int i=0; i<int(this->expr_snps.size()); i++)
        idx_to_fill.push_back(this->expr_snps.at(i)->first );
}


int QTL_Calculator::get_size_gene_idx(){
    // number of probesets to analyze
    return int(this->idx.size());
}


int QTL_Calculator::get_n_matched_samples(){
    // number of matched gene expression-SNP samples
    return (int)this->expr_snps.size();
}


void QTL_Calculator::set_robust(bool robust){
    this->is_robust=robust;
}

void QTL_Calculator::set_cis_interval(int cis_interval){
    this->cis_interval=cis_interval;
}

void QTL_Calculator::set_n_perms(int n_perms){
    this->n_perms = n_perms;
}


void QTL_Calculator::set_logfile( std::string fn_log ){
    this->fn_log = fn_log;
    this->fs_log.open( fn_log.c_str() );
    if( !this->fs_log.is_open() ){
        throw std::string( "Unable to open file for writing: " + fn_log );
    }
}


void QTL_Calculator::close_logfile_if_open(){
    this->write_to_log( std::string("Closed log file cleanly") );
    if( this->fs_log.is_open() )
        this->fs_log.close();
}

void QTL_Calculator::set_analysis_method( int method ){
    if(method != METHOD_REGRESSION && method != METHOD_STUDENT_T_TEST )
        throw(std::string("Unrecognized method value, should be 1 or 2"));
    this->analysis_method=method;
}


void QTL_Calculator::set_thread_count(int n_threads){
    this->n_threads=n_threads;
}

void QTL_Calculator::set_n_perms_max(int n_perms_max){
    this->n_perms_max = n_perms_max;
}

int QTL_Calculator::get_min_obs_per_genotype(){
    return this->min_obs_per_genotype;
}

void QTL_Calculator::set_min_obs_per_genotype(int min){
    this->min_obs_per_genotype=min;
}

void QTL_Calculator::set_verbose(bool verbose){
    this->verbose = verbose;
}


void QTL_Calculator::sort_QTLs(){
    std::sort(this->qtls.begin(), this->qtls.end(), sort_regressed_qtl);   
}

//***************************************************
// Setup and bookkeeping code
//***************************************************

void QTL_Calculator::check_probes_in_KV( KeyValuesParser& KV ){
    // ensure that all probesets specified in KV are in the dataset
    std::stringstream ss;
    std::vector<std::string> gene_probe_ids, snp_probes;
    
    KV.keys(gene_probe_ids);
    for(int i=0; i<int(gene_probe_ids.size()); i++){
        if( this->data_expr->raw_data->identifier2idx.find(gene_probe_ids.at(i)) == this->data_expr->raw_data->identifier2idx.end() ){
            ss << "Phenotype probe not found in dataset: " << gene_probe_ids.at(i);
            throw( std::string(ss.str() ) );
        }
        KV.get(gene_probe_ids.at(i), snp_probes);
        for(int s=0; s<int(snp_probes.size()); s++){
            if( this->data_snps->raw_data->identifier2idx.find( snp_probes.at(s) ) == this->data_snps->raw_data->identifier2idx.end() ){
                ss << "Genotype probe not found in dataset: " << snp_probes.at(s);
                throw( std::string(ss.str() ) );
            }
        }
    }
}


void QTL_Calculator::match_expr_and_snps(std::string shared_g_p_col){
	// populates this->expr_snps
    //
    // Use the column indicated by shared_g_p_col to match up identifiers in genotype with identifeirs
	// in phenotype.  A 0 indicates a sample that is not shared.
	// pairs of indices that are shared are pushed into expr_snps as (idx_expr, idx_snps)
	std::vector<int>* shared = new std::vector<int>();
	sa_snps->indices_without_property(shared_g_p_col, "0", shared);
	// indices in a_idx are correct for data_expr; need to translate them into indices in snps sa.
	HASH_S_I allowed_expr_ids;
	for(int i=0; i<(int)this->data_expr->a_idx->size(); i++){
		int idx_in_expr = this->data_expr->a_idx->at(i);
		allowed_expr_ids[ this->data_expr->sample_names->at( idx_in_expr ) ] = idx_in_expr;
	}
    HASH_S_I snp_ids2idx;
    for(int i=0; i<(int)this->data_snps->sample_names->size(); i++){
        snp_ids2idx[ this->data_snps->sample_names->at(i) ]= i;
    }
    
	string expr_id;
	string snp_id;
	std::pair<int, int>* expr_snp;
    boost::unordered_map<int, bool> seen_snps, seen_genes;
	for(int i=0; i<(int)shared->size(); i++){
        snp_id = this->sa_snps->identifiers.at( shared->at(i) );
		expr_id = this->sa_snps->prop_for_identifier(snp_id, shared_g_p_col);
        if( allowed_expr_ids.find(expr_id) != allowed_expr_ids.end() ){
			expr_snp = new std::pair<int, int>();
			expr_snp->first = allowed_expr_ids[ expr_id ];
            expr_snp->second = snp_ids2idx[ snp_id ];
            if( seen_genes.find( expr_snp->first ) != seen_genes.end() ){
                std::stringstream ss;
                ss << "phenotype identifier " << expr_id << " matched more than once between SNPs and phenotypes with index " << expr_snp->first << "\n";
                throw std::string(ss.str());
            }
            if( seen_snps.find( expr_snp->second ) != seen_snps.end() ){
                std::stringstream ss;
                ss << "genotype identifier " << snp_id << " matched more than once between SNPs and phenotypes with index " << expr_snp->second << "\n";
                throw std::string(ss.str());
            }
			// convert index in SNP Sample Attributes to index in SNP expr
			//snp_id = this->sa_snps->identifiers.at( shared->at(i) );
			//expr_snp->second = this->data_snps->raw_data->sample_name2idx[snp_id];
			this->expr_snps.push_back( expr_snp );
            seen_genes[ expr_snp->first ] = 1;
            seen_snps[ expr_snp->second ] = 1;
        }
	}
	delete shared;
}


void QTL_Calculator::calculate_snp_idx2chr(boost::unordered_map<int, int>& snp_idx2chr, std::string& chr_colname, int& n_chr){
    // create hash of snp idx -> chromosome number.
    // assumes that X chromosome, if present, is at the end of the list.  
    // assumes that chromosome that is not a number is the X chromosome.
    // Otherwise does not assume sorted order
    std::string c_val;
    int chr, max_chr_so_far = 0;
    for(int i=0; i<(int)this->ga_snps->identifiers.size(); i++){
        c_val = ga_snps->prop_for_identifier(this->ga_snps->identifiers.at(i),chr_colname);
        try{
            chr = boost::lexical_cast<int>(c_val);
            snp_idx2chr[i] = chr;
            if( chr > max_chr_so_far )
                max_chr_so_far = chr;
            n_chr = max_chr_so_far;
        }
        catch( bad_lexical_cast & ){
            // should be X, numbered as last seen Chr + 1
            snp_idx2chr[i] = max_chr_so_far + 1;
            n_chr = max_chr_so_far + 1;
        }
    }
}


void QTL_Calculator::calculate_sample_indexes(std::vector<int>* true_sample_idx, std::vector< std::vector<int>* >& perm_sample_idx, int n_samples){
    // utility function to calculate true and permutation sample indexes
    // Note that we create n_perms_max permutation indexes
    
    for( int s=0; s<n_samples; s++){
        true_sample_idx->push_back( expr_snps.at(s)->second ); // don't permute these
    }
    boost::mt19937 state;
    RNG generator(state);
    generator._state.seed(static_cast<unsigned int>(std::time(0)) );
    
    for( int n=0; n<this->n_perms_max; n++){
        std::vector<int>* v = new std::vector<int>();
        for( int s=0; s<n_samples; s++){
            v->push_back( expr_snps.at(s)->second ); // push snp idx
        }
        std::random_shuffle( v->begin(), v->end(), generator );
        perm_sample_idx.push_back(v);
    }
}


//***************************************************
// Simple statistical operations: mean, variance, etc
//***************************************************

double QTL_Calculator::find_var(std::vector<double>* C){
	int N = (int)C->size();
	if(N<=1)
		return 0;
	else{
		double sum=0;
		double mean=0;
		double sum_of_diff_squared=0.0;
		double diff;
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


void QTL_Calculator::find_mean(std::vector<double>* V, std::vector<int>* ab, double& mean_a, double& mean_b){
	int N = (int)V->size();
	mean_a = 0;
	mean_b = 0;
	double na=0, nb=0;
	double sum_a=0, sum_b=0;
	const int LBL_A = -1;
	const int LBL_B = 1;                                                                                                                    
	if(N==0)
		return;
	else{
		for(int i=0; i<N; i++){
			if( ab->at(i)==LBL_A ){
				sum_a += V->at(i);
				++na;
			}
			else if(ab->at(i) == LBL_B ){
				++nb;
				sum_b += V->at(i);
			}
		}
		if(na>0)
			mean_a = sum_a/na;
		if(nb>0)
			mean_b = sum_b/nb;
	}
}

void QTL_Calculator::find_mean_three_way(std::vector<double>* V, std::vector<double>* ab, double& mean_a, double& mean_b, double& mean_c){
    // Returns mean of values in V, segregated by genotypes {0,1,2} in ab
	// Used by calculate_genome_wide().  Assumption is that Y axis (ab) is genotype and therefore
    // actually an integer in {0,1,2}
    std::vector<int>* ab_int = new std::vector<int>();
    for(int i=0; i<(int)ab->size(); i++)
        ab_int->push_back( (int)ab->at(i) );
    find_mean_three_way(V, ab_int, mean_a, mean_b, mean_c);
    delete ab_int;
}


void QTL_Calculator::find_mean_three_way(std::vector<double>* V, std::vector<int>* ab, double& mean_a, double& mean_b, double& mean_c){
	// Returns mean of values in V, segregated by genotypes {0,1,2} in ab
	int N = (int)V->size();
	mean_a = 0; mean_b = 0; mean_c = 0;
	double na=0, nb=0, nc=0;
	double sum_a=0, sum_b=0, sum_c=0;
	const int LBL_A = 0; const int LBL_B = 1; const int LBL_C = 2;
	if(N==0)
		return;
	else{
		for(int i=0; i<N; i++){
			if( ab->at(i)==LBL_A ){
				sum_a += V->at(i);
				++na;
			}
			else if(ab->at(i) == LBL_B ){
				++nb;
				sum_b += V->at(i);
			}
            else if(ab->at(i) == LBL_C ){
				++nc;
				sum_c += V->at(i);
			}
		}
		if(na>0)
			mean_a = sum_a/na;
		if(nb>0)
			mean_b = sum_b/nb;
        if(nc>0)
			mean_c = sum_c/nc;
	}
}


void QTL_Calculator::find_mean_sd_three_way(std::vector<double>* V, std::vector<double>* ab, double& mean_a, double& mean_b, double& mean_c, double& sd_a, double& sd_b, double& sd_c){
	// Returns mean and SD of values in V, segregated by genotypes {0,1,2} in ab
	int N = (int)V->size();
	mean_a = 0; mean_b = 0; mean_c = 0;
    sd_a = 0; sd_b = 0; sd_c = 0;
	double na=0, nb=0, nc=0;
	double sum_a=0, sum_b=0, sum_c=0;
    std::vector<double> A, B, C;
	const int LBL_A = 0; const int LBL_B = 1; const int LBL_C = 2;
	if(N==0)
		return;
	else{
		for(int i=0; i<N; i++){
			if( int(ab->at(i))==LBL_A ){
                A.push_back(V->at(i));
				sum_a += V->at(i);
				++na;
			}
			else if(int(ab->at(i)) == LBL_B ){
				B.push_back(V->at(i));
                ++nb;
				sum_b += V->at(i);
			}
            else if(int(ab->at(i)) == LBL_C ){
				C.push_back(V->at(i));
                ++nc;
				sum_c += V->at(i);
			}
		}
		if(na>0){
			mean_a = sum_a/na;
            sd_a = sqrt(find_variance(mean_a, A));
        }
		if(nb>0){
			mean_b = sum_b/nb; 
            sd_b = sqrt(find_variance(mean_b, B));
            
        }
        if(nc>0){
			mean_c = sum_c/nc;
            sd_c = sqrt(find_variance(mean_c, C));
        }
	}
}


double QTL_Calculator::find_variance(double mean, const vector<double>& acc){
	// calculate sample variance
	double sum_of_diff_squared=0.0;
	double diff;
	int N = (int)acc.size();
	for(int i=0; i<N; i++){
		diff = acc.at(i) - mean;
		sum_of_diff_squared += diff * diff;
	}
	if( N-1>0 )
		return sum_of_diff_squared / (N - 1);
	else
		return 0;
}



//***************************************************
// Calculating which probesets to include
//***************************************************

void QTL_Calculator::limit_ids_by_var(double min_var){
    if(min_var==0){
        return;
    }
    std::vector<int> keep;

	std::vector<double>* C = new std::vector<double>();
	int n, c;
	double variance;
	for(int r=0; r<(int)this->data_expr->raw_data->identifiers.size(); r++){
		C->clear();
		for(n=0; n<(int)this->data_expr->a_idx->size(); n++){
			c = this->data_expr->a_idx->at(n);
			C->push_back( this->data_expr->raw_data->data->arr[r][c] );
		}
		for(n=0; n<(int)this->data_expr->b_idx->size(); n++){
			c = this->data_expr->b_idx->at(n);
			C->push_back( this->data_expr->raw_data->data->arr[r][c] );
		}
		variance = find_var(C);
		if( variance >= min_var )
			keep.push_back(r);
	}
	delete C;
    this->idx.clear();
    for(int i=0; i<int(keep.size()); i++)
        this->idx.push_back( keep.at(i) );
}


void QTL_Calculator::limit_ids_by_NA(double fraction_required){
    // restrict idx by elminating any rows where 
    // either class has fewer than fraction_required percent present
    
    if( this->data_expr->raw_data->data->has_missing_values == false )
        return;
    
    std::vector<int> keep;
    int n_in_a, a;

    int n_req_a = (int)((double)this->data_expr->a_idx->size() * fraction_required );
    for(int i=0; i<int(idx.size()); i++){
        n_in_a = 0;
        if(this->data_expr->raw_data->data->row_has_missing_value(i)){
            for(a=0;a<(int)this->data_expr->a_idx->size();a++){
                if(!this->data_expr->raw_data->data->is_missing(i,a)){
                    n_in_a++;
                }
            }
            if( n_in_a>=n_req_a ){
                keep.push_back(i);
            }
        }
        else{
            keep.push_back(i);
        }
    }
    this->idx.clear();
	for(int i=0; i<(int)keep.size(); i++)
		this->idx.push_back( keep.at(i) );
}


void QTL_Calculator::limit_ids_by_subset(int this_fraction, int n_fractions){
    // divide idx into n_fractions sections, keep only those in this_fraction
    // used to reproducibly divide up work without dividing raw data into physical files.
    std::vector<int> keep;
    int ctr=0;
    for(int i=0; i<int(idx.size()); i++){
        if( ctr+1 == this_fraction )
            keep.push_back(idx.at(i));
        ctr++;
        if(ctr==n_fractions)
            ctr=0;
    }
    this->idx.clear();
    for(int i=0; i<int(keep.size()); i++)
		this->idx.push_back( keep.at(i) );
}


void QTL_Calculator::limit_ids_by_NA_min_max(double fraction_required, double fraction_max){
    // further restrict idx by elminating any rows where 
    // either class has fewer than two samples present
    std::vector<int> keep;
    int n_in_a, a;
    int n_req_a_min = (int)((double)this->data_expr->a_idx->size() * fraction_required );
    int n_req_a_max = (int)((double)this->data_expr->a_idx->size() * fraction_max );
    for(int i=0; i<(int)idx.size(); i++){
        n_in_a = 0;
        for(a=0;a<(int)this->data_expr->a_idx->size();a++){
            if(!this->data_expr->raw_data->data->is_missing(i,a)){
                n_in_a++;
            }
        }
        if( n_in_a>=n_req_a_min && n_in_a < n_req_a_max ){
            keep.push_back(i);
        }
    }
    this->idx.clear();
	for(int i=0; i<(int)keep.size(); i++)
		this->idx.push_back( keep.at(i) );
}


void QTL_Calculator::limit_ids_by_probe_id(std::string probe_id){
    // restrict idx to a single probe ID
    this->idx.clear();
    if( this->data_expr->raw_data->identifier2idx.find(probe_id) != this->data_expr->raw_data->identifier2idx.end() ){
        this->idx.push_back( this->data_expr->raw_data->identifier2idx[ probe_id ] );
    }
    else{
        throw std::string("Probeid requested for limit was not found in the list of identifiers.");
    }
}



void QTL_Calculator::restrict_probesets_and_snps( std::vector<std::string>& probesets, std::vector<std::string>& snps, std::vector<int> & idx_snps){
    // Overrides any other limitations, and is therefore incompatible with other limitations.
    // It is an error to pass a probe ID that is not known present in data_expr
    std::vector<int> keep;
    for(int j=0; j<int(probesets.size()); j++){
        if( data_expr->raw_data->identifier2idx.find( probesets.at(j) ) != data_expr->raw_data->identifier2idx.end() )
            keep.push_back( data_expr->raw_data->identifier2idx[probesets.at(j)] );
        else{
            std::stringstream ss;
            ss << "unrecognized probeset: " << probesets.at(j);
            throw ss.str();
        }
    }
    this->idx.clear();
	for(int i=0; i<(int)keep.size(); i++)
		this->idx.push_back( keep.at(i) );
    keep.clear();
    for(int j=0; j<int(snps.size()); j++){
        if( data_snps->raw_data->identifier2idx.find( snps.at(j) ) != data_snps->raw_data->identifier2idx.end() )
            keep.push_back( data_snps->raw_data->identifier2idx[snps.at(j)] );
        else{
            std::stringstream ss;
            ss << "unrecognized snp: " << snps.at(j);
            throw ss.str();
        }
    }
    idx_snps.clear();
	for(int i=0; i<(int)keep.size(); i++)
		idx_snps.push_back( keep.at(i) );
}



//***************************************************
// Code to perform linear regressions
//***************************************************


void QTL_Calculator::fill_XY( std::vector<double>& X, std::vector<double>& Y, int snps_g_idx, int expr_g_idx, std::vector<int>* sample_idx){
    // given a SNP and gene, write the gene value into X and SNP value into Y for
    // all pairs where both are present in data_expr and data_snps respectively.
    // Only used for naive regession approach; too slow for frequent use
    bool missing_expr, missing_snps, expr_is_missing, snp_is_missing;
    int expr_s_idx, snps_s_idx;
    X.clear();
    Y.clear();
    missing_expr = this->data_expr->raw_data->data->row_has_missing_value(expr_g_idx);
	missing_snps = this->data_snps->raw_data->data->row_has_missing_value(snps_g_idx);
	if( missing_expr || missing_snps ){
		for(int p=0;p<(int)expr_snps.size();p++){
			expr_s_idx = this->expr_snps.at(p)->first;
			snps_s_idx = sample_idx->at(p);
			if( missing_expr )
				expr_is_missing = this->data_expr->raw_data->data->is_missing(expr_g_idx, expr_s_idx);
			else
				expr_is_missing = false;
			if( missing_snps )
				snp_is_missing = this->data_snps->raw_data->data->is_missing(snps_g_idx, snps_s_idx);
			else
				snp_is_missing = false;
			if( !expr_is_missing && !snp_is_missing ){
                snps_s_idx = sample_idx->at(p);
				Y.push_back( this->data_snps->raw_data->data->arr[snps_g_idx][snps_s_idx] );
				X.push_back( this->data_expr->raw_data->data->arr[expr_g_idx][expr_s_idx] );
			}
		}
	}
	else{
		for(int p=0;p<(int)expr_snps.size();p++){
			expr_s_idx = this->expr_snps.at(p)->first;
			snps_s_idx = sample_idx->at(p);
			Y.push_back( this->data_snps->raw_data->data->arr[snps_g_idx][snps_s_idx] );
            X.push_back( this->data_expr->raw_data->data->arr[expr_g_idx][expr_s_idx] );
		}
	}
}


void QTL_Calculator::modify_XY_to_be_robust( std::vector<double>& X, std::vector<double>& Y){
    // After calling fill_XY to get all gene-SNP pairs
    // eliminate all observations in a genotype if there are less than MIN_OBS observations for that genotype
    // if there are only two genotypes, require 20 observations (need to make these parameters eventually)
    
    int n_0=0, n_1=0, n_2=0, i, geno;
    std::vector<double> Y_rob, X_rob;
    
    for(i=0; i<int(Y.size()); i++){
        geno = int(Y.at(i));
        if(geno==0)
            n_0++;
        else if(geno==1)
            n_1++;
        else if(geno==2)
            n_2++;   
    }
    for(i=0; i<int(Y.size()); i++){
        if(Y.at(i)==0 && n_0>=this->min_obs_per_genotype ){
            X_rob.push_back(X.at(i)); 
            Y_rob.push_back(Y.at(i));
        }
        else if(Y.at(i)==1 && n_1> this->min_obs_per_genotype ){
            X_rob.push_back(X.at(i)); 
            Y_rob.push_back(Y.at(i));
        }
        else if(Y.at(i)==2 && n_2>=this->min_obs_per_genotype ){
            X_rob.push_back(X.at(i)); 
            Y_rob.push_back(Y.at(i));
        }
    }
    X.clear();
    Y.clear();
    for(i=0; i<int(Y_rob.size()); i++){
        X.push_back(X_rob.at(i));
        Y.push_back(Y_rob.at(i));
    }
}


bool QTL_Calculator::student_t_direct(std::vector<int>* sample_idx, int expr_g_idx, int snps_g_idx, double& t_stat){
	// find two-tailed p-value of a the Student's t-test
	// Assumes equal variance, so DF = n_a + n_b - 2  xxx
    return false;
    /*
    bool missing_expr = this->data_expr->raw_data->data->row_has_missing_value(expr_g_idx);
    bool missing_snps = this->data_snps->raw_data->data->row_has_missing_value(snps_g_idx);
    
    // two possible models: group A with AB, group B with AB
    // calculate both simultaneously
    double mean_a_AA_AB_vs_BB=0, mean_b_AA_AB_vs_BB=0, mean_a_AA_vs_AB_BB=0, mean_b_AA_VS_AB_BB=0; 
    double n_a_AA_AB_vs_BB=0, n_b_AA_AB_vs_BB=0, n_a_AA_vs_AB_BB=0, n_b_AA_VS_AB_BB=0; 
	double M2_a_AA_AB_vs_BB=0, M2_b_AA_AB_vs_BB=0, M2_a_AA_vs_AB_BB=0, M2_b_AA_VS_AB_BB=0;

    double x=0, delta;
    int genotype;
    int expr_s_idx, snps_s_idx;
    
    if( !missing_expr && !missing_snps ){
        n = int( expr_snps.size() );
        
        for(int p=0; p<n; p++){
			expr_s_idx = this->expr_snps.at(p)->first;
			snps_s_idx = sample_idx->at(p);
			genotype = this->data_snps->raw_data->data->arr[snps_g_idx][snps_s_idx];
			x = this->data_expr->raw_data->data->arr[expr_g_idx][expr_s_idx];
            
            // AA_AB vs BB model
            if( genotype==0 || genotype==1){
                n_a_AA_AB_vs_BB++;
                delta = x - mean_a_AA_AB_vs_BB;
                mean_a_AA_AB_vs_BB = mean_a_AA_AB_vs_BB + (delta/n_a_AA_AB_vs_BB);
                if(n_a_AA_AB_vs_BB>1)
                    M2_a_AA_AB_vs_BB += delta*(x-mean_a_AA_AB_vs_BB);
            }
        }
    }
    
    if( n_a_AA_AB_vs_BB > 1 )
        variance_a_AA_AB_vs_BB = M2_a_AA_AB_vs_BB / (n_a_AA_AB_vs_BB - 1);
    if( n_b_AA_AB_vs_BB > 1 )
        variance_b_AA_AB_vs_BB = M2_b_AA_AB_vs_BB / (n_b_AA_AB_vs_BB - 1);
    if( n_a_AA_vs_AB_BB > 1 )
        variance_a_AA_vs_AB_BB = M2_a_AA_vs_AB_BB / (n_a_AA_vs_AB_BB - 1);
    if( n_b_AA_vs_AB_BB > 1 )
        variance_b_AA_vs_AB_BB = M2_b_AA_vs_AB_BB / (n_b_AA_vs_AB_BB - 1);    

    

            if( genotype==0 ){
                sum_a_AA_AB_vs_BB += x;
                
            }
            if( genotype==1 ){
                sum_a_AA_AB_vs_BB += x;
                n_a_AA_AB_vs_BB++; // it's an A if we're grouping AB with AA
                sum_b_AA_vs_AB_BB += x;
                n_b_AA_vs_AB_BB++  // it's a B if we're grouping AB with AA
                
            }            
            if( genotype==2 ){
                sum_b_AA_AB_vs_BB += x;
                n_b_AA_AB_vs_BB++
            }
        }
    }
    bool allow_AA_AB_vs_BB = false;
    bool allow_AA_vs_AB_BB = false
    if( n_a_AA_AB_vs_BB > this->min_obs_per_genotype && n_b_AA_AB_vs_BB > this->min_obs_per_genotype ){
        allow_AA_AB_vs_BB = true;
        if( n_a_AA_AB_vs_BB>0  && n_b_AA_AB_vs_BB>0 ){
            mean_a_AA_AB_vs_BB = sum_a_AA_AB_vs_BB / n_a_AA_AB_vs_BB;
            mean_b_AA_AB_vs_BB = sum_b_AA_AB_vs_BB / n_b_AA_AB_vs_BB;
            allow_AA_AB_vs_BB = true;
        }
    }
    if( n_a_AA_vs_AB_BB > this->min_obs_per_genotype && n_b_AA_vs_AB_BB > this->min_obs_per_genotype ){
        if( n_a_AA_VS_AB_BB>0 && n_b_AA_VS_AB_BB>0 ){
            mean_a_AA_VS_AB_BB = sum_a_AA_VS_AB_BB / n_a_AA_VS_AB_BB;
            mean_b_AA_VS_AB_BB = sum_b_AA_VS_AB_BB / n_b_AA_VS_AB_BB;
            allow_AA_AB_vs_BB = true;
        }
    }
    if( !allow_AA_AB_vs_BB && !allow_AA_vs_AB_BB )
        return true; // NA
    
	for(int i=0; i<(int)ab->size(); i++){
		if(ab->at(i)==LBL_A){
			A.push_back( vals->at(i) );
			sum_a += vals->at(i);
		}
		else if(ab->at(i)==LBL_B){
			B.push_back( vals->at(i) );
			sum_b += vals->at(i);
		}
	}
	n_a = (int)A.size();
	n_b = (int)B.size();
	if( n_a>0 )
		mean_a = sum_a/n_a;
	if( n_b>0 )
		mean_b = sum_b/n_b;
    */
    /*
     double sum_of_diff_squared=0.0;
     double diff;
     int N = (int)acc.size();
     for(int i=0; i<N; i++){
        diff = acc.at(i) - mean;
        sum_of_diff_squared += diff * diff;
     }
     if( N-1>0 )
     return sum_of_diff_squared / (N - 1);
     else
     return 0;
     */
    /*
	double var_a = find_variance(mean_a, A);
	double var_b = find_variance(mean_b, B);
    
    double diff_means = mean_a - mean_b;
	double standard_error = (double)sqrt( (var_a/n_a) + (var_b/n_b) );
	t_stat = diff_means / standard_error;
 */   	
}



bool QTL_Calculator::regress_direct(std::vector<int>* sample_idx, int expr_g_idx, int snps_g_idx, double& m, double& b, double& t_stat){
    // This code performs the same operation as regress()
    // However, it does so without creating the arrays X and Y.
    double sumx=0, sumy=0, sumx2=0, sumxy=0;
    bool missing_expr = this->data_expr->raw_data->data->row_has_missing_value(expr_g_idx);
    bool missing_snps = this->data_snps->raw_data->data->row_has_missing_value(snps_g_idx);
    double n=0;
    double x, y;
    int expr_s_idx, snps_s_idx;
    if( !missing_expr && !missing_snps ){
        n = int( expr_snps.size() );
        for(int p=0; p<n; p++){
			expr_s_idx = this->expr_snps.at(p)->first;
			snps_s_idx = sample_idx->at(p);
			y = this->data_snps->raw_data->data->arr[snps_g_idx][snps_s_idx];
			x = this->data_expr->raw_data->data->arr[expr_g_idx][expr_s_idx];
            //std::cout << this->data_expr->raw_data->sample_names[expr_s_idx] << "\t" << x << "\t" << this->data_snps->raw_data->sample_names[snps_s_idx] << "\t" << y << "\n";
            sumy += y;
            sumx += x;
            sumx2 += (x*x);
            sumxy += (x*y);
        }
    }
    else{
        n=0;
        for(int p=0; p<int(expr_snps.size()); p++){
            expr_s_idx = this->expr_snps.at(p)->first;
			snps_s_idx = sample_idx->at(p);			
            if( this->data_expr->raw_data->data->is_missing(expr_g_idx, expr_s_idx) )
                continue;
            if( this->data_snps->raw_data->data->is_missing(snps_g_idx, snps_s_idx) )
                continue;
            n++;
			y = this->data_snps->raw_data->data->arr[snps_g_idx][snps_s_idx];
			x = this->data_expr->raw_data->data->arr[expr_g_idx][expr_s_idx];
            sumy += y;
            sumx += x;
            sumx2 += (x*x);
            sumxy += (x*y);
        }
    }
    if( n < 2 ){
		m = b = t_stat = 0;
		return true;
	}
	double divisor = (sumx2 - ((sumx * sumx) / n));
	if(divisor == 0){
		m = b = t_stat = 0;
		return true;
	}
	else{
		m = (sumxy - ((sumx * sumy) / n)) / divisor;
		b = (sumy - ((m) * sumx)) / n;
        
		// code for standard error of the residuals and slope
		double mean = sumx / n;
		double mean_sum=0, ye2_sum=0;
        double d, predict;
        if( !missing_expr && !missing_snps ){
            for(int p=0; p<n; p++){
                expr_s_idx = this->expr_snps.at(p)->first;
                snps_s_idx = sample_idx->at(p);
                y = this->data_snps->raw_data->data->arr[snps_g_idx][snps_s_idx];
                x = this->data_expr->raw_data->data->arr[expr_g_idx][expr_s_idx];
                d = x - mean;
                mean_sum += d*d;
                predict = m*x + b;
                ye2_sum += (y - predict) * (y - predict);
            }
        }
        else{
            for(int p=0; p<int(expr_snps.size()); p++){
                expr_s_idx = this->expr_snps.at(p)->first;
                snps_s_idx = sample_idx->at(p);
                if( this->data_expr->raw_data->data->is_missing(expr_g_idx, expr_s_idx) )
                    continue;
                if( this->data_snps->raw_data->data->is_missing(snps_g_idx, snps_s_idx) )
                    continue;
                y = this->data_snps->raw_data->data->arr[snps_g_idx][snps_s_idx];
                x = this->data_expr->raw_data->data->arr[expr_g_idx][expr_s_idx];
                d = x - mean;
                mean_sum += d*d;
                predict = m*x + b;
                ye2_sum += (y - predict) * (y - predict);
            }    
        }
		double sx2 = sqrt(mean_sum);
		double residual_se = sqrt(ye2_sum / (n - 2)) ;
		double slope_se = residual_se/sx2;
        if( slope_se==0 )
            return true;
		t_stat = m / slope_se; // div by zero -> Infinite t, which turns out to be fine.
	}
    return false;
    
}
    
bool QTL_Calculator::regress(vector<double>& X, vector<double>& Y, double& m, double& b, double& t_stat){
    // returns false if answer is NaN or otherwise compromised
	// http://tolstoy.newcastle.edu.au/R/e2/help/06/09/0985.html
	// http://david.swaim.com/cpp/linregc.htm
	// pvalue = pt(t_stat, df, 0, 0);
    // Naive operation; frequently called code must use regress_direct
	double sumx=0, sumy=0, sumx2=0, sumxy=0;
	double n = (double)X.size();
	if( X.size() != Y.size() )
		throw string("x and y must have same length");
	if( n < 2 ){
		m = b = t_stat = 0;
		return true;
	}
	double x,y;
	for(int i=0; i<(int)n; i++){
		x = X.at(i);
		y = Y.at(i);
		sumx += x;
		sumy += y;
		sumx2 += (x*x);
		sumxy += (x*y);
	}

	double divisor = (sumx2 - ((sumx * sumx) / n));
	if(divisor == 0){
		m = b = t_stat = 0;
		return true;
	}
	else{
		m = (sumxy - ((sumx * sumy) / n)) / divisor;
		b = (sumy - ((m) * sumx)) / n;

		// code for standard error of the residuals and slope
		double mean = sumx / n;
		double mean_sum=0, ye2_sum=0;
        double d, predict;
		for(int i=0; i<(int)X.size(); i++){
			d = X.at(i) - mean;
			mean_sum += d*d;
			predict = m*X.at(i) + b;
			ye2_sum += (Y.at(i) - predict) * (Y.at(i) - predict);
		}
		double sx2 = sqrt(mean_sum);
		double residual_se = sqrt(ye2_sum / (n - 2)) ;
		double slope_se = residual_se/sx2;
        if( slope_se==0 )
            return true;
		t_stat = m / slope_se; // div by zero -> Infinite t, which turns out to be fine.
	}
    return false;
}


bool QTL_Calculator::calculation_wrapper( int snps_g_idx, int expr_g_idx, std::vector<int>* sample_idx, double& t_stat, double value_to_exceed){
    // value_to_exceed is passed if we only want to bother checking robust results if the result appear to exceed some bounds
    std::vector<double> Y, X;
    double t_stat_robust, m, b;
    
    // calculate observed statistic
    bool is_nan;
    if( this->analysis_method==METHOD_REGRESSION )
       is_nan = regress_direct(sample_idx, expr_g_idx, snps_g_idx, m, b, t_stat);
    else
       is_nan = student_t_direct(sample_idx, expr_g_idx, snps_g_idx, t_stat);        
        
    // If user requests "robust" analysis and the raw T statistic exceeds some threshold, recalculate as robust
    if( this->is_robust && abs(t_stat) > value_to_exceed ){
        // be conservative; report "robust" tstat if it is WEAKER than original.   
        // If robust analysis causes us to return NAN, than there is a problem with the data 
        // (e.g. we've eliminated all but one genotype) and we NEED to return that NAN.
        // TODO: improve this by combining these steps
        fill_XY( X, Y, snps_g_idx, expr_g_idx, sample_idx);
        this->modify_XY_to_be_robust( X, Y );
        is_nan = regress( X, Y, m, b, t_stat_robust );
        if( !is_nan && abs(t_stat_robust) < abs(t_stat) )
            t_stat = t_stat_robust;
    }
    return is_nan;
}


//***************************************************
// Scan the full or partial genome for a given gene
//***************************************************

void QTL_Calculator::find_best_tstat_across_loci( int expr_g_idx, std::vector<int>* sample_idx, std::vector<int>* genotype_idx, int& idx_best_locus, double& t_obs ){
    // regress gene at idx expr_g_idx against all samples in sample_idx and all SNPs indexed by genotype_idx
    // used for cis-only analysis, but could be used with arbitrary SNP indices.
    int snps_g_idx;
    double t_stat;
    bool is_nan = false;
    t_obs=0;
    for(int i=0; i<int(genotype_idx->size()); i++){
		snps_g_idx = genotype_idx->at(i);
        is_nan = calculation_wrapper( snps_g_idx, expr_g_idx, sample_idx, t_stat, abs(t_obs) );
        if(!is_nan && abs(t_stat) > abs(t_obs) ){
            t_obs=t_stat;
            idx_best_locus = snps_g_idx;   
        }
    }   
}



bool QTL_Calculator::permutation_tstat_across_loci_exceeds( int expr_g_idx, std::vector<int>* expr_g_indexes, std::vector<int>* sample_idx, double tstat_observed, double& max_tstat_perm ){
    // Given a permuted ordering of samples, is there a single t statistic >= tstat_observed? 
    // same process as find_best_tstat_across_loci, but stops when a permutation tstat exceeds tstat_observed
    // since once we've exceeded the observed statistic with this permutation any further analysis is needless
	double t_stat;
    bool is_nan=false;
    double abs_tstat_observed = abs(tstat_observed);
    
    max_tstat_perm=0;
    int snps_g_idx;
    
    for(int i=0; i<int(expr_g_indexes->size()); i++){
        snps_g_idx = expr_g_indexes->at(i);
        is_nan = calculation_wrapper( snps_g_idx, expr_g_idx, sample_idx, t_stat, max_tstat_perm ); 
        if( !is_nan && abs(t_stat) >= abs_tstat_observed ){
            // a permutation Tstat exceeds the observed T statistic. Stop looking.
            max_tstat_perm=abs(t_stat); 
            return true;
        }
        else{
            if( !is_nan && abs(t_stat) > max_tstat_perm){
                // The permutation value did not exceed the observed Tstat, but it exceeds any previously seen value
                max_tstat_perm=abs(t_stat);    
            }
        }
    }
    return false;
}


void QTL_Calculator::find_best_tstat_by_chr( int expr_g_idx, std::vector<int>* sample_idx, std::vector<int>& idx_Tstat_maxes, std::vector<double>& Tstat_stopvals, std::vector<double>& Tstat_maxes, boost::unordered_map<int,int>& snp_idx2chr ){
    // Given a probeset at index expr_g_idx and an arrangment of samples sample_idx
    // for each chromosome, calculate all Tstat and store abs(Tstat_max) and index[ abs(Tstat_max) ] 
    // The length of idx_Tstat_maxes and Tstat_maxes should equal the number of chromosomes in snp_idx2chr
    // The length of Tstat_stopvals should be 0 or the same as Tstat_maxes.
    //
    // USE IN PERMUTATION ANALYSIS:
    // if values are present in Tstat_stopvals, stop calculating regressions for a given chromosome
    // if we see a Tstat that exceeds the Tstat_stopval. 
    //
    // If we are performing a robust analysis and we see a new "best observed t statistic", 
    // check whether that statistic is still the best after correcting for outliers by heuristic
    // If so, record the original observed statistic.
	
    double t_stat;
    int chr_idx;
    bool doing_permutations = false, is_nan=false;

    if( idx_Tstat_maxes.size() != Tstat_maxes.size() )
    	throw std::string("Tstat_maxes and idx_Tstat_maxes must have same length");
    if( Tstat_stopvals.size() != 0 && Tstat_stopvals.size() != Tstat_maxes.size() )
        throw std::string("Tstat_stopvals must have length 0 or same length as Tstat_maxes");
    
    if( Tstat_stopvals.size() > 0 )
        doing_permutations = true;
    
    for( int i=0; i<int(Tstat_maxes.size()); i++){
        Tstat_maxes.at(i)=0;
        idx_Tstat_maxes.at(i)=0;
    }
    
	for(int snps_g_idx=0; snps_g_idx<(int)data_snps->identifiers->size(); snps_g_idx++){
		chr_idx = snp_idx2chr[snps_g_idx] - 1;
        if( doing_permutations ){
            if(Tstat_maxes.at(chr_idx) >= Tstat_stopvals.at(chr_idx) ) // T_perm has exceeded T_obs for this chromosome
                continue;
        }
        is_nan = calculation_wrapper( snps_g_idx, expr_g_idx, sample_idx, t_stat, Tstat_maxes.at(chr_idx) );           
        if( !is_nan && abs(t_stat) > Tstat_maxes.at(chr_idx) ){
            Tstat_maxes.at(chr_idx) = abs(t_stat);
            idx_Tstat_maxes.at(chr_idx) = snps_g_idx;
        }
    }
}

double QTL_Calculator::permute_by_chromosome( int expr_g_idx, double T_obs, int chr_idx, std::vector< std::vector<int>* >& perm_sample_idx, boost::unordered_map<int,int>& snp_idx2chr ){
    // For one gene on one chromosome, perform n_perms_max permutations and return the resulting P value
    // There should already be the correct number of perm_sample_idx sample-sets created.
      
    double P=0;
    std::vector<int> idx_snps_in_this_chr;
    bool is_nan=false;
    double T_perm=0;
    int snps_g_idx;
    // identify the SNPs on chromosome chr_idx. Note that chr_idx equals real chromosome minus one.
    for(snps_g_idx=0; snps_g_idx<int(data_snps->identifiers->size()); snps_g_idx++){
        if( snp_idx2chr[snps_g_idx] - 1 == chr_idx ){
            idx_snps_in_this_chr.push_back(snps_g_idx);
        }
    }
    for( int itr=0; itr<this->n_perms_max; itr++){
        for(int i=0; i<int(idx_snps_in_this_chr.size()); i++){
            snps_g_idx = idx_snps_in_this_chr.at(i);
            is_nan = calculation_wrapper( snps_g_idx, expr_g_idx, perm_sample_idx.at(itr), T_perm, 0 ); 
            if(!is_nan && abs(T_perm)>=T_obs){
                P++;
                break;
            }
        }
    }
    return P / double( this->n_perms_max );
}


//***************************************************
// Thread worker functions
//***************************************************


int QTL_Calculator::request_next_index_to_process(){
    // If we're done, return -1
    // We're returning the index into idx, not the actual value of idx.
    // This is so individual threads can report where we are in the vector; 
    // otherwise we wouldn't be able to see how far we've gotten
    {
        // MUST be locked to continue
        boost::mutex::scoped_lock lock( this->thread_iter_mutex); 
        
        int index_of_gene_to_process = -1;
        
        if( this->index_for_next_thread<this->idx.size()){
            index_of_gene_to_process = this->index_for_next_thread;
            this->index_for_next_thread++;
        }
        return index_of_gene_to_process;
    }
}


void QTL_Calculator::process_genomewide_in_thread(int thread_id, std::vector<int>* true_sample_idx, std::vector< std::vector<int>* >& perm_sample_idx){
    int idx_best_locus, pval_final, idx_probe_id;
    double mu_A, mu_B, mu_C;
    double tstat_observed, pval_observed, max_tstat_perm;
    int n_samples = (int)this->expr_snps.size();
    
    rqtl* new_qtl;
    std::vector<double> X,Y;
    {
        boost::mutex::scoped_lock lock( this->io_mutex );
        std::cout << "MESSAGE: Initiated thread " << thread_id+1 << ".\n";
    }
    // In some calls to permutation_tstat_across_loci_exceeds we wish to restrict the set of loci under consideration
    // Here we use the whole available set of snps (normally meaning the whole genome)
    std::vector<int>* snp_indexes = new std::vector<int>();
    for(int i=0; i<int(this->data_snps->raw_data->identifiers.size()); i++)
        snp_indexes->push_back(i);
    
    std::string result_str;
    bool is_single_probe=false;
    if( this->idx.size()==1 )
        is_single_probe=true;
    
    int current_idx = request_next_index_to_process();
    while( current_idx != -1 ){
        idx_probe_id = this->idx.at(current_idx);
        pval_final = 0;
        find_best_tstat_across_loci( idx_probe_id, true_sample_idx, snp_indexes, idx_best_locus, tstat_observed );
        tstat_observed = abs(tstat_observed);
        for( int x=0; x<this->n_perms; x++){
            if( permutation_tstat_across_loci_exceeds( idx_probe_id, snp_indexes, perm_sample_idx.at(x), abs(tstat_observed), max_tstat_perm ) )
                pval_final++;
            if( is_single_probe && x % 50==0 ){
                report_thread_progress( thread_id, x+1, int(this->n_perms) );
            }
        }
        pval_observed = pt(tstat_observed, (double)n_samples-2, 0, 0)*2;
        new_qtl = new rqtl(idx_best_locus, idx_probe_id, pval_observed, (double) pval_final / this->n_perms, tstat_observed);
        fill_XY( X, Y, idx_best_locus, idx_probe_id, true_sample_idx);
        find_mean_three_way(&X, &Y, mu_A, mu_B, mu_C);
        new_qtl->mean_a = mu_A; new_qtl->mean_b = mu_B; new_qtl->mean_c = mu_C;
        {
            boost::mutex::scoped_lock lock( this->results_mutex );
            this->qtls.push_back( new_qtl );
        }
        std::stringstream ss;
        this->spear_result_string(new_qtl, result_str);
        ss << "[thread " << thread_id+1 << " gene " << current_idx+1 << " of " << this->idx.size() << "]\t" << result_str;
        this->write_to_log(ss.str());
        if( this->verbose ){
            report_thread_progress( thread_id, current_idx+1, int(this->idx.size()) );
        }
        current_idx = request_next_index_to_process();
    }
    
    delete snp_indexes;
}


void QTL_Calculator::process_by_chr_in_thread(int thread_id, std::vector<int>* true_sample_idx, std::vector< std::vector<int>* >& perm_sample_idx){
    
    // Identify the chromosome assignent for each SNP
    boost::unordered_map<int,int> snp_idx2chr;
    int n_chr=0;
    int n_samples = (int)this->expr_snps.size();
    std::string chr_column_name = this->ga_snps->get_chr_column();
    this->calculate_snp_idx2chr(snp_idx2chr, chr_column_name, n_chr);
    
    qtl_by_chrom* Q;
    std::vector<double> max_tstats_perm;
    std::vector<int> idx_max_tstats_perm;    
    for(int i=0; i<n_chr; i++){
        max_tstats_perm.push_back(0);
        idx_max_tstats_perm.push_back(0); // don't actually care about this, need for function signature
    }
    int idx_probe_id=0;
    std::vector<double> empty_vector;
    std::string result_str;
    
    double MAX_P_FOR_DEEP_PERM=0.1;
    if(this->n_perms >= 100)
        MAX_P_FOR_DEEP_PERM=0.01;
    
    int current_idx = request_next_index_to_process();
    while( current_idx != -1 ){
        idx_probe_id = this->idx.at(current_idx);
        Q = new qtl_by_chrom(idx_probe_id, n_chr);        
        find_best_tstat_by_chr( idx_probe_id, true_sample_idx, Q->snp_indexes, empty_vector, Q->T_obs, snp_idx2chr );
        for( int x=0; x<this->n_perms; x++){
            find_best_tstat_by_chr( idx_probe_id, perm_sample_idx.at(x), idx_max_tstats_perm, Q->T_obs, max_tstats_perm, snp_idx2chr );
            Q->update_perm_pvalue(max_tstats_perm); // if max_tstat_perm >= observed T statistic, increment pval_perm
        }
        Q->calculate_final_pvalues(n_samples, n_perms); 
        if( this->n_perms_max > this->n_perms ){
            for( int chr_idx=0; chr_idx<n_chr; chr_idx++){
                if( Q->pvalues_perm.at(chr_idx) <= MAX_P_FOR_DEEP_PERM ){
                    Q->pvalues_perm.at(chr_idx) = permute_by_chromosome( idx_probe_id, Q->T_obs.at(chr_idx), chr_idx, perm_sample_idx, snp_idx2chr );
                }
            }
        }
        this->spear_result_string(Q, result_str);
        this->write_to_log(result_str);
        {
            boost::mutex::scoped_lock lock( this->results_mutex );
            this->qtls_by_chrom.push_back( Q );
        }
        if( this->verbose ){
            report_thread_progress( thread_id, current_idx+1, int(this->idx.size()) );
        }
        current_idx = request_next_index_to_process();
    }
}



void QTL_Calculator::process_cis_only_in_thread(int thread_id, std::vector<int>* true_sample_idx, std::vector< std::vector<int>* >& perm_sample_idx, std::vector<int>* idx_all_snps, int cis_window_size ){
    // Given a defined window size around each gene locus (cis_window_size), find the single strongest eQTL
    // Originally, permutations were defined against all such windows, not just the local window
    // Now permutations are defined against the local window    
    {
        boost::mutex::scoped_lock lock( this->io_mutex );
        std::cout << "MESSAGE: Initiated thread " << thread_id+1 << ".\n";
    }
    std::vector<int>* genotype_idx = new std::vector<int>();
    int idx_best_locus, locus, idx_probe_id, loc_start, loc_end, pval_final;
    int half_window = cis_window_size/2;
    std::string chromosome, identifier;
    double tstat_observed, pval_observed;
    rqtl* Q;
    double mu_A, mu_B, mu_C, max_tstat_perm;
    std::vector<double> Y, X;
    
    int n_samples = (int)this->expr_snps.size();

    int current_idx = request_next_index_to_process();
    std::string result_str;
    while( current_idx != -1 ){
        idx_probe_id = this->idx.at(current_idx);
        identifier = this->data_expr->ga->identifiers.at(idx_probe_id);
        if( this->data_expr->ga->get_chromosome_and_locus_by_identifier( identifier, chromosome, locus) ){
            loc_start = locus-half_window;
            loc_end = locus+half_window;
            this->data_snps->ga->get_idx_by_genomic_range(chromosome, loc_start, loc_end, *genotype_idx);
            
            tstat_observed=0;
            if( int(genotype_idx->size()) == 0 ){
                pval_final=-1;
                pval_observed = -1;
            }
            else{
                pval_final=0;
                find_best_tstat_across_loci( idx_probe_id, true_sample_idx, genotype_idx, idx_best_locus, tstat_observed );
                tstat_observed = abs(tstat_observed);
                for( int x=0; x<this->n_perms; x++){
                    if( permutation_tstat_across_loci_exceeds( idx_probe_id, genotype_idx, perm_sample_idx.at(x), abs(tstat_observed), max_tstat_perm ) )
                        pval_final++;
                }
                pval_observed = pt(tstat_observed, (double)n_samples-2, 0, 0)*2;
            }
                        
            Q = new rqtl(idx_best_locus, idx_probe_id, pval_observed, (double) pval_final / this->n_perms, tstat_observed);
            fill_XY( X, Y, idx_best_locus, idx_probe_id, true_sample_idx);
            find_mean_three_way(&X, &Y, mu_A, mu_B, mu_C);
            Q->mean_a = mu_A;
            Q->mean_b = mu_B;
            Q->mean_c = mu_C;
            this->spear_result_string(Q, result_str);
            this->write_to_log(result_str);
            {
                boost::mutex::scoped_lock lock( this->results_mutex );
                this->qtls.push_back( Q );
            }
            if( this->verbose ){
                report_thread_progress( thread_id, current_idx+1, int(this->idx.size()) );
            }
        }
        current_idx = request_next_index_to_process();
    }
    delete genotype_idx;
}


//***************************************************
// Call scanning functions to process genes
//***************************************************


void QTL_Calculator::calculate_specific_SNPs(Matrix<double>& results, std::vector<int>& idx_snps){
    // Called when calculating "probe networks", i.e. when given a particular set of SNPs and probes to evaluate
    // Calculates all SNPs in idx_snps for the probes in this->idx
    // T statistics are store in Matrix results
    // NOT YET MULTITHREADED
    
    int snps_g_idx, expr_g_idx;
    double t_stat;
    bool is_nan;
    
    int n_samples = (int)this->expr_snps.size();
    std::vector<int>* sample_idx = new std::vector<int>();
    for( int s=0; s<n_samples; s++){
        sample_idx->push_back( expr_snps.at(s)->second ); // don't permute these
    }

    results.resize( int( idx_snps.size() ), int(this->idx.size() ) );
    for(int row=0; row<(int)idx_snps.size(); row++){
        snps_g_idx = idx_snps.at(row);
        for(int col=0; col<int(this->idx.size()); col++){
            expr_g_idx = this->idx.at(col);
            is_nan = calculation_wrapper(snps_g_idx, expr_g_idx, sample_idx, t_stat, 0 );
            results.arr[row][col] = t_stat;
        }
    }
    delete sample_idx;
}





void QTL_Calculator::calculate_permutations_for_gene_probe_pairs( KeyValuesParser& KV ){
    // NOT YET MULTITHREADED
    
    int idx_gene, idx_snp;
    rqtl* Q;
    bool is_nan;
    
    std::vector<int>* true_sample_idx = new std::vector<int>();
    std::vector< std::vector<int>* > perm_sample_idx;
    int n_samples = (int)this->expr_snps.size();
    calculate_sample_indexes(true_sample_idx, perm_sample_idx, n_samples);
    
    double tstat_observed, max_tstat_perm, pval_perm, pval_observed;
    std::vector<std::string> gene_probe_ids, snp_probes;
    KV.keys(gene_probe_ids);
    
    std::vector<int>* snp_indexes;
    boost::unordered_map<int,int> snp_idx2chr;
    int n_chr=0, snp_chr=0;
    std::string chr_column_name = this->ga_snps->get_chr_column();
    if( chr_column_name.size()>0 ){
        //restrict SNP idx to SNPs on the same chromosome as idx_snp
        this->calculate_snp_idx2chr(snp_idx2chr, chr_column_name, n_chr);
        snp_indexes = new std::vector<int>();
    }
    else{
        snp_indexes = data_snps->a_idx; // use whole genome 
    }
    
    std::vector<double> X, Y;
    double mu_A, mu_B, mu_C;
    
    for(int i=0; i<int(gene_probe_ids.size()); i++){
        idx_gene = this->data_expr->raw_data->identifier2idx[gene_probe_ids.at(i)];
        KV.get(gene_probe_ids.at(i), snp_probes);
        
        for(int s=0; s<int(snp_probes.size()); s++){
            idx_snp = this->data_snps->raw_data->identifier2idx[ snp_probes.at(s) ];
            std::cout.flush();            
            is_nan = calculation_wrapper( idx_snp, idx_gene, true_sample_idx, tstat_observed, max_tstat_perm );
            pval_observed = pt(tstat_observed, (double)n_samples-2, 0, 0)*2;
            if(chr_column_name.size()>0){
                snp_indexes->clear();
                snp_chr = snp_idx2chr[idx_snp];
                for(int i=0; i<int(data_snps->identifiers->size()); i++){
                    if( snp_idx2chr[i] == snp_chr )
                        snp_indexes->push_back(i);
                }
            }
            pval_perm=0.0;
            for(int i=0;i<n_perms;i++){
                if( permutation_tstat_across_loci_exceeds( idx_gene, snp_indexes, perm_sample_idx.at(i), tstat_observed, max_tstat_perm ) )
                    pval_perm += 1;
            }
            pval_perm = (pval_perm / n_perms);
            Q = new rqtl( idx_snp, idx_gene, pval_observed, pval_perm, tstat_observed);
            fill_XY( X, Y, idx_snp, idx_gene, true_sample_idx);
            find_mean_three_way(&X, &Y, mu_A, mu_B, mu_C);
            Q->mean_a = mu_A;
            Q->mean_b = mu_B;
            Q->mean_c = mu_C;
            
            this->qtls.push_back(Q);
        }
        report_progress(i+1, gene_probe_ids.size());
    }
    if( chr_column_name.size()>0)
        delete snp_indexes;

}


void QTL_Calculator::calculate_genome_wide( bool genomewide_by_chromosome ){
    // For each probe
    //     if cis_interval==0,
    //          find best pval across all loci using observed locus labels
    //     otherwise
    //          find best pval across all loci within cis-window using observed locus 
    
    // For each of N permutations:
    //     find best pval across all loci using permuted locus labels
    //     report pvalue as position in list of sorted "best permutation" pvals
    // MULTITHREADED
    
    std::stringstream ss;
    ss << "Beginning analysis of " << this->idx.size() << " phenotypes.";
    this->write_to_log(std::string(ss.str()));

    
    // Store the true sample order in true_sample_idx
    // calculate and store permutations of sample order in vector perm_sample_idx.
    // Use same permutation orders for all genes.  
    std::vector< std::vector<int>* > perm_sample_idx;
    std::vector<int>* true_sample_idx = new std::vector<int>();
    int n_samples = (int)this->expr_snps.size();    
    calculate_sample_indexes(true_sample_idx, perm_sample_idx, n_samples);
    
    // create n_threads threads to process data
    std::vector<thread*> threads;
    std::vector<int>* idx_all_snps = new std::vector<int>();
    
    if( this->cis_interval>0 ){
        // any permutation should be performed against ALL SNPs that will be considered, not just those
        // which are relevant to a particular cis-locus. Identify these:
        int half_window = this->cis_interval/2;
        boost::unordered_map<int,int> all_snps_for_perm;
        int idx_probe_id, locus, loc_start, loc_end;
        std::vector<int>* genotype_idx = new std::vector<int>();
        std::string identifier,chromosome;
        for(int n=0; n<int(this->idx.size()); n++){
            idx_probe_id = this->idx.at(n);
            identifier = this->data_expr->ga->identifiers.at(idx_probe_id);
            if( this->data_expr->ga->get_chromosome_and_locus_by_identifier( identifier, chromosome, locus) ){
                loc_start = locus-half_window;
                loc_end = locus+half_window;
                this->data_snps->ga->get_idx_by_genomic_range(chromosome, loc_start, loc_end, *genotype_idx);
                for(int i=0; i<int(genotype_idx->size()); i++){
                    all_snps_for_perm[genotype_idx->at(i)]=1;
                }
            }            
        }
        delete genotype_idx;
        for(boost::unordered_map<int,int >::iterator itr=all_snps_for_perm.begin(); itr != all_snps_for_perm.end(); itr++){
            idx_all_snps->push_back( (*itr).second );
        }
        std::cout << "MESSAGE: Calculating with " << idx_all_snps->size() << " SNPs selected from all CIS windows\n";
    }
    for( int thread_id=0; thread_id<this->n_threads; thread_id++){
        if( this->cis_interval==0 ){
            if(genomewide_by_chromosome){
                // genome-wide one chromosome at a time
                threads.push_back(new thread( boost::bind( &QTL_Calculator::process_by_chr_in_thread, this, thread_id, true_sample_idx, perm_sample_idx) ));
            }
            else{
                // genome-wide all loci
                threads.push_back(new thread( boost::bind( &QTL_Calculator::process_genomewide_in_thread, this, thread_id, true_sample_idx, perm_sample_idx) ));
            }
        }
        else{
            // CIS-only
            threads.push_back(new thread( boost::bind( &QTL_Calculator::process_cis_only_in_thread, this, thread_id, true_sample_idx, perm_sample_idx, idx_all_snps, this->cis_interval) ));
        }
    }
    for( int thread_id=0; thread_id<this->n_threads; thread_id++){
        threads.at(thread_id)->join();
    }
    // clean up
    delete true_sample_idx;
    for( int n=0; n<this->n_perms; n++)
        delete perm_sample_idx.at(n);
    for( int thread_id=0; thread_id<this->n_threads; thread_id++){
        delete threads.at(thread_id);
    }

    delete idx_all_snps;
    this->close_logfile_if_open();
}


void QTL_Calculator::process_all_regressions_in_thread(int thread_id, double max_pvalue_to_report ){
    
    std::vector<int>* true_sample_idx = new std::vector<int>();
    int n_samples = (int)this->expr_snps.size();    
    for( int s=0; s<n_samples; s++){
        true_sample_idx->push_back( expr_snps.at(s)->second );
    }
    
    bool is_nan;
    int idx_probe_id;
    int n_snps = int(this->data_snps->raw_data->identifiers.size());
    std::vector<int>* snp_indexes = new std::vector<int>();
    for(int i=0; i<n_snps; i++)
        snp_indexes->push_back(i);
    
    double m, b, t_stat, pvalue;
    std::string identifier_probe, identifier_snp;
    int current_idx = request_next_index_to_process();
    std::stringstream ss;
    while( current_idx != -1 ){
        ss.clear();
        idx_probe_id = this->idx.at(current_idx);
        identifier_probe = this->ga_expr->identifiers.at( idx_probe_id );        
        for(int snps_g_idx=0; snps_g_idx<n_snps; snps_g_idx++){
            is_nan = regress_direct(true_sample_idx, idx_probe_id, snps_g_idx, m, b, t_stat);
            if(!is_nan){
                pvalue = pt(abs(t_stat), (double)n_samples-2, 0, 0)*2;
                if( pvalue <= max_pvalue_to_report ){
                    {
                        boost::mutex::scoped_lock lock( this->results_mutex );
                        this->qtls.push_back(new rqtl( snps_g_idx, idx_probe_id, pvalue, 1, t_stat ));
                    }
                }
            }
        }
        report_thread_progress( thread_id, current_idx+1, int(this->idx.size()) );
        current_idx = request_next_index_to_process();
    }
    
    delete snp_indexes;
    delete true_sample_idx;
}



void QTL_Calculator::calculate_all_regressions(ClassMinerOptions* cmo_expr, ClassMinerOptions* cmo_snps, double fraction_required, double max_pvalue_to_report){   
    
    std::vector<thread*> threads;
    for( int thread_id=0; thread_id<this->n_threads; thread_id++){
        threads.push_back(new thread( boost::bind( &QTL_Calculator::process_all_regressions_in_thread, this, thread_id, max_pvalue_to_report) ));
    }
    for( int thread_id=0; thread_id<this->n_threads; thread_id++){
        threads.at(thread_id)->join();
    }
    // clean up
    for( int thread_id=0; thread_id<this->n_threads; thread_id++){
        delete threads.at(thread_id);
    }
}


//***************************************************
// Reporting functions
//***************************************************


void QTL_Calculator::reportFileList( std::string fn, std::vector<std::string>& keys, std::vector<std::string>& fn_out ){
    std::ofstream f_out( fn.c_str());
    if( !f_out.is_open() ){
        std::cout << "Unable to open file for writing: " <<  fn << std::endl;
        exit(-1);
    }
    f_out << "key\tfilename\n";
    for(int i=0; i<int(keys.size()); i++){
        f_out << keys.at(i) << "\t" << fn_out.at(i) << "\n";
    }
    f_out.close();
}


void QTL_Calculator::report_progress( int current, int total){
    // Called during QTL processing
    char timebuf[80];
    struct tm* newtime;
    time_t long_time;
    time( &long_time );
    newtime = localtime( &long_time );
    strftime(timebuf, 80, "%H_%M_%S", newtime);
    std::string timeout(timebuf);
    std::cout << "MESSAGE: " << timeout << " Calculations complete on gene " << current << " of " << total << "\n";
    std::cout.flush();   
}

void QTL_Calculator::report_thread_progress( int thread_id, int current, int total ){
    
    // Get a lock on IO to write output
    boost::mutex::scoped_lock lock( this->io_mutex );  
    // Called during multithreaded QTL processing
    char timebuf[80];
    struct tm* newtime;
    time_t long_time;
    time( &long_time );
    newtime = localtime( &long_time );
    strftime(timebuf, 80, "%H_%M_%S", newtime);
    std::string timeout(timebuf);
    std::cout << "MESSAGE: " << timeout << " Complete " << current << " of " << total << " [Thread " << thread_id+1 << "]\n";
    std::cout.flush();   
}


void QTL_Calculator::write_to_log( std::string msg ){
    char timebuf[80];
    struct tm* newtime;
    time_t long_time;
    time( &long_time );
    newtime = localtime( &long_time );
    strftime(timebuf, 80, "%c", newtime);

    std::string timeout(timebuf);

    boost::mutex::scoped_lock lock( this->io_mutex );  

    if( this->fs_log.is_open()){
        fs_log << timeout << "\t" << msg << "\n";
        fs_log.flush();
    }
}

std::string QTL_Calculator::create_output_header(ClassMinerOptions* cmo_expr, ClassMinerOptions* cmo_snps, double max_var,
                                                 double fraction_required){
    std::string id_expr, id_snps;
	char timebuf[80];
	struct tm* newtime;
	time_t long_time;
	time( &long_time );
    newtime = localtime( &long_time );
	strftime(timebuf, 80, "%Y_%m_%d_%H_%M_%S", newtime);
    
    std::stringstream ss;
    ss << "# Source " << cmo_expr->version << "\n";
    ss << "# Date " << timebuf << "\n";
    ss << "# Format long\n";
    ss << "# Phenotype_Data_File " << cmo_expr->file_name_dis << "\n";
    ss << "# Phenotype_Sample_File " << cmo_expr->file_name_sa << "\n";
    ss << "# Phenotype_Genes_File " << cmo_expr->file_name_ga << "\n";
    ss << "# Genotype_Data_File " << cmo_snps->file_name_dis << "\n";
    ss << "# Genotype_Sample_File " << cmo_snps->file_name_sa << "\n";
    ss << "# Genotype_Genes_File " << cmo_snps->file_name_ga << "\n";
    ss << "# cis_interval_0_indicates_genomewide " << this->cis_interval << "\n";
    ss << "# min_observations_per_genotype " << this->get_min_obs_per_genotype() << "\n";
    ss << "# min_var " << max_var << "\n";
    ss << "# n_perms " << this->n_perms << "\n";
    ss << "# test_type regression\n";
    ss << "# fraction_expr_samples_required " << fraction_required << "\n";
    ss << "# experiment_wide_0.05_threshold 1\n";
    ss << "# sample restriction " << cmo_expr->class_a << "\n";
    return ss.str();
}


void QTL_Calculator::reportNetResults(Matrix<double>& results, std::vector<int>& idx_snps, ClassMinerOptions* cmo_expr, ClassMinerOptions* cmo_snps ){
    std::ofstream f_out(cmo_expr->file_name_out.c_str());
    if( !f_out.is_open() ){
        std::cout << "Unable to open file for writing: " << cmo_expr->file_name_out << std::endl;
        exit(-1);
    }
    f_out << create_output_header(cmo_expr, cmo_snps, 0, 0);
    f_out.precision(5);        
    f_out << "IDENTIFIER";
    for(int col=0; col<int( this->idx.size() ); col++){
        f_out << "\t" << this->data_expr->raw_data->identifiers.at( this->idx.at(col) );
    }
    f_out << "\n";
    for(int row=0; row<int( idx_snps.size() ); row++){
        f_out << this->data_snps->raw_data->identifiers.at( idx_snps.at(row) );
        for(int col=0; col<int( this->idx.size() ); col++){
            f_out << "\t" << results.arr[row][col];
        }
        f_out << "\n";
    }
    f_out.close();
}


void QTL_Calculator::spear_result_string(rqtl* q, std::string& spear_str){
    // This isn't a method of rqtl because it needs ga_expr and ga_snps to figure out the symbol and SNP names
    std::stringstream s;
    std::string id_expr( this->ga_expr->identifiers.at(q->idx_probe) );
    std::string id_snps( this->ga_snps->identifiers.at(q->idx_genotype) );
    s << scientific << q->pvalue;
    std::string gene_name_col( this->ga_expr->get_gene_name_column() );
    s << "\t" << this->ga_expr->prop_for_identifier(id_expr, gene_name_col) << "\t";
    s << fixed << id_expr << "\t" << id_snps << "\t" << q->min_permuted_pval << "\t" << q->mean_a << "\t" << q->mean_b << "\t" << q->mean_c << "\t" << q->t_stat;
    spear_str = s.str();
}


void QTL_Calculator::spear_result_string(qtl_by_chrom* q, std::string& spear_str){
    std::stringstream s;    
    std::string id_expr = this->ga_expr->identifiers.at(q->idx_probe); 
    int n_chrom=q->T_obs.size();
    std::string gene_name_col( this->ga_expr->get_gene_name_column() );

    for(int i=0; i<n_chrom; i++){
        s <<  id_expr << "\t" << this->ga_expr->prop_for_identifier(id_expr, gene_name_col) << "\t";
        s << this->ga_snps->identifiers.at(q->snp_indexes.at(i)) << "\t" << i+1 << "\t" << scientific;
        s << q->pvalues_raw.at(i) << "\t" << fixed << q->pvalues_perm.at(i);
        if( i<n_chrom-1 )
            s << "\n";
    }
    spear_str = s.str();
}



void QTL_Calculator::print_all_regressions(ClassMinerOptions* cmo_expr, ClassMinerOptions* cmo_snps, double fraction_required ){
    
    std::ofstream f_out(cmo_expr->file_name_out.c_str());
    if( !f_out.is_open() ){
        std::cout << "Unable to open file for writing:" << cmo_expr->file_name_out << std::endl;
        exit(-1);
    }
    std::cout << "MESSAGE: Writing to " << cmo_expr->file_name_out << "\n";
    f_out << create_output_header(cmo_expr, cmo_snps, 0, fraction_required);
    f_out.precision(5);        
    f_out << "probe\tsnp\tp.value\tchrom\tloc\n";
    rqtl* q;
    std::string chromosome, identifier_snp;
    int locus;
    for(int i=0; i<(int)this->qtls.size(); i++){
        q = this->qtls.at(i);
        identifier_snp = this->ga_snps->identifiers.at( q->idx_genotype );
        this->data_snps->ga->get_chromosome_and_locus_by_identifier( identifier_snp, chromosome, locus);
        f_out << this->ga_expr->identifiers.at( q->idx_probe ) << "\t" << identifier_snp << "\t";
        f_out << q->pvalue << "\t" << chromosome << "\t" << locus << "\n";
    }
    
    f_out.close();
    
}


void QTL_Calculator::print_rqtl( ClassMinerOptions* cmo_expr, ClassMinerOptions* cmo_snps, double max_var, double fraction_required ){
	rqtl* q=NULL;
    std::string id_expr, id_snps, result_str;
	if( cmo_expr->file_name_out.length() == 0 ){
        std::cout << create_output_header(cmo_expr, cmo_snps, max_var, fraction_required);
        std::cout.precision(5);
        std::cout << "#raw.p.val\tsymbol\tprobe.gene\tprobe.snp\tperm.p.val\tMean_A\tMean_B\tMean_C\tt.stat\n";
        std::string gene_name_col = this->ga_expr->get_gene_name_column();
        for(int i=0; i<(int)this->qtls.size(); i++){
            q = this->qtls.at(i);
            if(q->pvalue<=1 & q->pvalue >=0 ){
                this->spear_result_string(q, result_str);
                std::cout << result_str << "\n";
            }
        }
	}
	else{
		std::ofstream f_out(cmo_expr->file_name_out.c_str());
		if( !f_out.is_open() ){
			std::cout << "Unable to open file for writing:" << cmo_expr->file_name_out << std::endl;
			exit(-1);
		}
        f_out << create_output_header(cmo_expr, cmo_snps, max_var, fraction_required);
        f_out.precision(5);        
        f_out << "#raw.p.val\tsymbol\tprobe.gene\tprobe.snp\tperm.p.val\tMean_A\tMean_B\tMean_C\tt.stat\n";
        std::string gene_name_col = this->ga_expr->get_gene_name_column();
        for(int i=0; i<(int)this->qtls.size(); i++){
            q = this->qtls.at(i);
            if(q->pvalue<1 & q->pvalue >=0){
                this->spear_result_string(q, result_str);
                f_out << result_str << "\n";    
            }
        }
		f_out.close();
	}
}


void QTL_Calculator::print_by_chromosome( ClassMinerOptions* cmo_expr, ClassMinerOptions* cmo_snps, double fraction_required ){
	int n_chrom=0; 
    if( qtls_by_chrom.size()>0 ){
        n_chrom = qtls_by_chrom.at(0)->T_obs.size();
    }
    std::string id_expr, result_str;
    if( cmo_expr->file_name_out.length()==0){
        std::cout << create_output_header(cmo_expr, cmo_snps, 0, fraction_required);
        std::cout.precision(5);        
        std::cout << "probe.gene\tsymbol\tprobe.snp\tchromosome\tpval.raw\tpval.perm\n";
        std::string gene_name_col = this->ga_expr->get_gene_name_column();
        for(int i=0; i<(int)qtls_by_chrom.size(); i++){
            this->spear_result_string( this->qtls_by_chrom.at(i), result_str);
            std::cout << result_str << "\n";
        }
	}
	else{
		std::ofstream f_out(cmo_expr->file_name_out.c_str());
		if( !f_out.is_open() ){
			std::cout << "Unable to open file for writing:" << cmo_expr->file_name_out << std::endl;
			exit(-1);
		}
        f_out << create_output_header(cmo_expr, cmo_snps, 0, fraction_required);
        f_out.precision(5);        
        f_out << "probe.gene\tsymbol\tprobe.snp\tchromosome\tpval.raw\tpval.perm\n";
        std::string gene_name_col = this->ga_expr->get_gene_name_column();
        for(int i=0; i<(int)qtls_by_chrom.size(); i++){
            this->spear_result_string( this->qtls_by_chrom.at(i), result_str);
            f_out << result_str << "\n";
        }
		f_out.close();
	}
}


