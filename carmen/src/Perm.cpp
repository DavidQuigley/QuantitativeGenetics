#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <time.h>
using namespace std;
#include <boost/tokenizer.hpp>
#include <boost/random.hpp>
#include <boost/lexical_cast.hpp>
using namespace boost;
#include "DataStructures.h"
#include "Attributes.h"
#include "ClassMinerOptions.h"
#include "ParseOptions.h"
#include "Rawdata.h"
#include "chi2.h"
#include "Discretize.h"
#include "Dataset.h"
#include "Perm.h"

Perm::Perm(){
}


bool comp_perm(const perm_result *a, const perm_result *b){ 
	return a->t_val < b->t_val; 
}

void Perm::limit_ids_by_NA(std::vector<int> &idx, ClassifierDataset* data, double fraction_required){
    // restrict idx by elminating any rows where 
    // either class has fewer than two samples present

    std::vector<int> keep;
    int n_in_a, n_in_b, a, b;
    if( data->raw_data->data->has_missing_values == false )
        return;
    int n_req_a = (int)((double)data->a_idx->size() * fraction_required );
    int n_req_b = (int)((double)data->b_idx->size() * fraction_required );
    for(int i=0; i<(int)idx.size(); i++){
        n_in_a = 0;
        n_in_b = 0;
        if(data->raw_data->data->row_has_missing_value(i)){
            for(a=0;a<(int)data->a_idx->size();a++){
                if(!data->raw_data->data->is_missing(i,a)){
                    n_in_a++;
                }
            }
            if(data->b_is_target){
                for(b=0;b<(int)data->b_idx->size();b++){
                    if(!data->raw_data->data->is_missing(i,b)){
                        n_in_b++;
                    }
                }
                if( n_in_a>=n_req_a && n_in_b>=n_req_b )
                    keep.push_back(i);
            }
            else{
                if( n_in_a>=n_req_a )
                    keep.push_back(i);
            }
        }
        else{
            keep.push_back(i);
        }
    }
    idx.clear();
	for(int i=0; i<(int)keep.size(); i++)
		idx.push_back( keep.at(i) );
}


void Perm::load_settings( std::string file_name, double max_p_value ){
	this->max_p_value = max_p_value;
	perm_result* r;
	double p_value;
	std::string line, cmd;
	std::vector<std::string>* arr = new std::vector<std::string>();
	char space[1], tab[1];
	space[0] = ' ';
	tab[0] = '\t';
	if( file_name.size()>0 ){	
		std::ifstream f_cv(file_name.c_str());	
		if( !f_cv.is_open() ){ 
			throw std::string("Unable to open difference file for reading.");
		}
		
		while( !f_cv.eof() ){
			getline(f_cv, line);
			if(line.length()==0)
				break;

			if( line.at(0) == '#' ){
				stringhasher::split(line, arr, space);
				cmd = arr->at(1);
				
				if( cmd.compare("Data_File")==0){
					this->file_name_dis = arr->at(2);
				}
				if( cmd.compare("Sample_File")==0){
					this->file_name_sa = arr->at(2);
				}
				if( cmd.compare("Genes_File")==0){
					this->file_name_ga = arr->at(2);
				}
				if( cmd.compare("Class_A")==0){
					this->class_a = arr->at(2);
					this->n_a = boost::lexical_cast<int>( arr->at(3) );
				}
				if( cmd.compare("Class_B")==0){
					this->class_b = arr->at(2);
					this->n_b = boost::lexical_cast<int>( arr->at(3));
				}
				if( cmd.compare("N_Permutations")==0){
					this->n_perm= boost::lexical_cast<int>( arr->at(2) );
				}
				if( cmd.compare("Mean_Trim")==0){
					this->mean_trim = boost::lexical_cast<int>( arr->at(2) );
				}
			}
			else{
				stringhasher::split(line, arr, tab);
				p_value = boost::lexical_cast<double>(arr->at(7));
				if(p_value <= this->max_p_value ){
					r = new perm_result();
					r->identifier = arr->at(0);
					r->name = arr->at(1);
					r->mu_a = boost::lexical_cast<float>( arr->at(2) );
					r->var_a = boost::lexical_cast<float>( arr->at(3) );
					r->mu_b = boost::lexical_cast<float>( arr->at(4) );
					r->var_b = boost::lexical_cast<float>( arr->at(5) );
					r->t_val = boost::lexical_cast<float>( arr->at(6) );
					r->p_value = (float)p_value;
					this->results.push_back(r);
				}
			}
		}
	}
	std::sort(this->results.begin(), this->results.end(), comp_perm);
}



void Perm::write_as_ruleset( ClassMinerOptions* cmo, Rawdata* rd, vector<int>* a_idx, vector<int>* b_idx ){
	char timebuf[80];
	struct tm* newtime;
	time_t long_time;
	time( &long_time );
	newtime = localtime( &long_time );
	strftime(timebuf, 80, "%Y_%m_%d_%H_%M_%S", newtime);
	perm_result* r;
	int idx, n_a, n_b;
	n_a = this->n_a;
	n_b = this->n_b;
	double sda, sdb, cut_a, cut_b;
	double chi2, A, B, conf;
	if(cmo->file_name_out.length()>0){
		std::ofstream f_out(cmo->file_name_out.c_str());
		if( !f_out.is_open() ){ 
			throw std::string( "Unable to open file for writing: " ) + cmo->file_name_out ;
		}
		f_out << "# Source CARMEN 1.0\n";
		f_out << "# Mine_Type Difference\n";
		f_out << "# Date " << timebuf << "\n";
		f_out << "# Data_File " << cmo->file_name_dis << "\n";
		f_out << "# Sample_File " << cmo->file_name_sa << "\n";
		f_out << "# Genes_File " << cmo->file_name_ga << "\n";
		f_out << "# Class_A " << this->class_a << " " <<  this->n_a << " members.\n";
		f_out << "# Class_B " << this->class_b << " " << this->n_b << " members.\n";
		f_out << "# N_Permutations " << this->n_perm << "\n";
		f_out << "# Mean_Trim " << this->mean_trim << "\n";
		f_out << "# Max_Depth 0\n";
		f_out << "# Discretization none\n";
		f_out << "# Generated " << results.size() << " rules.\n";
		f_out << "# chi2	%conf	supA	%supA	supB	%supB	imp	Rule	Rule_ID	Rule_Dis\n";
		//1.77969e-012	100	0	0	27	71.0526	1	ECM2_is_1	206101_at_is_1	206101_at>8.1441
		float* row;
		int x;
		for(int i=0; i< (int)results.size(); i++){
			r = results.at(i);
			sda = sqrt(r->var_a);
			sdb = sqrt(r->var_b);
			idx = rd->identifier2idx[ r->identifier ];
			row = rd->data->arr[idx];
			
			if( r->mu_a > r->mu_b ){
				A = 0;
				B = 0;
				cut_a = r->mu_a - sda;
				cut_b = r->mu_b + sdb;
				for(x=0;x<n_a;x++){ if( row[ a_idx->at(x) ] >= cut_a ){ ++A; } }
				for(x=0;x<n_b;x++){ if( row[ b_idx->at(x) ] >= cut_a ){ ++B; } }
				chi2 = pchisq(calc_chi2(A, n_a-A, B, n_b-B), 1);
				A>B ? conf = A/(A+B) : conf = B/(A+B);
				f_out << chi2 << '\t' << conf << "\t" << A << "\t" << A/n_a << "\t" << B << "\t" << B/n_b << "\t" << conf << '\t';
				f_out << r->name << "_is_1\t" << r->identifier << "_is_1\t" << r->identifier << ">" << cut_a << "\n";
				A = 0;
				B = 0;
				for(x=0;x<n_a;x++){ if( row[ a_idx->at(x) ] <= cut_b ){ ++A; } }
				for(x=0;x<n_b;x++){ if( row[ b_idx->at(x) ] <= cut_b ){ ++B; } }
				chi2 = pchisq(calc_chi2(A, n_a-A, B, n_b-B), 1);
				A>B ? conf = A/(A+B) : conf = B/(A+B);
				f_out << chi2 << "\t" << conf << "\t" << A << "\t" << A/n_a << "\t" << B << "\t" << B/n_b << "\t" << conf << "\t";
				f_out << r->name << "_is_-1\t" << r->identifier << "_is_-1\t" << r->identifier << "<" << cut_b << "\n";
			}
			else{
				A = 0;
				B = 0;
				cut_a = r->mu_a + sda;
				cut_b = r->mu_b - sdb;
				for(x=0;x<n_a;x++){ if( row[ a_idx->at(x) ] <= cut_a ){ ++A; } }
				for(x=0;x<n_b;x++){ if( row[ b_idx->at(x) ] <= cut_a ){ ++B; } }
				chi2 = pchisq(calc_chi2(A, n_a-A, B, n_b-B), 1);
				A>B ? conf = A/(A+B) : conf = B/(A+B);
				f_out << chi2 << "\t" << conf << "\t" << A << "\t" << A/n_a << "\t" << B << "\t" << B/n_b << "\t" << conf << "\t";
				f_out << r->name << "_is_-1\t" << r->identifier << "_is_-1\t" << r->identifier << ">" << cut_a << "\n";
				A = 0;
				B = 0;
				for(x=0;x<n_a;x++){ if( row[ a_idx->at(x) ] >= cut_b ){ ++A; } }
				for(x=0;x<n_b;x++){ if( row[ b_idx->at(x) ] >= cut_b ){ ++B; } }
				chi2 = pchisq(calc_chi2(A, n_a-A, B, n_b-B), 1);
				A>B ? conf = A/(A+B) : conf = B/(A+B);
				f_out << chi2 << "\t" << conf << "\t" << A << "\t" << A/n_a << "\t" << B << "\t" << B/n_b << "\t" << conf << "\t";
				f_out << r->name << "_is_1\t" << r->identifier << "_is_1\t" << r->identifier << "<" << cut_b << "\n";
			}
		}
	}
}


inline int Perm::round(float x){
	return int(x > 0.0 ? x + 0.5 : x - 0.5);
}


void Perm::write( ClassMinerOptions* cmo ){
	
	char timebuf[80];
	struct tm* newtime;
	time_t long_time;
	time( &long_time );
	newtime = localtime( &long_time );
	strftime(timebuf, 80, "%Y_%m_%d_%H_%M_%S", newtime);

	if(cmo->file_name_out.length()>0){
		std::ofstream f_out(cmo->file_name_out.c_str());
		if( !f_out.is_open() ){ 
			throw std::string( "Unable to open file for writing: " ) + cmo->file_name_out ;
		}
		f_out << "# Source CARMEN 1.0\n";
		f_out << "# Date " << timebuf << "\n";
		f_out << "# Data_File " << cmo->file_name_dis << "\n";
		f_out << "# Sample_File " << cmo->file_name_sa << "\n";
		f_out << "# Genes_File " << cmo->file_name_ga << "\n";
		f_out << "# Class_A " << cmo->class_a << " " <<  this->n_a << " members.\n";
		f_out << "# Class_B " << cmo->class_b << " " << this->n_b << " members.\n";
		f_out << "# N_Permutations " << this->n_perm << "\n";
		f_out << "# Mean_Trim " << this->mean_trim << "\n";
		f_out << "# Max_p_value " << this->max_p_value << "\n";
		f_out << "# ID\tNAME\tMEAN_A\tVAR_A\tMEAN_B\tVAR_B\tT_STAT\tP-VAL\tDIFF\n";
		for(int i=0; i<(int)results.size(); i++){
			f_out << results.at(i)->identifier << "\t";
			f_out << results.at(i)->name << "\t";
			f_out << results.at(i)->mu_a << "\t" << results.at(i)->var_a << "\t";
			f_out << results.at(i)->mu_b << "\t" << results.at(i)->var_b << "\t";
			f_out << results.at(i)->t_val << "\t" << results.at(i)->p_value << "\t";
			f_out << results.at(i)->mu_a - results.at(i)->mu_b << "\n";
		}
		f_out.close();
	}
	else{
		std::cout << "# Source CARMEN 1.0\n";
		std::cout << "# Date " << timebuf << "\n";
		std::cout << "# Data_File " << cmo->file_name_dis << "\n";
		std::cout << "# Sample_File " << cmo->file_name_sa << "\n";
		std::cout << "# Genes_File " << cmo->file_name_ga << "\n";
		std::cout << "# Class_A " << cmo->class_a << " " <<  this->n_a << " members.\n";
		std::cout << "# Class_B " << cmo->class_b << " " << this->n_b << " members.\n";
		std::cout << "# N_Permutations " << this->n_perm << "\n";
		std::cout << "# Mean_Trim " << this->mean_trim << "\n";
		std::cout << "# ID\tNAME\tMEAN_A\tVAR_A\tMEAN_B\tVAR_B\tT_STAT\tP-VAL\tDIFF\n";
		for(int i=0; i<(int)results.size(); i++){
			std::cout << results.at(i)->identifier << "\t";
			std::cout << results.at(i)->name << "\t";
			std::cout << results.at(i)->mu_a << "\t" << results.at(i)->var_a << "\t";
			std::cout << results.at(i)->mu_b << "\t" << results.at(i)->var_b << "\t";
			std::cout << results.at(i)->t_val << "\t" << results.at(i)->p_value << "\t";
			std::cout << results.at(i)->mu_a - results.at(i)->mu_b << "\n";
		}
	}
}


int Perm::rnd(int aRange){
	return int((double )rand()/(double(RAND_MAX) + 1)*aRange);
}


void Perm::shuffle(vector<int>* V){
	// wrote this hack because despite working for hours I could not get 
	// std::random_shuffle to seed with a new value
	int a, b, tmp;
	int len = (int)V->size();
	for(int i=0;i<10000;i++){
		a = rnd(len);
		b = rnd(len);
		tmp = V->at(b);
		V->at(b) = V->at(a);
		V->at(a) = tmp;
	}
}

void Perm::permute( vector<int>* idx_a, vector<int>* idx_b, vector<int>* perm_a, vector<int>* perm_b){
	int n_shorter, i;
	int n_a = (int)idx_a->size();
	int n_b = (int)idx_b->size();
	perm_a->clear();
	perm_b->clear();
	//vector<int>* combined = new vector<int>();
	//vector<int>* pool = new vector<int>();
	//vector<int>* shorter = new vector<int>();
	bool a_is_longer=true;
	if(n_a>n_b){
		for(i=0; i<n_a; i++){ perm_a->push_back(idx_a->at(i)); }
		for(i=0; i<n_b; i++){ perm_b->push_back(idx_b->at(i)); }
		n_shorter = n_b;
	}
	else{
		a_is_longer=false;
		for(i=0; i<n_b; i++){ perm_a->push_back(idx_b->at(i)); }
		for(i=0; i<n_a; i++){ perm_b->push_back(idx_a->at(i)); }
		n_shorter = n_a;
	}
	shuffle(perm_a); // shuffle in place, does not exchange classes
	shuffle(perm_b); // shuffle in place, does not exchange classes
	int tmp;
	for(i=0; i<n_shorter; i++){
		// choose items from the shorter to swap with the longer.  Since 
		// both were shuffled, this is random selection of a uniform 
		// distribution of swaps
		if( rnd(2)==1 ){	
			tmp = perm_b->at(i);
			perm_b->at(i) = perm_a->at(i);
			perm_a->at(i) = tmp;
		}
	}
}


void Perm::trimmed_stats(Matrix<float>* data, int row, vector<int>* cols, int trim_percent, float& mean, float& var){
	// calculates mean and variance
	mean=0.0;
	var=0.0;
	int n_cols = (int)cols->size();
	int trim_start = round( n_cols * ((float)trim_percent/(float)100) );
	int trim_end, n_trim;	
	vector<float> acc;
	
	if(data->row_has_missing_value(row)){
		for(int c=0; c<n_cols; c++ ){
			if( data->mask[row][ cols->at(c) ] == 0){
				acc.push_back( data->arr[row][ cols->at(c) ] );
			}
		}
		trim_end = (int)acc.size()-trim_start;
	}
	else{
		for(int c=0; c<n_cols; c++ ){
			acc.push_back(data->arr[row][ cols->at(c) ] );
		}
		trim_end = n_cols-trim_start;
	}
	n_trim = trim_end-trim_start;
	std::sort(acc.begin(), acc.end());
	if(n_trim==0)
		return;

	float sum=0.0;
	for(int c=trim_start; c<trim_end; c++ ){
		sum += acc.at(c);
	}
	mean = sum / (float)n_trim;
	var = this->find_variance(mean, acc, trim_start, trim_end);
}


float Perm::find_t_stat( float mu_a, float var_a, int n_a, float mu_b, float var_b, int n_b){
	// calculate trimmed t_statistic
	float diff_means = mu_a - mu_b;
	float standard_error = (float)sqrt( (var_a/n_a) + (var_b/n_b) );
	return diff_means / standard_error ;
}


float Perm::find_variance(float mean, const vector<float>& acc, int idx_start, int idx_end){
	// calculate trimmed sample variance
	int N = idx_end - idx_start;
	float sum_of_diff_squared=0.0;
	float diff;
	for(int i=idx_start; i<idx_end; i++){
		diff = acc.at(i) - mean;
		sum_of_diff_squared += diff * diff;
	}
	return sum_of_diff_squared / (N - 1) ;
}



float Perm::find_p_value( float r, vector< vector<int>*> * perms_a, vector< vector<int>*> * perms_b, Matrix<float>* data, int row, int trim_percent){
	int n_perms = (int)perms_a->size();
	float* vals = new float[n_perms];
	int n_a = (int)perms_a->at(0)->size();
	int n_b = (int)perms_b->at(0)->size();
	float mean_a, mean_b, var_a, var_b;
	for(int p=0; p<n_perms; p++){
		trimmed_stats(data, row, perms_a->at(p), trim_percent, mean_a, var_a);
		trimmed_stats(data, row, perms_b->at(p), trim_percent, mean_b, var_b);
		vals[p] = this->find_t_stat(mean_a, var_a, n_a, mean_b, var_b, n_b);
	}
	sort(vals, vals+n_perms);
	int p;
	float pval=0;
	float p_inc = (float)1 / n_perms;
	if( r>=0 ){
		for( p=n_perms-1; p>=0; p--){
			if( r >= vals[p] )
				break;
			pval+=p_inc;
		}
	}
	else{
		for( p=0; p<n_perms; p++){
			if( r <= vals[p] )
				break;
			pval+=p_inc;
		}
	}
	delete [] vals;
	return pval;
}


void Perm::permutations( Rawdata* rd, std::vector<int> &idx, vector<int>* idx_a, vector<int>* idx_b, Attributes* ga, bool verbose){
	int n_identifiers = (int)idx.size();
	int i, gene_idx;
	vector<float>* mean_a = new vector<float>(n_identifiers);
	vector<float>* mean_b = new vector<float>(n_identifiers);
	vector<float>* variance_a = new vector<float>(n_identifiers);
	vector<float>* variance_b = new vector<float>(n_identifiers);

	this->n_a = (int)idx_a->size();
	this->n_b = (int)idx_b->size();
	
	//vector<int>* perm_a = new vector<int>();
	//vector<int>* perm_b = new vector<int>();
	
	if( verbose ){
		std::cout << "MESSAGE: Calculating " << this->n_perm << " permutations on " << n_identifiers << " genes.\n";
		std::cout.flush();
	}
	vector< vector<int>* >* perms_a = new vector< vector<int>* >();
	vector< vector<int>* >* perms_b = new vector< vector<int>* >();
	
	float ma, mb, va, vb;

	for(i=0; i<n_identifiers; i++){
        gene_idx = idx.at(i);
		trimmed_stats(rd->data, gene_idx, idx_a, this->mean_trim, ma, va);
		trimmed_stats(rd->data, gene_idx, idx_b, this->mean_trim, mb, vb);
		mean_a->at(i)=ma;
		mean_b->at(i)=mb;
		variance_a->at(i)=va;
		variance_b->at(i)=vb;
	}
	
	for(int p=0; p<this->n_perm; p++){
		vector<int>* perm_a = new vector<int>();
		vector<int>* perm_b = new vector<int>();
		this->permute( idx_a, idx_b, perm_a, perm_b);
		perms_a->push_back(perm_a);
		perms_b->push_back(perm_b);
	}

	float r;
	float p_value;
	perm_result * result;
    if( this->n_perm == 1 )
        this->max_p_value = 1.0;
    std::string symbol_column = ga->get_gene_name_column();
	for(i=0; i<n_identifiers; i++){ 
        gene_idx = idx.at(i);
		if(verbose && i%1000==0){
			std::cout << "MESSAGE: Completed item " << i << " of " << n_identifiers << "...\n";
			std::cout.flush();
		}
		r = this->find_t_stat(mean_a->at(i), variance_a->at(i), this->n_a, mean_b->at(i), variance_b->at(i), this->n_b);
		p_value = 1;
        if( this->n_perm > 1 )
            p_value = find_p_value( r, perms_a, perms_b, rd->data, i, this->mean_trim);
		
		if( p_value <= (float)this->max_p_value ){
			result = new perm_result();
			result->identifier = rd->identifiers.at( gene_idx );
			result->name = ga->prop_for_identifier( rd->identifiers.at( gene_idx ), symbol_column);
			result->mu_a = mean_a->at(i);  
			result->var_a = variance_a->at(i);
			result->mu_b = mean_b->at(i);  
			result->var_b = variance_b->at(i);
			result->t_val = this->find_t_stat(result->mu_a, result->var_a, this->n_a, result->mu_b, result->var_b, this->n_b);
			result->p_value = p_value;
			result->idx = i;
			this->results.push_back(result);
		}
	}
	
	if( verbose ){ 
		std::cout << "MESSAGE: Done.\n";
		std::cout.flush();
	}
	delete mean_a;
	delete mean_b;
	delete variance_a;
	delete variance_b;
	std::sort(this->results.begin(), this->results.end(), comp_perm);
}
