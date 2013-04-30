#include <math.h>
#include <vector>
#include <string>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/random.hpp>
using namespace boost;
#include "DataStructures.h"
#include "Rawdata.h"
#include "Attributes.h"
#include "chi2.h"
#include "Discretize.h"

float Discretize::find_median(std::vector<float>* C){
	std::sort(C->begin(), C->end());
	int len = (int)C->size();
	if( len==0 )
        return 0;
    if( len % 2 ==0 ){
		float a = C->at( (len-1)/2 );
		float b = C->at( ((len-1)/2)+1 );
		return (a + b) / (float)2;
	}
	else
		return( C->at( (len-1) / 2 ) );
}

float Discretize::find_mean(std::vector<float>* C){
	int N = (int)C->size();
	if(N==0)
		return 0;
	float sum=0;
    for(int i=0; i<N; i++){
        sum += C->at(i);
    }
	return sum/N;
}


float Discretize::find_var(float mean, std::vector<float>* C){
	int N = (int)C->size();
	float sum_of_diff_squared=0.0;
	float diff;
	if(N==0)
		return 0;
	for(int i=0; i<N; i++){
		diff = C->at(i) - mean;
		sum_of_diff_squared += diff * diff;
	}
	return (float)( sum_of_diff_squared / (N-1) );
}


float Discretize::find_stdev(float mean, std::vector<float>* C){
	return (float)( sqrt( find_var(mean, C ) ));
}


float Discretize::find_mad(float median, std::vector<float>* C){
	int N = (int)C->size();
	for(int i=0; i<N; i++){
		C->at(i) = abs(C->at(i) - median);
	}
	return find_median(C);
}

int Discretize::disc(float f, float dn, float up){
	if( f <= dn )
        return -1;
	else if( f >= up )
        return 1; 
	return 0;
}


void Discretize::discretize_SD_Samples( float value_a, float value_b, std::vector<int>& idx, Rawdata* rd){
	// discretize_SD discretizes by taking the mean/SD across samples for a given gene
	// this function discretizes by taking the mean/SD across genes for a given sample
	// it is intended for use with CGH data, where this approach is more meaningful.
	std::vector<float> C;
	int r, c, i;
	float mid_up, mid_dn, mean=0, stdev=0;
	int n_idx = (int)idx.size();
	for( i=0; i<n_idx; i++ ){
		C.clear();
		c = idx.at(i);
        //float sum=0;
		for(r=0; r<rd->data->rows(); r++){
			if( !rd->data->is_missing(r,c) ){
                C.push_back( rd->data->arr[r][c] );
			}
		}
		mean = find_mean(&C);
		stdev = find_stdev(mean, &C);	
		mid_dn = mean - (stdev * value_a);
		mid_up = mean + (stdev * value_b);
		for(r=0; r<rd->data->rows(); r++){
			if( rd->data->is_missing(r, c) ){
				this->dis->arr[r][c]=0;
			}
			else{
				this->dis->arr[r][c] = disc( rd->data->arr[r][c], mid_dn, mid_up);
			}
		}
	}
	// Since we're discretizing by columns, the cutoff_up/down vectors are not meaningful
	for(r=0; r<rd->data->rows(); r++){
		this->cutoff_up.push_back(0);
		this->cutoff_dn.push_back(0);
	}
}


void Discretize::discretize_SD(float value_a, float value_b, std::vector<int>& idx, Rawdata* rd){
	std::vector<float>* C = new std::vector<float>();
	int r, c, i;
	float mid_up, mid_dn, mean=0, stdev=0;
	if( value_a < 0 )
		value_a = value_a * -1;
	if( value_b < 0 )
		value_b = value_b * -1;
	int n_idx = (int)idx.size();
	for(r=0; r<rd->data->rows(); r++){
		C->clear();
		if( rd->data->row_has_missing_value(r) ){
			for( i=0; i<n_idx; i++ ){
				c = idx.at(i);
				if( ! rd->data->is_missing(r,c) ){
					C->push_back( rd->data->arr[r][c] );
				}
			}
		}
		else{
			for( i=0; i<n_idx; i++ ){
				C->push_back( rd->data->arr[r][ idx.at(i) ] );
			}
		}	
		mean = find_mean(C);
		stdev = find_stdev(mean, C);	
		
        if( value_a==0 ){
            mid_up = mean + (stdev * value_b);
            mid_dn = mid_up; // special case: low vs. not low
        }
        else if( value_b==0 ){
            mid_dn = mean - (stdev * value_a);
            mid_up = mid_dn; // special case: high vs. not high
        }
        else{
            mid_dn = mean - (stdev * value_a);
		    mid_up = mean + (stdev * value_b);
        }
		if( rd->data->row_has_missing_value(r) ){
			for( i=0; i<n_idx; i++ ){
				c = idx.at(i);
				if( rd->data->is_missing(r, c) )
					this->dis->arr[r][c]=0;
				else
					this->dis->arr[r][c] = disc( rd->data->arr[r][c], mid_dn, mid_up);
			}
		}
		else{
			for( i=0; i<n_idx; i++ ){
				c = idx.at(i);
				this->dis->arr[r][c] = disc( rd->data->arr[r][c], mid_dn, mid_up);
			}
		}
		this->cutoff_up.push_back(mid_up);
		this->cutoff_dn.push_back(mid_dn);
	}
	delete C;
}


void Discretize::discretize_MAD(float value_a, float value_b, std::vector<int>& idx, Rawdata* rd){
	std::vector<float>* C = new std::vector<float>();
	int r, c, i;
	if( value_a < 0 )
		value_a = value_a * -1;
	if( value_b < 0 )
		value_b = value_b * -1;
	float mid_up, mid_dn, median=0, mad=0;
	int n_idx = (int)idx.size();
	for(r=0; r<rd->data->rows(); r++){
		C->clear();
		if( rd->data->row_has_missing_value(r) ){
			for( i=0; i<n_idx; i++ ){
				c = idx.at(i);
				if( ! rd->data->is_missing(r,c) )
					C->push_back( rd->data->arr[r][c] );
			}
		}
		else{
			for( i=0; i<n_idx; i++ )
				C->push_back( rd->data->arr[r][ idx.at(i) ] );
		}	
		median = find_median(C);
		mad = find_mad( median, C);
		mid_dn = median - (mad * value_a);
		mid_up = median + (mad * value_b);
			
		if( rd->data->row_has_missing_value(r) ){
			for( i=0; i<n_idx; i++ ){
				c = idx.at(i);
				if( rd->data->is_missing(r, c) )
					this->dis->arr[r][c]=0;
				else
					this->dis->arr[r][c] = disc( rd->data->arr[r][c], mid_dn, mid_up);
			}
		}
		else{
			for( i=0; i<n_idx; i++ ){
				c = idx.at(i);
				this->dis->arr[r][c] = disc( rd->data->arr[r][c], mid_dn, mid_up);
			}
		}
		this->cutoff_up.push_back(mid_up);
		this->cutoff_dn.push_back(mid_dn);
	}
	delete C;
}


void Discretize::discretize_ABS(float value_a, float value_b, std::vector<int> &idx, Rawdata* rd){
	int r, c, i;
	int n_idx = (int)idx.size();
	for(r=0; r<rd->data->rows(); r++){
		if( rd->data->row_has_missing_value(r) ){
			for( i=0; i<n_idx; i++ ){
				c = idx.at(i);
				if( rd->data->is_missing(r,c) )
					this->dis->arr[r][c]=0;
				else		
					this->dis->arr[r][c] = disc( rd->data->arr[r][c], value_a, value_b);
 			}
		}
		else{
			for( i=0; i<n_idx; i++ ){
				c = idx.at(i);
				this->dis->arr[r][c] = disc( rd->data->arr[r][c], value_a, value_b);
			}
		}
		this->cutoff_up.push_back(value_b);
		this->cutoff_dn.push_back(value_a);
	}
}

void Discretize::discretize_PER(float value_a, float value_b, std::vector<int> &idx, Rawdata* rd){
	std::vector<float>* C = new std::vector<float>();
	int r, c=-1, i, idx_a, idx_b;
	float max_a, min_b;
	int n_idx = (int)idx.size();
	std::vector<float> cp;
	for(r=0; r<rd->data->rows(); r++){
		C->clear();
		if( rd->data->row_has_missing_value(r) ){
			for( i=0; i<n_idx; i++ ){
				c = idx.at(i);
				if( ! rd->data->is_missing(r,c) )
					C->push_back( rd->data->arr[r][c] );
			}
		}
		else{
			for( i=0; i<n_idx; i++ )
				C->push_back( rd->data->arr[r][ idx.at(i) ] );
		}

		cp.clear();
		cp.insert(cp.begin(), C->begin(), C->end());
		std::sort(cp.begin(), cp.end());
		idx_a = (int)floor(value_a * (double)C->size());
		idx_b = (int)floor(value_b * (double)C->size());
		if( idx_a==0 || idx_b >= (int)cp.size() ){
			for( i=0; i<n_idx; i++ )	// pathological cases, almost all or all missing
				this->dis->arr[r][c]=0;
			this->cutoff_dn.push_back(0);
			this->cutoff_up.push_back(0);
			continue;
		}
		max_a = cp.at(idx_a-1);
		min_b = cp.at(idx_b);
		if(max_a==min_b)
			max_a -= (float)0.00001;
		if( rd->data->row_has_missing_value(r) ){
			for( i=0; i<n_idx; i++ ){
				c = idx.at(i);
				if( rd->data->is_missing(r,c) )
					this->dis->arr[r][c]=0;
				else
					this->dis->arr[r][c] = disc( rd->data->arr[r][c], max_a, min_b);
			}
		}
		else{
			for( i=0; i<n_idx; i++ ){
				c = idx.at(i);
				this->dis->arr[r][c] = disc( rd->data->arr[r][c], max_a, min_b);
			}
		}
		this->cutoff_dn.push_back(max_a);
		this->cutoff_up.push_back(min_b);
	}
	delete C;
}


void Discretize::discretize_NONE(std::vector<int> &idx, Rawdata* rd){
	int r, c, i;
	int n_idx = (int)idx.size();
	for(r=0; r<rd->data->rows(); r++){
		if( rd->data->row_has_missing_value(r) ){
			for( i=0; i<n_idx; i++ ){
				c = idx.at(i);
				if( rd->data->is_missing(r,c) )
					this->dis->arr[r][c]=0;
				else		
					this->dis->arr[r][c] = (int)rd->data->arr[r][c];
			}
		}
		else{
			for( i=0; i<n_idx; i++ ){
				c = idx.at(i);
				this->dis->arr[r][c] = (int)rd->data->arr[r][c];
			}
		}
		this->cutoff_up.push_back(0);
		this->cutoff_dn.push_back(0);
	}
}

Discretize::Discretize(Rawdata* rd, std::vector<int>& idx, std::string method_str, float value_a, float value_b){
	this->method = method_str;
    this->mult_lower = value_a;
    this->mult_upper = value_b;
	if((int)idx.size() > rd->data->cols() )
		throw std::string("Incorrect idx specified; too many columns.");
	for(int i=0; i<(int)idx.size(); i++)
		if( idx.at(i) >= rd->data->cols() )
			throw std::string("Incorrect idx specified; column out of bounds.");
	this->dis = new Matrix<int>(rd->data->rows(), rd->data->cols(), rd->data->has_missing_values);
	if( method_str.compare("MAD") == 0 || method_str.compare("mad") == 0)
		discretize_MAD( value_a, value_b, idx, rd );
	else if ( method_str.compare("SD") == 0 || method_str.compare("sd") == 0)
		discretize_SD( value_a, value_b, idx, rd );
	else if (method_str.compare("ABS") == 0 || method_str.compare("abs") == 0 )
		discretize_ABS( value_a, value_b, idx, rd );
	else if (method_str.compare("PER") == 0 || method_str.compare("per") == 0 )
		discretize_PER( value_a, value_b, idx, rd );
	else if (method_str.compare("NONE") == 0 || method_str.compare("none") == 0 )
		discretize_NONE(idx, rd);
	else if( method_str.compare("SD_SAMPLES")==0 || method_str.compare("sd_samples")==0 )
		discretize_SD_Samples(value_a, value_b, idx, rd);
	else
		throw std::string("Discretization method must be none, sd, mad, per, sd_samples, or abs.");
}

Discretize::~Discretize(){
	delete this->dis;
}
