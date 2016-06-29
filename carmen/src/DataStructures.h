#include <iostream>
#include <fstream>
#include <vector>
#include <time.h>
#include <string>
#include <algorithm>
#include "boost/unordered_map.hpp"
#include "boost/tokenizer.hpp"
#include "boost/lexical_cast.hpp"

#ifndef HASH_TYPEDEFS
typedef boost::unordered_map<std::string, int> HASH_S_I;
typedef boost::unordered_map<std::string, std::string> HASH_S_S;
typedef boost::unordered_map<int, std::vector<std::string>* > HASH_I_VECTOR_STR;
#define HASH_TYPEDEFS 1
#endif

class QTL{
public:
	QTL(std::string locus_id, std::string locus_name, std::string probe_id, std::string probe_name, double pval);
	std::string locus_id, locus_name, probe_id, probe_name;
	double pval;
};

class stringhasher{
public:
	static const int EQ=0;
	static const int LT=-1;
	static const int GT=1;
	static void split(std::string line, std::vector<std::string>* arr, const char* token){
		arr->clear();
		typedef boost::tokenizer< boost::char_separator<char> > tokenizer;
		boost::char_separator<char> sep(token);
		tokenizer tok(line, sep);
		for(tokenizer::iterator iter=tok.begin(); iter!=tok.end(); iter++){
			arr->push_back( (*iter) );
		}
	}
	static void split_comma(std::string line, std::vector<std::string>* arr){
			arr->clear();
			typedef boost::tokenizer< boost::char_separator<char> > tokenizer;
			boost::char_separator<char> sep(",");
			tokenizer tok(line, sep);
			for(tokenizer::iterator iter=tok.begin(); iter!=tok.end(); iter++){
				arr->push_back( (*iter) );
			}
		}
	static void join(std::string& result, std::vector<std::string>* arr, std::string token){
		result = "";
		int N =(int)arr->size();
		for(int i=0; i<N; i++){
			result += arr->at(i);
			if(i<N-1)
				result += token;
		}
	}

	static void get_idval(std::string rulepart, std::string& identifier, std::string& value){
		int loc_is = (int)rulepart.find("_is_");
		int len = (int)rulepart.length();
		if(loc_is>-1){
			identifier = rulepart.substr(0,loc_is);
			value = rulepart.substr(loc_is+4, len-loc_is);
		}
		else
			throw new std::string("Rule part does not contain _is_, cannot parse it.");
	}

	static void get_id_dir_val(std::string str, std::string& identifier, int& dir, float& value){
		int loc_LT = (int)str.find("<");
		int loc_GT = (int)str.find(">");
		int loc_EQ = (int)str.find("=");
		int len = (int)str.length();
		if( loc_EQ != -1 ){
			dir = stringhasher::EQ;
			identifier = str.substr(0,loc_EQ);
			value = boost::lexical_cast<float>( str.substr(loc_EQ+1, len-loc_EQ ) );
		}
		else if( loc_LT != -1 ){
			dir = stringhasher::LT;
			identifier = str.substr(0,loc_LT);
			value = boost::lexical_cast<float>( str.substr(loc_LT+1, len-loc_LT ) );
		}
		else if( loc_GT != -1 ){
			dir = stringhasher::GT;
			identifier = str.substr(0,loc_GT);
			value = boost::lexical_cast<float>( str.substr(loc_GT+1, len-loc_GT ) );
		}
		else
			throw new std::string("Rule part does not contain <, >, or =, cannot parse it.");
	}

	static void parse_limit(std::string limit, std::string& str_attrib, std::string& str_value, bool& is_equal){
		// limit should have the form foo=bar, foo!bar, or foo*bar. set str_attrib to foo and str_value to bar.
		std::vector<std::string> arr;
		if(limit.find(std::string("!")) != std::string::npos){
			is_equal=false;
			split(limit, &arr, std::string("!").c_str() );
			str_attrib = arr.at(0);
			str_value = arr.at(1);
		}
		else if(limit.find(std::string("*")) != std::string::npos){ // required for compatibility with bash shell
			is_equal=false;
			split(limit, &arr, std::string("*").c_str() );
			str_attrib = arr.at(0);
			str_value = arr.at(1);
		}
		else{
			is_equal=true;
			if( limit.find(std::string("=") ) == std::string::npos )
				throw std::string("Incorrect format for class limit: must include = or ! symbol.");
			split(limit, &arr, std::string("=").c_str() );
			str_attrib = arr.at(0);
			str_value = arr.at(1);
		}
	}
};

// RNG is a random number generator that uses mt19937 parameterization of the twister
// Used by shuffle_by_rows()
struct RNG : std::unary_function<int, int> {
    boost::mt19937 &_state;
    unsigned operator()(unsigned i) {
    	boost::uniform_int<> rng(0, i - 1);
        return rng(_state);
    }
    RNG(boost::mt19937 &state) : _state(state) {}
};

template <typename T> class Matrix {
public:
	Matrix();
	Matrix(int rows, int cols, bool has_missing_values=false, T val=0);
	Matrix(std::vector<T>* V);
	~Matrix();

	T ** arr;
	bool** mask;
	bool has_missing_values;
	bool* row_has_missing_values;

	void add(T val);
	void add(Matrix<T>* C);
	void cbind(Matrix<T>* C, Matrix<T>* result);
	void clone(Matrix<T>* C);
	int cols();
	void colSums(std::vector<T>* result);
	void div(T val);
	void div(Matrix<T> * B);
	void exponent();
    void eliminate_infrequent_rowvalues(std::vector<int>& idx, int min_examples );
	bool is_missing(int row, int col);
	void mult(T val);
	void mult(Matrix<T> * B);
	void sub(Matrix<T> * B);
	T mult_sum(Matrix<T> * B);
	T min(int& min_x, int& min_y);
	T min(std::vector<int>* min_x, std::vector<int>* min_y);
	void mm(Matrix<T> * B, Matrix<T> * C);
	void power(double p);
	void print(std::string spacer);
	void resize(int rows, int cols);
    void restrict_by_rows(std::vector<int> idx_keep );
    void restrict_by_cols(std::vector<int> idx_keep );
	bool row_has_missing_value(int row);
	int rows();
	void rowSums(std::vector<T>* result);
	void set_all_values_present();
	void set_missing_value(int row, int col);
	void slice_row(int row, Matrix<T>* X);
	void shuffle_by_rows( std::vector<int> valid_columns);
	T sum();
	void transpose(Matrix<T> * X);
	void where_LT(Matrix<T> * X, T val);
	void where_EQ(Matrix<T> * X, T val);
	void where_GT(Matrix<T> * X, T val);
    void write(std::string file_name_out, std::string spacer);
    void write_with_names(std::string file_name_out, std::string spacer, std::vector<std::string> col_names, std::vector<std::string> row_names);
	
private:
    boost::mt19937 state;
    RNG* generator;
	int r,c;
	void clean_up();
	void set_dimensions(int n_rows, int n_cols, bool has_missing_values, T val);
	int round(T x);
};


class Ant{
public:
	Ant(float conf, std::vector<int>* idx, int sup_A, int sup_B, float imp, double chi2);
	~Ant();
	std::string rule_as_str();
	float conf;
	int sup_A;
	int sup_B;
	float imp;
	double chi2;
	double t_stat;
	std::vector<int> * idx;
};


template <typename T> Matrix<T>::Matrix(){
	this->r = 0;
	this->c = 0;
	this->arr  = NULL;
	this->mask = NULL;
	this->has_missing_values = false;
	this->row_has_missing_values = NULL;
	this->generator = new RNG(this->state);
	this->generator->_state.seed((int)time(0));
}


template <typename T> Matrix<T>::Matrix(std::vector<T>* V){
	this->r = 1;
	int n_cols = (int)V->size();
	this->c = n_cols;
	int row=0;
	this->arr = new T*[ 1 ];
	T* new_t = new T[ this->c ];
	this->arr[row] = new_t;
	int col;
	
	for(col=0;col<n_cols;col++){
		this->arr[row][col]=V->at(col);
	}
	this->has_missing_values=false;
	this->mask = NULL;
	this->row_has_missing_values = NULL;
	this->generator = new RNG(this->state);
	this->generator->_state.seed((int)time(0));
}


template <typename T> Matrix<T>::Matrix(int n_rows, int n_cols, bool has_missing_values, T val){
	// Given the number of rows and cols, initialize 2-D array of ints
	this->set_dimensions(n_rows, n_cols, has_missing_values, val);
	this->generator = new RNG(this->state);
	this->generator->_state.seed((int)time(0));
}

template <typename T> void Matrix<T>::set_dimensions(int n_rows, int n_cols, bool has_missing_values, T val){
	this->r = n_rows;
	this->c = n_cols;
	this->arr = new T*[n_rows];
	int row,col;
	this->row_has_missing_values = new bool[n_rows];
	for(row=0;row<n_rows;row++){
		this->row_has_missing_values[row]=false;
		T* new_t = new T[n_cols];
		this->arr[row] = new_t;
		for(col=0;col<n_cols;col++){
			this->arr[row][col]=val;
		}
	}
	this->has_missing_values = has_missing_values;
	if(has_missing_values){
		this->mask = new bool*[n_rows];
		for(row=0;row<n_rows;row++){
			this->mask[row] = new bool[n_cols];
			for(col=0;col<n_cols;col++){
				this->mask[row][col]=false;
			}
		}
	}
	else
		this->mask = NULL;
}

template <typename T> void Matrix<T>::transpose(Matrix<T>* M){
	if( (this->r==0 && this->c==0) ) // empty 
		return;
	if( (this->r==1 && this->c==0) ) // empty vector case
		return;

	M->resize(this->c, this->r);
	int row, col;
	for(row=0; row<this->r; row++){
		for(col=0; col<this->c; col++){
			M->arr[col][row] = this->arr[row][col];
		}
	}
	if( this->has_missing_values ){
		for(row=0; row<this->r; row++){
			for(col=0; col<this->c; col++){
				if( this->is_missing(row, col) )
					M->set_missing_value(col, row);
			}
		}
	}
}

template <typename T> void Matrix<T>::clone(Matrix<T>* C){
    // copy this INTO C
	int n_rows = this->r;
	int n_cols = this->c;
	C->resize(n_rows, n_cols);
	int row, col;
	for(row=0; row<n_rows; row++){
		C->row_has_missing_values[row]=false;
		for(col=0; col<n_cols; col++){
			C->arr[row][col] = this->arr[row][col];
		}
	}
	if(this->has_missing_values){
		C->has_missing_values = true;
		for(row=0;row<n_rows;row++){
			for(col=0;col<n_cols;col++){
				if( this->mask[row][col] ){
					C->set_missing_value(row, col);
				}
			}
		}
	}
}


template <typename T> void Matrix<T>::cbind(Matrix<T>* C, Matrix<T>* R){
	// mimics R cbind(); matrix C is appended to the right side of this matrix, stored in R
	if( this->rows() != C->rows() )
		throw std::string( "Cannot bind matrices with different number of rows");
	R->resize(this->rows(), ( this->cols() + C->cols() ) );
	int row, col;
	int n_rows = this->r;
	int n_cols_A = this->c;
	int n_cols_C = C->cols();
	for(row=0; row<n_rows; row++){
		R->row_has_missing_values[row]=false;
		for(col=0; col<n_cols_A; col++){
			R->arr[row][col] = this->arr[row][col];
		}
		for(col=0; col<n_cols_C; col++){
			R->arr[row][ n_cols_A + col ] = C->arr[row][col];
		}
	}
	if(this->has_missing_values || C->has_missing_values){
		R->has_missing_values = true;
		for(row=0;row<n_rows;row++){
			for(col=0;col<n_cols_A;col++){
				if( this->mask[row][col] ){
					R->set_missing_value(row, col);
				}
			}
			for(col=0; col<n_cols_C; col++){
				if( C->is_missing(row,col) ){
					R->set_missing_value(row, n_cols_A + col );
				}
			}
		}
	}
}


template <typename T> void Matrix<T>::add(T val){
	// element-wise matrix addition of a scalar
	// Adding to missing values is harmless as their value has no meaning
	int row, col;
	for(row=0; row<this->r; row++){
		for(col=0; col<this->c; col++){
			(this->arr[row][col]) += val;
		}
	}
}


template <typename T> void Matrix<T>::add(Matrix<T>* B){
	// element-wise matrix addition of another matrix
	// Defined so NA values propigate (if Bi,j == NA, then THISi,j is set to NA)
	int row, col;
	if(B==NULL)
		throw std::string( "Must send initialized Matrix pointer into add");
	if( B->rows() != this->r && B->cols() != this->c )
		throw std::string( "Cannot add matrices with different dimensions");

	if( B->has_missing_values ){
		this->has_missing_values = true;
		for(row=0; row<this->r; row++){
			for(col=0; col<this->c; col++){
				if( B->is_missing(row, col) )
					this->set_missing_value(row, col);
				else
					this->arr[row][col] += B->arr[row][col];
			}
		}
	}
	else{
		for(row=0; row<this->r; row++){
			for(col=0; col<this->c; col++){
				this->arr[row][col] += B->arr[row][col];
			}
		}
	}
}


template <typename T> void Matrix<T>::div(T val){
	// Element-wise divide of this matrix by val.
	// Dividing missing values by val is harmless as the value still has no meaning
	int row, col;
	if( val==0 )
		throw "Divide by zero error in call to div";
	for(row=0; row<this->r; row++){
		for(col=0; col<this->c; col++){
			(this->arr[row][col]) /= val;
		}
	}
}

template <typename T> void Matrix<T>::div(Matrix<T>* B){
	// element-wise matrix division
	// Must check for and propigate missing values.
	int row, col;
	if(B==NULL)
		throw std::string( "Must send initialized Matrix pointer into div");
	if( B->rows() != this->r && B->cols() != this->c )
		throw std::string( "Cannot mult matrices with different dimensions");
	if( B->has_missing_values ){
		this->has_missing_values = true;
		for(row=0; row<this->r; row++){
			for(col=0; col<this->c; col++){
				if( B->is_missing(row, col) )
					this->set_missing_value(row, col);
				else{
					if( B->arr[row][col] == 0 )
						this->set_missing_value(row, col);
					else
						this->arr[row][col] /= B->arr[row][col];
				}
			}
		}
	}
	else{
		for(row=0; row<this->r; row++){
			for(col=0; col<this->c; col++){
				if( B->arr[row][col] == 0 )
					this->set_missing_value(row, col);
				else
					this->arr[row][col] /= B->arr[row][col];
			}
		}
	}
}

template <typename T> void Matrix<T>::exponent(){
	// Element-wise call exp( Mij )
	// exp(missing values) is harmless as the value still has no meaning
	int row, col;
	double e;
	for(row=0; row<this->r; row++){
		for(col=0; col<this->c; col++){
			e = exp((double)this->arr[row][col]);
			this->arr[row][col] = (T)e;
		}
	}
}


template <typename T> void Matrix<T>::eliminate_infrequent_rowvalues(std::vector<int>& idx, int min_examples ){
    // for each row, count the number of times we see each element
    // limit the analysis to samples with indexes in idx
    // if any element appears fewer than min_examples times, mark those appearances as 
    // missing data. Primary use case is for genotype data (integers 0,1,2) where we are 
    // attempting to eliminate outliers
    for( int i=0; i<int(idx.size()); i++){
        if(idx.at(i)>= this->cols()){
            throw std::string( "idx value is greater than largest column index");
        }
    }
    boost::unordered_map<int, int> seen, nuke;
    int key=-1, value=-1;
    int r=-1,c=-1;
    for(r=0; r< this->rows(); r++){
        seen.clear();
        nuke.clear();
        if( this->row_has_missing_value(r) ){
            for(int i=0; i<int( idx.size() ); i++ ){
                c = idx.at(i);
                if( !this->is_missing(r,c) ){
                    value=this->arr[r][c];
                    if(seen.find(value) == seen.end() ){
                        seen[value]=1;
                    }
                    else{
                        seen[ value ] = seen[value]+1;
                    }    
                }
            }        
        }
        else{
            for(int i=0; i<int( idx.size() ); i++ ){
                c = idx.at(i);
                value=this->arr[r][c];
                if(seen.find(value) == seen.end() ){
                    seen[ value ]=1;
                }
                else{
                    seen[ value ] = seen[value]+1;
                }
            }
        }
        // for each key, if value < min_examples mark for NA
        for( boost::unordered_map<int, int>::iterator itr = seen.begin(); itr != seen.end(); itr++ ){
            key = (*itr).first;
            value = seen[key];
            if( value < min_examples ){
                nuke[key]=1;
            }
        }
        if( nuke.size()>0 ){
            // run through row marking NA if 
            for(int i=0; i<int( idx.size() ); i++ ){
                c = idx.at(i);
                if( seen.find( this->arr[r][c] ) != seen.end() ){
                    this->set_missing_value(r,c);
                }
            }
        }
    }
}




template <typename T> void Matrix<T>::mult(T val){
	// element-wise matrix multiplication by a scalar
	// Multiplying missing values by val is harmless as the value still has no meaning
	int row, col;
	for(row=0; row<this->r; row++){
		for(col=0; col<this->c; col++){
			(this->arr[row][col]) *= val;
		}
	}
}


template <typename T> void Matrix<T>::mult(Matrix<T>* B){
	// element-wise matrix multiplication
	// Must check for and propigate missing values.
	int row, col;
	if(B==NULL)
		throw std::string( "Must send initialized Matrix pointer into mult");
	if( B->rows() != this->r && B->cols() != this->c )
		throw std::string( "Cannot mult matrices with different dimensions");
	if( B->has_missing_values ){
		this->has_missing_values = true;
		for(row=0; row<this->r; row++){
			for(col=0; col<this->c; col++){
				if( B->is_missing(row, col) )
					this->set_missing_value(row, col);
				else
					this->arr[row][col] *= B->arr[row][col];
			}
		}
	}
	else{
		for(row=0; row<this->r; row++){
			for(col=0; col<this->c; col++){
				this->arr[row][col] *= B->arr[row][col];
			}
		}
	}
}

template <typename T> void Matrix<T>::power(double p){
	// raise all elements in this to power p
	int row, col;
	for(row=0; row<this->r; row++){
		for(col=0; col<this->c; col++){
			this->arr[row][col] = (T) pow( this->arr[row][col], p);
		}
	}
}

template <typename T> void Matrix<T>::sub(Matrix<T>* B){
	// element-wise matrix subtraction
	// If B has missing values, propigate them
	int row, col;
	if(B==NULL)
		throw std::string( "Must send initialized Matrix pointer into sub");
	if( B->rows() != this->r && B->cols() != this->c )
		throw std::string( "Cannot subtract matrices with different dimensions");
	if( B->has_missing_values ){
		this->has_missing_values = true;
		for(row=0; row<this->r; row++){
			for(col=0; col<this->c; col++){
				if( B->is_missing(row, col) )
					this->set_missing_value(row, col);
				else
					this->arr[row][col] -= B->arr[row][col];
			}
		}
	}
	else{
		for(row=0; row<this->r; row++){
			for(col=0; col<this->c; col++){
				this->arr[row][col] -= B->arr[row][col];
			}
		}
	}
}


template <typename T> T Matrix<T>::mult_sum(Matrix<T>* B){
	// simultanious element-wise matrix multiplication and summation
	// If either matrix has a missing value, do not add it to sum.
	if(B==NULL)
		throw std::string( "Must send initialized Matrix pointer into mult_sum");
	if( B->rows() != this->r && B->cols() != this->c )
		throw std::string( "Cannot mult_sum matrices with different dimensions");
	T sum=0;
	int row, col;

	if( this->has_missing_values || B->has_missing_values ){
		for(row=0; row<this->r; row++){
			for(col=0; col<this->c; col++){
				if( !this->is_missing(row, col) && !B->is_missing(row, col) )
					sum += this->arr[row][col] * B->arr[row][col];
			}
		}
	}
	else{
		for(row=0; row<this->r; row++){
			for(col=0; col<this->c; col++){
				sum += this->arr[row][col] * B->arr[row][col];
			}
		}
	}
	return sum;
}

template <typename T> void Matrix<T>::mm(Matrix<T>* B, Matrix<T>* C){
	// naive implementation of matrix multiplication
	int m = this->r;
	int n = this->c;
	int p = B->cols();
	if(B==NULL)
		throw std::string( "Must send initialized Matrix pointer into mm at first parameter");
	if(C==NULL)
		throw std::string( "Must send initialized Matrix pointer into mm at second parameter");
	if( C->rows() != m || C->cols() != p)
		C->resize(m,p);
	if (B->rows() != n)
		throw "Matrix inner dimensions must agree.";
	
	T s=0;
	for (int j = p; --j >= 0; ) {
		for (int i = m; --i >= 0; ) {
			s = 0;
			for (int k = n; --k >= 0; ) {
				s += this->arr[i][k] * B->arr[k][j];
			}
			C->arr[i][j] = s;
		}
	}
}


template <typename T> T Matrix<T>::sum(){
	int row, col;
	T s=0;
	if( this->has_missing_values ){
		for(row=0; row<this->r; row++){
			for(col=0; col<this->c; col++){
				if( !this->is_missing(row, col) ){
					s+=this->arr[row][col];
				}
			}
		}
	}
	else{
		for(row=0; row<this->r; row++){
			for(col=0; col<this->c; col++){
				s+=this->arr[row][col];
			}
		}
	}
	return s;
}


template <typename T> T Matrix<T>::min(int& min_x, int& min_y){
	int row, col;
	T m;
	if(this->r==0 || this->c==0 )
		throw std::string( "Error: cannot call min on matrix of size 0,0");
	m = this->arr[0][0];
	min_x = 0; 
	min_y = 0;
	if( this->has_missing_values ){
		for(row=0; row<this->r; row++){
			for(col=0; col<this->c; col++){
				if( !this->is_missing(row, col) ){
					if( this->arr[row][col] < m ){
						m = this->arr[row][col];
						min_x = row;
						min_y = col;
					}
				}
			}
		}
	}
	else{
		for(row=0; row<this->r; row++){
			for(col=0; col<this->c; col++){
				if( this->arr[row][col] < m ){
					m = this->arr[row][col];
					min_x = row;
					min_y = col;
				}
			}
		}
	}
	return m;
}



template <typename T> T Matrix<T>::min(std::vector<int>* min_x, std::vector<int>* min_y){
	int row, col;
	T m;
	if(this->r==0 || this->c==0 )
		throw std::string( "Error: cannot call min on matrix of size 0,0");
	m = this->arr[0][0];
	min_x->clear(); 
	min_y->clear();
	if( this->has_missing_values ){
		for(row=0; row<this->r; row++){
			for(col=0; col<this->c; col++){
				if( !this->is_missing(row,col) ){
					if( this->arr[row][col] < m ){
						m = this->arr[row][col];
						min_x->clear();
						min_y->clear();
						min_x->push_back(row);
						min_y->push_back(col);
					}
					else if( this->arr[row][col] == m ){
						min_x->push_back(row);
						min_y->push_back(col);
					}
				}
			}
		}
	}
	else{
		for(row=0; row<this->r; row++){
			for(col=0; col<this->c; col++){
				if( this->arr[row][col] < m ){
					m = this->arr[row][col];
					min_x->clear();
					min_y->clear();
					min_x->push_back(row);
					min_y->push_back(col);
				}
				else if( this->arr[row][col] == m ){
					min_x->push_back(row);
					min_y->push_back(col);
				}
			}
		}
	}
	return m;
}


template <typename T> void Matrix<T>::where_LT(Matrix<T>* M, T val){
	int row, col;
	if(M==NULL)
		throw "Must send initialized Matrix pointer into where_LT";
	if(M->cols() != this->c || M->rows() != this->r )
		M->resize(this->r, this->c);
	if( this->has_missing_values ){
		for(row=0; row<this->r; row++){
			for(col=0; col<this->c; col++){
				if( this->is_missing(row, col) )
					M->set_missing_value(row, col);
				else
					this->arr[row][col] < val ? M->arr[row][col]=1 : M->arr[row][col]=0;
			}
		}
	}
	else{
		for(row=0; row<this->r; row++){
			for(col=0; col<this->c; col++){
				this->arr[row][col] < val ? M->arr[row][col]=1 : M->arr[row][col]=0;
			}
		}
	}
}


template <typename T> void Matrix<T>::where_EQ(Matrix<T>* M, T val){
	int row, col;
	if(M==NULL)
		throw std::string( "Must send initialized Matrix pointer into where_EQ");
	if(M->cols() != this->c || M->rows() != this->r )
		M->resize(this->r, this->c);
	if( this->has_missing_values ){
		for(row=0; row<this->r; row++){
			for(col=0; col<this->c; col++){
				if( this->is_missing(row, col) )
					M->set_missing_value(row, col);
				else
					this->arr[row][col] == val ? M->arr[row][col]=1 : M->arr[row][col]=0;
			}
		}
	}
	else{
		for(row=0; row<this->r; row++){
			for(col=0; col<this->c; col++){
				this->arr[row][col] == val ? M->arr[row][col]=1 : M->arr[row][col]=0;
			}
		}
	}
}


template <typename T> void Matrix<T>::where_GT(Matrix<T>* M, T val){
	int row, col;
	if(M==NULL)
		throw std::string( "Must send initialized Matrix pointer into where_GT");
	if(M->cols() != this->c || M->rows() != this->r )
		M->resize(this->r, this->c);
	if( this->has_missing_values ){
		for(row=0; row<this->r; row++){
			for(col=0; col<this->c; col++){
				if( this->is_missing(row, col) )
					M->set_missing_value(row, col);
				else
					this->arr[row][col] > val ? M->arr[row][col]=1 : M->arr[row][col]=0;
			}
		}
	}
	else{
		for(row=0; row<this->r; row++){
			for(col=0; col<this->c; col++){
				this->arr[row][col] > val ? M->arr[row][col]=1 : M->arr[row][col]=0;
			}
		}
	}
}



template <typename T> void Matrix<T>::colSums(std::vector<T>* results){
	int c,r;
	T sum;
	if(results==NULL)
		throw std::string( "Must send initialized vector pointer into colSums");
	results->clear();
	if( this->has_missing_values ){
		for(c=0; c<this->c; c++){
			sum=0;
			for(r=0; r<this->r; r++){
				if( !this->is_missing(r,c))
					sum += this->arr[r][c];
			}
			results->push_back(sum);
		}
	}
	else{
		for(c=0; c<this->c; c++){
			sum=0;
			for(r=0; r<this->r; r++){
				sum += this->arr[r][c];
			}
			results->push_back(sum);
		}
	}
}


template <typename T> void Matrix<T>::rowSums(std::vector<T>* results){
	int c,r;
	T sum;
	if(results==NULL)
		throw std::string( "Must send initialized vector pointer into colSums");
	results->clear();
	if( this->has_missing_values ){
		for(r=0; r<this->r; r++){
			sum=0;
			for(c=0; c<this->c; c++){
				if( !this->is_missing(r,c))
					sum += this->arr[r][c];
			}
			results->push_back(sum);
		}
	}
	else{
		for(r=0; r<this->r; r++){
			sum=0;
			for(c=0; c<this->c; c++){
				sum += this->arr[r][c];
			}
			results->push_back(sum);
		}
	}
}



template <typename T> int Matrix<T>::round(T x){
	return int(x > 0.0 ? x + 0.5 : x - 0.5);
}


template <typename T> bool Matrix<T>::is_missing(int row, int col){
	if( this->has_missing_values )
		return this->mask[row][col];
	else
		return false;
}


template <typename T> bool Matrix<T>::row_has_missing_value(int row){
	return this->row_has_missing_values[row];
}


template <typename T> void Matrix<T>::set_missing_value(int row, int col){
	if( row>=this->r || row<0 )
		throw std::string("Bounds error in set_missing_value: row");
	if( col>=this->c || col<0 )
		throw std::string("Bounds error in set_missing_value: col");
	if( this->mask == NULL ){
		this->has_missing_values = true;
		this->mask = new bool*[ this->r ];
		for(int row=0; row<this->r; row++){
			this->mask[row] = new bool[this->c];
			for(int col=0; col<this->c; col++){
				this->mask[row][col]=false;
			}
		}
	}
	this->mask[row][col]=true;
	this->row_has_missing_values[row] = true;
	this->has_missing_values = true;
}


template <typename T> void Matrix<T>::set_all_values_present(){
	// reset all mask (missing) values to read "present"
	if( this->mask==NULL ){
		return;
	}
	for(int row=0; row<this->r; row++){
		this->row_has_missing_values[row] = false;
		for(int col=0; col<this->c; col++){
			this->mask[row][col]=false;
		}
	}
	this->has_missing_values = false;
}


// Allow matrix to be trimmed to arbitrary sub-selection of rows.
template <typename T> void Matrix<T>::restrict_by_rows(std::vector<int> idx_keep ){
	// check that all values in idx_keep are valid indexes for Matrix, and no duplicates
    boost::unordered_map<int, int> seen;
    for(int i=0; i<int(idx_keep.size()); i++){
        if( idx_keep.at(i) >= this->r ){
            throw std::string("index in restrict_by_rows larger than number of rows");
        }
        if(seen.find( idx_keep.at(i) ) != seen.end() ){
            throw std::string("Duplicate index in restrict_by_rows");
        }
    }
    // copy existing matrix into tmp. resize this to new number of rows.
    // copy requested rows into this. Clean up temp.
    Matrix<T>* tmp = new Matrix<T>();
    this->clone(tmp);
    this->clean_up();
	this->set_dimensions(int(idx_keep.size()), this->c, false, 0);
    int old_row_idx;
    for(int new_row_idx=0; new_row_idx<int(idx_keep.size()); new_row_idx++){
        old_row_idx = idx_keep.at(new_row_idx);
        for( int col=0; col<tmp->c; col++){
            this->arr[new_row_idx][col] = tmp->arr[old_row_idx][col];
            if( tmp->is_missing(old_row_idx, col)){
                this->set_missing_value(new_row_idx, col);
            }
        }
    }
    tmp->clean_up();
}


// Allow matrix to be trimmed to arbitrary sub-selection of columns.
template <typename T> void Matrix<T>::restrict_by_cols(std::vector<int> idx_keep ){
	// check that all values in idx_keep are valid indexes for Matrix, and no duplicates
    boost::unordered_map<int, int> seen;
    for(int i=0; i<int(idx_keep.size()); i++){
        if( idx_keep.at(i) >= this->c ){
            throw std::string("index in restrict_by_cols larger than number of columns");
        }
        if(seen.find( idx_keep.at(i) ) != seen.end() ){
            throw std::string("Duplicate index in restrict_by_cols");
        }
    }
    // copy existing matrix into tmp. resize this to new number of rows.
    // copy requested rows into this. Clean up temp.
    Matrix<T>* tmp = new Matrix<T>();
    this->clone(tmp);
    this->clean_up();
	this->set_dimensions(this->r, int(idx_keep.size()), false, 0);
    int old_col_idx;
    for(int new_col_idx=0; new_col_idx<int(idx_keep.size()); new_col_idx++){
        old_col_idx = idx_keep.at(new_col_idx);
        for( int row=0; row<tmp->r; row++){
            this->arr[row][new_col_idx] = tmp->arr[row][old_col_idx];
            if( tmp->is_missing(row, old_col_idx)){
                this->set_missing_value(row, new_col_idx);
            }
        }
    }
    tmp->clean_up();
}


template <typename T> void Matrix<T>::resize(int n_rows, int n_cols){ 
	if(n_rows==0 && n_cols==0 && this->r==0 && this->c==0)
		return;
	if(n_rows<1 || n_cols<1){
		throw std::string( "Values for resize must be greater than zero") ;
	}
	this->clean_up();
	this->set_dimensions(n_rows, n_cols, this->has_missing_values, 0);
}


template <typename T> void Matrix<T>::slice_row(int row, Matrix<T>* X){
	if(row>=this->r || row<0)
		throw std::string( "Row requested in slice_row is out of range." );
	if(X==NULL)
		throw std::string( "Must send initialized Matrix pointer into slice_row" );
	if(X->rows() != 1 || X->cols() != this->c)
		X->resize(1, this->c);
	if( this->has_missing_values ){
		for(int col=0; col<this->c; col++){
			if( this->is_missing(row, col) )
				X->set_missing_value(0, col);
			else
				X->arr[0][col] = this->arr[row][col];
			}
	}
	else{
		for(int col=0; col<this->c; col++){
			X->arr[0][col] = this->arr[row][col];
		}
	}
}


template <typename T> int Matrix<T>::rows(){ 
	return this->r; 
}


template <typename T> int Matrix<T>::cols(){ 
	return this->c; 
}



template <typename T> void Matrix<T>::shuffle_by_rows( std::vector<int> valid_columns ){
	// For each row, it shuffles the values in that row.
	// valid_columns is required because in some cases we will only want to shuffle a subset of
	// all columns. In this case, we must exclude invalid columns from the shuffle, or they can
	// distort the distribution of values when we later consult a subset of rows (e.g. for
	// permutation analysis in Spearman code
	int n_valid_col = int(valid_columns.size());
	for(int i=0; i<n_valid_col; i++){
		if( valid_columns.at(i) >= this->c || valid_columns.at(i) < 0 ){
			std::stringstream ss;
			ss << "element " << valid_columns.at(i) << " out of bounds in valid_columns parameter for shuffle_rows";
			throw ss.str();
		}
	}
    std::vector<int> row_order(n_valid_col);
    std::vector<T> tmp(n_valid_col);
    std::vector<bool> tmp_mask(n_valid_col);
    int i,j;
    for(j=0; j<n_valid_col; j++){
        row_order.at(j)=valid_columns.at(j);
    }
    for(i=0; i<this->r; i++){
        random_shuffle( row_order.begin(), row_order.end(), *this->generator );
        for(j=0; j<n_valid_col; j++){
            tmp[j] = this->arr[i][row_order.at(j)];
        }
        for(j=0; j<n_valid_col; j++){
            this->arr[i][valid_columns[j]] = tmp[j]; // bugfix 12/10/10 was this->arr[i][j] = tmp[j];
        }
        if( this->row_has_missing_value(i) ){
            for(j=0; j<n_valid_col; j++){
                tmp_mask[j] = this->mask[i][row_order.at(j)];
            }
            for(j=0; j<n_valid_col; j++){
                this->mask[i][j] = tmp_mask[j];
            }
        }
    }
}


template <typename T> void Matrix<T>::print(std::string spacer){
	int row, col;
	for(row=0;row<this->r;row++){
		for(col=0;col<this->c;col++){
			if( this->has_missing_values && this->mask[row][col])
				std::cout << "NA";
			else
				std::cout << this->arr[row][col];
			if( col<this->c-1 )
				std::cout << spacer;
		}
		std::cout << std::endl;
	}
}


template <typename T> void Matrix<T>::write(std::string file_name_out, std::string spacer){
	std::ofstream f_out(file_name_out.c_str());
	if( !f_out.is_open() ){ 
		throw std::string( "Unable to open requested file for writing." );
	}
	int row, col;
	for(row=0;row<this->r;row++){
		for(col=0;col<this->c;col++){
			if( this->has_missing_values && this->mask[row][col])
				f_out << "NA";
			else
				f_out << this->arr[row][col];
			if(col<this->c-1)
				f_out << spacer;
		}
		f_out << std::endl;
	}
}


template <typename T> void Matrix<T>::write_with_names(std::string file_name_out, std::string spacer, std::vector<std::string> col_names, std::vector<std::string> row_names){
	std::ofstream f_out(file_name_out.c_str());
	if( !f_out.is_open() ){ 
		throw std::string( "Unable to open requested file for writing." );
	}
    if( this->r != (int)row_names.size() ){
        throw std::string( "Length of row_names does not match matrix row length." );
    }
    if( this->c != (int)col_names.size() ){
        throw std::string( "Length of col_names does not match matrix column length." );
    }
	f_out << "IDENTIFIER";
    for(int i=0; i<this->c; i++)
        f_out << spacer << col_names.at(i);
    
    f_out << "\n";
    int row, col;
	for(row=0;row<this->r;row++){
		f_out << row_names.at(row) << spacer;
        for(col=0;col<this->c;col++){
			if( this->has_missing_values && this->mask[row][col])
				f_out << "NA";
			else
				f_out << this->arr[row][col];
			if(col<this->c-1)
				f_out << spacer;
		}
		f_out << std::endl;
	}
}


template <typename T> void Matrix<T>::clean_up(){
	// Correctly free 2-D array of T.  
	for(int i=0; i<this->r; i++){
		delete [] this->arr[i];
	}
	delete [] this->arr;
	if( this->mask != NULL ){
		for(int i=0; i<this->r; i++){
			delete [] this->mask[i] ;
		}
		delete [] this->mask;
	}
	if( this->row_has_missing_values != NULL )
		delete [] this->row_has_missing_values;
}

template <typename T> Matrix<T>::~Matrix(){
	this->clean_up();
	delete this->generator;
}

