#include <iostream>
#include <stdio.h>
#include <string>
#include <fstream>
#include <vector>
#include <algorithm>
#include <math.h>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/random.hpp>
using namespace boost;
#include <time.h>
using namespace std;
using namespace boost;
#include "ClassMinerOptions.h"
#include "DataStructures.h"
#include "Rawdata.h"
#include "Attributes.h"
#include "Discretize.h"
#include "Dataset.h"
#include "graph.h"
#include "ParseOptions.h"
#include "ClassMiner.h"
#include "chi2.h"

#include "Perm.h"

// Author: David Quigley, UCSF Cancer Research Institute.
// Email:  dquigley@cc.ucsf.edu
// Date:   2007-2009

bool comp_desc(const Ant *a, const Ant *b){
	if( a->chi2 == b->chi2 )
		return a->t_stat > b->t_stat;
	else
		return a->chi2 < b->chi2;
}

void ClassMiner::print_ants_stdout(std::vector<Ant*>* ants, ClassifierDataset * data, ClassMinerOptions* cmo){

	std::vector<int*>* translate= data->translate;
	Matrix<int>* dis = data->features;

	std::vector<Ant*>::iterator iter = ants->begin();
	std::vector<int> * rows;
	int idx_in_dis;
	int i,c;
	bool* feature = new bool[ dis->cols() ];
	int n_cols = dis->cols();
	int ident_idx;
	int n_a = (int)data->a_idx->size();
	int n_b = (int)data->b_idx->size();

	char timebuf[80];
	struct tm* newtime;
	time_t long_time;
	time( &long_time );
	newtime = localtime( &long_time );
	strftime(timebuf, 80, "%Y_%m_%d_%H_%M_%S", newtime);
	std::cout << "# Source CARMEN 1.0\n";
	std::cout << "# Mine_Type " << cmo->mine_type << "\n";
	std::cout << "# Date " << timebuf << "\n";
	std::cout << "# Data_File " << cmo->file_name_dis << "\n";
	std::cout << "# Sample_File " << cmo->file_name_sa << "\n";
	std::cout << "# Genes_File " << cmo->file_name_ga << "\n";
	std::cout << "# Class_A " << cmo->class_a << " " << data->a_idx->size() << " members.\n";
	std::cout << "# Class_B " << cmo->class_b << " " << data->b_idx->size() << " members.\n";
	std::cout << "# S_C_I " << cmo->min_sup << "," << cmo->min_conf*100 << "," << cmo->min_imp*100 << "\n";
	std::cout << "# Max_Depth " << cmo->max_depth << "\n";
	std::cout << "# Which_Class " << cmo->which_class << "\n";
	if(cmo->gene_limit.length()==0)
		std::cout << "# Gene_Limit none\n";
	else
		std::cout << "# Gene_Limit " << cmo->gene_limit << "\n";
	if( cmo->discretization.compare("none")==0)
		std::cout << "# Discretization None Using raw values, no discretization.\n";
	else
		std::cout << "# Discretization " << cmo->discretization << " Lower,Upper multiples " << cmo->disc_lower << "," << cmo->disc_upper << "\n";
	std::cout << "# Maximum_Chi2 " << cmo->max_chi2 << "\n";
	double pval=1;
	if( cmo->n_tests>0)
		pval = 0.05 / cmo->n_tests;
	std::cout << "# Tested " << cmo->n_tests << " hypotheses.  Bonferroni 0.05 corrected p-value " << pval << "\n";
	std::cout << "# Generated " << ants->size() << " rules.\n";

	//boost::unordered_map<std::string, std::string> * id2name = new boost::unordered_map<std::string, std::string>();
	Attributes* ga = new Attributes("NA");
	ga->load(cmo->file_name_ga);

	float sup_a, sup_b, imp;
	std::cout << "# chi2\t%conf\tsupA\t%supA\tsupB\t%supB\timp\tRule\tRule_ID\tRule_Dis\tt_stat\n";;
	while(iter != ants->end() ){
		n_a>0 ? sup_a = (float)(*iter)->sup_A / (float)n_a * 100 : sup_a = 0.0;
		n_b>0 ? sup_b = (float)(*iter)->sup_B / (float)n_b * 100 : sup_b = 0.0;
		imp = (*iter)->imp;
		std::cout << (*iter)->chi2 << "\t" << (*iter)->conf*100 << "\t" << (*iter)->sup_A << "\t" << sup_a << "\t" << (*iter)->sup_B << "\t" << sup_b << "\t" << imp << "\t";
		rows = (*iter)->idx;
		for(c=0;c<n_cols;c++)
			feature[c]=true;
		for(i=0;i<(int)rows->size();i++){
			idx_in_dis = (*rows)[i];
			// [GENE_INDEX]_is_[DISCRETIZED_VALUE]
			ident_idx = data->translate->at(idx_in_dis)[0];
			std::cout << ga->prop_for_identifier( data->identifiers->at( ident_idx ), "Gene Name" ) << "_is_" << translate->at(idx_in_dis)[1];
			for(c=0;c<n_cols;c++)
				if( dis->arr[idx_in_dis][c]==0 )
					feature[c]=false;
			if(i<(int)rows->size()-1)
				std::cout << "|";
		}
		std::cout << "\t";
		for(i=0;i<(int)rows->size();i++){
			idx_in_dis = (*rows)[i];
			// IDENTIFIER_is_[DISCRETIZED_VALUE]
			ident_idx = data->translate->at(idx_in_dis)[0];
			std::cout << data->identifiers->at( ident_idx ) << "_is_" << translate->at(idx_in_dis)[1];
			if(i<(int)rows->size()-1)
				std::cout << "|";
		}
		std::cout << "\t";
		for(i=0;i<(int)rows->size();i++){
			idx_in_dis = (*rows)[i];
			// IDENTIFIER_>_cutoff
			ident_idx = data->translate->at(idx_in_dis)[0];

			std::cout << data->identifiers->at( ident_idx );
			if( cmo->discretization.compare("none")==0 ){
				std::cout << "=" << translate->at(idx_in_dis)[1];
			}
			else{
				if( translate->at(idx_in_dis)[1]>0 )
					std::cout << ">" << data->discrete->cutoff_up.at(ident_idx);
				else
					std::cout << "<" << data->discrete->cutoff_dn.at(ident_idx);
			}
			if(i<(int)rows->size()-1)
				std::cout << "|";
		}
		std::cout << "\t" << (*iter)->t_stat << "\n";
		iter++;
	}
	delete ga;
}

void ClassMiner::print_ants_file(std::vector<Ant*>* ants, ClassifierDataset* data, ClassMinerOptions* cmo){
	std::vector<int*>* translate = data->translate;
	Matrix<int>* dis = data->features;

	std::ofstream f_out(cmo->file_name_out.c_str());
	if( !f_out.is_open() ){
		std::cout << "Unable to open file for writing:" << cmo->file_name_out << std::endl;
		exit(-1);
	}

	char timebuf[80];
	struct tm* newtime;
	time_t long_time;
	time( &long_time );
	newtime = localtime( &long_time );
	strftime(timebuf, 80, "%Y_%m_%d_%H_%M_%S", newtime);
	f_out << "# Source CARMEN 1.0\n";
	f_out << "# Date " << timebuf << "\n";
	f_out << "# Mine_Type " << cmo->mine_type << "\n";
	f_out << "# Data_File " << cmo->file_name_dis << "\n";
	f_out << "# Sample_File " << cmo->file_name_sa << "\n";
	f_out << "# Genes_File " << cmo->file_name_ga << "\n";
	f_out << "# Class_A " << cmo->class_a << " " << data->a_idx->size() << " members.\n";
	f_out << "# Class_B " << cmo->class_b << " " << data->b_idx->size() << " members.\n";
	f_out << "# S_C_I " << cmo->min_sup << "," << cmo->min_conf*100 << "," << cmo->min_imp*100 << "\n";
	f_out << "# Which_Class " << cmo->which_class << "\n";
	f_out << "# Max_Depth " << cmo->max_depth << "\n";
	if( cmo->discretization.compare("none")==0)
		f_out << "# Discretization none Using raw values, no discretization.\n";
	else
		f_out << "# Discretization " << cmo->discretization << " Lower,Upper multiples " << cmo->disc_lower << "," << cmo->disc_upper << "\n";
	f_out << "# Maximum_Chi2 " << cmo->max_chi2 << "\n";

	double pval=1;
	if( cmo->n_tests>0)
		pval = 0.05 / cmo->n_tests;
	f_out << "# Tested " << cmo->n_tests << " hypotheses.  Bonferroni 0.05 corrected p-value " << pval << "\n";
	f_out << "# Generated " << ants->size() << " rules.\n";

	//boost::unordered_map<std::string, std::string> * id2name = new boost::unordered_map<std::string, std::string>();
	Attributes * ga = new Attributes("NA");
	ga->load(cmo->file_name_ga);

	std::vector<Ant*>::iterator iter = ants->begin();
	std::vector<int> * rows;
	int idx_in_dis;
	int i,c;
	bool* feature = new bool[ dis->cols() ];
	int n_cols = dis->cols();
	int n_a = (int)data->a_idx->size();
	int n_b = (int)data->b_idx->size();
	int ident_idx;
	if(cmo->out_style==3){
		f_out << "A_IDX\t";
		for(int x=0;x<(int)data->a_idx->size(); x++){ f_out << data->a_idx->at(x); if(x<(int)(data->a_idx->size())-1){f_out<<"\t";} } f_out << "\n";
		f_out << "B_IDX\t";
		for(int x=0;x<(int)data->b_idx->size(); x++){ f_out << data->b_idx->at(x); if(x<(int)(data->b_idx->size())-1){f_out<<"\t";} } f_out << "\n";
	}
	float sup_a, sup_b, imp;
	f_out << "# chi2\t%conf\tsupA\t%supA\tsupB\t%supB\timp\tRule\tRule_ID\tRule_Dis\tt_stat\n";

	while(iter != ants->end() ){
		n_a>0 ? sup_a = (float)(*iter)->sup_A / (float)n_a * 100 : sup_a = 0.0;
		n_b>0 ? sup_b = (float)(*iter)->sup_B / (float)n_b * 100 : sup_b = 0.0;
		imp = (*iter)->imp;
		f_out << (*iter)->chi2 << "\t" << (*iter)->conf*100 << "\t" << (*iter)->sup_A << "\t" << sup_a << "\t" << (*iter)->sup_B << "\t" << sup_b << "\t" << imp << "\t";
		rows = (*iter)->idx;
		for(c=0;c<n_cols;c++)
			feature[c]=true;

		for(i=0;i<(int)rows->size();i++){
			idx_in_dis = (*rows)[i];
			for(c=0;c<n_cols;c++)
				if( dis->arr[idx_in_dis][c]==0 )
					feature[c]=false;
			ident_idx = data->translate->at(idx_in_dis)[0];
			f_out << ga->prop_for_identifier( data->identifiers->at( ident_idx ), "Gene Name" ) << "_is_" << translate->at(idx_in_dis)[1];

			if(i<(int)rows->size()-1)
				f_out << "|";
		}
		f_out << "\t";
		for(i=0;i<(int)rows->size();i++){
			idx_in_dis = (*rows)[i];
			// IDENTIFIER_is_[DISCRETIZED_VALUE]
			ident_idx = data->translate->at(idx_in_dis)[0];
			f_out << data->identifiers->at( ident_idx ) << "_is_" << translate->at(idx_in_dis)[1];
			if(i<(int)rows->size()-1)
				f_out << "|";
		}
		f_out << "\t";
		for(i=0;i<(int)rows->size();i++){
			idx_in_dis = (*rows)[i];
			// IDENTIFIER_>_cutoff
			ident_idx = data->translate->at(idx_in_dis)[0];
			f_out << data->identifiers->at( ident_idx );
			if( cmo->discretization.compare("none")==0 ){
				f_out << "=" << translate->at(idx_in_dis)[1];
			}
			else{
				if( translate->at(idx_in_dis)[1]>0 )
					f_out << ">" << data->discrete->cutoff_up.at(ident_idx);
				else
					f_out << "<" << data->discrete->cutoff_dn.at(ident_idx);
			}
			if(i<(int)rows->size()-1)
				f_out << "|";

		}
		f_out << "\t" << (*iter)->t_stat << "\n";
		iter++;
	}
	delete ga;
}


inline int ClassMiner::round(float x) { return int(x > 0.0 ? x + 0.5 : x - 0.5); }
inline int ClassMiner::round_up(float x) {
	int rv = int(x > 0.0 ? x + 0.5 : x - 0.5);
	if( rv<x )
		return rv+1;
	else
		return rv;
}

void ClassMiner::append_L1_L2_antecedents(ClassifierDataset * data, std::vector<Ant*>* ants, ClassMinerOptions* cmo, int& n_tests){
	// find all single rows and pairs of rows that are of interest.
	// A single or pair is interesting if sum(dis IN class_A) - sum(dis IN class_B) >= min_dff, or vice versa.
	//
	// If cmo->which_class is not "both", then we're only interested in reporting on one class.
	// Improvement is confidence[N] - confidence[N-1]
	//float min_conf = cmo->min_conf;
	int min_sup = cmo->min_sup;
	//float min_imp = cmo->min_imp;
	int max_depth = cmo->max_depth;

	n_tests=0;
	Matrix<int>* dis = data->features;

	int row1, row2, x, r1, r2;
	float A, B;

	int n_a = (int)data->a_idx->size();
	int n_b = (int)data->b_idx->size();
	double n_a_double = (double)n_a;
    double n_b_double = (double)n_b;
    int min_sup_a = round((float)min_sup / 100 * n_a);
	int min_sup_b = round((float)min_sup / 100 * n_b);
	float imp_a, imp_b, conf_a, conf_b;
	double chi2, per_sup_A, per_sup_B;
	bool report_A = false;
	bool report_B = false;
	if( cmo->which_class.compare("both")==0 || cmo->which_class.compare("a")==0 )
		report_A = true;
	if( cmo->which_class.compare("both")==0 || cmo->which_class.compare("b")==0 )
		report_B = true;

	for(row1=0; row1<dis->rows(); row1++){
		A=0; B=0;
		for(x=0;x<n_a;x++){ A += dis->arr[row1][ data->a_idx->at(x) ]; }
		for(x=0;x<n_b;x++){ B += dis->arr[row1][ data->b_idx->at(x) ]; }
		++n_tests;
		if( (data->a_is_target && A>=min_sup_a ) || (data->b_is_target && B>=min_sup_b )){
			if(n_a==0 || n_b==0)
				chi2 = 1.0;
			else{
				chi2 = calc_chi2(A, n_a-A, B, n_b-B);
				chi2 = pchisq(chi2, 1);
			}
            A==0 ? per_sup_A=0 : per_sup_A = (double)A / n_a_double;
            B==0 ? per_sup_B=0 : per_sup_B = (double)B / n_b_double;
			if(report_A && (per_sup_A > per_sup_B)){
				std::vector<int>* v = new std::vector<int>();
				v->push_back(row1);
				ants->push_back( new Ant( A/(A+B), v, (int)A, (int)B, A/(A+B), chi2 ) );
			}
            else if(report_B && (per_sup_B > per_sup_A) ){
				std::vector<int>* v = new std::vector<int>();
				v->push_back(row1);
				ants->push_back( new Ant( B/(A+B), v, (int)A, (int)B, B/(A+B), chi2 ) );
			}
		}
	}

	int n_singles = (int)ants->size();
	if( cmo->verbose ){
		std::cout << "MESSAGE: Found " << n_singles << " one-gene rules.\n";
		if(cmo->which_class.compare("a")==0 )
			std::cout << "MESSAGE: Limiting reported class to class A\n";
		if(cmo->which_class.compare("b")==0 )
			std::cout << "MESSAGE: Limiting reported class to class B\n";
		std::cout.flush();
	}
	if(max_depth<2)
		return;
	for(r1=0; r1<n_singles; r1++){
		row1 = ants->at(r1)->idx->at(0);
		for(r2=r1+1; r2<n_singles; r2++){
			row2 = ants->at(r2)->idx->at(0);
			A=0; B=0;
			for(x=0;x<n_a;x++){
				if( dis->arr[row1][ data->a_idx->at(x) ] + dis->arr[row2][ data->a_idx->at(x) ] == 2 )
					A += 1;
			}
			for(x=0;x<n_b;x++){
				if( dis->arr[row1][ data->b_idx->at(x) ] + dis->arr[row2][ data->b_idx->at(x) ] == 2 )
					B += 1;
			}
			conf_a = A/(A+B);
			conf_b = B/(A+B);
			imp_a = conf_a - ants->at(r2)->conf;
			imp_b = conf_b - ants->at(r2)->conf;
			++n_tests;
			if( ( data->a_is_target && A>=min_sup_a )
				||
				( data->b_is_target && B>=min_sup_b ) ){
				chi2 = calc_chi2(A, n_a-A, B, n_b-B);
				chi2 = pchisq(chi2, 1);
				if( chi2 <= cmo->max_chi2 ){
					std::vector<int>* v = new std::vector<int>();
					v->push_back(row1);
					v->push_back(row2);
                    A==0 ? per_sup_A=0 : per_sup_A = (double)A / n_a_double;
                    B==0 ? per_sup_B=0 : per_sup_B = (double)B / n_b_double;
                    if(report_A && (per_sup_A > per_sup_B) )
						ants->push_back(new Ant( conf_a, v, (int)A, (int)B, imp_a, chi2));
					else if(report_B && (per_sup_B > per_sup_A) )
						ants->push_back(new Ant( conf_b, v, (int)A, (int)B, imp_b, chi2));
				}
			}
		}
	}
	if(cmo->verbose){
		std::cout << "MESSAGE: Found " << ants->size() - n_singles << " two-gene rules.\n";
		if(cmo->which_class.compare("a")==0 )
			std::cout << "MESSAGE: Limiting reported class to class A\n";
		if(cmo->which_class.compare("b")==0 )
			std::cout << "MESSAGE: Limiting reported class to class B\n";
		std::cout.flush();
	}
}


void ClassMiner::append_L_ants(ClassifierDataset * data, std::vector<Ant*>* ants, ClassMinerOptions* cmo, int L, int& n_tests){

	//float min_conf = cmo->min_conf;
	int min_sup = cmo->min_sup;
	//float min_imp = cmo->min_imp;
	int timeout = cmo->timeout;
	time_t start_time = cmo->start_time;
	float conf_a, conf_b, imp_a, imp_b;
	Matrix<int>* dis = data->features;

	bool report_A = false;
	bool report_B = false;
	if( cmo->which_class.compare("both")==0 || cmo->which_class.compare("a")==0 )
		report_A = true;
	if( cmo->which_class.compare("both")==0 || cmo->which_class.compare("b")==0 )
		report_B = true;

	std::vector<int> * single_idx = new std::vector<int>();
	std::vector<int> * single_val = new std::vector<int>();
	std::vector<int> * L_paths = new std::vector<int>();

	Node root = Node(Node::ROOT);
	for( std::vector<Ant*>::iterator iter_ant = ants->begin(); iter_ant != ants->end(); iter_ant++ ){
		root.add_path( (*iter_ant)->idx );
	}
	root.get_kids(single_val);

	std::vector<int>::iterator itr_single;

	Ant* ant;
	int n_ants = (int)ants->size();

	int i, sum, x;
	float A, B;
	int n_a = (int)data->a_idx->size();
	int n_b = (int)data->b_idx->size();
    double n_a_double = (double)n_a;
    double n_b_double = (double)n_b;

	int min_sup_a = round_up((float)min_sup / 100 * n_a);
	int min_sup_b = round_up((float)min_sup / 100 * n_b);

	vector<int>* rows = new vector<int>(L);
	int ant_ctr;
	int idx_size;
	bool duplicate;
	std::vector<int> * values = new std::vector<int>();
	double chi2, per_sup_A, per_sup_B;
	time_t last_time = time(NULL);
	for(ant_ctr=0; ant_ctr<n_ants; ant_ctr++ ){
		if(ant_ctr % 500 ==0 ){
			if(timeout>0 && time(NULL)-start_time>timeout)
				return;
			last_time = time(NULL);
		}
		ant = ants->at(ant_ctr);
		idx_size = (int)(ant->idx->size());
		if( idx_size==L-1){
			itr_single = single_val->begin();
			for(i=0; i<idx_size; i++ )
				rows->at(i) = ant->idx->at(i);
			while(itr_single != single_val->end()){
				A=0; B=0;
				rows->at(L-1) = (*itr_single); // write potential new value
				++itr_single;
				// if new value is already present in idx, skip to next value
 				duplicate=false;
				for(i=0;i<idx_size;i++)
					if( ant->idx->at(i)==rows->at(L-1))
						duplicate=true;
				if(duplicate){ continue; }

				// if path has already been discovered, pop off new item and skip to next value
				if( root.has_path(rows) ){
					continue;
				}

				// if sum( dis[r,c] ) = L, increment A or B
				for(x=0;x<n_a;x++){
					sum=0;
					for( i=0; i<L; i++){ sum += dis->arr[rows->at(i)][data->a_idx->at(x)]; }
					if(sum==L){ A += 1; }
				}
				for(x=0;x<n_b;x++){
					sum=0;
					for( i=0; i<L; i++){ sum += dis->arr[rows->at(i)][ data->b_idx->at(x) ]; }
					if(sum==L){ B += 1; }
				}
				conf_a = A/(A+B);
				conf_b = B/(A+B);
				imp_a = conf_a - ant->conf;
				imp_b = conf_b - ant->conf;
				++n_tests;
				if( (data->a_is_target && A>=min_sup_a )
					||
					(data->b_is_target && B>=min_sup_b ) ){
					chi2 = calc_chi2(A, n_a-A, B, n_b-B);
					chi2 = pchisq(chi2, 1);
					if( chi2 <= cmo->max_chi2 ){
						std::vector<int>* v = new std::vector<int>();
						for(i=0; i<L; i++)
							v->push_back(rows->at(i));
						root.add_path( rows );  // add to tree
                        A==0 ? per_sup_A=0 : per_sup_A = (double)A / n_a_double;
                        B==0 ? per_sup_B=0 : per_sup_B = (double)B / n_b_double;
						if(report_A && per_sup_A > per_sup_B )
							ants->push_back(new Ant( conf_a, v, (int)A, (int)B, imp_a, chi2 ));
						else if( report_B && per_sup_B > per_sup_A )
							ants->push_back(new Ant( conf_b, v, (int)A, (int)B, imp_b, chi2 ));
					}
				}
			}
		}
	}
	delete single_idx;
	delete single_val;
	delete L_paths;
	delete values;
	delete rows;
}


void ClassMiner::filter_ants(std::vector<Ant*>* ants, ClassifierDataset* data, float min_conf, int min_sup, float min_imp, double max_chi2){
	int n_a = (int)data->a_idx->size();
	int n_b = (int)data->b_idx->size();
	int min_sup_a = round_up((float)min_sup / 100 * n_a);
	int min_sup_b = round_up((float)min_sup / 100 * n_b);

	Ant* ant;
	std::vector<Ant*> ants_filtered = std::vector<Ant*>();
	for(int i=0;i<(int)ants->size();i++){
		ant = ants->at(i);
		if( ant->chi2 <= max_chi2 && 
			ant->conf >= min_conf && 
			(ant->sup_A>=min_sup_a || ant->sup_B>=min_sup_b) && 
			ant->imp >= min_imp)
		{
			ants_filtered.push_back( ant );
		}
	}
	ants->clear();
	ants->insert(ants->begin(), ants_filtered.begin(), ants_filtered.end());
	std::sort(ants->begin(), ants->end(), comp_desc);
}


void ClassMiner::assign_t_stats(std::vector<Ant*>* ants, ClassifierDataset* data){
	// for each ant, calculate the minimum value for abs(t_stat) and
	// store it in ant->t_stat.  Use this to rank rules that have
	// similar conf/support.

	if( data->a_idx->size()==0 || data->b_idx->size()==0)
		return;
	Ant* ant;
	std::vector<int>* rows;
	int idx_in_dis, a, i, r;
	float mu_a, var_a, mu_b, var_b;
	double t_stat, min_t_stat=0;
	int n_a = (int)data->a_idx->size();
	int n_b = (int)data->b_idx->size();
	Perm perm;
	boost::unordered_map<int, double> row2tstat;
	for(a=0;a<(int)ants->size(); a++){
		ant = ants->at(a);
		rows = ant->idx;
		for(i=0;i<(int)rows->size();i++){
			idx_in_dis = (*rows)[i];
			r = data->translate->at(idx_in_dis)[0];
			if( row2tstat.find(r) == row2tstat.end() ){
				perm.trimmed_stats(data->raw_data->data, r, data->a_idx, 0, mu_a, var_a);
				perm.trimmed_stats(data->raw_data->data, r, data->b_idx, 0, mu_b, var_b);
				t_stat = abs((double)perm.find_t_stat( mu_a, var_a, n_a, mu_b, var_b, n_b));
				if(i==0)
					min_t_stat = t_stat;
				else{
					if(t_stat < min_t_stat)
						min_t_stat = t_stat;
				}
				row2tstat[r] = min_t_stat;
			}
			else{
				if(i==0){
					min_t_stat = row2tstat[r];
				}
			}
			ant->t_stat = row2tstat[r];
		}
	}
}
