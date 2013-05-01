#include <iostream>
#include <stdio.h>
#include <string>
#include <fstream>
#include <vector>
#include <algorithm>
#include <math.h>
#include <time.h>
using namespace std;
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/random.hpp>
using namespace boost;
#include "DataStructures.h"
#include "ClassMinerOptions.h"
#include "Rawdata.h"
#include "Attributes.h"
#include "Discretize.h"
#include "Dataset.h"
#include "graph.h"
#include "ParseOptions.h"
#include "chi2.h"
#include "Rule.h"
#include "ClassMiner.h"


// Author: David Quigley, UCSF Cancer Research Institute.
// Email:  dquigley@cc.ucsf.edu
// Date:   2007-2010

void find_core_rules(std::vector< std::vector<Ant*>* >&  ant_farm, std::vector<Ant*>* full_ants, int tolerance){
	// For each vector of ants (ant_farm is a vector of vectors), count how many times we see
	// each rule.  If there are N vectors and we see the rule N - tolerance times, the rule is considered
	// a "core rule".  Core rules are deposited in full_ants.
	HASH_S_I rule2count;
	std::vector<Ant*>* ants;
	std::string rule;
	int n_folds = (int)ant_farm.size();

	for(int i=0; i<n_folds; i++){
		ants = ant_farm.at(i);
		for(int r=0; r< (int)ants->size(); r++){
			rule = ants->at(r)->rule_as_str();
			if( rule2count.find(rule) == rule2count.end() )
				rule2count[rule] = 1;
			else
				rule2count[rule] = rule2count[rule] + 1;
		}
	}
	Ant* ant;

	std::vector<Ant*> ants_filtered = std::vector<Ant*>();
	for(int i=0;i<(int)full_ants->size();i++){
		ant = full_ants->at(i);
		if( rule2count[ ant->rule_as_str() ] >= n_folds-tolerance ){
			ants_filtered.push_back( ant );
		}
	}
	full_ants->clear();
	full_ants->insert(full_ants->begin(), ants_filtered.begin(), ants_filtered.end());
}


bool filter_out_idx(std::vector<int>* V, int target){
	// Pop all instances of an arbitrary value "target" off of V.
	// in practice, used to remove a sample from a list of sample indices.
	// Returns whether an instance of "target" was found.
	std::vector<int> cp;
	bool found = false;
	for(int i=0; i< (int)V->size(); i++){
		if(V->at(i)==target)
			found = true;
		else
			cp.push_back( V->at(i) );
	}
	V->clear();
	for(int i=0; i< (int)cp.size(); i++)
		V->push_back( cp.at(i) );
	return found;
}


int main(int argc, char *argv[]){
	vector<Option*>* options = new vector<Option*>();
	options->push_back( new Option("ruleset_file", "r", "Path to file with ruleset.", "", OPT_REQUIRED));
	options->push_back( new Option("tolerance", "t", "# of rulesets where a rule may be missing, default to 0", "0", OPT_OPTIONAL));
	options->push_back( new Option("out_file", "o", "Output file, defaults to standard output", "", OPT_OPTIONAL));
	options->push_back( new Option("verbose", "v", "Print out progress to screen, T or F.  Defaults to F.", "F", OPT_OPTIONAL));

	std::stringstream ssp;
	ssp << "coreminer\n\nPerforms core rule discovery on an association rule mining dataset.\n";
	int retval=read_args(argc, argv, options, ssp.str());
	if( retval != 0 ){
		std::cout << "ERROR reading arguments, bailing out.\n";
		return(0);
	}
	ClassMinerOptions * cmo = new ClassMinerOptions();
	int r=0;
	RuleSet* rs = new RuleSet(options->at(r++)->value, 1);
	int tolerance = boost::lexical_cast<int>(options->at(r++)->value);
	cmo->file_name_dis = rs->file_name_dis;
	cmo->file_name_sa = rs->file_name_sa;
	cmo->file_name_ga = rs->file_name_ga;
	cmo->discretization = rs->disc;
	cmo->disc_lower = (float)rs->disc_lower;
	cmo->disc_upper = (float)rs->disc_upper;
	cmo->class_a = rs->class_a;
	cmo->class_b = rs->class_b;
	cmo->min_conf = (float)rs->minconf / (float)100.0;;
	cmo->min_imp = (float)rs->minimp / (float)100.0;;
	cmo->min_sup = (int)rs->minsup;
	cmo->max_chi2 = rs->max_chi2;
	cmo->max_depth = rs->max_depth;
	cmo->file_name_out = options->at(r++)->value;
	if( options->at(r++)->value.compare("T")==0 )
		cmo->verbose = true;
	Attributes* sa;
	Attributes* ga;
	ClassifierDataset* data = new ClassifierDataset();
	try{
		if(cmo->verbose){
			std::cout << "MESSAGE: Loading data set...\n";
			std::cout.flush();
		}
		sa = new Attributes("NA");
		ga = new Attributes("NA");
		sa->load(cmo->file_name_sa);
		ga->load(cmo->file_name_ga);
		data->load(sa, ga, cmo);
		if(cmo->verbose){
			std::cout << "MESSAGE: Loaded data set.\n";
			std::cout.flush();
		}
	}
	catch(std::string msg){
		std::cout << "ERROR: " << msg << "\n";
		std::cout.flush();
		exit(0);
	}

	std::vector<ClassMiner*> miners;
	std::vector<int> ids_to_pop;
	std::vector<Ant*>* ants;
	std::vector< std::vector<Ant*>* > ant_farm;
	for(int i=0; i<(int)data->a_idx->size(); i++){
		ids_to_pop.push_back( data->a_idx->at(i) );
	}
	for(int i=0; i<(int)data->b_idx->size(); i++){
		ids_to_pop.push_back( data->b_idx->at(i) );
	}
	std::string base_class_a = cmo->class_a;
	std::string base_class_b = cmo->class_b;
	std::string fold_class_a, fold_class_b;
	int n_tests=0, idx_to_pop;
	int ant_size_before=0;
	bool sample_in_class_a;

	if(cmo->verbose){
		std::cout << "MESSAGE: Beginning to build " << ids_to_pop.size() << " rule sets.\n";
		std::cout.flush();
	}
	for(int i=0; i<(int)ids_to_pop.size(); i++){
		idx_to_pop = ids_to_pop.at(i);
		sample_in_class_a = filter_out_idx(data->a_idx, idx_to_pop);
		if(!sample_in_class_a)
			filter_out_idx(data->b_idx, idx_to_pop);

		ClassMiner* cm = new ClassMiner();
		ants = new std::vector<Ant*>();
		n_tests=0;
		ant_size_before=0;
		cm->append_L1_L2_antecedents(data, ants, cmo, n_tests);
		int L=3;
		while(L<=cmo->max_depth && ant_size_before != (int) ants->size()){
			ant_size_before = (int) ants->size();
			cm->append_L_ants(data, ants, cmo, L, n_tests);
			++L;
		}
		cmo->n_tests = n_tests;
		cm->filter_ants(ants, data, cmo->min_conf, cmo->min_sup, cmo->min_imp, cmo->max_chi2 );
		ant_farm.push_back(ants);
		if(sample_in_class_a)
			data->a_idx->push_back(idx_to_pop);
		else
			data->b_idx->push_back(idx_to_pop);
		if(cmo->verbose){
			std::cout << "MESSAGE: " << ants->size() << " rules remain after filtering.\n";
			std::cout << "MESSAGE: Finished " << i+1 << " of " << ids_to_pop.size() << " rule sets.\n";
			std::cout.flush();
		}
	}

	//std::vector<std::string*>* core_rules = new std::vector<std::string*>();
	if(cmo->verbose){
		std::cout << "MESSAGE: Calculating core rules...\n";
		std::cout.flush();
	}
	cmo->class_a = base_class_a;
	cmo->class_b = base_class_b;

	ClassMiner* cm = new ClassMiner();
	ants = new std::vector<Ant*>();
	n_tests=0;
	ant_size_before=0;
	cm->append_L1_L2_antecedents(data, ants, cmo, n_tests);
	int L=3;
	while(L<=cmo->max_depth && ant_size_before != (int) ants->size()){
		ant_size_before = (int) ants->size();
		cm->append_L_ants(data, ants, cmo, L, n_tests);
		++L;
	}
	cmo->n_tests = n_tests;
	cm->assign_t_stats(ants, data);
	cm->filter_ants(ants, data, cmo->min_conf, cmo->min_sup, cmo->min_imp, cmo->max_chi2 );

	find_core_rules(ant_farm, ants, tolerance);
	std::stringstream ss;
	ss << "Core_Rules tolorance ";
	ss << tolerance;
	cmo->mine_type = ss.str();
	if( cmo->file_name_out.length()==0)
		cm->print_ants_stdout(ants, data, cmo);
	else
		cm->print_ants_file(ants, data, cmo);
}
