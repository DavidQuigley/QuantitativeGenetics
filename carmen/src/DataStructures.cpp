#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <sstream>
#include <boost/random.hpp>
#include "DataStructures.h"

QTL::QTL(std::string locus_id, std::string locus_name, std::string probe_id, std::string probe_name, double pval){
	this->locus_id = locus_id;
	this->locus_name = locus_name;
	this->probe_id = probe_id;
	this->probe_name = probe_name;
	this->pval = pval;
}

Ant::Ant(float conf, std::vector<int> * idx, int sup_A, int sup_B, float imp, double chi2=0){
	this->conf = conf;
	this->sup_A = sup_A;
	this->sup_B = sup_B;
	this->idx = idx;
	this->imp = imp;
	this->chi2 = chi2;
	this->t_stat = 0;
}

std::string Ant::rule_as_str(){
	std::vector<int> t;
	for(int i=0; i<(int)this->idx->size(); i++)
		t.push_back(this->idx->at(i));
	sort(t.begin(), t.end());
	
	std::ostringstream oss;
	oss << t.at(0);
	for(int i=1; i<(int)t.size(); i++){
		oss << "," << t.at(i);
	}
	return oss.str();
}

Ant::~Ant(){
	delete this->idx;
}
