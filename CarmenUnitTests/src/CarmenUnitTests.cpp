#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <set>
#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/unordered_map.hpp>
#include <boost/random.hpp>
#include <boost/thread/mutex.hpp>
namespace alg = boost::algorithm;
using namespace std;
using namespace boost;
#include "DataStructures.h"
#include "Parser.h"
#include "ClassMinerOptions.h"
#include "Rawdata.h"
#include "Attributes.h"
#include "Discretize.h"
#include "Dataset.h"
#include "graph.h"
#include "ParseOptions.h"
#include "ClassMiner.h"
#include "chi2.h"
#include "Perm.h"
#include "Rule.h"
#include "Evaluation.h"
#include "Classifier.h"
#include "adaboost.h"
#include "spear.h"

bool verbose = false;

void test(std::string testname, int test_number, int a, int b){
	if(a!=b)
		std::cout << "*** FAIL *** " << testname << " (" << test_number << ") : " << a << " != expected " << b << "\n";
	else if(verbose)
		std::cout << "PASS " << testname << " (" << test_number << ") : " << a << " = " << b << "\n";
}

void test(std::string testname, int test_number, float a, float b){
	if(a!=b)
		std::cout << "*** FAIL *** " << testname << " (" << test_number << ") : " << a << " != expected " << b << "\n";
	else if(verbose)
		std::cout << "PASS " << testname << " (" << test_number << ") : " << a << " = " << b << "\n";
}

void test(std::string testname, int test_number, int a, double b){
	if(a!=int(b))
		std::cout << "*** FAIL *** " << testname << " (" << test_number << ") : " << a << " != expected " << b << "\n";
	else if(verbose)
		std::cout << "PASS " << testname << " (" << test_number << ") : " << a << " = " << b << "\n";
}

void test(std::string testname, int test_number, bool a, bool b){
	if(a!=b)
		std::cout << "*** FAIL *** " << testname << " (" << test_number << ") : " << a << " != expected " << b << "\n";
	else if(verbose)
		std::cout << "PASS " << testname << " (" << test_number << ") : " << a << " = " << b << "\n";
}


void run_CacheHash(){
	std::cout << "Unit Test: CacheHash\n";
	Matrix<float>* data = new Matrix<float>(2,5);
	std::vector<int> valid_cols;
	valid_cols.push_back(0);
	valid_cols.push_back(1);
	valid_cols.push_back(2);
	valid_cols.push_back(3);
	valid_cols.push_back(4);
	data->set_missing_value(0,1);

	CacheHash* hsh= new CacheHash(data, valid_cols);

	std::set<int>* missing=NULL;
	missing = hsh->missing_by_idx(0);
	test(std::string("CacheHash missing one item"), 1, missing->find(1) != missing->end(), true );
	missing = hsh->missing_by_idx(1);
	test(std::string("CacheHash missing empty"), 2, missing->size(), 0);

	std::set<int> vals;
	std::set<int> vals_not_present;
	std::vector<int>* ranks1= new std::vector<int>();
	std::vector<int>* ranks2= new std::vector<int>();
	vals.insert(1);
	ranks1->push_back(0); ranks1->push_back(1); ranks1->push_back(2);
	ranks2->push_back(0); ranks2->push_back(1); ranks2->push_back(2);

	hsh->store(0, vals, ranks1);
	hsh->store(1, vals, ranks2);
	test(std::string("CacheHash has"), 3, hsh->has(0, vals), true);
	test(std::string("CacheHash has idx not present"), 4, hsh->has(2, vals), false);
	test(std::string("CacheHash has good idx bad values"), 5, hsh->has(0, vals_not_present), false);
	test(std::string("CacheHash has before clear"), 6, hsh->has(1, vals), true);
	hsh->clear(1);
	test(std::string("CacheHash has after clear"), 7, hsh->has(1, vals), false);
	std::vector<int>* ranks_back = hsh->retreive(0, vals );
	test(std::string("CacheHash retreive size"), 8, int(ranks_back->size()), 3);
	test(std::string("CacheHash retreive 0"), 9, ranks_back->at(0), 0);
	test(std::string("CacheHash retreive 1"), 10, ranks_back->at(1), 1);
	test(std::string("CacheHash retreive 2"), 11, ranks_back->at(2), 2);
	delete hsh;

}

void run_Perm(){
	std::cout << "Unit Test: Perm\n";
	string tn = "Perm";
	int r=1;
	test( tn, r++, int(round( float(2.6) )), int(3));
	test( tn, r++, int(round( float(-1.2) )), int(-1));

	Matrix<float>* data = new Matrix<float>(1, 6, true);
	data->arr[0][0]=-99;  data->arr[0][1]=2;  data->arr[0][2]=3;
	data->arr[0][3]=5;  data->arr[0][4]=6;  data->arr[0][5]=2000;
	int n_cols = 6;
	vector<int>* cols = new vector<int>();
	for(int i=0;i<n_cols;i++){ cols->push_back(i);}
	float mu, var;
	Perm p;
	p.trimmed_stats( data, 0, cols, 17, mu, var);
	test( tn, r++, mu, (float)4 );
	var *= 10;

	test( tn, r++, int(round(var)), int(33) );
	data->set_missing_value(0,2);
	data->set_missing_value(0,3);
	p.trimmed_stats( data, 0, cols, 17, mu, var);
	test( tn, r++, mu, (float)4 );
	test( tn, r++, var, (float)8 );

	float mu_a = (float)1.5;
	float var_a = (float)0.277777778;
	float mu_b = (float)4.2;
	float var_b = (float)3.066666667;
	int n_a = 10;
	int n_b = 10;
	float t = p.find_t_stat( mu_a, var_a, n_a, mu_b, var_b, n_b);
	test( tn, r++, int(round(t*100)), int(-467) );
	int x = p.rnd(2);
	test( tn, r++, x==0||x==1, true );

	vector<float> acc;
	acc.push_back(1); acc.push_back(2); acc.push_back(3); acc.push_back(4); acc.push_back(5);
 	test( tn, r++, p.find_variance((float)3, acc, 0, 5), (float)2.5);
	test( tn, r++, p.find_variance((float)3, acc, 1, 4), (float)1);

	vector<int>* sh = new vector<int>();
	sh->push_back(1); sh->push_back(2); sh->push_back(3);
	p.shuffle(sh);
	test( tn, r++, (int)sh->size(), 3);  // just to call this
	//TODO:
	// write tests:
	//permute( vector<int>* idx_a, vector<int>* idx_b, vector<int>* perm_a, vector<int>* perm_b)
	//find_p_value( float r, vector< vector<int>*> * perms_a, vector< vector<int>*> * perms_b, Matrix<float>* data, int row, int trim_percent)

}


void run_Rule(){
	std::cout << "Unit Test: Rule\n";
	Rule* r = new Rule();
	r->rules_g.push_back("Gene_1_is_1");
	r->rules_g.push_back("Gene_2_is_1");
	r->rules_i.push_back("1_at_is_1");
	r->rules_i.push_back("2_at_is_1");
	r->rules_r.push_back("1_at>4");
	r->rules_r.push_back("2_at>6");
	r->chi2=1e-10;
	r->conf_p=100;
	r->imp=100;
	r->rule_class=Rule::CLASS_A;
	r->sup_a_n=50;
	r->sup_a_p=95;
	r->sup_b_n=0;
	r->sup_b_p=0;
	r->t_stat=0;
	std::vector<Rule*> rule_bucket;
	r->clone_out_probes(rule_bucket);
	int x=1;
	test( "RULE clone_probes", x++, (int)rule_bucket.size(), 2);
	test( "RULE clone_probes", x++, (int)rule_bucket.at(0)->rules_i.at(0).compare("1_at_is_1"), 0);
	test( "RULE clone_probes", x++, (int)rule_bucket.at(0)->rules_i.size(), 1);
	test( "RULE clone_probes", x++, (int)rule_bucket.at(1)->rules_i.at(0).compare("2_at_is_1"), 0);
	test( "RULE clone_probes", x++, (int)rule_bucket.at(1)->rules_i.size(), 1);
	delete r;
}

void run_RuleSet(std::string basedir){
	std::cout << "Unit Test: RuleSet\n";
	RuleSet* rs = new RuleSet(basedir + "/test.ruleset", false);
	int r=1;
	test("RULESET n_rules", r++,  (int)rs->rules.size(), 4 );
	test("RULESET r0", r++,  rs->rules.at(0)->rules_g.at(0).compare("Gene1_is_1"), 0 );
	delete rs;
	rs = new RuleSet(basedir + "/test.ruleset", false);
	test("RULESET n_rules", r++,  (int)rs->rules.size(), 4 );
	test("RULESET r0", r++,  rs->rules.at(0)->rules_g.at(0).compare("Gene1_is_1"), 0 );
	test("RULESET r1", r++,  rs->rules.at(1)->rules_g.at(0).compare("Gene2_is_-1"), 0 );
	test("RULESET r2", r++,  rs->rules.at(2)->rules_g.at(0).compare("Gene1_is_-1"), 0 );
	test("RULESET r3", r++,  rs->rules.at(3)->rules_g.at(0).compare("Gene1_is_1"), 0 );
	test("RULESET r3 2", r++,  rs->rules.at(3)->rules_g.at(1).compare("Gene3_is_1"), 0 );
}



void run_Adaboost(){
	std::cout << "Unit Test: Adaboost\n";
	int r=1;

	Matrix<double>* features = new Matrix<double>(3,4);
	//  labels    - - + +
	//  features: 1 1 0 0    (perfect neg)
	//            0 0 1 1    (perfect pos)
	//            1 0 1 0    (useless)
	//  is train  Y N Y N
	features->arr[0][0]=1;  features->arr[0][1]=1;
	features->arr[1][2]=1;  features->arr[1][3]=1;
	features->arr[2][0]=1;  features->arr[2][2]=1;

	std::vector<double>* targets = new std::vector<double>();
	targets->push_back(-1); targets->push_back(-1); targets->push_back(1); targets->push_back(1);
	std::vector<int>* idx_train = new std::vector<int>();
	idx_train->push_back(0);  idx_train->push_back(2);
	std::vector<int>* idx_test = new std::vector<int>();
	idx_test->push_back(1);  idx_test->push_back(3);

	adaboost* ada = new adaboost(features, targets, idx_train, idx_test);
	ada->test_mode = true;

	test( "ada n_tst", r++, ada->n_tst, 2);
	test( "ada n_tr", r++, ada->n_tr, 2);
	test( "ada HO_tr", r++, ada->HO_tr->rows(), 4);
	test( "ada HO_tr", r++, (int)ada->HO_tr->arr[0][0], 1);
	test( "ada HO_tr", r++, (int)ada->HO_tr->arr[1][0], 0);
	test( "ada HO_tr", r++, (int)ada->HO_tr->arr[2][0], 1);
	test( "ada HO_tr", r++, (int)ada->HO_tr->arr[3][0], 0);
	test( "ada HO_tst", r++, ada->HO_tst->rows(), 4);
	test( "ada HO_tst", r++, (int)ada->HO_tst->arr[0][0], 0);
	test( "ada HO_tst", r++, (int)ada->HO_tst->arr[1][0], 1);
	test( "ada HO_tst", r++, (int)ada->HO_tst->arr[2][0], 0);
	test( "ada HO_tst", r++, (int)ada->HO_tst->arr[3][0], 1);
	test( "ada L", r++, ada->L->rows(), 4);
	test( "ada L", r++, (int)ada->L->arr[0][0], -1);
	test( "ada L", r++, (int)ada->L->arr[1][0], -1);
	test( "ada L", r++, (int)ada->L->arr[2][0], 1);
	test( "ada L", r++, (int)ada->L->arr[3][0], 1);
	test( "ada HO_trUp", r++, ada->HO_trUp->rows(), 4);
	test( "ada HO_trUp", r++, (int)ada->HO_trUp->arr[0][0], 0);
	test( "ada HO_trUp", r++, (int)ada->HO_trUp->arr[1][0], 0);
	test( "ada HO_trUp", r++, (int)ada->HO_trUp->arr[2][0], 1);
	test( "ada HO_trUp", r++, (int)ada->HO_trUp->arr[3][0], 0);
	test( "ada HO_trDn", r++, ada->HO_trDn->rows(), 4);
	test( "ada HO_trDn", r++, (int)ada->HO_trDn->arr[0][0], 1);
	test( "ada HO_trDn", r++, (int)ada->HO_trDn->arr[1][0], 0);
	test( "ada HO_trDn", r++, (int)ada->HO_trDn->arr[2][0], 0);
	test( "ada HO_trDn", r++, (int)ada->HO_trDn->arr[3][0], 0);
	ada->initialize();
	test( "ada init", r++, ada->W->rows(), 4);
	test( "ada init", r++,  abs(ada->W->arr[0][0]), (double)0.5);
	test( "ada init", r++,  abs(ada->W->arr[1][0]), (double)0);
	test( "ada init", r++,  abs(ada->W->arr[2][0]), (double)0.5);
	test( "ada init", r++,  abs(ada->W->arr[3][0]), (double)0);

	double sc = 0.005;

	test( "ada init", r++, ada->pred_tr->rows(), 4);
	test( "ada init", r++,  abs(ada->pred_tr->arr[0][0]), sc);
	test( "ada init", r++,  abs(ada->pred_tr->arr[1][0]), (double)0);
	test( "ada init", r++,  abs(ada->pred_tr->arr[2][0]), sc);
	test( "ada init", r++,  abs(ada->pred_tr->arr[3][0]), (double)0);
	test( "ada init", r++,  ada->pred_tst->rows(), 4);
	test( "ada init", r++,  abs(ada->pred_tst->arr[0][0]), (double)0);
	test( "ada init", r++,  abs(ada->pred_tst->arr[1][0]), sc);
	test( "ada init", r++,  abs(ada->pred_tst->arr[2][0]), (double)0);
	test( "ada init", r++,  abs(ada->pred_tst->arr[3][0]), sc);

	ada->update_weights();
	bool w0 = ( (int)(ada->W->arr[0][0]*10000) == 5024 || (int)(ada->W->arr[0][0]*10000)== 4975 );
	bool w2 = ( (int)(ada->W->arr[2][0]*10000) == 5024 || (int)(ada->W->arr[2][0]*10000)== 4975 );
	test( "ada update", r++, w0, true);
	test( "ada update", r++, (int)(ada->W->arr[1][0]*10000), 0);
	test( "ada update", r++, w2, true);
	test( "ada update", r++, (int)(ada->W->arr[3][0]*10000), 0);

	delete ada;

	ada = new adaboost(features, targets, idx_train, idx_test);

	ada->boost(2);
	weak_learner* wl0 = ada->weak_learners->at(0);
	weak_learner* wl1 = ada->weak_learners->at(1);
	test( "ada boost idx", r++,  (wl0->idx==0 && wl1->idx==1) || (wl0->idx==1 && wl1->idx==0), true);
	delete ada;

	ada = new adaboost(features, targets, idx_train, idx_test);
	ada->test_mode=true;
	ada->initialize();
	ada->update_weights();
	ada->W->arr[0][0] = 0.5024;
	ada->W->arr[2][0] = 0.5024;
	std::vector<weak_learner*> wls;
	ada->choose_weak_learners(wls, ada->HO_tr, ada->HO_tst);
	test( "ada update", r++, (int)wls.size(), 2);
	test( "ada update", r++,  wls.at(0)->idx != wls.at(1)->idx, true);
	test( "ada update", r++,  wls.at(0)->idx ==1 || wls.at(0)->idx==0, true);
	test( "ada update", r++,  wls.at(1)->idx ==1 || wls.at(1)->idx==0, true);
	delete ada;
	delete features;
	delete targets;
	delete idx_train;
	delete idx_test;
}

void run_Discretize(){
	std::cout << "Unit Test: Discretize\n";
	int r=1;
	std::vector<float>* f = new std::vector<float>();
	f->push_back(1);
	f->push_back(2);
	f->push_back(3);
	f->push_back(4);
	f->push_back(15);
	std::vector<int> idx;
	idx.push_back(0); idx.push_back(1); idx.push_back(2); idx.push_back(3); idx.push_back(4);
	Rawdata* rd = new Rawdata();
	Matrix<float>* d = new Matrix<float>(2,5, true);
	// 1   5  5  5  12
	// NA  5  5  5  12
	d->arr[0][0]=1; d->arr[0][1]=5; d->arr[0][2]=5; d->arr[0][3]=5; d->arr[0][4]=12;
	d->arr[1][0]=1; d->arr[1][1]=5; d->arr[1][2]=5; d->arr[1][3]=5; d->arr[1][4]=12;
	d->set_missing_value(1,0);
	rd->data = d;
	Discretize* disc = new Discretize(rd, idx, "none", 0, 0);
	test("dis find_mean", r++, int(disc->find_mean(f)), 5);		// 1
	test("dis find_median", r++, int(disc->find_median(f)), 3);		// 2
	test("dis find_stdev", r++, int(round(disc->find_stdev(disc->find_mean(f), f)*10)), 57);	// 3
	test("dis none", r++, (int)disc->dis->arr[0][0], 1); // 4
	test("dis none", r++, (int)disc->dis->arr[0][1], 5);	// 5
	test("dis none", r++, (int)disc->dis->arr[0][2], 5);	// 6
	test("dis none", r++, (int)disc->dis->arr[0][3], 5);	// 7
	test("dis none", r++, (int)disc->dis->arr[0][4], 12);	// 8
	test("dis none", r++, (int)disc->dis->arr[1][0], 0);	// 9
	test("dis none", r++, (int)disc->dis->arr[1][1], 5);	// 10
	test("dis none", r++, (int)disc->dis->arr[1][2], 5);	// 11
	test("dis none", r++, (int)disc->dis->arr[1][3], 5);	// 12
	test("dis none", r++, (int)disc->dis->arr[1][4], 12);	// 13
	delete disc;

	disc = new Discretize(rd, idx, "sd", 0.5, 0.5);
	test("dis SD", r++, (int)disc->dis->arr[0][0], -1);
	test("dis SD", r++, (int)disc->dis->arr[0][1], 0);
	test("dis SD", r++, (int)disc->dis->arr[0][2], 0);
	test("dis SD", r++, (int)disc->dis->arr[0][3], 0);
	test("dis SD", r++, (int)disc->dis->arr[0][4], 1);
	test("dis SD", r++, (int)disc->dis->arr[1][0], 0);
	test("dis SD", r++, (int)disc->dis->arr[1][1], -1);
	test("dis SD", r++, (int)disc->dis->arr[1][2], -1);
	test("dis SD", r++, (int)disc->dis->arr[1][3], -1);
	test("dis SD", r++, (int)disc->dis->arr[1][4], 1);
	delete disc;

	disc = new Discretize(rd, idx, "abs", .5, 10);
	test("dis abs", r++, (int)disc->dis->arr[0][0], 0); //14
	test("dis abs", r++, (int)disc->dis->arr[0][1], 0); //15
	test("dis abs", r++, (int)disc->dis->arr[0][2], 0); //16
	test("dis abs", r++, (int)disc->dis->arr[0][3], 0); //17
	test("dis abs", r++, (int)disc->dis->arr[0][4], 1); //18
	test("dis abs", r++, (int)disc->dis->arr[1][0], 0); //19
	test("dis abs", r++, (int)disc->dis->arr[1][1], 0); //20
	test("dis abs", r++, (int)disc->dis->arr[1][2], 0); //21
	test("dis abs", r++, (int)disc->dis->arr[1][3], 0); //22
	test("dis abs", r++, (int)disc->dis->arr[1][4], 1); //23
	delete disc;

	disc = new Discretize(rd, idx, "MAD", .5, .5);
	test("dis MAD", r++, (int)disc->dis->arr[0][0], -1); //24
	test("dis MAD", r++, (int)disc->dis->arr[0][1], -1); //25
	test("dis MAD", r++, (int)disc->dis->arr[0][2], -1);
	test("dis MAD", r++, (int)disc->dis->arr[0][3], -1);
	test("dis MAD", r++, (int)disc->dis->arr[0][4], 1);
	test("dis MAD", r++, (int)disc->dis->arr[1][0], 0);
	test("dis MAD", r++, (int)disc->dis->arr[1][1], -1);
	test("dis MAD", r++, (int)disc->dis->arr[1][2], -1);
	test("dis MAD", r++, (int)disc->dis->arr[1][3], -1);
	test("dis MAD", r++, (int)disc->dis->arr[1][4], 1);
	delete disc;

	// 1   4.5  5  5.2  12
	// NA  4.5  5  5.2  12
	d->arr[0][0]=1; d->arr[0][1]=(float)4.5; d->arr[0][2]=5; d->arr[0][3]=(float)5.2; d->arr[0][4]=12;
	d->arr[1][0]=1; d->arr[1][1]=(float)4.5; d->arr[1][2]=5; d->arr[1][3]=(float)5.2; d->arr[1][4]=12;
	disc = new Discretize(rd, idx, "per", (float).4, (float).6);
	test("dis per", r++, (int)disc->dis->arr[0][0], -1);
	test("dis per", r++, (int)disc->dis->arr[0][1], -1);
	test("dis per", r++, (int)disc->dis->arr[0][2], 0);
	test("dis per", r++, (int)disc->dis->arr[0][3], 1);
	test("dis per", r++, (int)disc->dis->arr[0][4], 1);
	test("dis per", r++, (int)disc->dis->arr[1][0], 0);
	test("dis per", r++, (int)disc->dis->arr[1][1], -1);
	test("dis per", r++, (int)disc->dis->arr[1][2], 0);
	test("dis per", r++, (int)disc->dis->arr[1][3], 1);
	test("dis per", r++, (int)disc->dis->arr[1][4], 1);
	delete disc;

	// 5 5 5 5 5
	d->arr[0][0]=5; d->arr[0][1]=5; d->arr[0][2]=5; d->arr[0][3]=(float)5; d->arr[0][4]=5;
	disc = new Discretize(rd, idx, "per", (float).4, (float).6);
	test("dis per", r++, (int)disc->dis->arr[0][0], 1);
	test("dis per", r++, (int)disc->dis->arr[0][1], 1);
	test("dis per", r++, (int)disc->dis->arr[0][2], 1);
	test("dis per", r++, (int)disc->dis->arr[0][3], 1);
	test("dis per", r++, (int)disc->dis->arr[0][4], 1);
	delete disc;
	delete rd;

	Rawdata* rd2 = new Rawdata();
	Matrix<float>* d2 = new Matrix<float>(5,2, true);
	// 1  NA
	// 5  5
	// 5  5
	// 5  5
	//12 12
	d2->arr[0][0]=1;  d2->arr[0][1]=1;
	d2->arr[1][0]=5;  d2->arr[1][1]=5;
	d2->arr[2][0]=5;  d2->arr[2][1]=5;
	d2->arr[3][0]=5;  d2->arr[3][1]=5;
	d2->arr[4][0]=12; d2->arr[4][1]=12;
	d2->set_missing_value(0,1);
	rd2->data = d2;
	idx.pop_back();
	idx.pop_back();
	idx.pop_back();
	Discretize* disc2 = new Discretize(rd2, idx, "sd_samples", (float).5, (float).5);
	test("dis sd_samples", r++, (int)disc2->dis->arr[0][0], -1);
	test("dis sd_samples", r++, (int)disc2->dis->arr[1][0], 0);
	test("dis sd_samples", r++, (int)disc2->dis->arr[2][0], 0);
	test("dis sd_samples", r++, (int)disc2->dis->arr[3][0], 0);
	test("dis sd_samples", r++, (int)disc2->dis->arr[4][0], 1);
	test("dis sd_samples", r++, (int)disc2->dis->arr[0][1], 0);
	test("dis sd_samples", r++, (int)disc2->dis->arr[1][1], -1);
	test("dis sd_samples", r++, (int)disc2->dis->arr[2][1], -1);
	test("dis sd_samples", r++, (int)disc2->dis->arr[3][1], -1);
	test("dis sd_samples", r++, (int)disc2->dis->arr[4][1], 1);

}

void run_ClassifierDataset( std::string fn_dir){
	std::cout << "Unit Test: ClassifierDataset\n";
	int r=1;
	std::string tn = "ClassifierDataset";

	ClassMinerOptions* cmo = new ClassMinerOptions();
	cmo->file_name_dis = fn_dir + "/test_dataset.txt";
	cmo->file_name_sa = fn_dir + "/test_sample_attributes.txt";
	cmo->file_name_ga = fn_dir + "/test_gene_attributes.txt";
	cmo->class_a = "PAVM=0";
	cmo->class_b = "PAVM=1";

	Attributes* sa = new Attributes("NA");
	Attributes* ga = new Attributes("NA");
	sa->load(cmo->file_name_sa);
	ga->load(cmo->file_name_ga);
	cmo->discretization="none";
	ClassifierDataset* D = new ClassifierDataset();
	try{ D->load(sa, ga, cmo); }
	catch(std::string msg){
		std::cout << "ERROR: " << msg << "\n";
		exit(0);
	}
	test( tn, r++, D->features->arr[2][0], 0);	//1
	test( tn, r++, D->features->rows(), 5);		//2
	test( tn, r++, D->features->cols(), 8);		//3
	std::vector<int>* idx = new std::vector<int>();
	D->samples_with_gene_value("GENE1", 2, idx);
	test( tn, r++, (int)idx->size(), 3);	//4
	test( tn, r++, idx->at(0), 5);		//5
	test( tn, r++, idx->at(1), 6);		//6
	test( tn, r++, idx->at(2), 7);		//7
	D->samples_with_gene_value("GENE3", 1, idx);
	test( tn, r++, (int)idx->size(), 0);	//8
	delete D;

	D = new ClassifierDataset();
	cmo->gene_limit = "Allowed=1";
	try{
		D->load(sa, ga, cmo);
	}
	catch(std::string msg){
		std::cout << "ERROR: " << msg << "\n";
		exit(0);
	}
	test( tn, r++, D->features->rows(), 4);		//9
	test( tn, r++, (int)D->identifiers->size(), 2);	//10
	test( tn, r++, D->identifiers->at(0).compare("GENE1"), 0); // 11
	test( tn, r++, D->identifiers->at(1).compare("GENE2"), 0); // 12
	delete D;

	D = new ClassifierDataset();
	cmo->gene_limit = "";
	cmo->class_a = "gene:GENE2=1";
	cmo->class_b = "gene:GENE3=2";
	try{
		D->load(sa, ga, cmo);
	}
	catch(std::string msg){
		std::cout << "ERROR: " << msg << "\n";
		exit(0);
	}
	test( tn, r++, D->features->rows(), 2);		// 13
	test( tn, r++, (int)D->a_idx->size(), 5); // 14
	test( tn, r++, (int)D->b_idx->size(), 3); // 15
	delete D;

	D = new ClassifierDataset();
	cmo->gene_limit = "";
	cmo->class_a = "PAVM=1,gene:GENE1=1";
	cmo->class_b = "PAVM=1,gene:GENE1=2";
	try{
		D->load(sa, ga, cmo);
	}
	catch(std::string msg){
		std::cout << "ERROR: " << msg << "\n";
		exit(0);
	}
	test( tn, r++, (int)D->a_idx->size(), 1); // 16
	test( tn, r++, (int)D->b_idx->size(), 3); // 17
}


void run_Dataset(std::string fn_data_dir){
	std::cout << "Unit Test: Dataset\n";
	int r=1;
	std::string tn = "Dataset";
	ClassMinerOptions* cmo = new ClassMinerOptions();
	cmo->file_name_dis = fn_data_dir + "/test_dataset.txt";
	cmo->file_name_sa = fn_data_dir + "/test_sample_attributes.txt";
	cmo->file_name_ga = fn_data_dir + "/test_gene_attributes.txt";
	cmo->class_a = "PAVM=0";
	cmo->class_b = "PAVM=1";
	Attributes* sa = new Attributes("NA");
	Attributes* ga = new Attributes("NA");
	sa->load(cmo->file_name_sa);
	ga->load(cmo->file_name_ga);
	cmo->discretization="none";
	Dataset* D = new Dataset();
	try{
		D->load(sa, ga, cmo);
	}
	catch(std::string msg){
		std::cout << "ERROR: " << msg << "\n";
		exit(0);
	}
	test( tn, r++, D->raw_data->data->rows(), 3);
	test( tn, r++, D->raw_data->data->cols(), 8);
	test( tn, r++, D->raw_data->data->has_missing_values, true);
	test( tn, r++, (int)D->identifiers->size(), 3);
	test( tn, r++, D->identifiers->at(0).compare("GENE1"), 0);
	test( tn, r++, D->identifiers->at(1).compare("GENE2"), 0);
	test( tn, r++, D->a_is_target, true);
	test( tn, r++, D->b_is_target, true);
	test( tn, r++, (int)D->a_idx->size(), 4);
	test( tn, r++, (int)D->b_idx->size(), 4);
	delete D;
	cmo->file_name_dis = fn_data_dir + "/test_dataset_for_sorting.txt";
	D = new Dataset();
	try{ D->load(sa, ga, cmo); }
	catch(std::string msg){
		std::cout << "ERROR: " << msg << "\n";
		exit(0);
	}
	std::vector< std::vector<double> *> values;
	std::vector<std::string> genes, labels, probes;
	std::string sample_limit("IDENTIFIER!NULL");
	std::string group_by("");
	bool sort_by_value = false;
	probes.push_back(std::string("GENE1"));
	std::string label_attribute("IDENTIFIER");
	D->extract(values, labels, genes, probes, sample_limit, group_by, sort_by_value, label_attribute );
	// should return -99999,7,6,5,4,3,2,1
	test( tn, r++, (int)values.size(), 1); // return one value
	test( tn, r++, (int)values.at(0)->size(), 8); // with eight items
	test( tn, r++, (int)values.at(0)->at(0), -99999);
	test( tn, r++, (int)values.at(0)->at(1), 7);
	test( tn, r++, (int)values.at(0)->at(2), 6);
	test( tn, r++, (int)values.at(0)->at(3), 5);
	test( tn, r++, (int)values.at(0)->at(4), 4);
	test( tn, r++, (int)values.at(0)->at(5), 3);
	test( tn, r++, (int)values.at(0)->at(6), 2);
	test( tn, r++, (int)values.at(0)->at(7), 1);
	test( tn, r++, (int)genes.size(), 1);
	test( tn, r++, std::string("Foo").compare(genes.at(0)), 0);
	sort_by_value=true;
	D->extract(values, labels, genes, probes, sample_limit, group_by, sort_by_value, label_attribute );
		// should return -99999,1,2,3,4,5,6,7
	test( tn, r++, (int)values.size(), 1); // return one value
	test( tn, r++, (int)values.at(0)->size(), 8); // with eight items
	test( tn, r++, (int)values.at(0)->at(0), -99999);
	test( tn, r++, (int)values.at(0)->at(1), 1);
	test( tn, r++, (int)values.at(0)->at(2), 2);
	test( tn, r++, (int)values.at(0)->at(3), 3);
	test( tn, r++, (int)values.at(0)->at(4), 4);
	test( tn, r++, (int)values.at(0)->at(5), 5);
	test( tn, r++, (int)values.at(0)->at(6), 6);
	test( tn, r++, (int)values.at(0)->at(7), 7);
	sort_by_value=false;
	group_by = std::string("ERP");
	D->extract(values, labels, genes, probes, sample_limit, group_by, sort_by_value, label_attribute );
	// should return -99999,6,4,2,7,5,3,1
	test( tn, r++, (int)values.size(), 1); // return one value
	test( tn, r++, (int)values.at(0)->size(), 8); // with eight items
	test( tn, r++, (int)values.at(0)->at(0), -99999);
	test( tn, r++, (int)values.at(0)->at(1), 6);
	test( tn, r++, (int)values.at(0)->at(2), 4);
	test( tn, r++, (int)values.at(0)->at(3), 2);
	test( tn, r++, (int)values.at(0)->at(4), 7);
	test( tn, r++, (int)values.at(0)->at(5), 5);
	test( tn, r++, (int)values.at(0)->at(6), 3);
	test( tn, r++, (int)values.at(0)->at(7), 1);

	sort_by_value=true;
	group_by = std::string("PAVM");
	try{D->extract(values, labels, genes, probes, sample_limit, group_by, sort_by_value, label_attribute ); }
	catch(std::string err){
		std::cout << "Threw " << err << "\n";
		return;
	}
	// should return -99999,5,6,7,1,2,3,4
	test( tn, r++, (int)values.size(), 1); // return one value
	test( tn, r++, (int)values.at(0)->size(), 8); // with eight items
	test( tn, r++, (int)values.at(0)->at(0), -99999);
	test( tn, r++, (int)values.at(0)->at(1), 5);
	test( tn, r++, (int)values.at(0)->at(2), 6);
	test( tn, r++, (int)values.at(0)->at(3), 7);
	test( tn, r++, (int)values.at(0)->at(4), 1);
	test( tn, r++, (int)values.at(0)->at(5), 2);
	test( tn, r++, (int)values.at(0)->at(6), 3);
	test( tn, r++, (int)values.at(0)->at(7), 4);
    
}

void run_Matrix(){
	std::cout << "Unit Test: Matrix\n";
	int r=1;

	Matrix<int> m = Matrix<int>();
	std::string tn = "Matrix";
	test( std::string("M rows"), r++, m.rows(), 0);
	test( std::string("M cols"), r++, m.cols(), 0);
	test( std::string("M arr"), r++, m.arr==NULL, true);

	Matrix<int> a = Matrix<int>(2,2);
	a.arr[0][0]=12;
	test( std::string("M rows"), r++, a.rows(), 2);
	test( std::string("M cols"), r++, a.cols(), 2);
	m.resize(2,2);
	test( std::string("M resize"), r++, m.rows(), 2);
	test( std::string("M resize"), r++, m.cols(), 2);

	//std::vector<float>*  trimmed_means = new std::vector<float>();
	std::vector<int>* cols = new std::vector<int>();
	cols->push_back(0);
	cols->push_back(1);
	cols->push_back(2);

	a.arr[0][0]=2;	a.arr[0][1]=3;
	a.arr[1][0]=4;	a.arr[1][1]=5;
	std::vector<int>* results = new std::vector<int>();
	a.colSums(results);
	test( std::string("M colSums"), r++, results->at(0), 6);
	test( std::string("M colSums"), r++, results->at(1), 8);
	a.rowSums(results);
	test( std::string("M rowSums"), r++, results->at(0), 5);
	test( std::string("M rowSums"), r++, results->at(1), 9);
	Matrix<int>* X = new Matrix<int>();
	a.transpose(X);
	test( std::string("M transpose"), r++, X->arr[0][0], 2); test( std::string("M transpose"), r++, X->arr[0][1], 4);
	test( std::string("M transpose"), r++, X->arr[1][0], 3);  test( std::string("M transpose"), r++, X->arr[1][1], 5);

	Matrix<int>* xx_miss = new Matrix<int>(1,3);
	xx_miss->arr[0][0]=1;      xx_miss->arr[0][1]=2;      xx_miss->arr[0][2]=3;
	xx_miss->set_missing_value(0,2);
	xx_miss->transpose(X);
	test( std::string("M transpose missing"), r++, X->rows(), 3);
	test( std::string("M transpose missing"), r++, X->cols(), 1);
	test( std::string("M transpose missing"), r++, X->is_missing(2,0), true);


	Matrix<int>* tt = new Matrix<int>(2,2);
	tt->arr[0][0] = 1;	tt->arr[0][1] = 2;
	tt->arr[1][0] = 3;	tt->arr[1][1] = 4;
	tt->add(5);
	test( std::string("M add"), r++, tt->arr[0][0], 6); test( std::string("M add"), r++, tt->arr[0][1], 7);
	test( std::string("M add"), r++, tt->arr[1][0], 8);	test( std::string("M add"), r++, tt->arr[1][1], 9);
	tt->add(tt);
	test( std::string("M add matrix"), r++, tt->arr[0][0], 12); test( std::string("M add matrix"), r++, tt->arr[0][1], 14);
	test( std::string("M add matrix"), r++, tt->arr[1][0], 16); test( std::string("M add matrix"), r++, tt->arr[1][1], 18);

	// test whether missing value inb tt_miss will propigate into tt after matrix add
	Matrix<int>* tt_miss = new Matrix<int>(2,2, true);
	tt_miss->arr[0][0] = 1;	tt_miss->arr[0][1] = 1; tt_miss->arr[1][0] = 1;
	tt_miss->set_missing_value(1,1);
	tt->add(tt_miss);
	test( std::string("M add matrix missing"), r++, tt->arr[0][0], 13); test( std::string("M add matrix missing"), r++, tt->arr[0][1], 15);
	test( std::string("M add matrix missing"), r++, tt->arr[1][0], 17); test( std::string("M add matrix missing"), r++, tt->is_missing(1,1),true);

	Matrix<int>* t2 = new Matrix<int>(2,2);
	t2->arr[0][0] = 1;	t2->arr[0][1] = 2;
	t2->arr[1][0] = 3;	t2->arr[1][1] = 4;
	Matrix<int>* gtlt = new Matrix<int>(2,2);
	t2->where_GT( gtlt, 2 );
	test( std::string("M where_GT"), r++, gtlt->arr[0][0], 0); test( std::string("M where_GT"), r++, gtlt->arr[0][1], 0);
	test( std::string("M where_GT"), r++, gtlt->arr[1][0], 1); test( std::string("M where_GT"), r++, gtlt->arr[1][1], 1);
	t2->where_LT( gtlt, 3 );
	test( std::string("M where_LT"), r++, gtlt->arr[0][0], 1); test( std::string("M where_LT"), r++, gtlt->arr[0][1], 1);
	test( std::string("M where_LT"), r++, gtlt->arr[1][0], 0); test( std::string("M where_LT"), r++, gtlt->arr[1][1], 0);
	t2->mult(2);
	test( std::string("M mult"), r++, t2->arr[0][0], 2); test( std::string("M mult"), r++, t2->arr[0][1], 4);
	test( std::string("M mult"), r++, t2->arr[1][0], 6); test( std::string("M mult"), r++, t2->arr[1][1], 8);
	t2->set_missing_value(0,0);
	t2->where_GT( gtlt, 2 );
	test( std::string("M where_GT missing"), r++, gtlt->is_missing(0,0), true);
	t2->where_LT( gtlt, 2 );
	test( std::string("M where_LT missing"), r++, gtlt->is_missing(0,0), true);
	t2->where_EQ( gtlt, 2 );
	test( std::string("M where_EQ missing"), r++, gtlt->is_missing(0,0), true);
	delete t2;
	delete tt;
	tt = new Matrix<int>(2,2);
	tt->arr[0][0]=1;	tt->arr[0][1]=2;
	tt->arr[1][0]=3;	tt->arr[1][1]=4;
	tt->where_EQ( gtlt, 1 );
	test( std::string("M where_EQ"), r++, gtlt->arr[0][0], 1); // 32
	test( std::string("M where_EQ"), r++, gtlt->arr[0][1], 0);
	test( std::string("M where_EQ"), r++, gtlt->arr[1][0], 0);
	test( std::string("M where_EQ"), r++, gtlt->arr[1][1], 0);
	test( std::string("M sum"), r++, (int)tt->sum(), 10);
	t2 = new Matrix<int>(2,2);
	t2->arr[0][0]=5;	t2->arr[0][1]=6;
	t2->arr[1][0]=7;	t2->arr[1][1]=8;
	tt->mm(t2, gtlt);
	test( std::string("M mm"), r++, gtlt->arr[0][0], 19);  test( std::string("M mm"), r++, gtlt->arr[0][1], 22);
	test( std::string("M mm"), r++, gtlt->arr[1][0], 43);  test( std::string("M mm"), r++, gtlt->arr[1][1], 50);
	t2->set_missing_value(0,0);
	test( std::string("M sum missing"), r++, (int)t2->sum(), 21);
	delete t2;
	t2 = new Matrix<int>(2,2);
	t2->arr[0][0]=5;	t2->arr[0][1]=6;
	t2->arr[1][0]=7;	t2->arr[1][1]=8;
	int x, y;
	int minval = t2->min(x,y);
	test( std::string("M min"), r++, minval, 5);
	test( std::string("M min"), r++, x, 0);
	test( std::string("M min"), r++, y, 0);
	Matrix<int>* fromV = new Matrix<int>(cols);
	test( std::string("M Matrix(vector)"), r++, fromV->rows(), 1);
	test( std::string("M Matrix(vector)"), r++, fromV->cols(), 3);
	test( std::string("M Matrix(vector)"), r++, fromV->arr[0][0], 0);
	test( std::string("M Matrix(vector)"), r++, fromV->arr[0][1], 1);
	test( std::string("M Matrix(vector)"), r++, fromV->arr[0][2], 2);
	t2->arr[0][0]=2;	t2->arr[0][1]=4;
	t2->arr[1][0]=6;	t2->arr[1][1]=8;
	t2->div(2);

	xx_miss = new Matrix<int>(2,2, true);
	xx_miss->arr[0][0]=1; xx_miss->arr[0][1]=2; xx_miss->arr[1][0]=3; xx_miss->arr[1][1]=-1;
	xx_miss->set_missing_value(1,1);
	minval = xx_miss->min(x,y);
	test( std::string("M min missing"), r++, minval, 1);
	test( std::string("M min missing"), r++, x, 0);
	test( std::string("M min missing"), r++, y, 0);
	delete xx_miss;

	test( std::string("M div"), r++, t2->arr[0][0], 1);	test( std::string("M div"), r++, t2->arr[0][1], 2);
	test( std::string("M div"), r++, t2->arr[1][0], 3);	test( std::string("M div"), r++, t2->arr[1][1], 4);
	Matrix<int>* xx = new Matrix<int>(2,2);
	xx_miss = new Matrix<int>(2,2, true);
	xx->arr[0][0]=1;      xx->arr[0][1]=2;      xx->arr[1][0]=3;      xx->arr[1][1]=4;
	xx_miss->arr[0][0]=1; xx_miss->arr[0][1]=2; xx_miss->arr[1][0]=3; xx_miss->set_missing_value(1,1);
	xx->div(xx_miss);

	test( std::string("M div missing"), r++, xx->arr[0][0], 1);  test( std::string("M div missing"), r++, xx->arr[0][1], 1);
	test( std::string("M div missing"), r++, xx->arr[1][0], 1);  test( std::string("M div missing"), r++, xx->is_missing(1,1), true);
	delete xx;
	delete xx_miss;
	t2->clone(tt);
	test( std::string("M clone"), r++, tt->arr[0][0], 1);  test( std::string("M clone"), r++, tt->arr[0][1], 2);
	test( std::string("M clone"), r++, tt->arr[1][0], 3);  test( std::string("M clone"), r++, tt->arr[1][1], 4);
	Matrix<double>* cloneme = new Matrix<double>(2,2,true,1);
	cloneme->set_missing_value(0,0);
	Matrix<double>* cloned = new Matrix<double>();
	cloneme->clone(cloned);
	test( std::string("M clone with missing data"), r++, cloned->has_missing_values, true);
	test( std::string("M clone with missing data"), r++, cloned->mask[0][0], true);
	test( std::string("M clone with missing data"), r++, cloned->mask[0][1], false);
	test( std::string("M clone with missing data"), r++, (bool)cloned->mask[1][0], false);
	test( std::string("M clone with missing data"), r++, cloned->mask[1][1], false);
	test( std::string("M clone with missing data"), r++, cloned->row_has_missing_value(0), true);

	t2->add(tt);
	test( std::string("M add"), r++, t2->arr[0][0], 2);  test( std::string("M add"), r++, t2->arr[0][1], 4);
	test( std::string("M add"), r++, t2->arr[1][0], 6);  test( std::string("M add"), r++, t2->arr[1][1], 8);
	t2->mult(tt);
	test( std::string("M mult"), r++, t2->arr[0][0], 2);  test( std::string("M mult"), r++, t2->arr[0][1], 8);
	test( std::string("M mult"), r++, t2->arr[1][0], 18); test( std::string("M mult"), r++, t2->arr[1][1], 32);

	xx = new Matrix<int>(2,2);
	xx_miss = new Matrix<int>(2,2, true);
	xx->arr[0][0]=1;      xx->arr[0][1]=2;      xx->arr[1][0]=3;      xx->arr[1][1]=4;
	xx_miss->arr[0][0]=1; xx_miss->arr[0][1]=2; xx_miss->arr[1][0]=3; xx_miss->set_missing_value(1,1);
	xx->mult(xx_miss);
	test( std::string("M mult missing"), r++, xx->arr[0][0], 1);  test( std::string("M mult missing"), r++, xx->arr[0][1], 4);
	test( std::string("M mult missing"), r++, xx->arr[1][0], 9);  test( std::string("M mult missing"), r++, xx->is_missing(1,1), true);
	delete xx;
	delete xx_miss;

	xx = new Matrix<int>(2,2);
	xx_miss = new Matrix<int>(2,2, true);
	xx->arr[0][0]=1;      xx->arr[0][1]=2;      xx->arr[1][0]=3;      xx->arr[1][1]=4;
	xx_miss->arr[0][0]=1; xx_miss->arr[0][1]=2; xx_miss->arr[1][0]=3; xx_miss->set_missing_value(1,1);
	test( std::string("M mult_sum"), r++, xx->mult_sum(xx), 30);
	test( std::string("M mult_sum missing"), r++, xx->mult_sum(xx_miss), 14);
	delete xx;
	delete xx_miss;

	t2->sub(tt);
	test( std::string("M sub"), r++, t2->arr[0][0], 1);	test( std::string("M sub"), r++, t2->arr[0][1], 6);
	test( std::string("M sub"), r++, t2->arr[1][0], 15); test( std::string("M sub"), r++, t2->arr[1][1], 28);
	xx = new Matrix<int>(2,2);
	xx_miss = new Matrix<int>(2,2, true);
	xx->arr[0][0]=1;      xx->arr[0][1]=2;      xx->arr[1][0]=3;      xx->arr[1][1]=4;
	xx_miss->arr[0][0]=1; xx_miss->arr[0][1]=2; xx_miss->arr[1][0]=3; xx_miss->set_missing_value(1,1);
	xx->sub(xx_miss);
	test( std::string("M sub missing"), r++, xx->arr[0][0], 0);  test( std::string("M sub missing"), r++, xx->arr[0][1], 0);
	test( std::string("M sub missing"), r++, xx->arr[1][0], 0);  test( std::string("M sub missing"), r++, xx->is_missing(1,1), true);
	delete xx;
	delete xx_miss;

	t2->arr[0][0]=2;	t2->arr[0][1]=4;
	t2->arr[1][0]=6;	t2->arr[1][1]=8;
	tt->arr[0][0]=3;	tt->arr[0][1]=3;
	tt->arr[1][0]=3;	tt->arr[1][1]=3;
	test( std::string("M mult_sum"), r++, (int)t2->mult_sum(tt), 6+12+18+24);
	t2->slice_row(0, tt);
	test( std::string("M slice_row"), r++, tt->rows(), 1);
	test( std::string("M slice_row"), r++, tt->cols(), 2);
	test( std::string("M slice_row"), r++, tt->arr[0][0], 2);	test( std::string("M slice_row"), r++, tt->arr[0][1], 4);
	tt->power(2);
	test( std::string("M power"), r++, tt->arr[0][0], 4);	test( std::string("M power"), r++, tt->arr[0][1], 16);
	Matrix<double>* ee = new Matrix<double>(1,1,false,1);
	ee->exponent();
	test( std::string("M exponent"), r++, (int)((ee->arr[0][0])*1000), 2718);

	xx = new Matrix<int>(2,5);
	xx->arr[0][0] = 1;
	xx->arr[0][1] = 2;
	xx->arr[0][2] = 3;
	xx->arr[0][3] = 4;
	xx->arr[0][4] = 5;
	xx->arr[1][0] = 6;
	xx->arr[1][1] = 7;
	xx->arr[1][2] = 8;
	xx->arr[1][3] = 9;
	xx->arr[1][4] = 10;
	std::vector<int> valid_cols;
	valid_cols.push_back(0);
	valid_cols.push_back(1);
	valid_cols.push_back(2);
	Matrix<int>* shuf = new Matrix<int>(xx->rows(), xx->cols(), xx->has_missing_values );
	xx->clone(shuf);
	shuf->shuffle_by_rows(valid_cols);
	test( std::string("shuffle_rows"), r++, shuf->rows(), 2);
	test( std::string("shuffle_rows"), r++, shuf->cols(), 5);
	test( std::string("shuffle_rows valid_cols"), r++, shuf->arr[0][3], 4);
	test( std::string("shuffle_rows valid_cols"), r++, shuf->arr[0][4], 5);
	test( std::string("shuffle_rows valid_cols"), r++, shuf->arr[1][3], 9);
	shuf->shuffle_by_rows(valid_cols);
	test( std::string("shuffle_rows valid_cols"), r++, shuf->arr[1][4], 10);

	Matrix<int>* LBind = new Matrix<int>(2,2);
	LBind->arr[0][0] = 1;
	LBind->arr[0][1] = 2;
	LBind->arr[1][0] = 3;
	LBind->arr[1][1] = 4;
	Matrix<int>* RBind = new Matrix<int>(2,2);
	RBind->arr[0][0] = 5;
	RBind->arr[0][1] = 6;
	RBind->arr[1][0] = 7;
	RBind->arr[1][1] = 8;
	Matrix<int>* CBind = new Matrix<int>();
	LBind->cbind(RBind, CBind);
	test( std::string("cbind ncols"), r++, CBind->cols(), 4);
	test( std::string("cbind nrows"), r++, CBind->rows(), 2);
	test( std::string("cbind vals 0 0"), r++, CBind->arr[0][0], 1);
	test( std::string("cbind vals 0 1"), r++, CBind->arr[0][1], 2);
	test( std::string("cbind vals 0 2"), r++, CBind->arr[0][2], 5);
	test( std::string("cbind vals 1 3"), r++, CBind->arr[0][3], 6);
	test( std::string("cbind vals 1 0"), r++, CBind->arr[1][0], 3);
	test( std::string("cbind vals 1 1"), r++, CBind->arr[1][1], 4);
	test( std::string("cbind vals 1 2"), r++, CBind->arr[1][2], 7);
	test( std::string("cbind vals 1 3"), r++, CBind->arr[1][3], 8);

	delete LBind;
	delete RBind;
	delete CBind;

	delete cols;
	delete results;
	delete X;
	delete tt;
	delete t2;
	delete gtlt;
	delete cloneme;
	delete cloned;
	delete ee;
	delete shuf;
	Matrix<int> foo(2,2);
	foo.arr[0][0]=1; foo.arr[0][1]=1; foo.arr[1][0]=1; foo.arr[1][1]=1;
	foo.set_missing_value(0,0);
	foo.set_all_values_present();
	test( std::string("set_all_values_present"), r++, foo.is_missing(0,0), false);
    Matrix<int> bart(2,3);
    bart.arr[0][0]=1; bart.arr[0][1]=1; bart.arr[0][2]=1;
    bart.arr[1][0]=0; bart.arr[1][1]=0; bart.arr[1][2]=1; 
    bart.set_missing_value(0,0);
    std::vector<int> idx;
    idx.push_back(0); idx.push_back(1); idx.push_back(2);
    bart.eliminate_infrequent_rowvalues(idx, 2);
    test( std::string("eliminate_infrequent"), r++, bart.is_missing(1, 2), true);
    
    
    Matrix<int>* x1 = new Matrix<int>(3,4);
    Matrix<int>* y1 = new Matrix<int>(3,4);
	x1->arr[0][0] = 1;
	x1->arr[0][1] = 2;
	x1->arr[0][2] = 3;
	x1->arr[0][3] = 4;    // 1  2  3  4
	x1->arr[1][0] = 5;    // 5  6  7  8
	x1->arr[1][1] = 6;    // 9 10 11 12
	x1->arr[1][2] = 7;
	x1->arr[1][3] = 8;
	x1->arr[2][0] = 9;
	x1->arr[2][1] = 10;
	x1->arr[2][2] = 11;
	x1->arr[2][3] = 12;
    x1->clone(y1);
    std::vector<int> to_keep;
    to_keep.push_back(2);
    to_keep.push_back(3);
    x1->restrict_by_cols(to_keep);
    //  3  4
    //  7  8
    // 11 12
    test( std::string("restrict_by_cols number of cols"), r++, x1->cols(), 2);
    test( std::string("restrict_by_cols number of rows"), r++, x1->rows(), 3);
    test( std::string("restrict_by_cols 0,0"), r++, x1->arr[0][0], 3);
    test( std::string("restrict_by_cols 0,1"), r++, x1->arr[0][1], 4);
    test( std::string("restrict_by_cols 1,0"), r++, x1->arr[1][0], 7);
    test( std::string("restrict_by_cols 1,1"), r++, x1->arr[1][1], 8);
    
    to_keep.pop_back();
    to_keep.push_back(0); // push first row into second position
    // 9 10 11 12
    // 1  2  3  4
    y1->restrict_by_rows(to_keep);
    test( std::string("restrict_by_rows number of cols"), r++, y1->cols(), 4);
    test( std::string("restrict_by_rows number of rows"), r++, y1->rows(), 2);
    test( std::string("restrict_by_rows 0,0"), r++, y1->arr[0][0], 9);
    test( std::string("restrict_by_rows 1,0"), r++, y1->arr[1][0], 1);
    test( std::string("restrict_by_rows 0,0"), r++, y1->arr[0][1], 10);
    test( std::string("restrict_by_rows 0,1"), r++, y1->arr[1][1], 2);
    
    delete x1;
    delete y1;
}


void run_ProbeSets(std::string base_dir){
	std::cout << "Unit Test: ProbeSets\n";
	ProbeSets probe_sets;
	Attributes* ga = new Attributes(std::string("NA"));
	std::string fn_ga = base_dir + "/test_correlation_gene_attributes.txt";
	ga->load(fn_ga);
	std::vector<std::string> probes;
	probes.push_back(std::string("GENE1"));
	probes.push_back(std::string("GENE2"));
	int min_probe_set_size=0;
	probe_sets.add_from_probelist(std::string("set1"), ga, probes, min_probe_set_size);

	ProbeSet* PS=NULL;
	PS = probe_sets.next_probeset();
	test(std::string("ProbeSets next_probeset"), 1, PS==NULL, false);
	test(std::string("ProbeSets next_probeset title"), 2, PS->title.compare(std::string("set1")), 0);
	test(std::string("ProbeSets next_probeset correct size"), 3, int(PS->probelist.size()), 2);
	test(std::string("ProbeSets next_probeset idx0"), 4, PS->probelist.at(0).compare(std::string("GENE1")) , 0);
	test(std::string("ProbeSets next_probeset idx1"), 5, PS->probelist.at(1).compare(std::string("GENE2")) , 0);


	PS = probe_sets.next_probeset();
	test(std::string("ProbeSets next_probeset off end"), 6, PS==NULL, true); // should now be out of probesets
	probe_sets.reset_iterator();
	PS = probe_sets.next_probeset();
	test(std::string("ProbeSets next_probeset after reset_iterator"), 7, PS==NULL, false);
	probe_sets.clear();
	PS = probe_sets.next_probeset();
	test(std::string("ProbeSets next_probeset after clear"), 8, PS==NULL, true);

	std::vector<std::string> genelist;
	genelist.push_back(std::string("Foo"));
	genelist.push_back(std::string("Bar"));
	genelist.push_back(std::string("Gene that is not there"));
	probe_sets.add_from_genelist(std::string("set2"), ga, genelist, min_probe_set_size);
	PS = probe_sets.next_probeset();
	test(std::string("ProbeSets next_probeset from genes"), 9, PS==NULL, false);
	test(std::string("ProbeSets next_probeset from genes title"), 10, PS->title.compare(std::string("set2")), 0);
	test(std::string("ProbeSets next_probeset from genes correct size"), 11, int(PS->probelist.size()), 2);
	test(std::string("ProbeSets next_probeset idx0"), 4, PS->probelist.at(0).compare(std::string("GENE1")) , 0);
	test(std::string("ProbeSets next_probeset idx1"), 5, PS->probelist.at(1).compare(std::string("GENE2")) , 0);
	probe_sets.clear();
	probe_sets.add_from_genelist(std::string("set2"), ga, genelist, 5);
	PS = probe_sets.next_probeset();
	test(std::string("ProbeSets next_probeset from list that was too small"), 14, PS==NULL, true);
}

void run_spear(std::string base_dir){


	std::cout << "Unit Test: Spear\n";
	std::string fn_expr = base_dir + "/test_correlation_dataset.txt";
	std::string fn_dc = base_dir + "/test_differential_correlation_dataset.txt";
	std::string fn_ga = base_dir + "/test_correlation_gene_attributes.txt";
	std::string fn_sa = base_dir + "/test_sample_attributes.txt";
	Spearman* sp = new Spearman();
	Matrix<float>* A = new Matrix<float>(1,6, true);
	Matrix<float>* B = new Matrix<float>(1,6, true);
	double corr_r, corr_rho;

	A->arr[0][0] = 1.0; A->arr[0][1] = 2.0; A->arr[0][2] = 3.0; A->arr[0][3] = 4.0; A->arr[0][4] = 5.0; A->arr[0][5] = 6.0;
	B->arr[0][0] = 1.0; B->arr[0][1] = 2.0; B->arr[0][2] = 3.0; B->arr[0][3] = 4.0; B->arr[0][4] = 5.0; B->arr[0][5] = 200.0;

	corr_r = sp->correlation(A, B, std::string("pearson") );
	corr_rho = sp->correlation(A, B, std::string("spearman") );

	test( std::string("pearson  1,2,3,4,5,6 vs 1,2,3,4,5,200"), 1, (int)(corr_r*1000), 667 );
	test( std::string("spearman 1,2,3,4,5,6 vs 1,2,3,4,5,200"), 2, (int)(corr_rho*1000), 1000 );

	B->arr[0][4] = 6.0;
	B->arr[0][5] = 5.0;
	corr_r = sp->correlation(A, B, std::string("pearson") );
	corr_rho = sp->correlation(A, B, std::string("spearman") );
	test( std::string("pearson  1,2,3,4,5,6 vs 1,2,3,4,6,5"), 3, (int)(corr_r*1000), 942 );
	test( std::string("spearman 1,2,3,4,5,6 vs 1,2,3,4,6,5"), 4, (int)(corr_rho*1000), 942);

	A->set_missing_value(0,4);
	corr_r = sp->correlation(A, B, std::string("pearson") );
	corr_rho = sp->correlation(A, B, std::string("spearman") );
	test( std::string("pearson  1,2,3,4,_,6 vs 1,2,3,4,6,5"), 5, (int)(corr_r*1000), 986 );
	test( std::string("spearman 1,2,3,4,_,6 vs 1,2,3,4,6,5"), 6, (int)(corr_rho*1000), 1000);

	B->set_missing_value(0,0);
	corr_r = sp->correlation(A, B, std::string("pearson") );
	corr_rho = sp->correlation(A, B, std::string("spearman") );
	test( std::string("pearson  1,2,3,4,_,6 vs _,2,3,4,6,5"), 7, (int)(corr_r*1000), 982 );
	test( std::string("spearman 1,2,3,4,_,6 vs _,2,3,4,6,5"), 8, (int)(corr_rho*1000), 1000);
	delete A;
	delete B;

	A = new Matrix<float>(1,8, true);
	B = new Matrix<float>(1,8, true);
	A->arr[0][0] = 1.0; A->arr[0][1] = 2.0; A->arr[0][2] = 3.0; A->arr[0][3] = 4.0; A->arr[0][4] = 5.0; A->arr[0][5] = 6.0; A->arr[0][6] = 7.0; A->arr[0][7] = 8.0;
	B->arr[0][0] = -99999; B->arr[0][1] = -99999; B->arr[0][2] = -99999; B->arr[0][3] = 4.0; B->arr[0][4] = 5.0; B->arr[0][5] = 6.0; B->arr[0][6] = 7.0; B->arr[0][7] = 8.0;
	B->set_missing_value(0,0);
	B->set_missing_value(0,1);
	B->set_missing_value(0,2);
	corr_r = sp->correlation(A, B, std::string("pearson") );
	corr_rho = sp->correlation(A, B, std::string("spearman") );
	test( std::string("pearson  1,2,3,4,5,6,7,8 vs _,_,_,4,5,6,7,8"), 9, (int)(corr_r*1000), 1000 );
	test( std::string("spearman 1,2,3,4,5,6,7,8 vs _,_,_,4,5,6,7,8"), 10, (int)(corr_rho*1000), 1000);
	delete A;
	delete B;
	Matrix<double>* RR = new Matrix<double>(1,5);
	RR->arr[0][0]=0; RR->arr[0][1]=1; RR->arr[0][2]=2; RR->arr[0][3]=3; RR->arr[0][4]=4;
	std::set<int> miss;
	delete RR;
	std::vector<Spear*> spears;
	sp->set_verbose(false);
	sp->set_input_files(fn_expr, fn_sa, fn_ga);
    sp->run();
    sp->results( spears );
	test( std::string("spearman testdata all n spears"), 11, (int)spears.size(), 6);
	test( std::string("spearman testdata all"), 12, spears.at(0)->row_1(), 0);
	test( std::string("spearman testdata all"), 13, spears.at(0)->row_2(), 3);
	test( std::string("spearman testdata all"), 14, int(spears.at(0)->rho_a()),  1);
	test( std::string("spearman testdata all"), 15, spears.at(1)->row_1(), 0);
	test( std::string("spearman testdata all"), 16, spears.at(1)->row_2(), 1);

	for(int i=0; i<(int)spears.size(); i++){
		delete spears.at(i);
	}
	sp->set_percent_required(0.9);
	sp->run();
	sp->results(spears);
	test( std::string("spearman testdata 90% required"), 17, (int)spears.size(), 3);
	test( std::string("spearman testdata 90% required"), 18, spears.at(0)->row_1(), 0);
	test( std::string("spearman testdata 90% required"), 19, spears.at(0)->row_2(), 1);
	for(int i=0; i<(int)spears.size(); i++){
		delete spears.at(i);
	}
	sp->set_corr_abs_a(0.9);
	sp->run();
	sp->results(spears);
	test( std::string("spearman testdata 90% req abs 0.9"), 20, (int)spears.size(), 1);
	test( std::string("spearman testdata 90% req abs 0.9"), 21, spears.at(0)->row_1(), 0);
	test( std::string("spearman testdata 90% req abs 0.9"), 22, spears.at(0)->row_2(), 1);
	for(int i=0; i<(int)spears.size(); i++){
		delete spears.at(i);
	}
	std::vector<std::string> seeds;
	seeds.push_back(std::string("GENE1"));
	sp->set_corr_abs_a(0);
	sp->set_percent_required(0);
	sp->set_include_seed_neighbor_correlations(false);
	sp->set_seeds(seeds);
	sp->set_limit_network_to_seeds(false);
	sp->run();
	sp->results(spears);
	test( std::string("spearman testdata seed"), 23, (int)spears.size(), 3);

	for(int i=0; i<(int)spears.size(); i++){
		delete spears.at(i);
	}
	sp->set_include_seed_neighbor_correlations(true); // fill in neighbor values
	sp->run();
	sp->results(spears);
	test( std::string("spearman testdata seed with neighbors"), 24, (int)spears.size(), 6); // 0-1, 0-2, 1-3, 1-2, 1-3, 2-
	Attributes* ga = new Attributes("NA");
	ga->load(base_dir + "/test_correlation_gene_attributes.txt");
	sp->set_fn_cytoscape(std::string("/cytoscape/cytoscape.exe"));
	sp->set_fn_cytoscape_props(std::string("/cytoscape/props.txt"));
	sp->write_to_cytoscape(0.4,  base_dir+ "/cytotest", ga);
	for(int i=0; i<(int)spears.size(); i++){
		delete spears.at(i);
	}
	seeds.push_back(std::string("GENE2"));
	sp->set_corr_abs_a(0);
	sp->set_percent_required(0);
	sp->set_include_seed_neighbor_correlations(false);
	sp->set_seeds(seeds);
	sp->run();
	sp->results(spears);

	for(int i=0; i<(int)spears.size(); i++){
		delete spears.at(i);
	}
	delete sp;

	/*
	Spearman* spdc = new Spearman();
	spdc->set_input_files(fn_expr, fn_sa, fn_ga);
	spdc->set_limit_a(std::string("PAVM=0"));
	spdc->set_limit_b(std::string("PAVM=1"));
	ProbeSets ps;
	std::vector<std::string> probe_idx;
	probe_idx.push_back(std::string("GENE1"));
	probe_idx.push_back(std::string("GENE2"));
	probe_idx.push_back(std::string("GENE3"));
	std::vector<DifferentalCorrelationResult*> dc_results;
	ps.add_from_probelist(std::string("set1"), ga, probe_idx, 0);
	spdc->set_n_permutations(0);
	spdc->calculate_differential_correlation_by_probesets(ps, dc_results);
	test( std::string("spearman differential_correlation n results"), 24, int(dc_results.size()), 1);
	test( std::string("spearman differential_correlation mean_diff"), 25, int(dc_results.at(0)->meanDiff*1000), -66);
	test( std::string("spearman differential_correlation p"), 26, int(dc_results.at(0)->pval), 1);
	delete spdc;
	spdc = new Spearman();
	spdc->set_input_files(fn_dc, fn_sa, fn_ga); // load different expression data
	spdc->set_limit_a(std::string("PAVM=0"));
	spdc->set_limit_b(std::string("PAVM=1"));
	spdc->set_n_permutations(1000);
	spdc->calculate_differential_correlation_by_probesets(ps, dc_results);
	test( std::string("spearman differential_correlation p"), 27, int(dc_results.at(0)->pval), 0);
	delete ga;
	delete spdc;
*/
	Spearman S;
	Matrix<double>* ranks1 = new Matrix<double>(1,5);
	Matrix<double>* ranks2 = new Matrix<double>(1,5);
	ranks1->arr[0][0] = 1; ranks2->arr[0][0] = 1;
	ranks1->arr[0][1] = 2; ranks2->arr[0][1] = 2;
	ranks1->arr[0][2] = 3; ranks2->arr[0][2] = 3;
	ranks1->arr[0][3] = 4; ranks2->arr[0][3] = 5;
	ranks1->arr[0][4] = 5; ranks2->arr[0][4] = 4;
	std::vector<int> valid_cols;
	valid_cols.push_back(0); valid_cols.push_back(1); valid_cols.push_back(2); valid_cols.push_back(3); valid_cols.push_back(4);
	Matrix<double> obs_a(2,2); obs_a.arr[0][0]=1; obs_a.arr[1][0]=1; obs_a.arr[0][1]=2; obs_a.arr[1][1]=2;
	Matrix<double> obs_b(2,2); obs_b.arr[0][0]=1; obs_b.arr[1][0]=2; obs_b.arr[0][1]=2; obs_b.arr[1][1]=1;
	Matrix<double> perm_a(2,2);
	Matrix<double> perm_b(2,2);
	//S.shuffle_ranks_between_labels(&obs_a, &obs_b, &perm_a, &perm_b );
	std::vector<double>* var = new std::vector<double>();
	var->push_back(0.25);
	var->push_back(0.5);
	var->push_back(0.75);
	var->push_back(1.00);
	var->push_back(1.25);
	double variance = sp->find_var(var);
	test( std::string("spearman variance"), 28, float(variance), float(0.15625));

}

void run_graph(){
	std::cout << "Unit Test: Graph\n";
	Graph* g = new Graph();
	g->add_edge(1,2,0.0);
    std::vector<std::vector<int>*>* results = new std::vector<std::vector<int>*>();

	std::vector<int> ids;
	std::vector<double> weights;
	g->neighbors(1, ids, weights);
	test( std::string("Graph num neighbors"), 1, (int)ids.size(), 1);
	test( std::string("Graph first neighbor"), 2, ids.at(0), 2 ) ;
	g->add_edge(1,3, 1.0);
	g->add_edge(2,3, 2.0);
    g->edges(results);
    test( std::string("Graph number edges"), 3, (int)results->size(), 3 ) ;
    std::vector<int> N;
    g->nodes(N);
    test( std::string("Graph nodes"), 4, (int)N.size(), 3 );
    std::sort(N.begin(), N.end() );
    test( std::string("Graph nodes"), 5, N.at(0), 1 );
    test( std::string("Graph nodes"), 6, N.at(1), 2 );
    test( std::string("Graph nodes"), 7, N.at(2), 3 );
    std::vector<int>* dq = new std::vector<int>();
    std::vector<int>* sq = new std::vector<int>();
    std::vector<int>* nudge = new std::vector<int>();
    dq->push_back(1); dq->push_back(2); dq->push_back(3); dq->push_back(4); dq->push_back(5);
    sq->push_back(1); sq->push_back(3); sq->push_back(4);
    int n_diff = g->n_diff(dq, sq);
    test( std::string("Graph n_diff"), 8, n_diff, 2 ) ;
    g->diff(sq, dq, nudge);
    test( std::string("Graph diff"), 9, (int)nudge->size(), 2 ) ;
    test( std::string("Graph diff"), 10, nudge->at(0), 2 ) ;
    test( std::string("Graph diff"), 11, nudge->at(1), 5 ) ;
    test( std::string("Graph n_diff"), 12, n_diff, 2 ) ;
    sq->clear();
    n_diff = g->n_diff(dq, sq);
    test( std::string("Graph n_diff"), 13, n_diff, 5 ) ;
    g->diff(sq, dq, nudge);
    test( std::string("Graph diff"), 14, (int)nudge->size(), 5 ) ;

    g->neighbors(1, ids, weights);
	test( std::string("Graph num neighbors"), 15, (int)ids.size(), 2);
	std::sort(ids.begin(), ids.end());
	test( std::string("Graph first neighbor"), 16, ids.at(0), 2);
	test( std::string("Graph second neighbor"), 17, ids.at(1), 3);
	double weight;
	test( std::string("Graph has_edge"), 18, g->has_edge(2,3, weight), true);
	test( std::string("Graph correct_weight"), 19, (int)weight, 2);
    g->add_edge(1,4, 1.0);
    g->add_edge(2,4, 1.0);
    g->add_edge(3,4, 1.0);
    test( std::string("Graph order"), 20, g->order(), 4);
    N.clear();
    std::vector<int>* C = new std::vector<int>();
    C->push_back(1);
    C->push_back(2);
    C->push_back(3);
    g->shared_neighbors(C, N);
    test( std::string("Graph shared neighbors"), 21, (int)N.size(), 1);
    test( std::string("Graph shared neighbors"), 22, N.at(0), 4);
    Graph* g1  = new Graph();
    double w;
    std::vector<int> sub_ids;
    sub_ids.push_back(1);
    sub_ids.push_back(2);
    g->subgraph( *g1, sub_ids );
    test( std::string("Graph subgraph"), 23, g1->order(), 2) ;
    test( std::string("Graph subgraph"), 24, g1->has_edge(1,2, w), true) ;

    std::vector< std::vector<int>* >* pool = new std::vector<std::vector<int>*>();
    std::vector<int>* not_redundant = new std::vector<int>();
    not_redundant->push_back(3);
    std::vector<int>* redundant = new std::vector<int>();
    redundant->push_back(1);
    redundant->push_back(2);
    std::vector<int>* target = new std::vector<int>();
    target->push_back(1);
    target->push_back(2);
    target->push_back(4);
    pool->push_back(not_redundant);
    pool->push_back(redundant);
    g->clear_redundant( pool, target);
    test( std::string("Graph clear redundant"), 25, (int)pool->size(), 2);
    test( std::string("Graph clear redundant"), 26, (int)pool->at(0)->size(), 1);
    test( std::string("Graph clear redundant"), 27, (int)pool->at(1)->size(), 0);
    g->remove_vectors_by_size(pool, 1, false);
    test( std::string("Graph remove_empty_vectors"), 28, (int)pool->size(), 1);
    test( std::string("Graph remove_empty_vectors"), 29, (int)pool->at(0)->at(0), 3);
    Graph* gn = new Graph();
    gn->add_edge(1,2);
    gn->add_edge(1,3);
    Graph gt1;
    std::vector<int> out;
    gn->reduce_to_triangles(gt1);
    gt1.nodes(out);
    test( std::string("Graph correct number nodes triangles"), 30, (int)out.size(), 0);
    out.clear();
    gn->add_edge(2,3); gn->add_edge(2,4); gn->add_edge(2,5); gn->add_edge(4,5);
    gn->add_edge(3,6); gn->add_edge(5,7);
    Graph gt;
    gn->reduce_to_triangles(gt);
    gt.nodes(out);
    std::sort(out.begin(), out.end());
    test( std::string("Graph correct number nodes triangles 2"), 31, (int)out.size(), 5);
    test( std::string("Graph correct triangle order 1"), 32, (int)out.at(0), 1 );
    test( std::string("Graph correct triangle order 2"), 33, (int)out.at(4), 5 );
    delete gn;
    gn = new Graph();
    //gn->clear();
    gn->add_edge(1,2);
    gn->add_edge(1,3);
    gn->add_edge(1,4);
    gn->add_edge(2,3);
    gn->add_edge(2,4);
    gn->add_edge(3,4); // quad 1,2,3,4
    gn->add_edge(1,5);
    gn->add_edge(1,6);
    gn->add_edge(5,6); // triangle 1,5,6
    gn->add_edge(5,7);
    gn->add_edge(2,7); // diamond (not quad) 1,2,7,5
    //gt.clear();
    Graph gx;
    gn->reduce_to_quads(gx);
    out.clear();
    gx.nodes(out);
    std::sort(out.begin(), out.end());
    test( std::string("Graph correct number nodes quads"), 34, (int)out.size(), 4);
    test( std::string("Graph correct quad order"), 35, (int)out.at(0), 1);
    gn->set_node_label(1, std::string("test label"));
    test( std::string("Graph get node by label"), 36, gn->get_node_index_by_label(std::string("test label")), 1);
    test( std::string("Graph get label by node"), 37, strcmp( gn->get_node_label_by_index(1).c_str(), "test label"), 0);
    std::vector<QTL*> qtls;
    qtls.push_back(new QTL(std::string("locus_id"), std::string("locus_name"), std::string("probe_id"), std::string("probe_name"), 0.0001));

    Graph G_EQTL;
    G_EQTL.add_qtls(qtls);
    std::vector<int> nodes;
    G_EQTL.nodes(nodes);
    test( std::string("Graph EQTL created"), 38, int(nodes.size()), 2);
    G_EQTL.add_edge(0,2,0.02, NODE_LOCUS, NODE_GENE);
    G_EQTL.nodes(nodes);
    test( std::string("Graph EQTL added edge"), 39, int(nodes.size()), 3);
    test( std::string("Graph EQTL has locus neighbor"), 40, G_EQTL.has_locus_neighbor(1), true);
    test( std::string("Graph EQTL has locus neighbor2"), 41, G_EQTL.has_locus_neighbor(0), false); // 0 is a locus
    test( std::string("Graph EQTL node_type"), 42, G_EQTL.node_type(0), NODE_LOCUS); // 0 is a locus
    test( std::string("Graph EQTL node_type2"), 43, G_EQTL.node_type(2), NODE_GENE); // 2 is a gene
    test( std::string("Graph EQTL get node with label"), 44, G_EQTL.get_node_index_by_label(std::string("locus_id") ), 0 );
    for(int i=0; i<int(qtls.size()); i++)
    	delete qtls.at(i);

    Graph GG;
    GG.add_edge(0,1);
    GG.add_edge(1,2);
    GG.add_edge(3,4);
    GG.remove_edge(0,1);

    test( std::string("Graph remove edge remove"), 45, GG.has_edge(0,1,w), false);
    test( std::string("Graph remove edge preserved adjacent"), 46, GG.has_edge(1,2,w), true); // should preserve this
    GG.set_node_label(3, std::string("My label"));
    GG.remove_edge(3,4);
    test( std::string("Graph remove edge label was removed"), 47, GG.has_node_label(3), false);
}

void run_node(){
	std::cout << "Unit Test: Node\n";
	std::string tn = "Tree";
	Node* root = new Node( Node::ROOT );
	root->add_child(1);
	test( tn, 1, root->has_child(1)!=NULL, true);
	std::vector<int>* path = new std::vector<int>();
	path->push_back(1);
	test( tn, 2, root->has_path(path), true);
	path->push_back(2);
	root->add_path(path);
	test( tn, 3, root->has_path(path), true );

	std::vector<int> * kids = new std::vector<int>();
	root->get_kids(kids);
	test( tn, 4, (int)kids->size(), 1);
	test( tn, 5, kids->at(0), 1);

	Node* x = new Node(Node::ROOT);
	Node* idx_1 = x->add_child(1);
	Node* idx_3 = idx_1->add_child(3);
	idx_3->add_child(5);

	std::vector<int> * set = new std::vector<int>();
	set->push_back(1);
	set->push_back(3);
	set->push_back(5);

	test( tn, 6, x->has_path(set), true);

	std::vector<int> * set2 = new std::vector<int>();
	set2->push_back(1);
	set2->push_back(5);
	set2->push_back(3); // should be same path as set1

	test( tn, 7, x->has_path(set2), true);

	Node* rt = new Node(Node::ROOT);
	std::vector<int>* p1 = new std::vector<int>; p1->push_back(1); p1->push_back(4);
	rt->add_path(p1);
	p1->push_back(2);
	rt->add_path(p1);
	std::vector<int>* p2 = new std::vector<int>;
	p2->push_back(1);
	p2->push_back(2);
	p2->push_back(4);
	test( tn, 8, rt->has_path(p2), true);
	delete rt;
	delete root;
	delete path;
	delete x;
	delete set;
	delete set2;
}

void run_GOAnnotationParser(std::string path_test_ontology){
	std::cout << "Unit Test: GOAnnotationParser\n";
	GOAnnotationParser* GOparser = new GOAnnotationParser();
	//std::cout << "Loading ontology at " << path_test_ontology << "\n";
	GOparser->load(path_test_ontology);
	GOAnnotation * go;

	test(std::string("GOAnnotationParser n_annotations"), 1, GOparser->size(), 31716 );
	test(std::string("GOAnnotationParser first_annotation"), 2, GOparser->get_idx(std::string("GO:0000001")), 0);
	go = GOparser->get_annotation(0);
	test(std::string("GOAnnotationParser description"), 3, strcmp( go->description.c_str(), "mitochondrion inheritance"), 0 );
	delete GOparser;
}

void run_GeneAnnotationParser(std::string path_test_ontology, std::string path_test_annotation){
	std::cout << "Unit Test: GeneAnnotationParser\n";
	GOAnnotationParser* GOparser = new GOAnnotationParser();
	GeneAnnotationParser* geneparser = new GeneAnnotationParser();
	//std::cout << "Loading ontology at " << path_test_ontology << "\n";
	//std::cout << "Loading gene annotation at " << path_test_annotation << "\n";
	GOparser->load(path_test_ontology);
	geneparser->load(path_test_annotation, *GOparser);
	GeneAnnotation* ga;
	int caught=0;
	try{
		ga = geneparser->get_annotation(std::string("foo"));
	}
	catch(std::string err){
		caught=1;
	}
	test(std::string("GeneAnnotationParser bogus_symbol_test"), 1, caught,1 );
	ga = geneparser->get_annotation(std::string("LGR5"));
	test(std::string("GeneAnnotationParser full_name"), 2, strcmp( ga->full_name.c_str(), "Homo sapiens G protein-coupled receptor LGR5 (LGR5) mRNA, complete cds."), 0 );
	std::vector<std::string> symbols;
	int idx = GOparser->get_idx(std::string("GO:0005198"));
	geneparser->get_symbols_with_GO_idx(idx, symbols);
	test(std::string("GeneAnnotationParser get_symbols_with_GO_idx1"), 3, (int)symbols.size(), 1);
	test(std::string("GeneAnnotationParser get_symbols_with_GO_idx2"), 4, strcmp( symbols.at(0).c_str(), "krt71"), 0 );
	delete GOparser;
	delete geneparser;
}


void run_KeyValuesParser(std::string path_keyvalue){
    std::cout << "Unit Test: KeyValueParser\n";
    KeyValuesParser kp(path_keyvalue);
    std::vector<std::string> k, vals;
    kp.keys(k);
    test(std::string("KeyValuesParser number keys"), 1, int(k.size()), 3);
    kp.get( "key1", vals);
    test(std::string("KeyValuesParser has_key true"), 2, kp.has_key(std::string("key1")), true);
    test(std::string("KeyValuesParser has_key false"), 3, kp.has_key(std::string("bong")), false);
    test(std::string("KeyValuesParser key one has three values"), 4, int(vals.size()), 3);
    test(std::string("KeyValuesParser key one value one is foo"), 5, std::string("foo").compare(vals[0]), 0);
}

void run_CarmenParser_EQTL(std::string path_test_eqtl){
	std::cout << "Unit Test: CarmenParser_EQTL\n";
	CarmenParser cp(path_test_eqtl);
	cp.PrepareToReadValues();
	std::vector<QTL*> qtls;
	cp.ExtractEQTL(qtls, 0.01);
	test(std::string("CarmenParser_EQTL pval"), 1, qtls.at(0)->pval, 0);
	test(std::string("CarmenParser_EQTL probe"), 2, strcmp( qtls.at(0)->probe_id.c_str(), "1422027_a_at"), 0) ;
	test(std::string("CarmenParser_EQTL n"), 3, int(qtls.size()), 4); // two eqtl should be excluded by pval
	for(int i=0; i<int(qtls.size()); i++)
		delete qtls.at(i);
}

void run_Attributes(std::string basedir){
	std::cout << "Unit Test: Attributes\n";
	std::string testname = "Attributes";

	Attributes* a = new Attributes("NA");
	a->load(basedir + "/test_sample_attributes.txt");
	test(std::string("Attributes load"), 1, (int)a->identifiers.size(), 8);
	test(std::string("Attributes prop_for_identifier"), 2, a->prop_for_identifier("S5", "PAVM").compare("1"), 0);
	test(std::string("Attributes prop_for_identifier2"), 3, a->prop_for_identifier("S5", "ERP").compare("0"), 0);
	vector<int>* rows = new vector<int>();
	a->indices_with_property(std::string("PAVM"), std::string("1"), rows);
	test(std::string("Attributes indices_with_property"), 4, (int)rows->size() , 4);
	rows->clear();
	a->indices_with_property(std::string("PAVM"), std::string("1"), rows, false);
	test(std::string("Attributes indices_with_property case insensitive"), 5, (int)rows->size() , 4);
	delete rows;
	vector<std::string>* A = new vector<std::string>();
	vector<std::string>* B = new vector<std::string>();
	a->find_identifiers_in_class("ERP=1,PAVM=0", A);
	a->find_identifiers_in_class("ERP=1,PAVM=1", B);
	sort(A->begin(), A->end());
	sort(B->begin(), B->end());
	test(std::string("Attributes find_identifiers_in_class"), 6, (int)A->size(), 2);
	test(std::string("Attributes find_identifiers_in_class2"), 7, (int)B->size(), 2);
	test(std::string("Attributes find_identifiers_in_class3"), 8, A->at(0).compare("S2"), 0);
	test(std::string("Attributes find_identifiers_in_class4"), 9, A->at(1).compare("S4"), 0);

	vector<int>* source = new vector<int>();
	vector<int>* diff = new vector<int>();
	source->push_back(1); source->push_back(2); source->push_back(3); source->push_back(4); source->push_back(5);
	diff->push_back(1); diff->push_back(2); diff->push_back(3);
	a->intersection(source, diff);
	test(std::string("Attributes intersection"), 10, (int)diff->size(), 3);
	std::sort(diff->begin(), diff->end() );
	test(std::string("Attributes intersection2"), 11, diff->at(0), 1);
	test(std::string("Attributes intersection3"), 12, diff->at(1), 2);
	test(std::string("Attributes intersection4"), 13, diff->at(2), 3);
    Attributes* G = new Attributes("NA");
	G->load(basedir + "/gene_attributes.txt");
    G->set_chromosome_and_locus(std::string("CHR"), std::string("loc_start") );
    std::vector<int> matching_loci;
    G->get_idx_by_genomic_range(std::string("1"), 10, 20, matching_loci);
    test(std::string("Attributes chromosome_and_locus"), 14, int(matching_loci.size()), 1);
    G->get_idx_by_genomic_range(std::string("1"), 10, 35, matching_loci);
    test(std::string("Attributes chromosome_and_locus"), 15, int(matching_loci.size()), 2);
    test(std::string("Attributes chromosome_and_locus"), 16, matching_loci.at(0), 0);
    test(std::string("Attributes chromosome_and_locus"), 17, matching_loci.at(1), 2);
    
    std::string chromosome;
    int locus;
    bool is_present;
    is_present = G->get_chromosome_and_locus_by_identifier( std::string("GENE1"), chromosome, locus);
    test(std::string("Attributes chromosome_and_locus"), 18, is_present, true);
    test(std::string("Attributes chromosome_and_locus"), 19, locus, 20);
    test(std::string("Attributes chromosome_and_locus"), 20, chromosome.compare("1"), 0);
    is_present = G->get_chromosome_and_locus_by_identifier( std::string("GENE2"), chromosome, locus);
    test(std::string("Attributes chromosome_and_locus"), 18, is_present, false); // has missing value    
    
    delete G;
    delete A;
	delete B;
	delete a;
	delete source;
	delete diff;
    
}


void run_Rawdata(std::string fn_data_dir){
	std::cout << "Unit Test: Rawdata\n";

	Rawdata* rd = new Rawdata();
	std::string fn(fn_data_dir + "/test_dataset_realvalued.txt");
	rd->load(fn);
	int r=1;
	test(std::string("Rawdata rows"), r++, rd->data->rows(), 3);
	test(std::string("Rawdata cols"), r++, rd->data->cols(), 8);
	test(std::string("Rawdata is_missing"), r++, rd->data->is_missing(0,0), true);
	test(std::string("Rawdata is_missing"), r++, rd->data->is_missing(0,1), false);
	test(std::string("Rawdata is_missing"), r++, rd->data->is_missing(2,3), true);
	test(std::string("Rawdata value"), r++, rd->data->arr[2][7], (float)1);
	test(std::string("Rawdata value"), r++, rd->data->arr[1][1], (float)1.5);
	test(std::string("Rawdata index"), r++, rd->index_of_identifier("GENE1"), 0 );
	test(std::string("Rawdata index not present"), r++, rd->index_of_identifier("BOGUS"), -1 );
	test(std::string("Rawdata sample names"), r++, rd->sample_names.at(0).compare("S1"), 0 );
	test(std::string("Rawdata file name"), r++, rd->file_name.compare(fn), 0 );
	test(std::string("Rawdata identifier"), r++, rd->identifiers.at(0).compare("GENE1"), 0);
	test(std::string("Rawdata identifier N"), r++, rd->identifiers.size(), 3);
	test(std::string("Rawdata sample names N"), r++, rd->sample_names.size(), 8);
	test(std::string("Rawdata sample name"), r++, rd->sample_names.at(0).compare("S1"), 0);

	delete rd;
}



int main(int argc, char *argv[]){

	vector<Option*>* options = new vector<Option*>();
	options->push_back( new Option("dir", "d", "Path to directory with test data.  Defaults to ../../test_data", "../../test_data", OPT_OPTIONAL));
	options->push_back( new Option("verbose", "v", "Verbose, if T report all tests, if F report only failures, defaults F", "F", OPT_OPTIONAL));
	int retval=read_args(argc, argv, options, std::string("Tests for various classes"));
	if( retval==-2 ){
		delete options->at(0);
		delete options;
		return(0);
	}
	else if( retval != 0 ){
		std::cout << "ERROR: Could not read arguments, bailing out.\n";
		delete options->at(0);
		delete options;
		return(0);
	}
	std::string dir_data = options->at(0)->value;
    std::cout << "Data directory: " << dir_data << "\n";
    if( options->at(1)->value.compare("T")==0 ){
		verbose=true;
		std::cout << "Running in verbose mode\n";
	}
	else
		std::cout << "Runing silently: only noting test failures\n";


	run_Attributes(dir_data);
	run_Matrix();
	run_spear(dir_data);

	//run_CacheHash();
	run_graph();
	run_ProbeSets(dir_data);

	run_GOAnnotationParser(dir_data + "/test_ontology.txt");
	run_GeneAnnotationParser(dir_data + "/test_ontology.txt", dir_data + "/test_annotation.txt");
	run_CarmenParser_EQTL(dir_data + "/test_eqtl.txt");
    run_KeyValuesParser(dir_data + "/keyvalue.txt");
	run_node();
	run_Rawdata(dir_data);
	run_Dataset(dir_data);
	run_ClassifierDataset(dir_data);
	run_Discretize();
	run_Adaboost();
	run_Rule();
	run_RuleSet(dir_data);
	run_Perm();

	delete options->at(0);
	delete options->at(1);
	delete options;

}
