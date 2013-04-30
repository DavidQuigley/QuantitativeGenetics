#include <time.h>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/random.hpp>
using namespace boost;
#include "DataStructures.h"
#include "ClassMinerOptions.h"
#include "Attributes.h"
#include "Rawdata.h"
#include "Discretize.h"
#include "Dataset.h"
#include "chi2.h"
#include "Rule.h"

Rule::Rule(){
	this->chi2=-1;
	this->conf_p=-1;
	this->imp=-1;
	this->rule_class=-1;
    this->rule_idx=-1;
	this->sup_a_n=-1;
	this->sup_a_p=-1;
	this->sup_b_n=-1;
	this->sup_b_p=-1;
	this->t_stat=0;
	this->weight=0;
	this->parent_idx=0;
}


Rule::Rule(Rule* rule_to_clone){
	this->clone(rule_to_clone);
}

void Rule::load_line(Ant *ant, ClassifierDataset* data, Attributes* ga, ClassMinerOptions* cmo){
	this->chi2 = ant->chi2;
	this->conf_p = (double)ant->conf;
	this->sup_a_n = ant->sup_A;
	this->sup_b_n = ant->sup_B;
	int n_a = (int)data->a_idx->size();
	int n_b = (int)data->b_idx->size();
	n_a>0 ? this->sup_a_p = (float)this->sup_a_n / (float)n_a * 100 : this->sup_a_p = 0.0;
	n_b>0 ? this->sup_b_p = (float)this->sup_b_n / (float)n_b * 100 : this->sup_b_p = 0.0;

	this->imp = (double)ant->imp;
	int idx_in_dis=0, ident_idx=0;
	for(int i=0; i<(int)ant->idx->size(); i++){
		idx_in_dis = ant->idx->at(i);
		ident_idx = data->translate->at(idx_in_dis)[0];
		std::stringstream ss_g, ss_i, ss_r;
		ss_g << ga->prop_for_identifier( data->identifiers->at( ident_idx ), "Gene Name" ) << "_is_" << data->translate->at(idx_in_dis)[1];
		ss_i << data->identifiers->at( ident_idx ) << "_is_" << data->translate->at(idx_in_dis)[1];
		this->rules_g.push_back( ss_g.str() );
		this->rules_i.push_back( ss_i.str() );
		
		if( cmo->discretization.compare("none")==0 ){
			ss_r << data->identifiers->at( ident_idx ) << "=" << data->translate->at(idx_in_dis)[1];
		}
		else{
			if( data->translate->at(idx_in_dis)[1]>0 )
				ss_r << data->identifiers->at( ident_idx ) << ">" << data->discrete->cutoff_up.at(ident_idx) ;
			else
				ss_r << data->identifiers->at( ident_idx ) << "<" << data->discrete->cutoff_dn.at(ident_idx) ;
		}
		this->rules_r.push_back( ss_r.str());
	}
}


void Rule::load_line(std::string line){

	std::vector<std::string>* arr = new std::vector<std::string>();
	stringhasher::split(line, arr, std::string("\t").c_str());
	this->chi2 = boost::lexical_cast<double>( arr->at(0) ); 
	this->conf_p = boost::lexical_cast<double>( arr->at(1) ); 
	this->sup_a_n = boost::lexical_cast<int>( arr->at(2) ); 
	this->sup_a_p = boost::lexical_cast<double>( arr->at(3) ); 
	this->sup_b_n = boost::lexical_cast<int>( arr->at(4) ); 
	this->sup_b_p = boost::lexical_cast<double>( arr->at(5) ); 
	this->imp = boost::lexical_cast<double>( arr->at(6) ) * 100;
	stringhasher::split( arr->at(7), &this->rules_g, std::string("|").c_str());
	stringhasher::split( arr->at(8), &this->rules_i, std::string("|").c_str());
	stringhasher::split( arr->at(9), &this->rules_r, std::string("|").c_str());
	try{
		this->t_stat = boost::lexical_cast<double>( arr->at(10) ); 
	}
	catch(bad_lexical_cast &){
		this->t_stat = 0;
	}
	if( this->sup_a_p > this->sup_b_p )
		this->rule_class = Rule::CLASS_A;
	else
		this->rule_class = Rule::CLASS_B;
	delete arr;
}


void Rule::evaluate(ClassifierDataset* data, bool disc_is_none){
	float A=0, B=0;
	std::string id, value_str;
	stringhasher::get_idval(this->rules_i.at(0), id, value_str);
	int val_rule = -1;
	if( value_str.compare("1")==0 )
		val_rule = 1;
	int row = data->raw_data->identifier2idx[id];
	int x, val_dis;
	int n_a = (int)data->a_idx->size();
	int n_b = (int)data->b_idx->size();
	for(x=0;x<n_a;x++){ 
		val_dis = data->discrete->dis->arr[row][ data->a_idx->at(x) ]; 
		if( val_dis == val_rule)
			A += 1;
	}
	for(x=0;x<n_b;x++){ 
		val_dis = data->discrete->dis->arr[row][ data->b_idx->at(x) ]; 
		if( val_dis == val_rule)
			B += 1;
	}
	double chi2 = calc_chi2(A, n_a-A, B, n_b-B);
	
	if( this->rule_class==Rule::CLASS_A ){
		std::stringstream ss;
		this->conf_p = A/(A+B);
		if(!disc_is_none){
			ss << id << "<" << data->discrete->cutoff_dn.at(row);
			this->rules_r.push_back( ss.str() );
		}
	}
	else{
		std::stringstream ss;
		this->conf_p = B/(A+B);
		if(!disc_is_none){
			ss << id << ">" << data->discrete->cutoff_up.at(row);
			this->rules_r.push_back( ss.str() );
		}
	}
	this->sup_a_n = (int)A;
	this->sup_a_p = A / (float)n_a;
	this->sup_b_n = (int)B;
	this->sup_b_p = B / (float)n_b;
	this->chi2 = pchisq(chi2, 1);
	this->imp = this->conf_p;
}


void Rule::clone_out_probes(std::vector<Rule*>& rules_from_split ){
	rules_from_split.clear();
	for(int n=0; n<(int)this->rules_i.size(); n++){
		Rule* s = new Rule();
		s->chi2 = this->chi2;
		s->conf_p = this->conf_p;
		s->imp = this->imp;
		s->parent_idx = this->parent_idx;
		s->rule_class = this->rule_class;
		s->rule_idx = this->rule_idx;
		s->sup_a_n = this->sup_a_n;
		s->sup_b_n = this->sup_b_n;
		s->sup_a_p = this->sup_a_p;
		s->sup_b_p = this->sup_b_p;
		s->t_stat = this->t_stat;
		s->weight = this->weight;
		s->rules_g.push_back(this->rules_g.at(n));
		s->rules_i.push_back(this->rules_i.at(n));
		s->rules_r.push_back(this->rules_r.at(n));
		rules_from_split.push_back( s );
	}
}

void Rule::clone(Rule* R){
	this->chi2 = R->chi2;
	this->conf_p = R->conf_p;
	this->imp = R->imp;
	this->parent_idx = R->parent_idx;
	this->rule_class = R->rule_class;
	this->rule_idx = R->rule_idx;
	this->sup_a_n = R->sup_a_n;
	this->sup_b_n = R->sup_b_n;
	this->sup_a_p = R->sup_a_p;
	this->sup_b_p = R->sup_b_p;
	this->t_stat = R->t_stat;
	this->weight = R->weight;
	this->rules_g.clear(); 
	this->rules_i.clear(); 
	this->rules_r.clear();
	for(int i=0; i<(int)R->rules_g.size(); i++)
		this->rules_g.push_back( R->rules_g.at(i) );
	for(int i=0; i<(int)R->rules_i.size(); i++)	
		this->rules_i.push_back( R->rules_i.at(i) );
	for(int i=0; i<(int)R->rules_r.size(); i++)	
		this->rules_r.push_back( R->rules_r.at(i) );
}


std::string Rule::get_friendly_label(){
	std::stringstream ss;
	std::string label = "";
	for(int i=0; i<(int)this->rules_g.size(); i++){
		ss << this->rules_g.at(i);
		if( i<(int)this->rules_g.size()-1 )
			ss << "|";
	}
	label = ss.str();
	return label;
}


std::string Rule::get_identifier_label(){
	std::stringstream ss;
	std::string label = "";
	for(int i=0; i<(int)this->rules_i.size(); i++){
		ss << this->rules_i.at(i);
		if( i<(int)this->rules_i.size()-1 )
			ss << "|";
	}
	label = ss.str();
	return label;
}


std::string Rule::get_identifier_cutoff_label(){
	std::stringstream ss;
	std::string label = "";
	for(int i=0; i<(int)this->rules_r.size(); i++){
		ss << this->rules_r.at(i);
		if( i<(int)this->rules_r.size()-1 )
			ss << "|";
	}
	label = ss.str();
	return label;
}


Edge::Edge(std::string n_1, std::string v_1, std::string n_2, std::string v_2, int edge_class, double sup, double conf, double p_value){
	// always list alphabetically first name as name1
	if(name1.compare(name2)>0){
		this->name1 = n_2;
		this->name2 = n_1;
		this->value1 = v_2;
		this->value2 = v_1;
	}
	else{
		this->name1 = n_1;
		this->name2 = n_2;
		this->value1 = v_1;
		this->value2 = v_2;
	}
	this->edge_class = edge_class;
	this->sup = sup;
	this->conf = conf;
	this->p_value = p_value;
}


RuleSet::RuleSet(){
	this->edges = new std::vector<Edge*>();
}


RuleSet::RuleSet(std::string file_name, bool strip_edges){
	this->edges = new std::vector<Edge*>();

	std::string line;
	std::vector<std::string>* arr = new std::vector<std::string>();
	std::string cmd;
	Rule* rule;
	std::vector<std::string>* sci = new std::vector<std::string>();
	HASH_S_I rules_seen;                 // used if stripping edges
	std::vector<Rule*> rules_from_split; // used if stripping edges
	char space[1], comma[1];
	space[0] = ' ';
	comma[0] = ',';
	if( file_name.size()>0 ){
		
		std::ifstream f_cv_testcrlf(file_name.c_str());
		if( !f_cv_testcrlf.is_open() ){
			throw std::string("ERROR: Unable to open rule file for reading: " + file_name);
		}
		bool has_CRLF = false;
		getline(f_cv_testcrlf, line);
		if( line.at( line.size()-1 ) == '\r' )
			has_CRLF=true; //Windows-originated file
		f_cv_testcrlf.close();
		std::ifstream f_cv(file_name.c_str());
		while( !f_cv.eof() ){
			getline(f_cv, line);
			if(line.length()==0)
				break;
			if(has_CRLF)
				line = line.substr(0,line.size()-1);
			try{
				if( line.at(0) == '#' ){
					stringhasher::split(line, arr, space);
					cmd = arr->at(1);
					if( cmd.compare("Data_File")==0){
						if( arr->size()>2 )
							this->file_name_dis = arr->at(2);
					}
					if( cmd.compare("Sample_File")==0){
						if( arr->size()>2 )
							this->file_name_sa = arr->at(2);
					}
					if( cmd.compare("Genes_File")==0){
						if( arr->size()>2 )
							this->file_name_ga = arr->at(2);
					}
					if( cmd.compare("Class_A")==0){
						if( arr->size()>3 ){
							this->class_a = arr->at(2);
							this->class_a_size = boost::lexical_cast<int>( arr->at(3) );
						}
					}
					if( cmd.compare("Class_B")==0){
						if( arr->size()>3 ){
							this->class_b = arr->at(2);
							if( this->class_b.compare("0")==0 ){
								this->class_b="";
								this->class_b_size=0;
							}
							else{
								this->class_b_size = boost::lexical_cast<int>( arr->at(3));
							}
						}
					}
					if( cmd.compare("Gene_Limit")==0){
						if( arr->size()>2 )
							this->gene_limit = arr->at(2);
					}
					if( cmd.compare("S_C_I")==0){
						if( arr->size()>2 ){
							stringhasher::split(arr->at(2), sci, comma);
							this->minsup = boost::lexical_cast<double>( sci->at(0) );
							this->minconf = boost::lexical_cast<double>( sci->at(1) ) / 100;
							this->minimp = boost::lexical_cast<double>( sci->at(2) );
						}
					}
					if( cmd.compare("Max_Depth")==0){
						if( arr->size()>2 )
							this->max_depth = boost::lexical_cast<int>( arr->at(2) );
					}
					if( cmd.compare("Discretization")==0){
						this->disc = arr->at(2);
						if(this->disc.compare("none")==0 || this->disc.compare("None")==0){
							this->disc_lower = -1;
							this->disc_upper = -1;
						}
						else{
							stringhasher::split(arr->at(5), sci, comma);
							this->disc_lower = boost::lexical_cast<double>(sci->at(0) );
							this->disc_upper = boost::lexical_cast<double>(sci->at(1) );
						}
					}
					if( cmd.compare("Maximum_Chi2")==0){
						if( arr->size()>2 )
							this->max_chi2 = boost::lexical_cast<double>(arr->at(2));
					}
					if( cmd.compare("Tested")==0){
						if( arr->size()>8 ){
							this->n_tested = boost::lexical_cast<int>(arr->at(2));
							this->bon_pval = boost::lexical_cast<double>(arr->at(8));
						}
					}
				}
			}
			catch( ... ){
				std::stringstream ss;
				ss << "error parsing ruleset header line: " << line;
				throw std::string( ss.str() );
			}
			if( cmd.compare("Generated")==0){
				try{
					int n_rules = boost::lexical_cast<int>(arr->at(2));
					getline(f_cv, line); // eat the header line

					for(int i=0; i<n_rules; i++){
						getline(f_cv, line);
						if(line.length()==0)
							break;
						if(has_CRLF)
							line = line.substr(0,line.size()-1);
						rule = new Rule();
						rule->load_line(line);
						if(strip_edges){
							int n_items = (int)rule->rules_g.size();
							if( n_items==1 ){
								if( rules_seen.find( rule->rules_i.at(0) ) == rules_seen.end() ){
									rules_seen[rule->rules_i.at(0)] = 1;
									this->rules.push_back(rule);
								}
							}
							else{
								//int n_added = 0;
								rule->clone_out_probes(rules_from_split);
								for(int n=0; n<n_items; n++){
									rule = rules_from_split.at(n);
									if( rules_seen.find( rule->rules_i.at(0) ) == rules_seen.end() ){
										rules_seen[rule->rules_i.at(0)] = 1;
										this->rules.push_back(rule);
									}
								}
							}
						}
						else{
							this->rules.push_back(rule);
						}
					}
				}
				catch( ... ){
					throw std::string("Error reading lines ruleset");
				}
			}
		}
		f_cv.close();
	}
	for(int i=0; i<(int)rules_from_split.size(); i++)
		delete rules_from_split.at(i);
	delete sci;
	delete arr;
}


RuleSet::RuleSet(ClassMinerOptions* cmo, std::vector<Ant*>* ants, ClassifierDataset* data, Attributes* ga){
	this->edges = new std::vector<Edge*>();

	const int UNKNOWN = -1;
	this->file_name_dis = cmo->file_name_dis;
	this->file_name_sa = cmo->file_name_sa;
	this->file_name_ga = cmo->file_name_ga;
	this->class_a = cmo->class_a;
	this->class_b = cmo->class_b;
	
	this->gene_limit = cmo->gene_limit;
	this->minsup = cmo->min_sup;
	this->minconf = cmo->min_conf;
	this->minimp = cmo->min_imp;
	this->max_depth = cmo->max_depth;
	this->disc = cmo->discretization;
	this->disc_lower = cmo->disc_lower;
	this->disc_upper = cmo->disc_upper;
	this->max_chi2 = cmo->max_chi2;
	this->class_a_size = UNKNOWN;
	this->class_b_size = UNKNOWN;
	this->n_tested = UNKNOWN;
	this->bon_pval = (double)UNKNOWN;
	
	Rule* rule;
	for(int i=0; i<(int)ants->size(); i++){
		rule = new Rule();
		rule->load_line(ants->at(i), data, ga, cmo);
		this->rules.push_back(rule);
	}
}


void RuleSet::write_header(std::ofstream& f_out, int n_rules=-1, double merge_threshold=-1){
	// can pass n_rules and merge_threshold to allow case when we have merged redundant rules.
	// default (-1) for these values prints the true number of rules
	char timebuf[80];
	struct tm* newtime;
	time_t long_time;
	time( &long_time );
	newtime = localtime( &long_time );
	strftime(timebuf, 80, "%Y_%m_%d_%H_%M_%S", newtime);
	f_out << "# Source CARMEN 1.0\n";
	f_out << "# Date " << timebuf << "\n";
	f_out << "# Data_File " << this->file_name_dis << "\n";
	f_out << "# Sample_File " << this->file_name_sa << "\n";
	f_out << "# Genes_File " << this->file_name_ga << "\n";
	f_out << "# Class_A " << this->class_a << " " << this->class_a_size << " members.\n";
	f_out << "# Class_B " << this->class_b << " " << this->class_b_size << " members.\n";
	f_out << "# S_C_I " << this->minsup << "," << this->minconf*100 << "," << this->minimp << "\n";
	f_out << "# Max_Depth " << this->max_depth << "\n";
	if(this->gene_limit.length()>0)
		f_out << "# Gene_Limit " << this->gene_limit << "\n";
	if( this->disc.compare("none")==0 || this->disc.compare("None")==0 )
		f_out << "# Discretization none Using raw values, no discretization.\n";
	else
		f_out << "# Discretization " << this->disc << " Lower,Upper multiples " << this->disc_lower << "," << this->disc_upper << "\n";
	f_out << "# Maximum_Chi2 " << this->max_chi2 << "\n";
	double pval=1;
	if( this->n_tested>0)
		pval = 0.05 / this->n_tested;
	f_out << "# Tested " << this->n_tested << " hypotheses.  Bonferroni 0.05 corrected p-value " << pval << "\n";
	if( merge_threshold > -1 )
		f_out << "# Merged " << merge_threshold << " percent similar rules\n";
	if( n_rules == -1 )
		f_out << "# Generated " << this->rules.size() << " rules.\n";
	else
		f_out << "# Generated " << n_rules << " rules.\n";
	f_out << "# chi2\t%conf\tsupA\t%supA\tsupB\t%supB\timp\tRule\tRule_ID\tRule_Dis\tt_stat\n";;
}


void RuleSet::write_rule(std::ofstream& f_out, int rule_idx){
	Rule* r = this->rules.at(rule_idx);
	f_out << r->chi2 << "\t" << r->conf_p << "\t" << r->sup_a_n << "\t" << r->sup_a_p << "\t" << r->sup_b_n << "\t" << r->sup_b_p << "\t" << r->imp << "\t";
	f_out << r->get_friendly_label() << "\t" << r->get_identifier_label() << "\t" << r->get_identifier_cutoff_label() << "\t";
	f_out << r->t_stat << "\n";
}


void RuleSet::write_expression(std::string file_name, bool verbose){
	
	// load data
	// load sample_attributes
	// calculate class_a, class_b members
	ClassifierDataset* data = new ClassifierDataset();
	Attributes* sa;
	Attributes* ga;
	sa = new Attributes("NA");
	ga = new Attributes("NA");
	sa->load(this->file_name_sa);
	ga->load(this->file_name_ga);
	ClassMinerOptions* cmo = new ClassMinerOptions();
	cmo->verbose = verbose;
	cmo->class_a = this->class_a;
	cmo->class_b = this->class_b;
	cmo->file_name_dis = this->file_name_dis;
	cmo->discretization = this->disc;
	cmo->disc_lower = (float)this->disc_lower;
	cmo->disc_upper = (float)this->disc_upper;
	if(cmo->verbose){
		std::cout << "MESSAGE: Loading data set...\n";
		std::cout.flush();
	}
	data->load(sa, ga, cmo);
	if(cmo->verbose){
		std::cout << "MESSAGE: Loaded data set...\n";
		std::cout.flush();
	}
	std::ofstream f_out(file_name.c_str());
	if( !f_out.is_open() ) 
		throw std::string("Unable to open for writing: " + file_name);
	this->write_header(f_out);

	HASH_S_S id2name;
	std::vector<std::string> unique_ids;
	Rule* rule;
	std::vector<std::string> rules_i;
	std::vector<std::string> rules_g;
	std::string id, id_1, id_2;

	for(int i=0; i<(int)this->rules.size(); i++){
		rule = this->rules.at(i);
		rules_i = rule->rules_i;
		rules_g = rule->rules_g;
		for(int j=0; j<(int)rules_i.size(); j++){
			id = rules_i.at(j);
			if( id.find("_is_-1") != std::string::npos )
				id = id.substr(0, id.length() - 6);
			else
				id = id.substr(0, id.length() - 5);	

			if( id2name.find(id) == id2name.end() ){
				unique_ids.push_back( id );
			}
			id2name[id] = rules_g.at(j);
		}
	}
	int row_idx, c;
	for(int i=0; i<(int)unique_ids.size(); i++){
		f_out << id2name[ unique_ids.at(i) ] << "\t" << unique_ids.at(i);
		row_idx = data->raw_data->identifier2idx[ unique_ids.at(i) ];
		for(int j=0; j<this->class_a_size; j++){
			c = data->a_idx->at(j);
			f_out << "\t" << data->raw_data->data->arr[row_idx][c];
		}
		for(int j=0; j<this->class_b_size; j++){
			c = data->b_idx->at(j);
			f_out << "\t" << data->raw_data->data->arr[row_idx][c];
		}
		f_out << "\n";
	}
	f_out.close();
}


bool comp_counts(const std::pair<std::string, int> *a, const std::pair<std::string, int> *b){
	return a->second > b->second;
}

void RuleSet::write_frequency_list(std::string fn_out, std::string fn_source, bool verbose){
	HASH_S_I gene2count;
	Rule* rule;
	std::string gene, value;
	for(int i=0; i<(int)this->rules.size(); i++){
		rule = this->rules.at(i);
		for(int g=0; g<(int)rule->rules_g.size(); g++){
			stringhasher::get_idval(rule->rules_g.at(g), gene, value);
			if( gene2count.find( gene ) == gene2count.end() )
				gene2count[gene] = 1;
			else
				gene2count[gene] = gene2count[gene]+1;
		}
	}
	std::vector< std::pair<std::string, int>* >  counts;
	for( HASH_S_I::iterator itr=gene2count.begin(); itr != gene2count.end(); itr++ ){
		std::pair<std::string, int>* p = new std::pair<std::string, int>();
		p->first = (*itr).first;
		p->second = (*itr).second;
		counts.push_back(p);
	}
	std::sort(counts.begin(), counts.end(), comp_counts);
	std::ofstream f_out(fn_out.c_str());
	if( !f_out.is_open() ) 
		throw std::string("Unable to open for writing: " + fn_out);
	f_out << "# Gene frequency list from source file " << fn_source << "\n";
	for(int i=0; i<(int)counts.size(); i++){
		f_out << counts.at(i)->first << "\t" << counts.at(i)->second << "\n";
		delete counts.at(i);
	}
	f_out.close();
}


void RuleSet::write_gene_name_list(std::string file_name, bool verbose){
	HASH_S_I gene2;
	Rule* rule;
	std::string rulepart, identifier, value;
	for(int i=0; i<(int)this->rules.size(); i++){
		rule = this->rules.at(i);
		for(int g=0; g<(int)rule->rules_g.size(); g++){
			stringhasher::get_idval(rule->rules_g.at(g), identifier, value);
			if( gene2.find(identifier) == gene2.end() )
				gene2[identifier] = 1;
			else
				gene2[identifier] = gene2[identifier] + 1;
		}
	}
	std::ofstream f_out(file_name.c_str());
	if( !f_out.is_open() ) 
		throw std::string("Unable to open for writing: " + file_name);
	std::vector<std::string> gene_names;
	for( HASH_S_I::iterator itr=gene2.begin(); itr != gene2.end(); itr++ ){
		gene_names.push_back( (*itr).first );
	}
	std::sort(gene_names.begin(), gene_names.end());
	for(int i=0; i<(int)gene_names.size(); i++)
		f_out << gene_names.at(i) << "\n";
	f_out << "# Gene name list from source file " << file_name << "\n";
	f_out.close();
}


void RuleSet::write_crosstab(std::string file_name){
	this->find_edges();
	
	HASH_S_I unique_edges;
	std::vector<std::string> edge_names;
	int n_edges=0;
	Edge* edge;
	std::string e1, e2;
	for(int i=0; i<(int)this->edges->size(); i++ ){
		edge = this->edges->at(i);
		e1 = edge->name1 + "_is_" + edge->value1;
		e2 = edge->name2 + "_is_" + edge->value2;
		if( unique_edges.find(e1) == unique_edges.end() ){
			unique_edges[e1] = n_edges++;
			edge_names.push_back(e1);
		}
		if( unique_edges.find(e2) == unique_edges.end() ){
			unique_edges[e2] = n_edges++;
			edge_names.push_back(e2);
		}
	}
	Matrix<int> xtab(n_edges+1, n_edges+1);
	int r,c;
	for(int i=0; i<(int)this->edges->size(); i++ ){
		edge = this->edges->at(i);
		e1 = edge->name1 + "_is_" + edge->value1;
		e2 = edge->name2 + "_is_" + edge->value2;
		r = unique_edges[e1];
		c = unique_edges[e2];
		if( r == c ){
			(xtab.arr[r][c])++;
		}
		else{
			(xtab.arr[r][c])++;
			(xtab.arr[c][r])++;
		}
	}
	std::ofstream f_out(file_name.c_str());
	if( !f_out.is_open() ) 
		throw std::string("Unable to open for writing: " + file_name);
	
	// first row, starts with tab
	for(int i=0; i<n_edges; i++ )
		f_out << "\t" << edge_names.at(i);
	f_out << "\n";

	for(r=0; r<n_edges; r++){
		f_out << edge_names.at(r);
		for(c=0; c<n_edges; c++){
			f_out << "\t" << xtab.arr[r][c];
		}
		f_out << "\n";
	}

	f_out.close();
}


void RuleSet::write_network(std::string file_name_base){
	std::string file_name_ea_class = file_name_base + "_class.eda";
	std::string file_name_ea_sup = file_name_base + "_sup.eda";
	std::string file_name_ea_conf = file_name_base + "_conf.eda";
	std::string file_name_ea_pval = file_name_base + "_pval.eda";
	std::string file_name_sif = file_name_base + ".sif";
	std::string file_name_names = file_name_base + "_names.noa";
	std::string file_name_disc = file_name_base + ".noa";

	std::ofstream f_net(file_name_sif.c_str());
	std::ofstream f_disc(file_name_disc.c_str());
	std::ofstream f_names(file_name_names.c_str());
	std::ofstream f_eda_class(file_name_ea_class.c_str());
	std::ofstream f_eda_sup(file_name_ea_sup.c_str());
	std::ofstream f_eda_conf(file_name_ea_conf.c_str());
	std::ofstream f_eda_pval(file_name_ea_pval.c_str());
	if( !f_net.is_open() ) 
		throw std::string("Unable to open for writing: " + file_name_sif);
	if( !f_names.is_open() ) 
		throw std::string("Unable to open for writing: " + file_name_names);
	if( !f_disc.is_open() ) 
		throw std::string("Unable to open for writing: " + file_name_disc);
	if( !f_eda_class.is_open() ) 
		throw std::string("Unable to open file for writing: " + file_name_ea_class);
	if( !f_eda_sup.is_open() ) 
		throw std::string("Unable to open file for writing: " + file_name_ea_sup);
	if( !f_eda_conf.is_open() ) 
		throw std::string("Unable to open file for writing: " + file_name_ea_conf);
	if( !f_eda_pval.is_open() ) 
		throw std::string("Unable to open file for writing: " + file_name_ea_pval);
	f_disc << "Discretization (java.lang.Integer)\n";
	f_names << "Names (java.lang.Integer)\n";
	f_eda_class << "EdgeClass (java.lang.String)\n";
	f_eda_sup << "Support (java.lang.Double)\ncytoscape (gg) bugfix = 0.0\n";
	f_eda_conf << "Confidence (java.lang.Double)\ncytoscape (gg) bugfix = 0.0\n";
	f_eda_pval << "PValue (java.lang.Double)\ncytoscape (gg) bugfix = 0.0\n";

	this->find_edges();
	Edge* edge;
	std::string e1, e2;
	for(int i=0; i<(int)this->edges->size(); i++ ){
		edge = this->edges->at(i);
		e1 = edge->name1 + "_is_" + edge->value1;
		e2 = edge->name2 + "_is_" + edge->value2;
		f_net << e1 << "\tgg\t" << e2 << "\n";
		f_disc << e1 << " = " << edge->value1 << "\n";
		f_disc << e2 << " = " << edge->value2 << "\n";
		f_names << e1 << " = " << edge->name1 << "\n";
		f_names << e2 << " = " << edge->name2 << "\n";
		if( edge->edge_class == Rule::CLASS_A )
			f_eda_class << e1 << " (gg) " << e2 << " = CLASS_A" << "\n";
		else
			f_eda_class << e1 << " (gg) " << e2 << " = CLASS_B" << "\n";
		f_eda_sup << e1 << " (gg) " << e2 << " = " << edge->sup << "\n";
		f_eda_conf << e1 << " (gg) " << e2 << " = " << edge->conf << "\n";
		f_eda_pval << e1 << " (gg) " << e2 << " = " << edge->p_value << "\n";
	}
	f_net.close();
	f_disc.close();
	f_names.close();
	f_eda_class.close();
	f_eda_sup.close();
	f_eda_conf.close();
	f_eda_pval.close();
}


void RuleSet::find_edges(){
	std::vector<std::string> nodes;
	int x,y, N;
	Rule* rule;
	std::string gene_name_1, value_1, gene_name_2, value_2;
	for(int i=0; i<(int)this->rules.size(); i++){
		rule = this->rules.at(i);
		nodes.clear();
		for(int j=0;j<(int)rule->rules_g.size();j++)
			nodes.push_back(rule->rules_g.at(j));
		N = (int)nodes.size(); 
		if(N==1){
			stringhasher::get_idval(nodes.at(0), gene_name_1, value_1);
			Edge* e;
			if(rule->sup_a_p>rule->sup_b_p)
				e = new Edge(gene_name_1, value_1, gene_name_1, value_1, Rule::CLASS_A, rule->sup_a_p, rule->conf_p, rule->chi2);
			else
				e = new Edge(gene_name_1, value_1, gene_name_1, value_1, Rule::CLASS_B, rule->sup_a_p, rule->conf_p, rule->chi2);
			this->edges->push_back(e);
		}
		else{
			x = 0;
            y = 1;
			while( x<N-1){
				stringhasher::get_idval(nodes.at(x), gene_name_1, value_1);
				stringhasher::get_idval(nodes.at(y), gene_name_2, value_2);
				Edge* e;
				if(rule->sup_a_p>rule->sup_b_p)
					e = new Edge(gene_name_1, value_1, gene_name_2, value_2, Rule::CLASS_A, rule->sup_a_p, rule->conf_p, rule->chi2);
				else
					e = new Edge(gene_name_1, value_1, gene_name_2, value_2, Rule::CLASS_B, rule->sup_a_p, rule->conf_p, rule->chi2);
				this->edges->push_back(e);
                ++y;
				if( y == N){
					++x;
					y = x+1;
				}
			}
		}
	}
}


void RuleSet::write_rules_by_samples(std::string filename, Matrix<double>* r_by_s, Attributes* sa){
	std::ofstream f(filename.c_str());
	if( !f.is_open() ) 
		throw std::string("Unable to open for writing: " + filename);
	f << "Probe_ID\tName\t";
	for(int i=0; i<(int)sa->identifiers.size(); i++){
		f << sa->identifiers.at(i);
		if(i<(int)sa->identifiers.size()-1)
			f << "\t";
	}
	f << "\n";
	for(int a=0; a<(int)sa->attrib.size(); a++){
		f << "\t" << sa->attrib.at(a) << "\t";
		for(int i=0; i<(int)sa->identifiers.size(); i++){
			f << sa->prop_for_identifier(sa->identifiers.at(i), sa->attrib.at(a));
			if(i<(int)sa->identifiers.size()-1)
				f << "\t";
		}
		f << "\n";
	}
	int n_cols = r_by_s->cols();
	for(int r=0; r<r_by_s->rows(); r++){
		f << this->rules.at(r)->get_friendly_label() << "\t" << this->rules.at(r)->get_identifier_label() << "\t";
		for(int c=0; c<n_cols; c++){
			if( r_by_s->has_missing_values && r_by_s->is_missing(r, c) )
				f << "NA";
			else
                f << r_by_s->arr[r][c];
			if( c < n_cols-1 )
				f << "\t";
		}
		f << "\n";
	}
	f.close();
}


void RuleSet::build_rules_by_samples(Matrix<double>* r_by_s, ClassifierDataset* data){
	Rule* rule;
	bool has_all;
	int n, c, r, idx;
	
	std::string id;
	std::string value;
	int val;
	
	r_by_s->resize((int)this->rules.size(), (int)data->raw_data->sample_names.size());
    
	for( c=0; c< (int)data->raw_data->sample_names.size(); c++){
		for( r=0; r<(int)rules.size(); r++ ){
			rule = rules.at(r);
			has_all = true;
			for( n=0; n<(int)rule->rules_r.size(); n++ ){
				stringhasher::get_idval(rule->rules_i.at(n), id, value);
				idx = data->raw_data->identifier2idx[id];
				val = boost::lexical_cast<int>(value);
				if( data->discrete->dis->arr[idx][c] != val )
					has_all=false;
			}
			if(has_all){
				if( rule->sup_a_p > rule->sup_b_p)
                    r_by_s->arr[r][c] = Rule::CLASS_A;
				else
                    r_by_s->arr[r][c] = Rule::CLASS_B;
			}
			else
				r_by_s->arr[r][c]=0;
		}
	}
}



RuleSet::~RuleSet(){
	for(int i=0; i<(int)this->rules.size(); i++)
		delete this->rules.at(i);
	for(int i=0; i<(int)this->edges->size(); i++)
		delete this->edges->at(i);
	delete this->edges;
}

CrossValidation::CrossValidation(std::string file_name_cv){
	this->rulesets = new std::vector<RuleSet*>();
	
	std::ifstream f_cv(file_name_cv.c_str());
	if( !f_cv.is_open() ){ 
		throw new std::string( "Unable to open file for reading:" + file_name_cv );
	}
	std::string line;
	std::vector<std::string>* arr = new std::vector<std::string>();
	std::string cmd;
	RuleSet* rs;
	Rule* rule;
	std::vector<std::string>* sci = new std::vector<std::string>();
	char space[1], comma[1];
	space[0] = ' ';
	comma[0] = ',';
	while( !f_cv.eof() ){
		getline(f_cv, line);
		if(line.length()==0)
			break;
		
		if( line.at(0) == '#' ){
			stringhasher::split(line, arr, space);
			cmd = arr->at(1);

			if( cmd.compare("Data_File")==0)
                rs = new RuleSet();
			if( cmd.compare("Class_A")==0){
                rs->class_a = arr->at(2);
                rs->class_a_size = boost::lexical_cast<int>( arr->at(3) );
			}
			if( cmd.compare("Class_B")==0){
                rs->class_b = arr->at(2);
                rs->class_b_size = boost::lexical_cast<int>( arr->at(3));
			}
			if( cmd.compare("S_C_I")==0){
                stringhasher::split(arr->at(2), sci, comma);
				rs->minsup = boost::lexical_cast<double>( sci->at(0) );
				rs->minconf = boost::lexical_cast<double>( sci->at(1) );
				rs->minimp = boost::lexical_cast<double>( sci->at(2) );
			}
			if( cmd.compare("Discretization")==0){
                rs->disc = arr->at(2);
				stringhasher::split(arr->at(5), sci, comma);
                rs->disc_lower = boost::lexical_cast<double>(sci->at(0) );
				rs->disc_upper = boost::lexical_cast<double>(sci->at(1) );
			}
            if( cmd.compare("Maximum_Chi2")==0){
                rs->max_chi2 = boost::lexical_cast<double>(arr->at(2));
			}
			if( cmd.compare("Tested")==0){
                rs->n_tested = boost::lexical_cast<int>(arr->at(2));
                rs->bon_pval = boost::lexical_cast<double>(arr->at(8));
			}
			if( cmd.compare("Generated")==0){
                int n_rules = boost::lexical_cast<int>(arr->at(2));
                getline(f_cv, line); // eat the header line
				
				for(int i=0; i<n_rules; i++){
                    getline(f_cv, line);
                    rule = new Rule();
                    rule->load_line(line);
                    rs->rules.push_back(rule);
				}
				this->rulesets->push_back(rs);
			}
		}
	}
}



void CrossValidation::find_core_rules(RuleSet* rs_core, Attributes* ga, double p_override){
    HASH_S_I core; 
	std::vector<std::string> core_keys;
    
	RuleSet* ruleset=NULL;
	Rule* r = NULL;
	std::string rule_sorted;
	int i,j,n;
	std::string out;
	double max_p;
	for( i=0;i<(int)this->rulesets->size();i++){
		ruleset = this->rulesets->at(i);
		if( p_override == -1)
			max_p = ruleset->bon_pval;
		else
			max_p = p_override;
        for( j=0; j<(int)ruleset->rules.size(); j++){
			r = ruleset->rules.at(j);	
			
			if( r->chi2 <= max_p){
				std::sort(r->rules_i.begin(), r->rules_i.end());
				stringhasher::join( rule_sorted, &r->rules_i, "|" );
				if( core.find(rule_sorted)==core.end() ){
					core_keys.push_back(rule_sorted);
					core[rule_sorted] = 1;
				}
				else
					core[rule_sorted] = core[rule_sorted] + 1;
			}
		}
	}
	int n_folds = (int)this->rulesets->size();
	std::string identifier;
	std::string value;
	char pipe[1];
	pipe[0] = '|';
	for(i=0; i<(int)core_keys.size(); i++){
		if(core[core_keys.at(i)]==n_folds){
			Rule* r_new = new Rule();
            stringhasher::split( core_keys.at(i), &r_new->rules_i, pipe);
			stringhasher::split( core_keys.at(i), &r_new->rules_g, pipe);
			for(n=0; n<(int)r_new->rules_i.size(); n++){
				stringhasher::get_idval(r_new->rules_i.at(n), identifier, value);
				r_new->rules_g.at(n) = ga->prop_for_identifier(identifier, "Gene Name") + "_is_" + value;
			}
			rs_core->rules.push_back(r_new);
		}
	}
}


CrossValidation::~CrossValidation(){
	for( int i=0;i<(int)this->rulesets->size();i++)
		delete this->rulesets->at(i);
	delete this->rulesets;
}
