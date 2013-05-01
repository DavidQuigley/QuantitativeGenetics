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
#include "ClassMiner.h"
#include "chi2.h"

// Author: David Quigley, UCSF Helen Diller Family Comprehensive Cancer Center.
// Email:  dquigley@cc.ucsf.edu
// Date:   2007-2010

int safe_parse(ClassMinerOptions* cmo, std::vector<Option*>* options){
    int r=0;
    try{ cmo->file_name_dis = options->at(r++)->value; }
    catch( boost::bad_lexical_cast err ){ std::cout << "ERROR: " << err.what() << " parsing option -d\n"; return -1; }

    try{ cmo->file_name_sa = options->at(r++)->value; }
    catch( boost::bad_lexical_cast err ){ std::cout << "ERROR: " << err.what() << " parsing option -f\n"; return -1; }

    try{ cmo->file_name_ga = options->at(r++)->value; }
    catch( boost::bad_lexical_cast err ){ std::cout << "ERROR: " << err.what() << " parsing option -g\n"; return -1; }

    try{ cmo->class_a = options->at(r++)->value; }
    catch( boost::bad_lexical_cast err ){ std::cout << "ERROR: " << err.what() << " parsing option -a\n"; return -1; }

    try{ cmo->class_b = options->at(r++)->value; }
    catch( boost::bad_lexical_cast err ){ std::cout << "ERROR: " << err.what() << " parsing option -b\n"; return -1; }

    try{ cmo->min_sup = boost::lexical_cast<int>( options->at(r++)->value); }
    catch( boost::bad_lexical_cast err ){ std::cout << "ERROR: " << err.what() << " parsing option -s\n"; return -1; }

    try{ cmo->min_conf = boost::lexical_cast<float>( options->at(r++)->value) / (float)100.0; }
    catch( boost::bad_lexical_cast err ){ std::cout << "ERROR: " << err.what() << " parsing option -c\n"; return -1; }

    try{ cmo->min_imp = boost::lexical_cast<float>( options->at(r++)->value) / (float)100.0; }
    catch( boost::bad_lexical_cast err ){ std::cout << "ERROR: " << err.what() << " parsing option -i\n"; return -1; }

    try{ cmo->max_chi2 = boost::lexical_cast<double>( options->at(r++)->value); }
    catch( boost::bad_lexical_cast err ){ std::cout << "ERROR: " << err.what() << " parsing option -p\n"; return -1; }

    try{ cmo->discretization = options->at(r++)->value; }
    catch( boost::bad_lexical_cast err ){ std::cout << "ERROR: " << err.what() << " parsing option -m\n"; return -1; }

    try{ cmo->disc_lower = boost::lexical_cast<float>( options->at(r++)->value); }
    catch( boost::bad_lexical_cast err ){ std::cout << "ERROR: " << err.what() << " parsing option -l\n"; return -1; }

    try{ cmo->disc_upper = boost::lexical_cast<float>( options->at(r++)->value); }
    catch( boost::bad_lexical_cast err ){ std::cout << "ERROR: " << err.what() << " parsing option -u\n"; return -1; }

    try{ cmo->gene_limit = options->at(r++)->value; }
    catch( boost::bad_lexical_cast err ){ std::cout << "ERROR: " << err.what() << " parsing option -r\n"; return -1; }

    try{ cmo->max_depth = boost::lexical_cast<int>( options->at(r++)->value); }
    catch( boost::bad_lexical_cast err ){ std::cout << "ERROR: " << err.what() << " parsing option -n\n"; return -1; }

    try{ cmo->timeout = boost::lexical_cast<int>( options->at(r++)->value); }
    catch( boost::bad_lexical_cast err ){ std::cout << "ERROR: " << err.what() << " parsing option -k\n"; return -1; }

    try{ cmo->which_class = options->at(r++)->value; }
    catch( boost::bad_lexical_cast err ){ std::cout << "ERROR: " << err.what() << " parsing option -w\n"; return -1; }

    try{ cmo->out_style = boost::lexical_cast<int>( options->at(r++)->value); }
    catch( boost::bad_lexical_cast err ){ std::cout << "ERROR: " << err.what() << " parsing option -x\n"; return -1; }

    try{ cmo->file_name_out = options->at(r++)->value; }
    catch( boost::bad_lexical_cast err ){ std::cout << "ERROR: " << err.what() << " parsing option -o\n"; return -1; }

    try{
    	if(options->at(r++)->value.compare("T")==0)
	    	cmo->verbose=true;
	    if( cmo->class_a.length()==0 && cmo->class_b.length()==0){
		    std::cout << "ERROR: Must supply at least one of: {--class_a, --class_b}.  Bailing out...\n";
    		std::cout.flush();
            return -1;
	    }
    }
    catch( boost::bad_lexical_cast err ){ std::cout << "ERROR: " << err.what() << " parsing option -v\n"; return 0; }
	if(cmo->which_class.compare("a")!=0 && cmo->which_class.compare("b")!=0 && cmo->which_class.compare("both")!=0){
		std::cout << "ERROR: which_class must be one of {a, b, both}.  Bailing out...\n";
		std::cout.flush();
        return -1;
	}
    return 0;
}

int main(int argc, char *argv[]){

	vector<Option*>* options = new vector<Option*>();
	options->push_back( new Option("data_file", "d", "Path to file with raw data.  Defaults to expr.txt", "expr.txt", OPT_OPTIONAL));
	options->push_back( new Option("sample_file", "f", "Path to sample attributes file. Default: sample_attributes.txt", "sample_attributes.txt", OPT_OPTIONAL));
	options->push_back( new Option("gene_file", "g", "Path to gene attributes file. Default: gene_attributes.txt", "gene_attributes.txt", OPT_OPTIONAL));
	options->push_back( new Option("class_a_limit", "a", "Comma-delimited list of sample attrib=value restrictions", "", OPT_OPTIONAL));
	options->push_back( new Option("class_b_limit", "b", "Comma-delimited list of sample attrib=value restrictions", "", OPT_OPTIONAL));
	options->push_back( new Option("min_support", "s", "Minimum percent of samples to which a rule must apply, defaults to 5", "5", OPT_OPTIONAL));
	options->push_back( new Option("min_conf", "c", "Minimum confidence of rule", "", OPT_REQUIRED));
	options->push_back( new Option("min_imp", "i", "Minimum percent improvement over simpler rule, defaults to 1", "1", OPT_OPTIONAL));
	options->push_back( new Option("max_chi2", "p", "Maximum chi2 p-value for rules, defaults to 1.", "1", OPT_OPTIONAL));

	options->push_back( new Option("discretization", "m", "Discretization Method (SD, MAD, abs, per, none).  Default: none", "none", OPT_OPTIONAL));
	options->push_back( new Option("lower_disc", "l", "Lower limit on discretization", "0", OPT_OPTIONAL));
	options->push_back( new Option("upper_disc", "u", "Upper limit on discretization", "0", OPT_OPTIONAL));

	options->push_back( new Option("gene_limit", "r", "Comma-delimited list of gene attrib=value restrictions", "", OPT_OPTIONAL));
	options->push_back( new Option("max_depth", "n", "Maximum number of genes to combine for rules, default 4", "4", OPT_OPTIONAL));
	options->push_back( new Option("timeout", "k", "Time out in seconds, default indefinite", "0", OPT_OPTIONAL));
	options->push_back( new Option("which_class", "w", "Which class to evaluate, {a, b, both}.  Defaults to both.", "both", OPT_OPTIONAL));

	options->push_back( new Option("out_style", "x", "Output style (default is 1, 2 to print class members)", "1", OPT_OPTIONAL));
	options->push_back( new Option("out_file", "o", "Output file, default to standard output", "", OPT_OPTIONAL));
	options->push_back( new Option("verbose", "v", "Print out progress to screen, T or F.  Defaults to F.", "F", OPT_OPTIONAL));

	std::stringstream ss;
	ss << "miner\nDavid Quigley, Balmain Lab, UCSF\n\n";
	ss << "Performs Association Rule Mining\n\n";
	ss << "GENERAL USE\n";
	ss << "Pass --data_file, --sample_file, --gene_file to indicate the data set to analyze.\n";
	ss << "By default all samples are included; limit the samples by passing limits to\n";
	ss << "--class_a_limit and --class_b_limit with the format FOO=BAR or FOO!BAR to include only samples in\n";
	ss << "--sample_file where the column FOO has (doesn't have) the value BAR. Multiple \n";
	ss << "constraints are combined with the logical AND, using the syntax \"FOO=BAR,BAZ=BIM\".\n\n";
	ss << "Minimal levels of support --min_support, confidence --min_conf, and improvement --min_imp\n";
	ss << "are set, and maximal chi-squared p-value for results can be set with --max_chi2.\n\n";
	ss << "The number of rules that can be combined is controls via --max_depth (defaults to 4); this\n";
	ss << "has a profound impact on execution time. By default both classes are evaluated; if only one limit\n";
	ss << "is passed, then only one class is evaluated.\n";

	int retval=read_args(argc, argv, options, ss.str());

	if( retval != 0 ){
		if(retval == -2)
			std::cout << "ERROR: Could not read arguments, bailing out.\n";
		return(0);
	}

	ClassMinerOptions * cmo = new ClassMinerOptions();

    if( safe_parse(cmo, options) != 0 )
        exit(0);

	ClassifierDataset* data = new ClassifierDataset();
	Attributes* sa;
	Attributes* ga;
	std::vector<Ant*>* ants = new std::vector<Ant*>();
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
	int ant_size_before=0;
	int n_tests=0;
	int L=3;
    if(cmo->verbose){
        std::cout << "MESSAGE: S / C / I: " << cmo->min_sup << " / " << cmo->min_conf*100 << " / " << cmo->min_imp*100 << "\n";
        std::cout << "MESSAGE: Discretization: " << cmo->discretization << " +" << cmo->disc_upper  << " / -" << cmo->disc_lower << ".\n";
		std::cout.flush();
    }
	cmo->start_time = time(NULL);
	ClassMiner cm = ClassMiner();

	cm.append_L1_L2_antecedents(data, ants, cmo, n_tests );

	time_t tm = time(NULL);
	if(cmo->verbose){
		std::cout << "MESSAGE: Completed level 1 and 2 in " << tm - cmo->start_time << " seconds, finding " << ants->size() << " candidates.\n";
		std::cout.flush();
	}
	while(L<=cmo->max_depth && ant_size_before != (int) ants->size()){
		ant_size_before = (int) ants->size();
		cm.append_L_ants(data, ants, cmo, L, n_tests);
		tm = time(NULL);
		if(cmo->verbose){
			std::cout << "MESSAGE: Completed level " << L << " after " << tm - cmo->start_time << " seconds, " << ants->size() << " total candidates.\n";
			std::cout.flush();
		}
		++L;

	}
	cmo->n_tests = n_tests;
	cm.assign_t_stats(ants, data);
	cm.filter_ants(ants, data, cmo->min_conf, cmo->min_sup, cmo->min_imp, cmo->max_chi2 );
	tm = time(NULL);
	if(cmo->verbose){
		std::cout << "MESSAGE: Completed mining, " << tm - cmo->start_time << " seconds, " << ants->size() << " total candidates.\n";
		std::cout.flush();
	}
	if( cmo->file_name_out.size() > 0 )
		cm.print_ants_file(ants, data, cmo);
	else
		cm.print_ants_stdout(ants, data, cmo);

	return(0);
}
