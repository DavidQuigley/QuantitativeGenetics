#define _SCL_SECURE_NO_DEPRECATE 1
#include <iostream>
#include <time.h>
#include <vector>
#include <set>
#include <string>
#include <algorithm>
#include "math.h"
using namespace std;
#include "boost/tokenizer.hpp"
#include "boost/lexical_cast.hpp"
#include "boost/algorithm/string.hpp"
#include <boost/random.hpp>
namespace alg = boost::algorithm;
using namespace boost;
#include "DataStructures.h"
#include "ClassMinerOptions.h"
#include "Rawdata.h"
#include "Attributes.h"
#include "Discretize.h"
#include "Dataset.h"
#include "ParseOptions.h"
#include "graph.h"
#include "Parser.h"
#include "spear.h"

/*
 * David Quigley
 * UCSF Cancer Research Institute, 2007-2013
 *
 * Spear is used to measure correlation between rows of elements.  Spearman
 * rank correlation is calculated; ties are ranked using the mean of the ranks
 * that are spanned, so (1,2,2) is ranked (1, 2.5, 2.5).
 * By passing in only one class, you can find genes that are
 * correlated in a single group.  By passing in two non-intersecting classes,
 * you can find genes whose correlation changes between two groups.
 */

int main(int argc, char *argv[]){

    vector<Option*>* options = new vector<Option*>();
	options->push_back( new Option("data_file", "d", "Path to file with raw data.  Defaults to expr.txt", "expr.txt", OPT_OPTIONAL));
	options->push_back( new Option("sample_file", "f", "Path to sample attributes file. Default: sample_attributes.txt", "sample_attributes.txt", OPT_REQUIRED));
	options->push_back( new Option("gene_file", "g", "Path to gene attributes file. Default: gene_attributes.txt", "gene_attributes.txt", OPT_REQUIRED));
	options->push_back( new Option("symbol_column", "y", "Column in gene attributes file indicating symbol to display, defaults 'Gene Name'", "", OPT_REQUIRED));
	options->push_back( new Option("class_a_limit", "a", "Sample definition for class A, comma sep. list attrib=value, defaults to all", "", OPT_OPTIONAL));
	options->push_back( new Option("class_b_limit", "b", "Sample definition for class A, comma sep. list attrib=value, defaults to none", "", OPT_OPTIONAL));
    options->push_back( new Option("corr_type", "t", "Method of correlation, {spearman,pearson}, default spearman", "spearman", OPT_OPTIONAL));
	options->push_back( new Option("probe", "p", "Restrict calculation to first degree neighbors of these seed probes, comma delimited", "", OPT_OPTIONAL));
	options->push_back( new Option("include_seed_neighbor_corr", "r", "If T, calculate correlation between seed neighbors, default F", "F", OPT_OPTIONAL));
	options->push_back( new Option("limit_to_seeds", "l", "If T, limit correlation to those between seeds default F", "F", OPT_OPTIONAL));
	options->push_back( new Option("min_var", "m", "Minimum variance across all samples", "0", OPT_OPTIONAL));
	options->push_back( new Option("percent_present", "n", "Require this fraction present in each group, default 0.9", "0.9", OPT_OPTIONAL));
    options->push_back( new Option("corr_diff", "c", "Minimum change in correlation between classes", "0", OPT_OPTIONAL));
	options->push_back( new Option("corr_abs_a", "s", "Minimum absolute correlation of class A", "0", OPT_OPTIONAL));
	options->push_back( new Option("corr_abs_b", "h", "Minimum absolute correlation of class B", "0", OPT_OPTIONAL));
	options->push_back( new Option("min_zscore", "x", "Minimum absolute Z score. Default 0", "0", OPT_OPTIONAL));
    options->push_back( new Option("perms", "z", "Number of permutations, defaults to 0. Return max rho per perm.", "0", OPT_OPTIONAL));
    options->push_back( new Option("distribution", "e", "If T, compute distribution of DC values. Defaults F.", "F", OPT_OPTIONAL));
    options->push_back( new Option("fn_out", "o", "Output file, default to standard output", "", OPT_OPTIONAL));
	options->push_back( new Option("verbose", "v", "Print out progress to screen, T or F.  Defaults to F.", "F", OPT_OPTIONAL));

	std::stringstream ss;
	ss << "SPEAR\nDavid Quigley, Balmain Lab, UCSF\n\n";

	ss << "Calculates correlation and returns a file in the .SPEAR format\n\n";
	ss << "GENERAL USAGE\n";
	ss << "Pass --data_path, --sample_path, --gene_path to indicate the data set to analyze.\n";
	ss << "By default all samples are included; limit the samples by passing limits to\n";
	ss << "--class_a_limit with the format FOO=BAR or FOO!BAR to include only samples in\n";
	ss << "--sample_file where the column FOO has (doesn't have) the value BAR. Multiple \n";
	ss << "constraints are combined with the logical AND, using the syntax \"FOO=BAR,BAZ=BIM\".\n";
	ss << "If a second limit is passed in --class_b_limit, correlation is calculated in both \n";
	ss << "classes. Results can be filtered for differential correlation using --corr_diff\n";
	ss << "--class_a_limit and --class_b_limit must be disjoint.\n";
	ss << "Limit analysis to probes with at variance >= using --min_var.\n";
	ss << "Limit analysis to samples present in at least --percent_present samples\n";
	ss << "Limit absolute value of correlation coefficient by passing --min_abs_corr.\n\n";

	ss << "LIMITING PROBES USING SEED_PROBES\n";
	ss << "If --seed_probes is not passed, all probes not excluded by--min_var are included.\n";
	ss << "If --seed_probes is passed with a comma-delimited list of seed identifiers,\n";
	ss << "analysis is limited to seeds and their immediate neighbors.\n";
	ss << "By default, only correlations between seeds and their neighbors are reported (excluding \n";
	ss << "correlations between neighbors). To include correlations between neighbors, pass\n";
	ss << "both --seed_probes and --include_seed_neighbor_corr=T\n\n";

	ss << "ANALYSIS OF ONLY A SUBSET OF PROBES\n";
	ss << "To limit analysis to a subset of probes (and not include neighbors) pass --seed_probes\n";
	ss << "to specify the probes and --limit_to_seeds=T to restrict analysis to only those probes.\n\n";

	ss << "PERMUTATION ANALYSIS\n";
	ss << "To perform a permutation analysis (to generate the experiment-wide GWER) pass a number\n";
	ss << "to --perms. The resulting spear file will contain bogus entries for genes; the only value\n";
	ss << "of use is the correlation coefficient for class A; this represents the single strongest\n";
	ss << "correlation result in all probe comparisons for one permutation of the true data. To\n";
	ss << "identify the 5% GWER, report the 50th-highest value when --perms=1000\n";

	int retval=read_args(argc, argv, options, ss.str());
	if( retval==-2 )
		return(0);
	else if( retval != 0 ){
		std::cout << "ERROR: Could not read arguments, bailing out.\n";
		return(0);
	}

	Spearman sp;
	int r=0;

	std::string fn_data = options->at(r++)->value;
	std::string fn_sa = options->at(r++)->value;
	std::string fn_ga = options->at(r++)->value;
	std::string symbol_column = options->at(r++)->value;
	sp.set_input_files(fn_data, fn_sa, fn_ga);
	if( symbol_column.size() > 0 )
		sp.set_gene_name_column(symbol_column);
	std::string limit_a( options->at(r++)->value);
	std::string limit_b( options->at(r++)->value);
	boost::algorithm::replace_all(limit_a, std::string("*"), std::string("!"));
	boost::algorithm::replace_all(limit_b, std::string("*"), std::string("!"));
	sp.set_limit_a( limit_a );
	sp.set_limit_b( limit_b );

	sp.set_method( options->at(r++)->value );
	std::vector<std::string> seeds;
	std::string seed_str = options->at(r++)->value;
	if( seed_str.size()>0 )
		boost::algorithm::split(seeds, seed_str, boost::algorithm::is_any_of(",") );
	if( seeds.size()>0 )
		sp.set_seeds( seeds );

	bool include_seed_neighbors=false;
	bool limit_network_to_seeds=false;
	if(options->at(r++)->value.compare("T")==0)
		include_seed_neighbors = true;
	if(options->at(r++)->value.compare("T")==0)
		limit_network_to_seeds = true;
	if( (include_seed_neighbors || limit_network_to_seeds) && seeds.size()==0 ){
		std::cout << "ERROR: must pass seed probes in order to limit network to seeds or include seed neighbors.";
		return 0;
	}
	sp.set_include_seed_neighbor_correlations(include_seed_neighbors);
	sp.set_limit_network_to_seeds(limit_network_to_seeds);

	double percent_required, min_var, delta_corr, abs_corr_a, abs_corr_b, min_zscore;
    try{ min_var = boost::lexical_cast<double>( options->at(r++)->value ); }
    catch( boost::bad_lexical_cast &){ std::cout << "ERROR: invalid parameter for -m, min_var"; return 0;}

    try{ percent_required = boost::lexical_cast<double>( options->at(r++)->value ); }
    catch( boost::bad_lexical_cast &){ std::cout << "ERROR: invalid parameter for -n, percent_present"; return 0;}

    try{ delta_corr = boost::lexical_cast<double>( options->at(r++)->value ); }
    catch( boost::bad_lexical_cast &){ std::cout << "ERROR: invalid parameter for -c, corr_diff"; return 0;}

    try{ abs_corr_a = boost::lexical_cast<double>( options->at(r++)->value ); }
    catch( boost::bad_lexical_cast &){ std::cout << "ERROR: invalid parameter for -s, corr_abs_a"; return 0;}

    try{ abs_corr_b = boost::lexical_cast<double>( options->at(r++)->value ); }
    catch( boost::bad_lexical_cast &){ std::cout << "ERROR: invalid parameter for -h, corr_abs_b"; return 0;}

    try{ min_zscore = boost::lexical_cast<double>( options->at(r++)->value ); }
    catch( boost::bad_lexical_cast &){ std::cout << "ERROR: invalid parameter for -x, min_zscore"; return 0;}

    int n_perms;
    try{ n_perms = boost::lexical_cast<int>( options->at(r++)->value ); }
    catch( boost::bad_lexical_cast &){ std::cout << "ERROR: invalid parameter for -z, perms"; return 0;}
    bool do_distribution = false;
    if(options->at(r++)->value.compare("T")==0)
    	do_distribution=true;
    sp.set_do_distribution(do_distribution);
    //sp.set_GWER_method(GWER_method);
    try{
		sp.set_min_var(min_var);
		sp.set_percent_required(percent_required);
		sp.set_corr_diff(delta_corr);
		sp.set_corr_abs_a(abs_corr_a);
		sp.set_corr_abs_b(abs_corr_b);
		sp.set_min_zscore(min_zscore);
		std::string fn_out = options->at(r++)->get_value();
		sp.set_fn_out( fn_out );
		sp.set_n_permutations( n_perms );
		if(options->at(r++)->value.compare("T")==0)
			sp.set_verbose(true);
		sp.run();
		if( sp.get_success() ){
			if( do_distribution ){
				if( fn_out.length() > 0 )
					sp.write_distribution();
				else
					sp.print_distribution();
			}
			else{
				if( fn_out.length() > 0 )
					sp.write_spears();
				else
					sp.print_spears();
			}
		}
    }

	catch(std::string err){
		std::cout << "ERROR: " << err << "\n";
	}
	catch(std::runtime_error rte){
		std::cout << "ERROR: " << rte.what() << "\n";
	}
    return 0;
}
