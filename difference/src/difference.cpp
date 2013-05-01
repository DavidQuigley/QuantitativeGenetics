#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
using namespace std;
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/random.hpp>
#include <boost/algorithm/string.hpp>
using namespace boost;
#include "DataStructures.h"
#include "Attributes.h"
#include "ClassMinerOptions.h"
#include "ParseOptions.h"
#include "Rawdata.h"
#include "Discretize.h"
#include "Dataset.h"
#include "Perm.h"


int main(int argc, char *argv[]){

	vector<Option*>* options = new vector<Option*>();
	options->push_back( new Option("data_file", "d", "Path to file with raw data, default: expr.txt", "expr.txt", OPT_OPTIONAL));
	options->push_back( new Option("sample_file", "f", "Path to sample attributes file, default: sample_attributes.txt", "sample_attributes.txt", OPT_OPTIONAL));
	options->push_back( new Option("gene_file", "g", "Path to gene attributes file, default: gene_attributes.txt", "gene_attributes.txt", OPT_OPTIONAL));
	options->push_back( new Option("symbol_column", "y", "Column in gene attributes file indicating symbol to display, defaults 'Gene Name'", "Gene Name", OPT_OPTIONAL));
	options->push_back( new Option("class_a", "a", "Comma-delimited list of attrib=value restrictions", "", OPT_OPTIONAL) );
	options->push_back( new Option("class_b", "b", "Comma-delimited list of attrib=value restrictions", "", OPT_OPTIONAL));

	// Discretization is only required if the user is passing in a gene: limit
	options->push_back( new Option("discretization", "m", "Discretization Method (SD, MAD, abs, none).  Default: none", "none", OPT_OPTIONAL));
	options->push_back( new Option("lower_disc", "l", "Lower limit on discretization", "0", OPT_OPTIONAL));
	options->push_back( new Option("upper_disc", "u", "Upper limit on discretization", "0", OPT_OPTIONAL));

	options->push_back( new Option("n_perm", "n", "Number of permutations, default:1000", "1000", OPT_OPTIONAL));
	options->push_back( new Option("percent_present", "r", "Require this fraction present in each group, default 0.9", "0.9", OPT_OPTIONAL));

    options->push_back( new Option("mean_trim", "t", "Percent to trim off of each end of mean, default:5", "5", OPT_OPTIONAL));
	options->push_back( new Option("difference_file", "i", "Path to file difference file to load for conversion", "", OPT_OPTIONAL));
	options->push_back( new Option("p_val", "p", "Maximum p-value to write out.  Defaults to 1.", "1.0", OPT_OPTIONAL));
	options->push_back( new Option("output_file", "o", "Path to output file", "", OPT_OPTIONAL));
	options->push_back( new Option("verbose", "v", "Print out progress to screen, T or F.  Defaults to F.", "F", OPT_OPTIONAL));

	std::stringstream ss;
	ss << "difference\nDavid Quigley, Balmain Lab, UCSF\n\n";
	ss << "Calculates t-tests for a dataset, comparing two classes\n\n";
	ss << "GENERAL USE\n";
	ss << "GENERAL USE\n";
	ss << "Pass --data_file, --sample_file, --gene_file to indicate the data set to analyze.\n";
	ss << "By default all samples are included; limit the samples by passing limits to\n";
	ss << "--class_a and --class_b with the format FOO=BAR or FOO!BAR to include only samples in\n";
	ss << "--sample_file where the column FOO has (doesn't have) the value BAR. Multiple \n";
	ss << "constraints are combined with the logical AND, using the syntax \"FOO=BAR,BAZ=BIM\".\n\n";
	ss << "By default results are mean-trimmed at 5% with --mean_trim. If --n_perm is passed, perform\n";
	ss << "permutation testing. If permutations are tun, limit reported results with --p_val. Note that\n";
	ss << "this only controls the comparison-wise error-rate, not the experiment-wise error rate.\n";

	int retval=read_args(argc, argv, options, ss.str());
	if( retval != 0 )
		return(0);

	int r=0;
	ClassMinerOptions* cmo = new ClassMinerOptions();
	cmo->file_name_dis = options->at(r++)->value;
	cmo->file_name_sa = options->at(r++)->value;
	cmo->file_name_ga = options->at(r++)->value;
	std::string symbol_column = options->at(r++)->value;
	cmo->class_a = options->at(r++)->value;
	cmo->class_b = options->at(r++)->value;
	boost::algorithm::replace_all(cmo->class_a, std::string("*"), std::string("!"));
	boost::algorithm::replace_all(cmo->class_b, std::string("*"), std::string("!"));
	cmo->discretization = options->at(r++)->value;
	cmo->disc_lower = boost::lexical_cast<float>(options->at(r++)->value);
	cmo->disc_upper = boost::lexical_cast<float>(options->at(r++)->value);
	int n_perm = boost::lexical_cast<int>( options->at(r++)->value );
	double fraction_required = boost::lexical_cast<double>( options->at(r++)->value );
    int mean_trim = boost::lexical_cast<int>( options->at(r++)->value );
	std::string difference_file = options->at(r++)->value;
	double max_p_value = boost::lexical_cast<double>( options->at(r++)->value );
	cmo->file_name_out = options->at(r++)->value;
	bool verbose = false;
	if( options->at(r++)->value.compare("T")==0 )
		verbose = true;
	ClassifierDataset* data = new ClassifierDataset();
	Attributes* sa;
	Attributes* ga;
	Perm pick;
	pick.n_perm = n_perm;
	pick.mean_trim = mean_trim;
	pick.max_p_value = max_p_value;
	if(difference_file.length()>0){
		pick.load_settings(difference_file, max_p_value);
		cmo->class_a = pick.class_a;
		cmo->class_b = pick.class_b;
		cmo->file_name_dis = pick.file_name_dis;
		cmo->file_name_ga = pick.file_name_ga;
		cmo->file_name_sa = pick.file_name_sa;
	}
	try{
		if(verbose){
			std::cout << "MESSAGE: Loading data set...\n";
			std::cout.flush();
		}
		sa = new Attributes("NA");
		ga = new Attributes("NA");
		sa->load(cmo->file_name_sa);
		ga->load(cmo->file_name_ga);
		if( symbol_column.size() > 0 )
			ga->set_gene_name_column(symbol_column);
		data->load(sa, ga, cmo);
		if(verbose){
			std::cout << "MESSAGE: Loaded data set.\n";
			std::cout << "MESSAGE: Loaded data set.  Class A has " << data->a_idx->size() << " members, Class B has " << data->b_idx->size() << " members\n";
			std::cout.flush();
		}
	}
	catch(string msg){
		std::cout << "ERROR: " << msg << "\n";
		exit(0);
	}
	if( (int)data->a_idx->size()<2 ){
		std::cout << "ERROR: Class A has only " << data->a_idx->size() << " elements, cannot calculate statistic.\n";
		exit(0);
	}
	if( (int)data->b_idx->size()<2 ){
		std::cout << "ERROR: Class B has only " << data->b_idx->size() << " elements, cannot calculate statistic.\n";
		exit(0);
	}
    if( ga->identifiers.size() != data->raw_data->identifiers.size() ){
		std::cout << "ERROR: Gene attributes and raw data file do not have same set of identifiers.\n";
		exit(0);
	}
    std::vector<int> idx;
    for(int i=0; i<(int)ga->identifiers.size(); i++)
        idx.push_back(i);
    int n_before_NA = (int)idx.size();
    pick.limit_ids_by_NA(idx, data, fraction_required);
    int n_after_NA = (int)idx.size();

    if(verbose){
        std::cout << "MESSAGE: Requiring " << fraction_required*100 << " percent of samples present for each group.\n";
        std::cout << "MESSAGE: Removed " << n_before_NA - n_after_NA << " identifiers for insufficient number of data points.\n";
        std::cout.flush();
    }

	if(difference_file.length()>0){
		pick.write_as_ruleset( cmo, data->raw_data, data->a_idx, data->b_idx);
	}
	else{
		pick.permutations( data->raw_data, idx, data->a_idx, data->b_idx, ga, verbose);
		pick.write( cmo );
	}
	delete cmo;
	delete sa;
	delete ga;
	delete data;
	return 0;
}
