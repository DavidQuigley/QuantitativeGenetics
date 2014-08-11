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
#include "boost/thread/mutex.hpp"
//#include "boost/filesystem/operations.hpp"
//#include "boost/filesystem/path.hpp"
#include "boost/random.hpp"
namespace alg = boost::algorithm;
using namespace boost;
#include "DataStructures.h"
#include "ClassMinerOptions.h"
#include "Rawdata.h"
#include "Attributes.h"
#include "Discretize.h"
#include "Dataset.h"
#include "ParseOptions.h"
#include "Parser.h"
#include "graph.h"
#include "spear.h"



// Compiling this requires you to build the boost::filesystem library:
// bjam --with-filesystem threading=multi link=static runtime-link=static stage

int main(int argc, char *argv[]){
	vector<Option*>* options = new vector<Option*>();
	options->push_back( new Option("data_path", "d", "Path to file with raw data", "", OPT_REQUIRED));
	options->push_back( new Option("sample_path", "f", "Path to sample attributes file", "", OPT_REQUIRED));
	options->push_back( new Option("gene_path", "g", "Path to gene attributes file", "", OPT_REQUIRED));
	options->push_back( new Option("symbol_column", "y", "Column in gene attributes file indicating symbol to display, defaults 'Gene Name'", "", OPT_OPTIONAL));
	options->push_back( new Option("spear", "r", "Path to existing Spear file", "", OPT_OPTIONAL));
	options->push_back( new Option("fn_eqtl", "q", "Path to existing eQTL file", "", OPT_OPTIONAL));
	options->push_back( new Option("max_perm_pval", "k", "Maximum permutation p-value for eQTL, default 1", "1", OPT_OPTIONAL));
	options->push_back( new Option("snp_required", "j", "Restrict network to probes with significant eQTL, default F", "F", OPT_OPTIONAL));
	options->push_back( new Option("allow_uncorrelated_loci_with_eQTL", "u", "Allow genes with eQTL but no correlation in network, default T", "T", OPT_OPTIONAL));
	options->push_back( new Option("sample_limit", "a", "Sample definition for class A, comma sep. list attrib=value, defaults to all", "", OPT_OPTIONAL));
	options->push_back( new Option("seed_probes", "p", "Restrict calculation to first degree neighbors of these seed probes, comma delimited", "", OPT_OPTIONAL));
	options->push_back( new Option("include_seed_neighbor_corr", "e", "If T, calculate correlation between seed neighbors, default F", "F", OPT_OPTIONAL));
	options->push_back( new Option("limit_to_seeds", "l", "Restrict analysis to seed probes. Defaults to F", "F", OPT_OPTIONAL));
	options->push_back( new Option("min_var", "m", "Minimum variance across all samples. Defaults to 0", "0", OPT_OPTIONAL));
	options->push_back( new Option("percent_present", "n", "Require this fraction present in each group, default 0", "0", OPT_OPTIONAL));
    options->push_back( new Option("min_abs_corr", "s", "Minimum absolute correlation, default 0", "0", OPT_OPTIONAL));
    options->push_back( new Option("min_clique_size", "c", "Minimum size for cliques, default 1 {1,2,3,4}", "1", OPT_OPTIONAL));
    options->push_back( new Option("cytoscape_path", "t", "Path to cytoscape script", "", OPT_REQUIRED));
    options->push_back( new Option("vizmap_path", "b", "Path to vizmap properties file", "", OPT_REQUIRED));
    options->push_back( new Option("extra_gene_attribute", "x", "Name of numeric gene attribute to append to network", "", OPT_OPTIONAL));
    options->push_back( new Option("ontology_path", "h", "Path to ontology file", "", OPT_OPTIONAL));
    options->push_back( new Option("annotation_path", "i", "Path to gene annotation file", "", OPT_OPTIONAL));

    options->push_back( new Option("base_filename", "o", "Base name of output files", "", OPT_REQUIRED));
	options->push_back( new Option("verbose", "v", "Print out progress to screen, T or F.  Defaults to F.", "F", OPT_OPTIONAL));

	std::stringstream ss;
	ss << "CORRELATION NETWORK\nDavid Quigley, Balmain Lab, UCSF\n\n";

	ss << "Calculates correlation networks and generates input files for Cytoscape. The path\n";
	ss << "to Cytoscape and a Cytoscape vizmap file must be provided by passing --cytoscape_path\n";
	ss << "and vizmap_path. An existing .SPEAR file (the output of the spear program) can be used by\n";
	ss << "passing --spear. Numerous files are written as --base_filename plus various extensions.\n\n";

	ss << "CALCULATING NEW CORRELATIONS\n";
	ss << "To calculate new correlations instead of using an existing SPEAR file:\n";
	ss << "Pass --data_path, --sample_path, --gene_path to indicate the data set to analyze.\n";
	ss << "By default all samples are included; limit the samples by passing limits to\n";
	ss << "--sample_limit with the format FOO=BAR or FOO!BAR to include only samples in\n";
	ss << "--sample_file where the column FOO has (doesn't have) the value BAR. Multiple \n";
	ss << "constraints are combined with the logical AND, using the syntax \"FOO=BAR,BAZ=BIM\".\n";
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

	ss << "INCLUSION OF EQTL DATA\n";
	ss << "eQTL data can be added to the network by passing --fn_eqtl.\n";
	ss << "By default all eQTL are included; set a maximal permutation p-value with --max_perm_pval\n";
	ss << "If --snp_required=T is passed, only genes with a significant eQTL are plotted.\n";
	ss << "By default genes with eQTL but no correlations at the specified minimal correlation strength\n";
	ss << "will be plotted; these can be excluded by passing --allow_uncorrelated_loci_with_eQTL=F\n\n";

	ss << "EXTRA ATTRIBUTES: GENE ONTOLOGY AND GENE ATTRIBUTES\n";
	ss << "If --ontology_path and --annotation_path are passed, using the CARMEN annotation formats,\n";
	ss << "genes in the network will have attributes corresponding to their gene ontology annotations\n";
	ss << "To add a numeric column from the --gene_path file as an attribute (e.g. chromosome, if numeric)\n";
	ss << "pass --extra_gene_attribute\n";

	int retval=read_args(argc, argv, options, ss.str());
	if( retval==-2 )
		return(0);
	else if( retval != 0 ){
		std::cout << "ERROR: Could not read arguments, bailing out.\n";
		return(0);
	}

	Spearman* sp;
	int r=0;
	// LOAD DATASET
	std::string fn_data = options->at(r++)->value;
	std::string fn_sa = options->at(r++)->value;
	std::string fn_ga = options->at(r++)->value;
	std::string symbol_column = options->at(r++)->value;
	std::string fn_spear = options->at(r++)->value;
	std::string fn_eQTL = options->at(r++)->value;
	double max_perm_pval;
	try{ max_perm_pval = boost::lexical_cast<double>( options->at(r++)->value ); }
	catch( boost::bad_lexical_cast &){ std::cout << "ERROR: invalid parameter for -k, --max_perm_pval\n"; return 0;}
	bool require_eQTL = false;
	bool allow_uncorrelated_loci_with_eQTL = true;
	if(options->at(r++)->value.compare("T")==0)
		require_eQTL=true;
	if(options->at(r++)->value.compare("F")==0)
		allow_uncorrelated_loci_with_eQTL=false;
	std::string sample_limits = options->at(r++)->value;
	std::string seed_str = options->at(r++)->value;
	bool include_seed_neighbors=false;
	if(options->at(r++)->value.compare("T")==0)
		include_seed_neighbors=true;
	bool limit_to_seeds = false;
	if(options->at(r++)->value.compare("T")==0)
		limit_to_seeds = true;
	double min_percent_required, min_var, min_abs_corr;
	try{ min_var = boost::lexical_cast<double>( options->at(r++)->value ); }
	catch( boost::bad_lexical_cast &){ std::cout << "ERROR: invalid parameter for -m, min_var\n"; return 0;}
	try{ min_percent_required = boost::lexical_cast<double>( options->at(r++)->value ); }
	catch( boost::bad_lexical_cast &){ std::cout << "ERROR: invalid parameter for -n, percent_present\n"; return 0;}
	try{ min_abs_corr = boost::lexical_cast<double>( options->at(r++)->value ); }
	catch( boost::bad_lexical_cast &){ std::cout << "ERROR: invalid parameter for -s, min_abs_corr\n"; return 0;}
	int min_clique_size;
	try{ min_clique_size = boost::lexical_cast<int>( options->at(r++)->value ); }
	catch( boost::bad_lexical_cast &){ std::cout << "ERROR: invalid parameter for -c, min_clique_size\n"; return 0;}
	std::string fn_cytoscape= options->at(r++)->value;
	std::string fn_vizmap = options->at(r++)->value;
	std::string extra_ga_column = options->at(r++)->value;
	std::string fn_ontology = options->at(r++)->value;
	std::string fn_annotation = options->at(r++)->value;
	std::string fn_base = options->at(r++)->value;
	bool verbose = false;
	if(options->at(r++)->value.compare("T")==0)
		verbose = true;
	bool generate_new_correlations = true;
	if( fn_annotation.size()>0 && fn_ontology.size()==0){
		std::cout << "ERROR: if annotation_path is passed, must also pass ontology_path\n"; std::cout.flush();
		exit(0);
	}
	if(fn_annotation.size()==0 && fn_ontology.size()>0){
		std::cout << "ERROR: if ontology_path is passed, must also pass annotation_path\n"; std::cout.flush();
		exit(0);
	}
	Attributes* ga;
	try{
		ga = new Attributes("NA");
		ga->load(fn_ga);
	}
	catch(std::string msg){
		std::cout << "ERROR: could not load gene attribute file " << fn_ga << "\n";
		std::cout.flush();
		exit(0);
	}
	if( fn_spear.size()>0 ){
		generate_new_correlations = false;
		try{
            sp = new Spearman(fn_spear, ga);
        }
        catch(std::string msg){
            std::cout << "ERROR: Could not load spear file " << fn_spear << "\n";
            std::cout.flush();
            exit(0);
        }
		if( symbol_column.size() > 0 )
			ga->set_gene_name_column(symbol_column); // required because we're not going to call run()
	}
	else{
		try{
			sp = new Spearman();
			sp->set_input_files(fn_data, fn_sa, fn_ga);
			if( symbol_column.size() > 0 )
				sp->set_gene_name_column(symbol_column);
        }
		catch(std::string msg){ 
			std::cout << "ERROR: could not load dataset\n";
			std::cout.flush();
			exit(0);
		}
	}
	sp->set_verbose(verbose);
	sp->set_require_eqtl(require_eQTL); // defaults to false
	sp->set_fn_eqtl(fn_eQTL); // defaults to empty
	sp->set_max_eqtl_pval(max_perm_pval); 
    sp->set_allow_uncorrelated_loci_with_eQTL(allow_uncorrelated_loci_with_eQTL);
	// SAMPLE LIMITS
	boost::algorithm::replace_all(sample_limits, std::string("*"), std::string("!"));
	sp->set_limit_a(sample_limits);

	// SEED PROBES
	std::vector<std::string> seeds;
	if( seed_str.size()>0 )
		boost::algorithm::split(seeds, seed_str, boost::algorithm::is_any_of(",") );
	if( seeds.size()>0 )
		sp->set_seeds( seeds );
	if(include_seed_neighbors)
		sp->set_include_seed_neighbor_correlations(true);
	//sp->set_include_seed_neighbor_correlations(true);
	sp->set_limit_network_to_seeds(limit_to_seeds);

	// CORRELATION SETTINGS
	sp->set_min_var(min_var);
	sp->set_percent_required(min_percent_required);
	sp->set_corr_abs_a(min_abs_corr);

	// NETWORK SETTINGS 
	if( min_clique_size>0 && min_clique_size<5 )
		sp->set_min_clique_size(min_clique_size);
	else{
		std::cout << "ERROR: invalid minimum clique size, must be in {1,2,3,4}\n";
		std::cout.flush();
		exit(0);
	}
	if(generate_new_correlations){
		try{
			sp->run();
		}
		catch( std::string errmsg){
			std::cout << "ERROR: " << errmsg << "\n";
			std::cout.flush();
			exit(0);
		}
		catch(std::runtime_error re){
			std::cout << "ERROR: " << re.what() << "\n";
			std::cout.flush();
			exit(0);
		}
	}
    if(verbose){
        std::cout << "MESSAGE: Correlation calculations are complete.\n";
        std::cout.flush(); 
    }
#ifdef WIN32
	sp->set_batch_extension(std::string("bat"));
#else
	sp->set_batch_extension(std::string("sh"));
#endif
	sp->set_fn_cytoscape(fn_cytoscape);
	sp->set_fn_cytoscape_props(fn_vizmap);
	if( extra_ga_column.size()>0 )
		sp->set_column_extra_attributes( extra_ga_column );
	try{
		GeneAnnotationParser* gene_parser=NULL;
		GOAnnotationParser* GO_parser=NULL;
		if( fn_annotation.size()>0 && fn_ontology.size()>0){
			GO_parser = new GOAnnotationParser(fn_ontology);
			gene_parser = new GeneAnnotationParser(fn_annotation, *GO_parser);
		}
		sp->write_to_cytoscape(min_abs_corr, fn_base, ga, GO_parser, gene_parser, fn_eQTL, max_perm_pval, require_eQTL);
		if( GO_parser != NULL ){
			delete GO_parser;
			delete gene_parser;
		}
	} 
	catch( std::string* err ){
		std::cout << "ERROR: " << err->c_str() << "\n";
		std::cout.flush();
		exit(0);
	}
	catch( std::string err ){
		std::cout << "ERROR: " << err << "\n";
		std::cout.flush();
		exit(0);
	}
	delete sp;
	if( fn_spear.size()==0 ){
		// if ga was passed into Spearman, it will be deleted and set to NULL by spearman when we call "delete sp"
		delete ga;
	}
}

