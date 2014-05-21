#include <iostream>
#include <sstream>
#include <vector>
#include <list>
#include <string>
#include <algorithm>
#include <time.h>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include "boost/algorithm/string.hpp"
#include <boost/unordered_map.hpp>
#include <boost/random.hpp>
#include <boost/thread.hpp>
#include <boost/bind.hpp>
#include <boost/thread/mutex.hpp>

namespace alg = boost::algorithm;
using namespace boost;
using namespace std;
#include "DataStructures.h"
#include "ClassMinerOptions.h"
#include "Rawdata.h"
#include "Attributes.h"
#include "Discretize.h"
#include "Dataset.h"
#include "ParseOptions.h"
#define MATHLIB_STANDALONE 1
#include "Parser.h"
#include "QTL_Calculator.h"
#include "Rmath.h"

void crash(std::string msg){
    std::cout << "ERROR: " << msg << "\n";
    std::cout.flush();
    exit(0);
}

void say_message(std::string msg, bool verbose ){
    if(verbose){
        std::cout << "MESSAGE: " << msg << "\n";
        std::cout.flush();
    }
}


void try_cast( double& val, std::string input, std::string error_msg ){
    try{ 
        val = boost::lexical_cast<double>(input); 
    }
    catch(boost::bad_lexical_cast &){
        crash(error_msg);
    }
}

void try_cast( int& val, std::string input, std::string error_msg ){
    try{ 
        val = boost::lexical_cast<int>(input); 
    }
    catch(boost::bad_lexical_cast &){
        crash(error_msg);
    }
}

void try_extract_columns( std::string& locus_columns, Attributes* ga_expr, Attributes* ga_snps, std::string& chr_expr, std::string& chr_snps, std::string& loc_expr, std::string& loc_snps ){
    // move this code out of main() just to make main() easier to read
    std::vector<std::string> columns;
    alg::split(columns, locus_columns, alg::is_any_of(",") );
    if( int(columns.size()) != 4 ){
        crash(std::string("ERROR: Value for --cis_columns must be four column names separated by commas"));
    }
    if( ga_expr->attribute2idx.find(columns.at(0)) == ga_expr->attribute2idx.end() ){
        crash(std::string("ERROR: Value for phenotype chromosome (first element in --cis_columns) not found in phenotype gene attributes file"));
    }
    if( ga_expr->attribute2idx.find(columns.at(1)) == ga_expr->attribute2idx.end() ){
        crash(std::string("ERROR: Value for phenotype locus (second element in --cis_columns) not found in phenotype gene attributes file"));
    }            
    if( ga_snps->attribute2idx.find(columns.at(2)) == ga_snps->attribute2idx.end() ){
        crash(std::string("ERROR: Value for genotype chromosome (third element in --cis_columns) not found in genotype gene attributes file"));
    }
    if( ga_snps->attribute2idx.find(columns.at(3)) == ga_snps->attribute2idx.end() ){
        crash(std::string("ERROR: Value for genotype locus (fourth element in --cis_columns) not found in genotype gene attributes file"));
        exit(0);
    }
    chr_expr = columns.at(0);
    loc_expr = columns.at(1);
    chr_snps = columns.at(2);
    loc_snps = columns.at(3);
}


int main(int argc, char *argv[]){

    vector<Option*>* options = new vector<Option*>();
	options->push_back( new Option("pheno_data_file", "d", "Path to file with phenotype data", "", OPT_REQUIRED));
	options->push_back( new Option("pheno_sample_file", "f", "Path to phenotype sample attributes file.", "", OPT_REQUIRED));
	options->push_back( new Option("pheno_attr_file", "g", "Path to phenotype attributes file.", "", OPT_REQUIRED));
	options->push_back( new Option("expr_limit", "a", "Comma-delimited list of sample attrib=value restrictions, default to use all", "", OPT_OPTIONAL));
	options->push_back( new Option("geno_file", "e", "Path to file with genotype data", "", OPT_REQUIRED));
	options->push_back( new Option("geno_sample_file", "h", "Path to genotype sample attributes file.", "", OPT_REQUIRED));
	options->push_back( new Option("geno_attr_file", "s", "Path to genotype attributes file.", "", OPT_REQUIRED));
    options->push_back( new Option("n_perms_max", "p", "Total number of permutations to run if P=0 after n_perms permutations, default 0", "0", OPT_OPTIONAL));
    options->push_back( new Option("fraction_required", "b", "Fraction of phenotype samples that must be present, default 0.9", "0.9", OPT_OPTIONAL));
    options->push_back( new Option("gene_fract", "t", "Limit analysis to subset N of X chunks, pass N,X, default 1,1", "1,1", OPT_OPTIONAL));
    options->push_back( new Option("cis_interval", "i", "Calculate cis-eQTL only, interval around pheno_locus for SNPs, default 0 (genome-wide search)", "0", OPT_OPTIONAL));
    options->push_back( new Option("by_chr", "k", "Calculate per-chromosome results, {T,F}, default F", "F", OPT_OPTIONAL));
    options->push_back( new Option("locus_columns", "l", "Column names for gene and SNP locations: pheno_chromosome,pheno_locus,geno_chromosome,geno_locus", "", OPT_OPTIONAL));
    options->push_back( new Option("n_perms", "n", "Number of permutations, default 0", "0", OPT_OPTIONAL));
    options->push_back( new Option("networks", "j", "Path to probe network key-value file, tab delimited, probes line 1 SNPs line 2", "", OPT_OPTIONAL));
    options->push_back( new Option("gene_snp", "m", "Path to gene-snp file, tab delimited, gene probe col 1 snp probe col 2 ", "", OPT_OPTIONAL));    
    options->push_back( new Option("min_per_geno", "q", "Minimum number of samples in each genotype for robust analysis, default 0", "0", OPT_OPTIONAL));
    options->push_back( new Option("probe", "r", "phenotype probe to use", "", OPT_OPTIONAL));
    options->push_back( new Option("id_match_col", "c", "Column name in geno_sample_file with matching pheno_id. Defaults to IDENTIFIER.", "IDENTIFIER", OPT_OPTIONAL));
	options->push_back( new Option("symbol_column", "y", "Column in pheno_attr_file with gene symbol, defaults symbol", "symbol", OPT_OPTIONAL));
	options->push_back( new Option("genetic_model", "u", "Use regression (additive) or ttest (recessive/dominant), default regression,", "regression", OPT_OPTIONAL));    
    options->push_back( new Option("n_threads", "w", "Number of threads on which to execute, default 1", "1", OPT_OPTIONAL));
    options->push_back( new Option("out_file", "o", "Output file, default to standard output", "", OPT_OPTIONAL));
    
	std::stringstream ss;
    
	ss << "eqtl\nDavid Quigley, Balmain Lab, UCSF\n\n";
	ss << "Performs eQTL analysis\n\n";
	ss << "GENERAL USE\n";
	ss << "Pass --pheno_data_file, --pheno_sample_file, --pheno_gene_file to indicate the expression\n";
	ss << "data set to analyze. Pass --geno_data_file, --geno_sample_file, --geno_gene_file to indicate\n";
	ss << "the genotype data set. To indicate which genotype in geno_sample_file matches the phenotype\n";
	ss << "in pheno_sample_file, pass --id_match_col. This is the column name in geno_sample_file where\n";
	ss << "the values will be found in pheno_id. Defaults to IDENTIFIER.\n\n";
	ss << "By default all samples are included; limit the samples by passing a limit to\n";
	ss << "--expr_limit with the format FOO=BAR or FOO!BAR to include only samples in\n";
	ss << "--pheno_sample_file where the column FOO has (doesn't have) the value BAR. Multiple \n";
	ss << "constraints are combined with the logical AND, using the syntax \"FOO=BAR,BAZ=BIM\".\n\n";
	ss << "Set the number of permutations with --n_perms.\n";

	int retval=read_args( argc, argv, options, ss.str());
	if( retval==-2 )
		return(0);
	else if( retval != 0 ){
        crash( std::string( "Could not read arguments, bailing out.") );
	}
	ClassMinerOptions * cmo_expr = new ClassMinerOptions();
	ClassMinerOptions * cmo_snps = new ClassMinerOptions();
	int r=0;
	cmo_expr->file_name_dis = options->at(r++)->value;
	cmo_expr->file_name_sa = options->at(r++)->value;
	cmo_expr->file_name_ga = options->at(r++)->value;
	cmo_expr->class_a = options->at(r++)->value;
    if(cmo_expr->class_a.size()==0)
        cmo_expr->class_a = std::string("IDENTIFIER!NULL");
	cmo_snps->file_name_dis = options->at(r++)->value;
	cmo_snps->file_name_sa = options->at(r++)->value;
	cmo_snps->file_name_ga = options->at(r++)->value;
    int n_perms_max = false;
    try_cast( n_perms_max, options->at(r++)->value, std::string( "Value for --n_perms_max must be an integer" ) );
    
    double fraction_required;
    try_cast(fraction_required, options->at(r++)->value, std::string("Value for --fraction_required must be an number"));
    
    std::string gene_fract( options->at(r++)->value );
    std::vector<std::string> fractions;
    alg::split(fractions, gene_fract, alg::is_any_of(",") );
    if( int(fractions.size()) != 2 ){
        crash( std::string("Value for --gene_fract must be two integers separated by commas."));
    }
    int this_fraction, n_fractions;
    try_cast(this_fraction, fractions.at(0), std::string("Value for first element of --gene_fract must be an integer"));
    try_cast(n_fractions, fractions.at(1), std::string("Value for second element of --gene_fract must be an integer"));
    
    int cis_interval;
    try_cast( cis_interval,options->at(r++)->value, std::string("Value for --cis_interval must be an integer") );
    
    bool per_chromosome=false;
    if( (options->at(r++)->value).compare("T") == 0 )
        per_chromosome = true;
    
    std::string locus_columns = options->at(r++)->value;
    int n_perms;
    try_cast( n_perms, options->at(r++)->value, std::string("Value for --n_perms must be an integer") );

    std::string fn_probe_networks = options->at(r++)->value;
    std::string fn_gene_snp = options->at(r++)->value;
    
    int min_per_geno=0;
    try_cast( min_per_geno, options->at(r++)->value, std::string("Value for --min_per_geno must be an integer"));
    
    std::string probe_id = options->at(r++)->value;
	std::string shared_geno_pheno_col = options->at(r++)->value;
    std::string symbol_col = options->at(r++)->value;
    std::string model = options->at(r++)->value;
    if( model.compare(std::string("regression") )!=0 && model.compare(std::string("ttest") )!=0 )
        crash(std::string("genetic_model must be regression or ttest") );
    
    int n_threads=1;
    try_cast( n_threads, options->at(r++)->value, std::string("Value for --n_threads must be an integer") );
    
	cmo_snps->file_name_out = options->at(r++)->get_value();
	cmo_expr->file_name_out = cmo_snps->file_name_out;
	cmo_expr->verbose=true;
	cmo_snps->verbose=true;
	cmo_expr->discretization="none";
	cmo_snps->discretization="none";

    Dataset* data_expr = new Dataset();
    Dataset* data_snps = new Dataset();
    Attributes* sa_expr, * ga_expr, * sa_snps, *ga_snps;

	try{
        say_message(std::string( "*** Loading data ***"), cmo_expr->verbose);
        say_message(std::string( "Reading genotypes from " + cmo_snps->file_name_dis ), cmo_expr->verbose);
        say_message(std::string( "Reading phenotypes from " + cmo_expr->file_name_dis ), cmo_expr->verbose);
		sa_expr = new Attributes("NA");
		ga_expr = new Attributes("NA");
        sa_expr->load(cmo_expr->file_name_sa);
		ga_expr->load(cmo_expr->file_name_ga);       
		ga_expr->set_gene_name_column(symbol_col);
        sa_snps = new Attributes("NA");
		ga_snps = new Attributes("NA");
		sa_snps->load(cmo_snps->file_name_sa);
		ga_snps->load(cmo_snps->file_name_ga);
        
        std::vector<std::string> columns;
        if( locus_columns.size()>0 ){
            std::string chr_expr, chr_snps, loc_expr, loc_snps;
            try_extract_columns( locus_columns, ga_expr, ga_snps, chr_expr, chr_snps, loc_expr, loc_snps );
            ga_expr->set_chromosome_and_locus(chr_expr, loc_expr);
            ga_snps->set_chromosome_and_locus(chr_snps, loc_snps);
            
        }
        if( ( per_chromosome || cis_interval != 0) && locus_columns.size()==0 ){
            crash(std::string( "Cannot specify cis-only analysis or per-chromosome analysis without passing --cis_columns parameter."));
        }
        if( per_chromosome && cis_interval != 0 ){
            crash( std::string("Cannot specify both cis-only analysis AND per-chromosome analysis simultaniously") );
        }
		data_expr->load(sa_expr, ga_expr, cmo_expr);
        say_message( std::string("Loaded phenotype data. Loading genotypes..."), cmo_expr->verbose );
		data_snps->load(sa_snps, ga_snps, cmo_snps);
        say_message( std::string("Loaded genotype data"), cmo_expr->verbose );

		if( data_expr->raw_data->identifiers.size() != ga_expr->identifiers.size() ) {
			crash(std::string("Raw data and gene attributes files for quantitiative data have a different number of identifiers"));
		}
		if( data_snps->raw_data->identifiers.size() != ga_snps->identifiers.size() ) {
			crash(std::string("Raw data and gene attributes files for genotype data have a different number of identifiers"));
		}
        if( data_expr->a_idx->size()==0 ){
            crash(std::string("No expression samples selected by class limit " + cmo_expr->class_a ) );
        }
        say_message( std::string("*** Data load complete ***"), cmo_snps->verbose );
	}
	catch(std::string msg){
		crash( msg );
	}
     
    QTL_Calculator* QC;
    try{
        QC = new QTL_Calculator(data_expr, sa_expr, ga_expr, data_snps, sa_snps, ga_snps, shared_geno_pheno_col);
    }
    catch(std::string msg){
        crash(msg);
    }
    ss.str("");
    ss << "Starting " << n_threads << " worker threads";
    say_message( ss.str(), cmo_snps->verbose );
    QC->set_thread_count(n_threads);
    if( cmo_snps->file_name_out.size()>0 ){
        ss.str("");
        ss << cmo_snps->file_name_out << ".log";
        try{
            QC->set_logfile( ss.str() );
        }
        catch( std::string msg){
            crash( msg );
        }
        say_message( std::string("LOGFILE will be written to " + ss.str() ), cmo_snps->verbose);
    }
    
    if( min_per_geno>0 ){
        QC->set_robust(true);
        QC->set_min_obs_per_genotype(min_per_geno);
        ss.str("");
        ss << "ROBUST ANALYSIS: For all T statistics, recalculating if there are not at least\n";
        ss << "           " << min_per_geno << " samples per genotype. Using the LOWER T statistic of the two results";
        say_message( ss.str(), cmo_snps->verbose );
    }
    if( QC->get_n_matched_samples() == 0 ){
        crash( std::string("No samples were matched between genotypes and quantitative traits.  Check id_match_col parameter.") );
    }
    if(n_perms_max>0 && n_perms==0){
        crash(std::string("Requested stepwise permutations but did not specify number of permutations") );
    }
    
    QC->set_n_perms(n_perms);
    if(n_perms_max==0)
        n_perms_max=n_perms;
    QC->set_n_perms_max(n_perms_max);
    
	// select subset of transcripts (either one probe_id or limited by variance and fraction present)
	// calculate genome-wide eQTL for all transcripts
	// store result on filesystem or stdio
    
	if( probe_id.size() > 0 ){
		// Restrict analysis to a single probe
        if( data_expr->raw_data->identifier2idx.find( probe_id ) == data_expr->raw_data->identifier2idx.end() ){
            crash(std::string( "Specified probe " + probe_id + " not found in phenotype identifiers") );
        }
        QC->limit_ids_by_probe_id(probe_id);
        std::string gene_name_col = ga_expr->get_gene_name_column();
        say_message(std::string("ANALYSIS TYPE: Restrict to single probe"), cmo_snps->verbose);
        ss.str("");
        ss << "Limiting to probe " << probe_id << ", symbol ";
        ss << ga_expr->prop_for_identifier( probe_id, gene_name_col );
        say_message( ss.str(), cmo_snps->verbose);
	}
	else{
		say_message(std::string("*** Sample and gene restrictions ***"), cmo_snps->verbose);

        // Apply other user-specified limits        
		int len_before_NA = QC->get_size_gene_idx();
		QC->limit_ids_by_NA(fraction_required); // only allow transcripts present at least fraction_required of the time
		int len_after_NA = QC->get_size_gene_idx();
        ss.str("");
        ss << "Removed " << len_before_NA - len_after_NA << " phenotypes because of " << fraction_required*100 << " % NA stringency";
        say_message( ss.str(), cmo_snps->verbose );
        
        // use only genes in this_fraction of n_fractions. Allows user to break up processing into equal chunks.
        QC->limit_ids_by_subset(this_fraction, n_fractions);
        int len_after_gene_fractions = QC->get_size_gene_idx();
        if( len_after_gene_fractions != len_after_NA ){
            ss.str("");
            ss << "RESTRICTING analysis to group " << this_fraction << " of " << n_fractions << " divisions of genes";
            say_message( ss.str(), cmo_snps->verbose );
        }
        ss.str("");
        ss << "PERMUTATION: Performing " << n_perms << " iterations of each SNP to assess thresholds.";
        say_message( ss.str(), cmo_snps->verbose );
        
        if( n_perms_max>n_perms ){
            ss.str("");
            ss << "Performing a maximum of " << n_perms_max << " iterations on SNPs in the top 1% after " << n_perms << " iterations";    
            say_message(ss.str(), cmo_snps->verbose);
        }
        ss.str("");
        ss << "Calculating with " << len_after_gene_fractions << " expression transcripts and " << ga_snps->identifiers.size() << " SNPs";
        say_message(ss.str(), cmo_snps->verbose);
	}
    ss.str("");
    ss << "Matched " << QC->get_n_matched_samples() << " samples between genotypes and expression";
    say_message( ss.str(), cmo_snps->verbose );
    say_message( std::string("*** Restrictions complete ***"), cmo_snps->verbose );
        
    if( fn_gene_snp.size()>0 ){
        say_message(std::string("ANALYSIS TYPE: exact lists of gene-snp pairs"), cmo_snps->verbose);
        try{
            KeyValuesParser KV(fn_gene_snp);

            if(cmo_snps->verbose){
                std::vector<std::string> keys;
                KV.keys(keys);
                ss.str("");
                ss << "Specifying " << keys.size() << " gene probes and exact matching SNPs for permutation analysis";
                say_message( ss.str(), cmo_snps->verbose );
                if(locus_columns.size()>0)
                    say_message( std::string("Limiting permutations to chromosome where locus is found"), cmo_snps->verbose );                
                else
                    say_message( std::string("performing genome-wide permutation analysis"), cmo_snps->verbose );
            }
            QC->check_probes_in_KV( KV ); // confirm that all probes requested actually exist in the dataset
            QC->calculate_permutations_for_gene_probe_pairs( KV );
            QC->sort_QTLs();
            cmo_expr->version = "eQTL 1.0";
            QC->print_rqtl( cmo_expr, cmo_snps, 0, 0);
        }
        catch(std::string msg){
            crash(msg);
        }        
    }
    else if( fn_probe_networks.size()>0 ){
        say_message(std::string("ANALYSIS TYPE: All combinations of set lists of genes and SNPs"), cmo_snps->verbose);
        say_message( std::string("Specifying probes and SNPs\n"), cmo_snps->verbose );
        try{
            KeyValuesParser KV(fn_probe_networks);
            Matrix<double> results;
            std::vector<std::string> keys, keys_written, files_written;
            KV.keys(keys);
            std::vector<int> idx_snps;
            std::vector<std::string> values;
            std::string fn_base = cmo_expr->file_name_out;
            std::vector<std::string> probesets, snps;
            ss.str("");
            ss << "Performing analysis on " << keys.size() << " networks";
            say_message(ss.str(), cmo_snps->verbose);
            for( int i=0; i<int(keys.size()); i++){
                KV.get(keys.at(i), values);
                std::stringstream s1, s2;
                s1 << "Network " << keys.at(i) << " " << i+1 << " of " <<  keys.size();
                say_message( s1.str(), cmo_snps->verbose);
                if( int( values.size() != 2 ) ){
                    s2 << "Network " << keys.at(i) << " contains zero probes or zero SNPs; skipping";
                    say_message( s2.str(), cmo_snps->verbose);
                }
                else{
                    boost::algorithm::split(probesets, values.at(0), boost::algorithm::is_any_of(" ") );
                    boost::algorithm::split(snps, values.at(1), boost::algorithm::is_any_of(" ") );
                    QC->restrict_probesets_and_snps(probesets, snps, idx_snps); 
                    QC->calculate_specific_SNPs(results, idx_snps);
                    cmo_expr->file_name_out = fn_base + "_" + keys.at(i);
                    keys_written.push_back(keys.at(i));
                    files_written.push_back( fn_base + "_" + keys.at(i) );
                    QC->reportNetResults(results, idx_snps, cmo_expr, cmo_snps);
                }
            }
            QC->reportFileList( fn_base + "_filelist", keys_written, files_written );
            ss.str("");
            ss << "Wrote network file list to " << fn_base << "_filelist";
            say_message(ss.str(), cmo_snps->verbose);
        }
        catch(std::string msg){
            crash(msg);
        }
        
    }
    else if(per_chromosome){
        say_message(std::string("ANALYSIS TYPE: genome-wide by chromosome"), cmo_snps->verbose);
        say_message(std::string("Identifying best SNP for each chromosome for each probeset"), cmo_snps->verbose );
        QC->set_verbose(cmo_expr->verbose);
        QC->calculate_genome_wide(0, true);
        QC->print_by_chromosome( cmo_expr, cmo_snps, fraction_required);
    }
    else{
        if(cis_interval>0){
            say_message(std::string("ANALYSIS TYPE: CIS-eQTL only analysis"), cmo_snps->verbose);
            ss.str("");
            ss << "Performing CIS-only analysis with window size " << cis_interval;
            say_message( ss.str(), cmo_snps->verbose );
        }
        else{
            say_message(std::string("ANALYSIS TYPE: genome-wide"), cmo_snps->verbose);
            say_message(std::string("Identifying best SNP for each probeset"), cmo_snps->verbose );
        }
        QC->calculate_genome_wide(cis_interval, false); 
        QC->sort_QTLs();
        cmo_expr->version = "eQTL 1.0";
        QC->print_rqtl( cmo_expr, cmo_snps, 0, fraction_required);
    }
    delete QC;
    return 0;
}
