#include <iostream>
#include <vector>
#include <string>

#include "ParseOptions.h"

Option::Option(std::string long_name, std::string short_name, std::string description, std::string default_value, bool required){
	this->long_name = long_name;
	this->short_name = short_name;
	this->description = description;
	this->required = required;
	this->value = default_value;
	this->default_value = default_value;
	this->is_set=false;
}

void Option::set_value(std::string value){
	this->value = value;
	this->is_set = true;
}

std::string Option::get_value(){
	if( this->is_set)
		return this->value;
	else
		return this->default_value;
}

void print_usage(std::vector<Option*>* options, std::string preamble){
	if( preamble.size()>0)
		std::cout << preamble << "\n";
	std::cout << "Usage:\n";
	std::vector<Option*>::iterator iter;
	for( iter=options->begin(); iter != options->end(); iter++){
		std::cout << " " << (*iter)->short_name << ", " << (*iter)->long_name << "   " << (*iter)->description << "\n";
	}
}


int read_args(int argc, char* argv[], std::vector<Option*>* options, std::string preamble){
	int i;
	std::vector<Option*>::iterator iter;
	for(i=1; i<argc; i++){
		std::string arg = std::string(argv[i]);
		if(arg.length() < 2){
			std::cout << "ERROR: Argument number " << i << ", " << arg << ", cannot be parsed.\n";
			return(-1);
		}
		int idx_equals;
		bool found = false;
        
		if(arg[0]=='-'){
			if(arg[1]=='-'){
				idx_equals = (int)arg.find("=");
				std::string passed;
				passed = arg.substr(2,idx_equals-2);
				std::string value;
				if(idx_equals==-1)
					value = "";
				else
					value = arg.substr(idx_equals+1, (int)arg.length()-idx_equals);
				if( passed.compare("help")==0){
					print_usage(options, preamble);
					return(-2);
				}
				 
				found=false;
				for( iter=options->begin(); iter != options->end(); iter++){
					if( passed.compare( (*iter)->long_name) == 0 ){
						(*iter)->set_value(value);
						found = true;
					}
				}
				if( !found ){
					std::cout << "Unrecognizd option: " << passed << ".  Try using --help for more information.\n";
					return -1;
				}
			}
			else{ 
				std::string passed = arg.substr(1,1);
				std::string value =  arg.substr(2, (int)arg.length()-2);
                if( value.substr(0,1).compare("=")==0){
					std::cout << "ERROR: Cannot read = sign in flag when using short parameter name: " << passed << "\n";
					return -1;
				}

				found=false;
				for( iter=options->begin(); iter != options->end(); iter++){
					if( passed.compare( (*iter)->short_name) == 0 ){
						(*iter)->set_value(value);
						found=true;
					}
				}
				if( !found ){
					std::cout << "Unrecognized option: " << passed << ".  Try using --help for more information.\n";
					return -1;
				}
			}
		}
		else{
			std::cout << "ERROR: Unrecognized flag: " << arg << "\n";
			return -1;
		}
	}
	for( iter=options->begin(); iter != options->end(); iter++){
		if( (*iter)->required && !(*iter)->is_set){
			std::cout << "ERROR: Missing required option: " << (*iter)->long_name << ".  Try using --help for more information.\n";
			return -1;
		}
	}
	return 0;
}
	
