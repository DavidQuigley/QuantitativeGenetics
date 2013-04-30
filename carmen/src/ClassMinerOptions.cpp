#include <string>
#include "ClassMinerOptions.h"
using namespace std;

ClassMinerOptions::ClassMinerOptions(){
	this->out_style=1;
	this->timeout=0;
	this->start_time=0;
	this->min_imp=0;
	this->mine_type = "All_Data";
	this->verbose = false;
	this->which_class="both";
}
