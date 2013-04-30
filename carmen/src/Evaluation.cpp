#include "Evaluation.h"

Evaluation::Evaluation(){
	// This is just a structure to hold the performance of a classifier
	this->n_train_a = 0;
	this->n_train_b = 0;
	this->n_test_a = 0;
	this->n_test_b = 0;
	this->n_train = 0;
	this->n_test = 0;
	this->c_train = 0;
	this->c_test = 0;
	this->c_train_a = 0;
	this->c_train_b = 0;
	this->c_test_a = 0;
	this->c_test_b = 0;
	this->w_test = 0;
	this->w_test_a = 0;
	this->w_test_b = 0;
	this->w_train = 0;
	this->w_train_a  = 0;
	this->w_train_b = 0;
	this->a_tr_all = 0;
	this->a_tr_a = 0;
	this->a_tr_b = 0;
	this->a_tst_all = 0;
	this->a_tst_a = 0;
	this->a_tst_b = 0;	

	this->sens_train = 0;
	this->sens_test = 0;
	this->spec_train = 0;
	this->spec_test = 0;
	this->ppv_train = 0;
	this->ppv_test = 0;
	this->npv_train = 0;
	this->npv_test = 0;
}


void Evaluation::sum_with(Evaluation* R){	

	this->n_train_a += R->n_train_a;
	this->n_train_b += R->n_train_b;
	this->n_test_a += R->n_test_a;
	this->n_test_b += R->n_test_b;
	this->n_train += R->n_train;
	this->n_test += R->n_test;
	this->c_train += R->c_train;
	this->c_test += R->c_test;
	this->c_train_a += R->c_train_a;
	this->c_train_b += R->c_train_b;
	this->c_test_a += R->c_test_a;
	this->c_test_b += R->c_test_b;
	this->w_test += R->w_test;
	this->w_test_a += R->w_test_a;
	this->w_test_b += R->w_test_b;
	this->w_train += R->w_train;
	this->w_train_a  += R->w_train_a;
	this->w_train_b += R->w_train_b;
	this->a_tr_all += R->a_tr_all;
	this->a_tr_a += R->a_tr_a;
	this->a_tr_b += R->a_tr_b;
	this->a_tst_all += R->a_tst_all;
	this->a_tst_a += R->a_tst_a;
	this->a_tst_b += R->a_tst_b;	

	this->sens_train += R->sens_train;
	this->sens_test += R->sens_test;
	this->spec_train += R->spec_train;
	this->spec_test += R->spec_test;
	this->ppv_train += R->ppv_train;
	this->ppv_test += R->ppv_test;
	this->npv_train += R->npv_train;
	this->npv_test += R->npv_test;
}

void Evaluation::divide_by(int n_iter){
	int n_tr = n_iter-1;
	if( n_tr>0 ){
		this->n_train_a /= n_tr;
		this->n_train_b /= n_tr;
		this->n_train /= n_tr;
		this->c_train /= n_tr;
		this->c_train_a /= n_tr;
		this->c_train_b /= n_tr;
		this->w_train /= n_tr;
		this->w_train_a  /= n_tr;
		this->w_train_b /= n_tr;
	}
	if( n_iter>0 ){
		this->n_test_a /= n_iter;
		this->n_test_b /= n_iter;
		this->n_test /= n_iter;
		this->c_test /= n_iter;
		this->c_test_a /= n_iter;
		this->c_test_b /= n_iter;
		this->w_test /= n_iter;
		this->w_test_a /= n_iter;
		this->w_test_b /= n_iter;
		this->a_tr_all /= n_tr;
		this->a_tr_a /= n_iter;
		this->a_tr_b /= n_iter;
		this->a_tst_all /= n_iter;
		this->a_tst_a = this->c_test_a / this->n_test_a * 100;
		this->a_tst_b = this->c_test_b / this->n_test_b * 100;
		// sens, spec, ppt are tested every time, not every-1 time.
		this->sens_test = this->c_test_b / (this->c_test_b + this->w_test_a) * 100;
		this->spec_test = this->c_test_a / (this->c_test_a + this->w_test_b) * 100;
		this->ppv_test = this->c_test_b / (this->c_test_b + this->w_test_b) * 100;
		this->npv_test = this->c_test_a / (this->c_test_a + this->w_test_a) * 100;
		this->sens_train /= n_iter;
		this->spec_train /= n_iter;
		this->ppv_train /= n_iter;
		this->npv_train /= n_iter;	
	}
}
