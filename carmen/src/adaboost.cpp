#include <vector>
#include <math.h>
#include <time.h>
#include <boost/random.hpp>
#include "DataStructures.h"
#include "adaboost.h"

weak_learner::weak_learner(){
	this->affected_tr = new Matrix<double>();
	this->affected_tst = new Matrix<double>();
}

weak_learner::~weak_learner(){
	delete this->affected_tr;
	delete this->affected_tst;
}

adaboost::adaboost(Matrix<double>* features, std::vector<double>* targets, std::vector<int>* idx_train, std::vector<int>* idx_test){
	// F: {feature x target} (boolean matrix of which features are present in the targets)
    // L: {target x 1}       (integer matrix of sample labels)
    // HO_: {target x 1}     (boolean matrix of which samples are held out)
    // W: {target x 1}       (weight matrix for targets)

    this->F = features;
	
	Matrix<double>* L_horiz = new Matrix<double>(targets);
	this->L = new Matrix<double>();
	L_horiz->transpose( this->L );
	delete L_horiz;
	this->N = L->rows();

	this->HO_tr = new Matrix<double>(N,1);
	this->HO_tst = new Matrix<double>(N,1);
	this->HO_trUp = new Matrix<double>(N,1);
	this->HO_trDn = new Matrix<double>(N,1);
	this->pred_tr = new Matrix<double>(N,1);
	this->pred_tst = new Matrix<double>(N,1);
	this->W = new Matrix<double>(N,1);
	int row=0;
	this->n_tr = (int)idx_train->size();
	this->n_tst = (int)idx_test->size();
	for(int i=0; i<this->n_tr; i++){
		row = idx_train->at(i);
		this->HO_tr->arr[row][0] = 1;
		if( this->L->arr[row][0] == 1 )
			this->HO_trUp->arr[row][0] = 1;
		else
			this->HO_trDn->arr[row][0] = 1;
	}
	for(int i=0; i<this->n_tst; i++){
		row = idx_test->at(i);
		this->HO_tst->arr[row][0] = 1;
	}
	this->predictions = new std::vector<double>();
	this->weak_learners = new std::vector<weak_learner*>();
	this->verbose = false;
	this->test_mode = false;
    this->stumps_only_mode = true;
}

inline int adaboost::rnd(int aRange){
	return int((double )rand()/(double(RAND_MAX) + 1)*aRange);
}

void adaboost::initialize(){
	// PY> W = 1.0/this->n_tr * self.HO.tr;    
	// PY> us = (W*self.HO.trUp).sum()
    // PY> ds = (W*self.HO.trDn).sum() + (10**-20)
	// PY> scroot = 0.5*log( us / ds )
    // PY> if scroot==0:
    // PY>    scroot = 0.005
	// PY> pred_tr = scroot * self.HO.tr
    // PY> pred_tst = scroot * self.HO.tst
	this->HO_tr->clone( this->W );
	this->W->mult( 1.0 / this->n_tr );
	double us = this->W->mult_sum(this->HO_trUp);
	double ds = this->W->mult_sum(this->HO_trDn) + pow((double)10, (double)-20);
	double scroot = 0.5*log( us / ds );
	if(scroot==0)
		scroot = 0.005;
	this->HO_tr->clone( this->pred_tr );
	this->pred_tr->mult(scroot);
	this->pred_tr->mult( this->HO_tr );
	this->HO_tst->clone( this->pred_tst );
    this->pred_tst->mult(scroot);
	this->pred_tst->mult( this->HO_tst );
	if( ! this->test_mode ){
		/*
		srand((unsigned int)time(0));
		for(int r=0; r<this->pred_tst->rows(); r++){
			if( rnd(2) == 1 ){
				this->pred_tst->arr[r][0] *= -1;
				this->pred_tr->arr[r][0] *= -1;
			}
		}
		*/
		// deliberately mis-label all scores to start
		for(int r=0; r<this->pred_tst->rows(); r++){
			if( (this->L->arr[r][0] * this->pred_tst->arr[r][0]) > 0) {this->pred_tst->arr[r][0]*= -1;}
            if( (this->L->arr[r][0] * this->pred_tr->arr[r][0]) > 0) {this->pred_tr->arr[r][0] *= -1;}
		}
	}
}


void adaboost::report_r_by_s(Matrix<double>* r_by_s){
	// update the rules by samples matrix to reflect the results of boosting.
	if( r_by_s == NULL )
		throw std::string("ERROR: r_by_s is NULL");

	// contain r_by_s to the rules we used
	int n_rules = (int)this->weak_learners->size();
	r_by_s->resize( n_rules, r_by_s->cols() );
	weak_learner* wl;
	for(int r=0; r<n_rules; r++){
		for(int c = 0; c<r_by_s->cols(); c++){
			wl = weak_learners->at(r);
			r_by_s->arr[r][c]= ( wl->affected_tr->arr[c][0] + wl->affected_tst->arr[c][0] );
		}
	}
}

void adaboost::update_weights(){
	
	// PY> W = ( e**(-1*L*pred_tr) ) * hoTr
	// PY> W = W/W.sum()
	this->pred_tr->clone(this->W);
	this->W->mult(this->L);
	this->W->mult(-1);	
	this->W->exponent();
	this->W->mult(this->HO_tr);
	double w_sum = this->W->sum();
	this->W->div( w_sum );
}


void adaboost::boost(int n_rounds){
	Matrix<double>* L_LTZ = new Matrix<double>(this->N, 1);
	Matrix<double>* L_GTZ = new Matrix<double>(this->N, 1);
	Matrix<double>* affected_all_t = new Matrix<double>();
	Matrix<double>* affected_all = new Matrix<double>();
	Matrix<double>* w_tr_pos = new Matrix<double>();
	Matrix<double>* w_tr_neg = new Matrix<double>();
	Matrix<double>* score_tr = new Matrix<double>();
	Matrix<double>* score_tst = new Matrix<double>();
	Matrix<double>* pred_tr_LTZ = new Matrix<double>();
	Matrix<double>* pred_tr_GTZ = new Matrix<double>();
	Matrix<double>* pred_tst_LTZ = new Matrix<double>();
	Matrix<double>* pred_tst_GTZ = new Matrix<double>();
	Matrix<double>* margin_tr = new Matrix<double>();
	Matrix<double>* margin_tst = new Matrix<double>();
	Matrix<double>* wrong_tr = new Matrix<double>();
	Matrix<double>* wrong_tst = new Matrix<double>();
	this->initialize();
	//this->update_weights();
	this->L->where_LT(L_LTZ, 0);
	this->L->where_GT(L_GTZ, 0);
	double eps = pow((double)10, (double)-20);
	if( this->verbose ){
		std::cout << "Boosting " << n_rounds << " rounds using " << this->F->rows() << " features.\n"; 
		std::cout.flush();
	}

	for(int r=0; r<n_rounds; r++){
		weak_learner* wl = new weak_learner();

		this->choose_weak_learner(wl, this->HO_tr, this->HO_tst); 		
        wl->id = r;
		wl->parent_idx = adaboost::ROOT;
        
        if( !this->stumps_only_mode ){
            for(int w=0; w<(int)this->weak_learners->size(); w++){
                weak_learner* wl_child = new weak_learner();
                this->choose_weak_learner(wl_child, this->weak_learners->at(w)->affected_tr, this->weak_learners->at(w)->affected_tst);
                if( wl_child->loss < wl->loss ){
                    wl->parent_idx = w;
                    wl_child->affected_tr->clone( wl->affected_tr );
                    wl_child->affected_tst->clone( wl->affected_tst );
                    wl->loss = wl_child->loss;
                    wl->idx = wl_child->idx;
                }
                delete wl_child;
            }
        }
		// find score for new feature
		// PY> w_tr_pos = W * wl.affected_tr * self.HO.trUp
        // PY> w_tr_neg = W * wl.affected_tr * self.HO.trDn 
		W->clone(w_tr_pos);
		W->clone(w_tr_neg);
		w_tr_pos->mult(wl->affected_tr);
		w_tr_neg->mult(wl->affected_tr);
		w_tr_pos->mult(this->HO_trUp);
		w_tr_neg->mult(this->HO_trDn);
		// PY> sum_pos = w_tr_pos.sum()
        // PY> sum_neg = w_tr_neg.sum()
		// PY> wl.score = 0.5*log( (sum_pos+eps) / ( sum_neg+eps))
		double sum_pos = w_tr_pos->sum();
		double sum_neg = w_tr_neg->sum();
		wl->score = 0.5 * log( (sum_pos+eps) / ( sum_neg+eps));
		
		// update prediction

		// PY> pred_tr = pred_tr + array( wl.affected_tr*wl.score )
		wl->affected_tr->clone(score_tr);
		score_tr->mult(wl->score);
		this->pred_tr->add( score_tr );
        
		// PY> pred_tst = pred_tst + affected_tst * wl.score
		wl->affected_tst->clone(score_tst);
        score_tst->mult(wl->score);
		this->pred_tst->add(score_tst);
		pred_tr->where_LT(pred_tr_LTZ, 0);
		pred_tr->where_GT(pred_tr_GTZ, 0);
		pred_tst->where_LT(pred_tst_LTZ, 0);
		pred_tst->where_GT(pred_tst_GTZ, 0);
        
		this->L->clone(margin_tr);  // want margin to be as high as possible.
		this->L->clone(margin_tst);
		margin_tr->mult(pred_tr);
		margin_tst->mult(pred_tst);
		margin_tr->where_LT(wrong_tr, 0);
		margin_tst->where_LT(wrong_tst, 0);
		
		double n_wrong_tr = wrong_tr->sum();
		double n_wrong_tst = wrong_tst->sum();
        wl->loss_tr = (double)n_wrong_tr / ((double)this->n_tr + eps);
        wl->loss_tst = (double)n_wrong_tst / ((double)this->n_tst + eps);

		wl->CN_tr = (int)L_LTZ->mult_sum(pred_tr_LTZ);
        wl->CP_tr = (int)L_GTZ->mult_sum(pred_tr_GTZ);
        wl->WN_tr = (int)L_LTZ->mult_sum(pred_tr_GTZ);
        wl->WP_tr = (int)L_GTZ->mult_sum(pred_tr_LTZ);
        wl->CN_tst = (int)L_LTZ->mult_sum(pred_tst_LTZ);
        wl->CP_tst = (int)L_GTZ->mult_sum(pred_tst_GTZ);
        wl->WN_tst = (int)L_LTZ->mult_sum(pred_tst_GTZ);
        wl->WP_tst = (int)L_GTZ->mult_sum(pred_tst_LTZ);
        
		if(this->verbose){
			if(this->feature_names.size()>0){
				std::cout << "MESSAGE: Round " << r+1 << " " << this->feature_names.at(wl->idx) << "\t";
				std::cout << wl->CN_tr << " " << wl->WN_tr << " " << wl->CP_tr << " " << wl->WP_tr << " | ";
				std::cout << wl->CN_tst << " " << wl->WN_tst << " " << wl->CP_tst << " " << wl->WP_tst << "\n";
				std::cout.flush();
			}
		}
		this->weak_learners->push_back(wl);
		this->update_weights();
	}
	// write predictions.  Key step!
	Matrix<double> * p = new Matrix<double>();
	this->pred_tr->clone(p);
	p->add(this->pred_tst);
	p->rowSums(this->predictions);

	delete p;
	delete L_LTZ;
	delete L_GTZ; 
	delete affected_all_t;
	delete affected_all;
	delete w_tr_pos;
	delete w_tr_neg; 
	delete score_tr; 
	delete score_tst;
	delete pred_tr_LTZ;
	delete pred_tr_GTZ;
	delete pred_tst_LTZ;
	delete pred_tst_GTZ; 
	delete margin_tr ;
	delete margin_tst;
	delete wrong_tr;
	delete wrong_tst;
}


void adaboost::boost_robust(int n_rounds){
	// The difference between boost and boost_robust is behavior when more than one 
	// feature can produce the lowest new training loss.  In boost, a single feature
	// is chosen at random.  This means that multiple classifiers built on the same 
	// features can be equally powerful but made of different features.
	//
	// In the multiple-feature situation, boost_robust adds all of the features 
	// with weight normalized by the number of features chosen.  This may result in 
	// superior interpretability and generalization performance.  However, it is no longer 
	// known exactly how many features will be produced by a given number of rounds.
	Matrix<double>* L_LTZ = new Matrix<double>(this->N, 1);
	Matrix<double>* L_GTZ = new Matrix<double>(this->N, 1);
	Matrix<double>* affected_all_t = new Matrix<double>();
	Matrix<double>* affected_all = new Matrix<double>();
	Matrix<double>* w_tr_pos = new Matrix<double>();
	Matrix<double>* w_tr_neg = new Matrix<double>();
	Matrix<double>* score_tr = new Matrix<double>();
	Matrix<double>* score_tst = new Matrix<double>();
	Matrix<double>* pred_tr_LTZ = new Matrix<double>();
	Matrix<double>* pred_tr_GTZ = new Matrix<double>();
	Matrix<double>* pred_tst_LTZ = new Matrix<double>();
	Matrix<double>* pred_tst_GTZ = new Matrix<double>();
	Matrix<double>* margin_tr = new Matrix<double>();
	Matrix<double>* margin_tst = new Matrix<double>();
	Matrix<double>* wrong_tr = new Matrix<double>();
	Matrix<double>* wrong_tst = new Matrix<double>();
	this->initialize();
	this->update_weights();
	this->L->where_LT(L_LTZ, 0);
	this->L->where_GT(L_GTZ, 0);
	double eps = pow((double)10, (double)-20);
	weak_learner* wl;
	for(int r=0; r<n_rounds; r++){
		std::vector<weak_learner*> wls;
		this->choose_weak_learners(wls, this->HO_tr, this->HO_tst); 
		for(int i=0; i<(int)wls.size(); i++){
			// find score
			wl = wls.at(i);
			W->clone(w_tr_pos);
			W->clone(w_tr_neg);
			w_tr_pos->mult(wl->affected_tr);
			w_tr_neg->mult(wl->affected_tr);
			w_tr_pos->mult(this->HO_trUp);
			w_tr_neg->mult(this->HO_trDn);
			wl->score = 0.5 * log( (w_tr_pos->sum()+eps) / ( w_tr_neg->sum()+eps));
			// KEY STEP: normalize score by number of features chosen at once
			wl->score = wl->score / (double)wls.size();

			wl->affected_tr->clone(score_tr);
			score_tr->mult(wl->score);
			this->pred_tr->add( score_tr );
        
			wl->affected_tst->clone(score_tst);
			score_tst->mult(wl->score);
			this->pred_tst->add(score_tst);
			pred_tr->where_LT(pred_tr_LTZ, 0);
			pred_tr->where_GT(pred_tr_GTZ, 0);
			pred_tst->where_LT(pred_tst_LTZ, 0);
			pred_tst->where_GT(pred_tst_GTZ, 0);
        
			this->L->clone(margin_tr);  // want margin to be as high as possible.
			this->L->clone(margin_tst);
			margin_tr->mult(pred_tr);
			margin_tst->mult(pred_tst);
			margin_tr->where_LT(wrong_tr, 0);
			margin_tst->where_LT(wrong_tst, 0);
		
			double n_wrong_tr = wrong_tr->sum();
			double n_wrong_tst = wrong_tst->sum();
			wl->loss_tr = (double)n_wrong_tr / ((double)this->n_tr + eps);
			wl->loss_tst = (double)n_wrong_tst / ((double)this->n_tst + eps);

			wl->CN_tr = (int)L_LTZ->mult_sum(pred_tr_LTZ);
			wl->CP_tr = (int)L_GTZ->mult_sum(pred_tr_GTZ);
			wl->WN_tr = (int)L_LTZ->mult_sum(pred_tr_GTZ);
			wl->WP_tr = (int)L_GTZ->mult_sum(pred_tr_LTZ);
			wl->CN_tst = (int)L_LTZ->mult_sum(pred_tst_LTZ);
			wl->CP_tst = (int)L_GTZ->mult_sum(pred_tst_GTZ);
			wl->WN_tst = (int)L_LTZ->mult_sum(pred_tst_GTZ);
			wl->WP_tst = (int)L_GTZ->mult_sum(pred_tst_LTZ);
			this->weak_learners->push_back(wl);
		}
		this->update_weights();

	}
	// write predictions.  Key step!
	Matrix<double> * p = new Matrix<double>();
	this->pred_tr->clone(p);
	p->add(this->pred_tst);
	p->rowSums(this->predictions);

	delete p;
	delete L_LTZ;
	delete L_GTZ; 
	delete affected_all_t;
	delete affected_all;
	delete w_tr_pos;
	delete w_tr_neg; 
	delete score_tr; 
	delete score_tst;
	delete pred_tr_LTZ;
	delete pred_tr_GTZ;
	delete pred_tst_LTZ;
	delete pred_tst_GTZ; 
	delete margin_tr ;
	delete margin_tst;
	delete wrong_tr;
	delete wrong_tst;
}


void adaboost::choose_weak_learner(weak_learner* wl, Matrix<double>* valid_tr, Matrix<double>* valid_tst){

	Matrix<double>* wpos = new Matrix<double>(); 
	Matrix<double>* wneg = new Matrix<double>(); 
	Matrix<double>* wpos_f = new Matrix<double>(); 
	Matrix<double>* wneg_f = new Matrix<double>(); 
	
	// PY> wpos = W * valid_tr * self.HO.trUp
	// PY> wneg = W * valid_tr * self.HO.trDn
	this->W->clone(wpos);
	this->W->clone(wneg);
	wpos->mult(valid_tr);
	wneg->mult(valid_tr);
	wpos->mult( this->HO_trUp );
	wneg->mult( this->HO_trDn );
    
	// PY> wpos = asarray(self.F * wpos) # {features, outcome}
    // PY> wneg = asarray(self.F * wneg)
	this->F->mm(wpos, wpos_f);
	this->F->mm(wneg, wneg_f);
    int r = wpos_f->rows();
	int c = 1;

	// PY> other = ones((r,c)) - wpos - wneg;
    Matrix<double>* other = new Matrix<double>(r, c, false, 1); 
    other->sub(wpos_f);
	other->sub(wneg_f);
	
	// PY> sqrt_term = (wpos*wneg)
    // PY> sqrt_term = sqrt_term**(0.5)
	Matrix<double>* sqrt_term = new Matrix<double>();
	wpos_f->clone(sqrt_term);
	sqrt_term->mult(wneg_f);
	sqrt_term->power(0.5);
    
    // PY> loss = (2*(sqrt_term)) + other;
	Matrix<double>* loss = new Matrix<double>();
    sqrt_term->clone(loss);
	loss->mult((double)2);
	loss->add(other);
	std::vector<int>* loss_rows = new std::vector<int>();
	std::vector<int>* loss_cols = new std::vector<int>();

	wl->loss = loss->min(loss_rows, loss_cols);
	if( loss_rows->size()>1 ){
		std::random_shuffle( loss_rows->begin(), loss_rows->end() );
		if(this->verbose){
			//std::cout << "MESSAGE: Choosing among";
			//if( this->feature_names.size()>0 ){
			//	for(int x=0;x<(int)loss_rows->size(); x++)
			//		std::cout << " " << this->feature_names.at(loss_rows->at(x));
			//	std::cout << "\n";
			//}
			//else{
			//	std::cout << " " << loss_rows->size() << " weak learners.\n"; 
			//}
			//std::cout.flush();
		}
	}
	wl->idx = loss_rows->at(0);

	valid_tr->clone(wl->affected_tr);
	valid_tst->clone(wl->affected_tst);
	// New_affected = currently_affected AND affected_by_feature
	for( int c=0; c<this->F->cols(); c++){
		if( this->F->arr[wl->idx][c]==0 ){
            wl->affected_tr->arr[c][0]=0;
            wl->affected_tst->arr[c][0]=0;
		}
	}
	delete loss_rows;
	delete loss_cols;
	delete wpos; 
	delete wneg; 
	delete wpos_f; 
	delete wneg_f;
	delete other;
	delete sqrt_term;
	delete loss;
}


void adaboost::choose_weak_learners(std::vector<weak_learner*>& wls, Matrix<double>* valid_tr, Matrix<double>* valid_tst){
	Matrix<double>* wpos = new Matrix<double>(); 
	Matrix<double>* wneg = new Matrix<double>(); 
	Matrix<double>* wpos_f = new Matrix<double>(); 
	Matrix<double>* wneg_f = new Matrix<double>(); 

	this->W->clone(wpos);
	this->W->clone(wneg);
	wpos->mult(valid_tr);
	wneg->mult(valid_tr);
	wpos->mult( this->HO_trUp );
	wneg->mult( this->HO_trDn );
	this->F->mm(wpos, wpos_f);
	this->F->mm(wneg, wneg_f);
    int r = wpos_f->rows();
	int c = 1;
    Matrix<double>* other = new Matrix<double>(r, c, false, 1); 
    other->sub(wpos_f);
	other->sub(wneg_f);
	Matrix<double>* sqrt_term = new Matrix<double>();
	wpos_f->clone(sqrt_term);
	sqrt_term->mult(wneg_f);
	sqrt_term->power(0.5);
	Matrix<double>* loss = new Matrix<double>();
    sqrt_term->clone(loss);
	loss->mult((double)2);
	loss->add(other);
	std::vector<int>* loss_rows = new std::vector<int>();
	std::vector<int>* loss_cols = new std::vector<int>();
	weak_learner* wl;
	double learner_loss = loss->min(loss_rows, loss_cols);
	wls.clear();
	for(int i=0; i<(int)loss_rows->size();i++){
		wl = new weak_learner();
		wl->loss = learner_loss;
		wl->idx = loss_rows->at(i);
		valid_tr->clone(wl->affected_tr);
		valid_tst->clone(wl->affected_tst);
		for( int c=0; c<this->F->cols(); c++){
			if( this->F->arr[wl->idx][c]==0 ){
				wl->affected_tr->arr[c][0]=0;
				wl->affected_tst->arr[c][0]=0;
			}
		}
		wls.push_back(wl);
	}
	delete loss_rows;
	delete loss_cols;
	delete wpos; 
	delete wneg; 
	delete wpos_f; 
	delete wneg_f;
	delete other;
	delete sqrt_term;
	delete loss;
}


void adaboost::print(){
	std::cout << "Round\tTrain\t(" << this->n_tr << ")\t\t\t\tTest\t(" << this->n_tst << ")\n";
	std::cout << "\t%Loss\tCP\tWP\tCN\tWN\t%Loss\tCP\tWP\tCN\tWN\tRule_IDX\n";

	char L[80];
	weak_learner* wl;
	for(int i=0; i<(int)this->weak_learners->size(); i++){
		wl = this->weak_learners->at(i);
		sprintf(L, "%d)\t%2.3f\t%d\t%d\t%d\t%d\t%2.3f\t%d\t%d\t%d\t%d\t%d\n", i+1, wl->loss_tr, wl->CP_tr, wl->WP_tr, wl->CN_tr, wl->WN_tr, wl->loss_tst, wl->CP_tst, wl->WP_tst, wl->CN_tst, wl->WN_tst, wl->idx);
		std::cout << L;
	}
}

adaboost::~adaboost(){
	delete this->HO_tr;
	delete this->HO_tst;
	delete this->HO_trUp;
	delete this->HO_trDn;
	delete this->pred_tr;
	delete this->pred_tst;
	delete this->L;
	delete this->W; 
	delete this->predictions;
	for(int i=0; i< (int)this->weak_learners->size(); i++)
		delete this->weak_learners->at(i);
	delete this->weak_learners;
}



