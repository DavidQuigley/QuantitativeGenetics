
class weak_learner{
public:
	weak_learner();
	~weak_learner();
	double loss;
	int id;
	int idx;
	int parent_idx;
	double score;
	double loss_tr, loss_tst;
	Matrix<double>* affected_tr;
	Matrix<double>* affected_tst;
	int CN_tr, CP_tr, WN_tr, WP_tr;
	int CN_tst, CP_tst, WN_tst, WP_tst;
};

class adaboost{
public:
	adaboost(Matrix<double>* features, std::vector<double>* targets, std::vector<int>* idx_train, std::vector<int>* idx_test);
	~adaboost();
	const static int ROOT = -1;
	Matrix<double>* F;
	Matrix<double>* L;
	int N;
	int n_tr;
	int n_tst;
	Matrix<double>* HO_tr;
	Matrix<double>* HO_tst;
	Matrix<double>* HO_trUp;
	Matrix<double>* HO_trDn;
	Matrix<double>* pred_tr;
	Matrix<double>* pred_tst;
	std::vector<double>* predictions;
	Matrix<double>* W;
	std::vector<weak_learner*>* weak_learners;
	std::vector<std::string> feature_names;
	bool verbose;
	bool test_mode;
	bool stumps_only_mode;
	void boost(int rounds);
	void boost_robust(int rounds);
	void report_r_by_s(Matrix<double>* r_by_s);

	void initialize();
	void update_weights();
	void choose_weak_learners(std::vector<weak_learner*>& wls, Matrix<double>* valid_tr, Matrix<double>* valid_tst);
private:
	void choose_weak_learner(weak_learner* wl, Matrix<double>* valid_tr, Matrix<double>* valid_tst);
	void print();
	int rnd(int aRange);
};

