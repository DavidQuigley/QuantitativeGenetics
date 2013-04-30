class Evaluation{
	/// Essentially a struct used to classifier performance data
public:
	Evaluation();

	void sum_with(Evaluation* recipient); // used by loocv
	void divide_by(int n_iter); // used by loocv
	double n_train_a, n_train_b, n_test_a, n_test_b;
	double n_train, n_test;
	double c_train, c_test, c_train_a, c_train_b, c_test_a, c_test_b;
	double w_test, w_test_a, w_test_b, w_train, w_train_a, w_train_b;
	double a_tr_all, a_tr_a, a_tr_b, a_tst_all, a_tst_a, a_tst_b;	

	double sens_train, sens_test;
	double spec_train, spec_test;
	double ppv_train, ppv_test;
	double npv_train, npv_test;
};
