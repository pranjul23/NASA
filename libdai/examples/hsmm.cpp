#include <iostream>
#include <vector>
#include <utility>
#include <fstream>
#include <iomanip>
#include <string>
#include <stdlib.h>

#include <dai/alldai.h>
#include <dai/trainHSMM.h>
#include <dai/testHSMM.h>

using namespace dai;
using namespace std;


int main( int argc, char *argv[] ){

	size_t ID = 0, max_num_iter = 0;

	if ( argc != 3) {
		cout << "Usage: " << argv[0] << " <ID> <max_num_iter>" << endl << endl;
		return 1;
	}
	else {
		ID = atoi(argv[1]);
		max_num_iter = atoi(argv[2]);
	}

	stringstream trainData, initFactor, testData, learntFactor;

	trainData << string("data/HSMMtraining_") << ID << string(".txt");
	testData << string("data/HSMMtesting_") << ID << string(".txt");
	initFactor << string("data/hsmm_factor_graph_init_") << ID << string(".fg");
	learntFactor << string("data/hsmm_factor_graph_learnt_") << ID << string(".fg");

	//===========================================================
	//train model using normal sequences

	TrainHSMM model(trainData.str().c_str());

    model.train(initFactor.str().c_str(), ID, max_num_iter);

	//===========================================================
	//now evaluate test sequences for anomaly

	TestHSMM evaluation(testData.str().c_str());

	evaluation.test_loglik(learntFactor.str().c_str(), ID);
	evaluation.test_marginal(learntFactor.str().c_str(), ID);

};












