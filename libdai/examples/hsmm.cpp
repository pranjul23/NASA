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

	stringstream trainData, initFactor, testData, learntFactor, trueFactor;

	trainData << string("data/HSMMtraining_") << ID << string(".txt");
	testData << string("data/HSMMtesting_") << ID << string(".txt");

	initFactor << string("data/hsmm_factor_graph_init_") << ID << string(".fg");
	learntFactor << string("data/hsmm_factor_graph_learnt_") << ID << string(".fg");
	trueFactor << string("data/hsmm_factor_graph_true_") << ID << string(".fg");

	//===========================================================
	//train model using normal sequences.

    TrainHSMM model(trainData.str().c_str());

	//NOTE:  in the code, if the variable "int dummy" is present
	//it means we use a code for different HSMM model

    //model.train(initFactor.str().c_str(), ID, max_num_iter);
	model.train(initFactor.str().c_str(), ID, max_num_iter, 0);

	//===========================================================
	//now evaluate test sequences for anomaly

//	TestHSMM evalTrain(trainData.str().c_str());
//	evalTrain.test_loglik(learntFactor.str().c_str(), ID, "train");

	//TestHSMM evaluation(testData.str().c_str());
	TestHSMM evaluation(testData.str().c_str());

	//evaluation.test_loglik(learntFactor.str().c_str(), ID, "test");

	evaluation.test_loglik(learntFactor.str().c_str(), ID, "test", 0);
	evaluation.test_loglik(trueFactor.str().c_str(), ID, "true", 0);

	//evaluation.test_loglik(initFactor.str().c_str(), ID, "init", 0);

	//evaluation.test_marginal(learntFactor.str().c_str(), ID);
	//evaluation.test_marginal_cut(learntFactor.str().c_str(), ID);

};












