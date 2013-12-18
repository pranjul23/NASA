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

//this file is to test hsmm EM algorithm and compare it to spectral algorithm in
//simulations

int main( int argc, char *argv[] ){

	size_t ID = 0, max_num_iter = 0, Ntrain;

	if ( argc != 4) {
		cout << "Usage: " << argv[0] << " <ID> <max_num_iter> <Ntrain>" << endl << endl;
		return 1;
	}
	else {
		ID = atoi(argv[1]);
		max_num_iter = atoi(argv[2]);
		Ntrain = atoi(argv[3]);
	}

	stringstream trainData, initFactor, testData, learntFactor, trueFactor;

	trainData << string("data/HSMMtraining_") << ID << string(".txt");
	testData << string("data/HSMMtesting_") << ID << string(".txt");

	initFactor << string("data/hsmm_factor_graph_init_") << ID << string(".fg");
	learntFactor << string("data/hsmm_factor_graph_learnt_") << ID << string(".fg");
	trueFactor << string("data/hsmm_factor_graph_true_") << ID << string(".fg");

	//===========================================================
	//true resutls

	TestHSMM evaluation(testData.str().c_str());
	evaluation.test_loglik(trueFactor.str().c_str(), ID, "true", 0);


	//===========================================================
	//dummy initialization

	//save hsmm_factor_graph_init in hsmm_factor_graph_learnt
	TrainHSMM model(trainData.str().c_str());

	//NOTE:  in the code, if the variable "int dummy" is present
	//it means we use a code for different HSMM model
	//last arg tells how many training data use in training
	model.train(initFactor.str().c_str(), ID, 0, 0, Ntrain);
	evaluation.test_loglik(learntFactor.str().c_str(), ID, "test", 0, 0);


	for(int num_iter = 0, iter = 5, c = 1; num_iter <= max_num_iter;  num_iter += iter, c++){

		//===========================================================
		//train model

		model.train(learntFactor.str().c_str(), ID, iter, 0, Ntrain);

		//===========================================================
		//test model

		evaluation.test_loglik(learntFactor.str().c_str(), ID, "test", 0, c);
	}
};












