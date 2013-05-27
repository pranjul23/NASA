#include <vector>
#include <iostream>
#include <fstream>

#include <dai/alldai.h>
#include <dai/testHMM.h>
#include <dai/observation.h>
#include <dai/HMMparam.h>


namespace dai{

using namespace std;

TestHMM::TestHMM(const char* filename){
	// Read test data
	Observation test( filename );
	test_data = test.getData();
}


void TestHMM::test(){

	HMMparam param("hmm_factor_graph_learnt.fg");

	vector<Real> likelihood_test;
	FactorGraph *graph;
	JTree *jt;

	// Set some constants
	size_t maxiter = 10000;
	Real   tol = 1e-9;
	size_t verb = 0;

	// Store the constants in a PropertySet object
	PropertySet opts;
	opts.set("maxiter",maxiter);  // Maximum number of iterations
	opts.set("tol",tol);          // Tolerance for convergence
	opts.set("verbose",verb);     // Verbosity (amount of output generated)


	cout << "Now we do sequence testing...\n";

	for(size_t i=0; i<test_data.size(); i++) {

		//initialize HSMM of size equal the number of observations
		graph = new FactorGraph();
		graph->createHMMFactorGraph(param.init, param.dist, test_data[i].size());

		jt = new JTree(*graph, opts("updates",string("HUGIN"))("heuristic",string("MINWEIGHT")) );

		//clamp the observation variables to their observed values
		for(size_t j = 0; j < test_data[i].size(); j++ ){
			jt->clamp(test_data[i][j].first, test_data[i][j].second);
		}

		jt->init();
		jt->run();

		likelihood_test.push_back(jt->logZ());

		delete jt;
		delete graph;

		//cout << "Tested point " << i << " out of " << test_data.size() <<"\n";
	}

	cout << "done.\n";

	ofstream os;
	os.open("HMMlikelihood_test.txt", ios::trunc);

	for(size_t i=0; i<likelihood_test.size(); i++){
		os << likelihood_test.at(i)<<"\n";
	}
}
}
