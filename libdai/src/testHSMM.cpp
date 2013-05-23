
#include <vector>
#include <iostream>

#include <dai/alldai.h>
#include <dai/testHSMM.h>
#include <dai/observation.h>
#include <dai/HSMMparam.h>



namespace dai{

using namespace std;

void TrainHSMM::train(){

}
//HSMMparam param("hsmm_factor_graph_learnt.fg");

	// Read test data
	Observation test( "/Users/igor/Documents/Projects/anomaly/code/HSMM/libdai/examples/HSMMtesting.txt" );
	vector<vector<pair<size_t, size_t> > > test_data = test.getData();

	vector<Real> likelihood_test;

	cout << "Now we do sequence testing...\n";

	for(size_t i=0; i<test_data.size(); i++) {

		//initialize HSMM of size equal the number of observations
		graph = new FactorGraph();
		graph->createHSMMFactorGraph(param.init, param.dist, test_data[i].size());

		jt = new JTree(*graph, opts("updates",string("HUGIN"))("heuristic",string("MINWEIGHT")) );

		//clamp the observation variables to their observed values
		for(size_t j = 0; j < test_data[i].size(); j++ ){
			//cout << "clamping var" << test_data[i][j].first << " to value " << test_data[i][j].second << "\n";
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

	//cout << "test likelihood: "<<"\n";
	//for(size_t i=0; i<likelihood_test.size(); i++){
	//	cout << likelihood_test.at(i)<<"\n";
	//}

	ofstream os;
	os.open("HSMMlikelihood_test.txt", ios::trunc);

	for(size_t i=0; i<likelihood_test.size(); i++){
		os << likelihood_test.at(i)<<"\n";
	}
