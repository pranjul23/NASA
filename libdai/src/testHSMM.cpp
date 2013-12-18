
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>

#include <dai/alldai.h>
#include <dai/testHSMM.h>
#include <dai/observation.h>
#include <dai/HSMMparam.h>


namespace dai{

using namespace std;

TestHSMM::TestHSMM(const char* filename){
	// Read test data
	Observation test( filename );

	test_data = test.getData();
}


void TestHSMM::test_loglik(const char* filename, size_t ID, string type){

	HSMMparam param(filename);

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


	cout << "Now we do testing...\n";

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

		//compute normalized loglikelyhood
		likelihood_test.push_back(jt->logZ()/test_data[i].size());

		delete jt;
		delete graph;

		cout << "Tested point " << i << " out of " << test_data.size() <<"\n";
	}

	cout << "done.\n";

	ofstream os;
	stringstream result;
	result << string("data/HSMMlikelihood_") << type << string("_") << ID << string(".txt");
	os.open(result.str().c_str(), ios::trunc);

	for(size_t i=0; i<likelihood_test.size(); i++){
		os << likelihood_test.at(i)<<"\n";
	}
}



void TestHSMM::test_loglik(const char* filename, size_t ID, string type, int dummy, int ID2){

	//"type" specifies the suffix of the output file "test" or "true"

	//read in HSMM parameters from textfile
	HSMMparam param(filename, 0);

	vector<Real> likelihood_test;
	FactorGraph *graph;
	JTree *jt;

	VarSet X;

	// Set some constants
	size_t maxiter = 10000;
	Real   tol = 1e-9;
	size_t verb = 0;

	// Store the constants in a PropertySet object
	PropertySet opts;
	opts.set("maxiter",maxiter);  // Maximum number of iterations
	opts.set("tol",tol);          // Tolerance for convergence
	opts.set("verbose",verb);     // Verbosity (amount of output generated)


	cout << "Now we do testing...\n";

	for(size_t i=0; i<test_data.size(); i++) {

		//initialize HSMM of size equal the number of observations
		graph = new FactorGraph();
		graph->createHSMMFactorGraph(param.init, param.dist, test_data[i].size(), 0);

		jt = new JTree(*graph, opts("updates",string("HUGIN"))("heuristic",string("MINWEIGHT")) );


		//clamp the observation variables to their observed values
		for(size_t j = 0; j < test_data[i].size(); j++ ){
			//cout << "clamping var" << test_data[i][j].first << " to value " << test_data[i][j].second << "\n";
			jt->clamp(test_data[i][j].first, test_data[i][j].second);
		}

		jt->init();
		jt->run();

		//compute normalized loglikelihood
		//likelihood_test.push_back(jt->logZ()/test_data[i].size());
		likelihood_test.push_back(jt->logZ());

		delete jt;
		delete graph;

		cout << "Tested point " << i << " out of " << test_data.size() <<"\n";

		//cout << "jt->logZ() = " << likelihood_test.at(0) << "\n";
		//exit(1);
	}

	cout << "done.\n";

	ofstream os;
	stringstream result;

	//ID2 is used to wirte test results with fixed number of training iterations
	if(ID2 >= 0){
		result << string("data/HSMMlikelihood_") << type << string("_") << ID << string("-") << ID2 << string(".txt");
	}
	else{
		result << string("data/HSMMlikelihood_") << type << string("_") << ID << string(".txt");
	}
	os.open(result.str().c_str(), ios::trunc);


	os.unsetf ( std::ios::floatfield );
	os.precision(18);
	for(size_t i=0; i<likelihood_test.size(); i++){
		os << likelihood_test.at(i)<<"\n";
	}
}




void TestHSMM::test_marginal(const char* filename, size_t ID){

	HSMMparam param(filename);

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

	Factor O_last;
	vector< vector<Real> > all_marginal;
	vector<Real> sequence_marginal;
	all_marginal.reserve(test_data.size());


	cout << "Now we do testing...\n";

	for(size_t i=0; i<test_data.size(); i++) {

		//allocate memory
		sequence_marginal.reserve(test_data[i].size());

		//initialize HSMM of ever increasing size up to test_data[i].size()
		graph = new FactorGraph();
		graph->createHSMMFactorGraph(param.init, param.dist, test_data[i].size());

		jt = new JTree(*graph, opts("updates",string("HUGIN"))("heuristic",string("MINWEIGHT")) );

		jt->init();
		jt->run();

		O_last = jt->calcMarginal(graph->var(4));
		sequence_marginal.push_back( log(O_last.p().get(test_data[i][0].second)) );

		for(size_t k=1; k < test_data[i].size(); k++){

			//clamp all the observation variables to their observed values, except last variable which is not clamped
			jt->clamp(test_data[i][k-1].first, test_data[i][k-1].second);

			jt->run();

			//compute p(o_last=c | o_1...o_{last-1})
			//this will give us a distribution: {o_last=1, o_last=2, ... o_last=M}
			O_last = jt->calcMarginal(graph->var(3*k+4));

			//since we have a specific observation at last time step: o_last=c, get its probability:
			sequence_marginal.push_back( log(O_last.p().get(test_data[i][k].second)) );
		}

		cout << "Tested point " << i << " out of " << test_data.size() <<"\n";

		all_marginal.push_back(sequence_marginal);
		sequence_marginal.clear();

		delete jt;
		delete graph;
	}

	cout << "Testing done.\n";

	ofstream os;
	stringstream result;
	result << string("data/HSMMmarginal_test_") << ID << string(".txt");
	os.open(result.str().c_str(), ios::trunc);

	for(size_t i=0; i<all_marginal.size(); i++){
		for(size_t j=0; j<all_marginal[i].size(); j++){
			os << all_marginal[i][j]<<" ";
		}
		os << "\n";
	}
}


void TestHSMM::test_marginal_cut(const char* filename, size_t ID){

	HSMMparam param(filename);

	FactorGraph *graph;
	JTree *jt;

	// Set some constants
	size_t maxiter = 10000;
	Real   tol = 1e-9;
	size_t verb = 0;

	//window size
    size_t W = 20;
    size_t start = 0;

	// Store the constants in a PropertySet object
	PropertySet opts;
	opts.set("maxiter",maxiter);  // Maximum number of iterations
	opts.set("tol",tol);          // Tolerance for convergence
	opts.set("verbose",verb);     // Verbosity (amount of output generated)

	Factor O_last;
	vector< vector<Real> > all_marginal;
	vector<Real> sequence_marginal;
	all_marginal.reserve(test_data.size());


	cout << "Now we do testing...\n";

	for(size_t i=0; i<test_data.size(); i++) {

		//allocate memory
		sequence_marginal.reserve(test_data[i].size());

		//initialize HSMM of ever increasing size up to test_data[i].size()
		graph = new FactorGraph();
		graph->createHSMMFactorGraph(param.init, param.dist, test_data[i].size());

		jt = new JTree(*graph, opts("updates",string("HUGIN"))("heuristic",string("MINWEIGHT")) );

		jt->init();
		jt->run();

		O_last = jt->calcMarginal(graph->var(4));
		sequence_marginal.push_back( log(O_last.p().get(test_data[i][0].second)) );

		delete jt;

		for(size_t k=1; k < test_data[i].size(); k++){

			jt = new JTree(*graph, opts("updates",string("HUGIN"))("heuristic",string("MINWEIGHT")) );

			//clamp a window of observable variables to their values, except last variable which is not clamped
			start = k-W;
			if(start < 0) start = 0;

			for(size_t j = start; j <= k-1; j++){
				jt->clamp(test_data[i][j].first, test_data[i][j].second);
			}

			jt->init();
			jt->run();

			//compute p(o_last=c | o_1...o_{last-1})
			//this will give us a distribution: {o_last=1, o_last=2, ... o_last=M}
			O_last = jt->calcMarginal(graph->var(3*k+4));

			//since we have a specific observation at last time step: o_last=c, get its probability:
			sequence_marginal.push_back( log(O_last.p().get(test_data[i][k].second)) );

			delete jt;
		}

		cout << "Tested point " << i << " out of " << test_data.size() <<"\n";

		all_marginal.push_back(sequence_marginal);
		sequence_marginal.clear();

		delete graph;
	}

	cout << "Testing done.\n";

	ofstream os;
	stringstream result;
	result << string("data/HSMMmarginal_test_") << ID << string(".txt");
	os.open(result.str().c_str(), ios::trunc);

	for(size_t i=0; i<all_marginal.size(); i++){
		for(size_t j=0; j<all_marginal[i].size(); j++){
			os << all_marginal[i][j]<<" ";
		}
		os << "\n";
	}
}


}
