#include <vector>
#include <iostream>

#include <dai/alldai.h>
#include <dai/trainHMM.h>
#include <dai/observation.h>
#include <dai/HMMparam.h>


namespace dai{

using namespace std;

TrainHMM::TrainHMM(const char* filename){
	// Read test data
		Observation train( filename );
		data = train.getData();
}

void TrainHMM::train(const char* filename, size_t ID, size_t max_num_iter){

	// Read FactorGraph from the file specified by the first command line argument
	HMMparam param(filename);

	stringstream hmmParam;
	hmmParam << string("data/hmm_param_init_") << ID << string(".txt");

	param.printHMMparam(hmmParam.str().c_str());

	// Set some constants
	size_t maxiter = 10000;
	Real   tol = 1e-9;
	size_t verb = 0;

	// Store the constants in a PropertySet object
	PropertySet opts;
	opts.set("maxiter",maxiter);  // Maximum number of iterations
	opts.set("tol",tol);          // Tolerance for convergence
	opts.set("verbose",verb);     // Verbosity (amount of output generated)

	FactorGraph *graph;
	JTree *jt;

	VarSet X, tmp;
	std::vector<VarSet> vs;
	Factor initA0(0), transit(0),  observ(0);

	Real likelihood_prev = -1.0;
	Real likelihood_curr = 0;
	vector<Real> likelihood;

	cout << "Model training...\n";

	for(size_t iter = 0; iter < max_num_iter; iter++){

		for(size_t i=0; i<data.size(); i++) {

			//initialize HMM of size equal the number of observations
			graph = new FactorGraph();
			graph->createHMMFactorGraph(param.init, param.dist, data[i].size());

			jt = new JTree(*graph, opts("updates",string("HUGIN"))("heuristic",string("MINWEIGHT")) );

			//clamp the observation variables to their observed values
			for(size_t j = 0; j < data[i].size(); j++ ){
				//cout << "clamping var" << data[i][j].first << " to value " << data[i][j].second << "\n";
				jt->clamp(data[i][j].first, data[i][j].second);
			}

			jt->init();
			jt->run();

			likelihood_curr += jt->logZ();

			//do inference

			//================= calculate initial distributions ===============
			initA0 += jt->calcMarginal(graph->var(0));


			//================= prepare variables to calculate transition distribution ===============
			vs.clear();

			X |= graph->var(0);
			X |= graph->var(1);

			vs.push_back(X);
			tmp = X;
			X /= tmp;

			for(size_t j = 2; j <= data[i].size(); j++){
				X |= graph->var(2*j-3);
				X |= graph->var(2*j-1);
				vs.push_back(X);

				tmp = X;
				X /= tmp;
			}

			//calculate state transition distribution
			//res = jt->calcDistrib(vs, 2);
			transit += jt->calcDistrib(vs);
//			cout << "transit = " << transit << "\n";


			//================= prepare variables to calculate observation distribution ===============
			vs.clear();

			for(size_t j = 1; j <= data[i].size(); j++){
				X |= graph->var(2*j);
				X |= graph->var(2*j-1);
				vs.push_back(X);

				tmp = X;
				X /= tmp;
			}

			//calculate state duration distribution
			//res = jt->calcDistrib(vs, 1);
			observ += jt->calcDistrib(vs);

			delete jt;
			delete graph;

			cout << "Processed data point " << i << " out of " << data.size() << "\n";

//			cout << "likelihood_curr "<<likelihood_curr << "\n";
		}

		//normalize the result and put back into init and dist structure
		initA0 /= data.size();

		param.init[0] = initA0;

		//transit contains {at-1, at}, marginalize at out to get {at-1}; then {at|at-1} is obtianed by dividing {at-1, at}/{at-1}
		//marginalize(remove) variable 1==at-1
		jt->normDistrib(transit, 1);
		param.dist[0] = transit;

		//{at, ot}, marginalize ot out to get {at}, and devide {at, ot}/{at}
		jt->normDistrib(observ, 1);
		param.dist[1] = observ;

		//clear factors
		initA0.fill(0);
		transit.fill(0);
		observ.fill(0);

		likelihood.push_back(likelihood_curr);

		cout << "Iteration # " << iter << ". LogLikelihood: " << likelihood_curr <<
				", diff: "<<  std::abs(likelihood_curr-likelihood_prev) <<"\n";

		if( std::abs(likelihood_curr-likelihood_prev) < 0.1 ){
			break;
		}

		likelihood_prev = likelihood_curr;
		likelihood_curr = 0;
	}

	cout << "Training done.\n";

	stringstream hmmParam_learnt;
	hmmParam_learnt << string("data/hmm_param_learnt_") << ID << string(".txt");

	param.printHMMparam(hmmParam_learnt.str().c_str());


	stringstream learntFactorGraph;
	learntFactorGraph << string("data/hmm_factor_graph_learnt_") << ID << string(".fg");

	param.saveHMMparam(learntFactorGraph.str().c_str());
}
}
