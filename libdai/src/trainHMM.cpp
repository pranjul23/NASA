#include <vector>
#include <iostream>

#include <dai/alldai.h>
#include <dai/trainHMM.h>
#include <dai/observation.h>
#include <dai/HMMparam.h>


namespace dai{

using namespace std;

void TrainHMM::train(){

	// Read FactorGraph from the file specified by the first command line argument
	HMMparam param("hmm_factor_graph_init.fg");

	// Set some constants
	size_t maxiter = 10000;
	Real   tol = 1e-9;
	size_t verb = 0;

	// Store the constants in a PropertySet object
	PropertySet opts;
	opts.set("maxiter",maxiter);  // Maximum number of iterations
	opts.set("tol",tol);          // Tolerance for convergence
	opts.set("verbose",verb);     // Verbosity (amount of output generated)

	// Read sample from file
	Observation obs( "HMMtraining.txt" );
	vector<vector<pair<size_t, size_t> > > data = obs.getData();

	FactorGraph *graph;
	JTree *jt;

	VarSet X, tmp;
	std::vector<VarSet> vs;
	Factor initA0(0), transit(0),  observ(0);

	Real likelihood_prev = -1.0;
	Real likelihood_curr = 0;
	vector<Real> likelihood;

	param.printHMMparam("hmm_param_init.txt");

	for(size_t iter = 0; iter<10; iter++){

		for(size_t i=0; i<data.size(); i++) {

			//initialize HSMM of size equal the number of observations
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

			//cout << "Processed data point " << i << " out of " << data.size() << "\n";
		}

		//normalize the result and put back into init and dist structure
		initA0 /= data.size();

		param.init[0] = initA0;

		//{at-1, at}, marginalize at-1 out: {at}, and devide {at-1, at}/{at}
		jt->normDistrib(transit, 0);
		param.dist[0] = transit;

		//{at, ot}, marginalize ot out: {at}, and devide {at, ot}/{at}
		jt->normDistrib(observ, 1);
		param.dist[1] = observ;

		//clear factors
		initA0.fill(0);
		transit.fill(0);
		observ.fill(0);

		likelihood.push_back(likelihood_curr);
		likelihood_prev = likelihood_curr;
		likelihood_curr = 0;

		cout << "Iteration # " << iter << ". LogLikelihood: " << likelihood_prev <<"\n";
	}

	param.printHMMparam("hmm_param_learnt.txt");
	param.saveHMMparam("hmm_factor_graph_learnt.fg");
}
}
