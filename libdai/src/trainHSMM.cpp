#include <vector>
#include <iostream>

#include <dai/alldai.h>
#include <dai/trainHSMM.h>
#include <dai/observation.h>
#include <dai/HSMMparam.h>

//#include "valgrind/callgrind.h"


namespace dai{

using namespace std;


TrainHSMM::TrainHSMM(const char* filename){
	// Read sample from file
	Observation obs(filename);
	data = obs.getData();

}

void TrainHSMM::train(const char* filename){

	// Read FactorGraph from the file specified by the first command line argument
	HSMMparam param(filename);

	param.printHSMMparam("hsmm_param_init.txt");

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
	Factor initA0(0), initD0(0), transit(0), durat(0), observ(0);

	Real likelihood_prev = -1.0;
	Real likelihood_curr = 0;
	vector<Real> likelihood;

	cout << "Model training...\n";


	for(size_t iter = 0; iter<50; iter++){

		for(size_t i=0; i<data.size(); i++) {

			//initialize HSMM of size equal the number of observations
			graph = new FactorGraph();

			graph->createHSMMFactorGraph(param.init, param.dist, data[i].size());

			jt = new JTree(*graph, opts("updates",string("HUGIN"))("heuristic",string("MINWEIGHT")) );

			//clamp the observation variables to their observed values
			for(size_t j = 0; j < data[i].size(); j++ ){
				//cout << "clamping var" << data[i][j].first << " to value " << data[i][j].second << "\n";
				jt->clamp(data[i][j].first, data[i][j].second);
			}

			jt->init();

			//CALLGRIND_START_INSTRUMENTATION;
			//CALLGRIND_TOGGLE_COLLECT;

			jt->run();

			//CALLGRIND_TOGGLE_COLLECT;
			//CALLGRIND_STOP_INSTRUMENTATION;

			likelihood_curr += jt->logZ();

			//do inference

			//================= calculate initial distributions ===============
			initD0 += jt->calcMarginal(graph->var(0));
			initA0 += jt->calcMarginal(graph->var(1));

			//================= prepare variables to calculate transition distribution ===============
			vs.clear();

			X |= graph->var(0);
			X |= graph->var(1);
			X |= graph->var(3);

			vs.push_back(X);

			tmp = X;
			X /= tmp;

			for(size_t j = 2; j <= data[i].size(); j++){
				X |= graph->var(3*j-4);
				X |= graph->var(3*j-3);
				X |= graph->var(3*j);
				vs.push_back(X);

				tmp = X;
				X /= tmp;
			}


			//calculate state transition distribution
			//res = jt->calcDistrib(vs, 2);
			transit += jt->calcDistrib(vs);

			//================= prepare variables to calculate duration distribution ===============
			vs.clear();

			X |= graph->var(0);
			X |= graph->var(2);
			X |= graph->var(3);

			vs.push_back(X);
			tmp = X;
			X /= tmp;


			for(size_t j = 2; j <= data[i].size(); j++){
				X |= graph->var(3*j-4);
				X |= graph->var(3*j-1);
				X |= graph->var(3*j);
				vs.push_back(X);

				tmp = X;
				X /= tmp;
			}


			//calculate state duration distribution
			//res = jt->calcDistrib(vs, 1);
			durat += jt->calcDistrib(vs);



			//================= prepare variables to calculate observation distribution ===============
			vs.clear();

			X |= graph->var(3);
			X |= graph->var(4);

			vs.push_back(X);
			tmp = X;
			X /= tmp;

			for(size_t j = 2; j <= data[i].size(); j++){
				X |= graph->var(3*j);
				X |= graph->var(3*j+1);
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
//			exit(1);

		}

		//normalize the result and put back into init and dist structure
		initD0 /= data.size();
		initA0 /= data.size();

		param.init[0] = initD0;
		param.init[1] = initA0;

		//{dt-1, dt, at}, marginalize dt out: {dt-1, at}, and devide {dt-1, dt, at}/{dt-1, at}
		jt->normDistrib(durat, 1);
		param.dist[0] = durat;

		//{at-1, dt, at}, marginalize at out: {at-1, dt}, and devide {at-1, dt, at}/{at-1, dt}
		jt->normDistrib(transit, 2);
		param.dist[1] = transit;

		//{at, ot}, marginalize ot out: {at}, and devide {at, ot}/{at}
		jt->normDistrib(observ, 1);
		param.dist[2] = observ;

		//clear factors
		initA0.fill(0);
		initD0.fill(0);
		transit.fill(0);
		durat.fill(0);
		observ.fill(0);

		cout << "Iteration # " << iter << ". LogLikelihood: " << likelihood_curr <<
				", diff: "<<  likelihood_curr-likelihood_prev <<"\n";

		likelihood.push_back(likelihood_curr);
		likelihood_prev = likelihood_curr;
		likelihood_curr = 0;
	}

	cout << "done.\n";

	param.printHSMMparam("hsmm_param_learnt.txt");

	param.saveHSMMparam("hsmm_factor_graph_learnt.fg");
}

}
