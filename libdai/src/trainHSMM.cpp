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

void TrainHSMM::train(const char* filename, size_t ID, size_t max_num_iter){

	// Read FactorGraph from the file specified by the first command line argument
	HSMMparam param(filename);

	stringstream hsmmParam;
	hsmmParam << string("data/hsmm_param_init_") << ID << string(".txt");

	param.printHSMMparam(hsmmParam.str().c_str());

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

//	std::vector<size_t> v;
//	for(size_t i=0; i<param.init.size(); i++){
//		param.init[i].getAllNonZeros(v);
//		for(size_t j=0; j<v.size(); j++){
//			std::cout << "init[i] ind = " << v[j] << "\n";
//		}
//		std::cout << "========================\n";
//	}
//
//	for(size_t i=0; i<param.dist.size(); i++){
//		param.dist[i].getAllNonZeros(v);
//		for(size_t j=0; j<v.size(); j++){
//			std::cout << "dist[i] ind = " << v[j] << "\n";
//		}
//		std::cout << "========================\n";
//	}


	for(size_t iter = 0; iter < max_num_iter; iter++){

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
			//transit += jt->calcDistrib(vs);
			transit += jt->calcDistrib(vs, param.dist[1]);

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
			//durat += jt->calcDistrib(vs);
			durat += jt->calcDistrib(vs, param.dist[0]);



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
			//observ += jt->calcDistrib(vs);
			observ += jt->calcDistrib(vs, param.dist[2]);


			delete jt;
			delete graph;

			cout << "Processed data point " << i << " out of " << data.size() << "\n";

//			cout << "likelihood_curr "<<likelihood_curr << "\n";
		}

		//normalize the result and put back into init and dist structure
		initD0 /= data.size();
		initA0 /= data.size();

		param.init[0] = initD0;
		param.init[1] = initA0;

		//{dt-1, dt, at}, marginalize dt out to get {dt-1, at}, and devide {dt-1, dt, at}/{dt-1, at}
		jt->normDistrib(durat, 1);
		param.dist[0] = durat;

		//{at-1, dt, at}, marginalize at out to get {at-1, dt}, and devide {at-1, dt, at}/{at-1, dt}
		jt->normDistrib(transit, 2);
		param.dist[1] = transit;

		//{at, ot}, marginalize ot out to get {at}, and devide {at, ot}/{at}
		jt->normDistrib(observ, 1);
		param.dist[2] = observ;

//		for(size_t i=0; i<param.init.size(); i++){
//			std::cout << "init[i] = " << param.init[i] << "\n";
//		}
//
//		for(size_t i=0; i<param.dist.size(); i++){
//			std::cout << "dist[i] = " << param.dist[i] << "\n";
//		}

		//clear factors
		initA0.fill(0);
		initD0.fill(0);
		transit.fill(0);
		durat.fill(0);
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

	stringstream hsmmParam_learnt;
	hsmmParam_learnt << string("data/hsmm_param_learnt_") << ID << string(".txt");

	param.printHSMMparam(hsmmParam_learnt.str().c_str());


	stringstream learntFactorGraph;
	learntFactorGraph << string("data/hsmm_factor_graph_learnt_") << ID << string(".fg");

	param.saveHSMMparam(learntFactorGraph.str().c_str());
}



//the EM training algorithm applied to updated version of HSMM

void TrainHSMM::train(const char* filename, size_t ID, size_t max_num_iter, int dummy, int Ntrain, int Ltrain){

	// Read FactorGraph from the file specified by the first command line argument
	HSMMparam param(filename, 0);

	stringstream hsmmParam;
	hsmmParam << string("data/hsmm_param_init_") << ID << string(".txt");

	param.printHSMMparam(hsmmParam.str().c_str(), 0);

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
	Factor initA0(0), transit(0), initdur(0), durat(0), observ(0);

	Real likelihood_prev = -1.0;
	Real likelihood_curr = 0;
	vector<Real> likelihood;

	cout << "Model training...\n";

//	std::vector<size_t> v;
//	for(size_t i=0; i<param.init.size(); i++){
//		param.init[i].getAllNonZeros(v);
//		for(size_t j=0; j<v.size(); j++){
//			std::cout << "init[i] ind = " << v[j] << "\n";
//		}
//		std::cout << "========================\n";
//	}
//
//	for(size_t i=0; i<param.dist.size(); i++){
//		param.dist[i].getAllNonZeros(v);
//		for(size_t j=0; j<v.size(); j++){
//			std::cout << "dist[i] ind = " << v[j] << "\n";
//		}
//		std::cout << "========================\n";
//	}

	if(Ntrain == -1 or Ntrain > data.size()){
		Ntrain = data.size();
	}

	size_t len = 0;

	for(size_t iter = 0; iter < max_num_iter; iter++){

		for(size_t i=0; i<Ntrain; i++) {


			if(Ltrain == -1 or Ltrain > data[i].size()){
				len = data[i].size();
			}
			else{
				len = Ltrain;
			}

			//initialize HSMM of size equal the number of observations
			graph = new FactorGraph();

			graph->createHSMMFactorGraph(param.init, param.dist, len, 0);

			jt = new JTree(*graph, opts("updates",string("HUGIN"))("heuristic",string("MINWEIGHT")) );

			//clamp the observation variables to their observed values
			for(size_t j = 0; j < len; j++ ){

				//cout << "clamping var" << data[i][j].first << " to value " << data[i][j].second << "\n";

				//decrement index of last measurement
				//this is needed because in new HSMM model at the last step there is no d_t
				//but if we simply truncate training data, then there will be d_t
				if(j==len-1 and len < data[i].size()){
					jt->clamp(data[i][j].first-1, data[i][j].second);
				}
				else{
					jt->clamp(data[i][j].first, data[i][j].second);
				}
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
			initA0 += jt->calcMarginal(graph->var(1));


			//================= prepare variables to calculate initial duration distribution =========
			vs.clear();

			X |= graph->var(0);
			X |= graph->var(1);

			vs.push_back(X);

			initdur += jt->calcDistrib(vs, param.dist[0]);

			//================= prepare variables to calculate transition distribution ===============
			vs.clear();

			for(size_t j = 1; j < len-1; j++){
				X |= graph->var(3*j-3);
				X |= graph->var(3*j-2);
				X |= graph->var(3*j+1);
				vs.push_back(X);

				tmp = X;
				X /= tmp;
			}

			X |= graph->var(3*(len-1) - 3);
			X |= graph->var(3*(len-1) - 2);
			X |= graph->var(3*(len-1));
			vs.push_back(X);

			tmp = X;
			X /= tmp;

			transit += jt->calcDistrib(vs, param.dist[2]);

			//================= prepare variables to calculate duration distribution ===============
			vs.clear();


			for(size_t j = 1; j < len-1; j++){
				X |= graph->var(3*j-3);
				X |= graph->var(3*j);
				X |= graph->var(3*j+1);
				vs.push_back(X);

				tmp = X;
				X /= tmp;
			}

			durat += jt->calcDistrib(vs, param.dist[1]);



			//================= prepare variables to calculate observation distribution ===============
			vs.clear();

			X |= graph->var(1);
			X |= graph->var(2);

			vs.push_back(X);
			tmp = X;
			X /= tmp;

			for(size_t j = 1; j < len-1; j++){
				X |= graph->var(3*j+1);
				X |= graph->var(3*j+2);
				vs.push_back(X);

				tmp = X;
				X /= tmp;
			}

			X |= graph->var(3*(len-1));
			X |= graph->var(3*(len-1) + 1);
			vs.push_back(X);

			tmp = X;
			X /= tmp;

			//calculate state duration distribution
			//res = jt->calcDistrib(vs, 1);
			//observ += jt->calcDistrib(vs);
			observ += jt->calcDistrib(vs, param.dist[3]);


			delete jt;
			delete graph;

			//cout << "Processed data point " << i << " out of " << data.size() << "\n";

//			cout << "likelihood_curr "<<likelihood_curr << "\n";
		}

		//normalize the result and put back into init and dist structure
		initA0 /= Ntrain;
		param.init[0] = initA0;

		//{dt, at}, marginalize out dt to get {at}, and devide {dt, at}/{at}
		jt->normDistrib(initdur, 0);
		param.dist[0] = initdur;

		//{dt-1, dt, at}, marginalize dt out to get {dt-1, at}, and devide {dt-1, dt, at}/{dt-1, at}
		jt->normDistrib(durat, 1);
		param.dist[1] = durat;

		//{at-1, dt, at}, marginalize at out to get {at-1, dt}, and devide {at-1, dt, at}/{at-1, dt}
		jt->normDistrib(transit, 2);
		param.dist[2] = transit;

		//{at, ot}, marginalize ot out to get {at}, and devide {at, ot}/{at}
		jt->normDistrib(observ, 1);
		param.dist[3] = observ;

//		for(size_t i=0; i<param.init.size(); i++){
//			std::cout << "init[i] = " << param.init[i] << "\n";
//		}
//
//		for(size_t i=0; i<param.dist.size(); i++){
//			std::cout << "dist[i] = " << param.dist[i] << "\n";
//		}

		//clear factors
		initA0.fill(0);
		initdur.fill(0);
		durat.fill(0);
		transit.fill(0);
		observ.fill(0);

		likelihood.push_back(likelihood_curr);

		cout << "Iteration # " << iter << ". LogLikelihood: " << likelihood_curr <<
				", diff: "<<  std::abs(likelihood_curr-likelihood_prev) <<"\n";

		//if( std::abs(likelihood_curr-likelihood_prev) < 2.0 ){
		//	break;
		//}

		likelihood_prev = likelihood_curr;
		likelihood_curr = 0;

	}

	cout << "Training done.\n";

	stringstream hsmmParam_learnt;
	hsmmParam_learnt << string("data/hsmm_param_learnt_") << ID << string(".txt");

	param.printHSMMparam(hsmmParam_learnt.str().c_str(), 0);

	stringstream learntFactorGraph;
	learntFactorGraph << string("data/hsmm_factor_graph_learnt_") << ID << string(".fg");

	param.saveHSMMparam(learntFactorGraph.str().c_str(), 0);
}



}









