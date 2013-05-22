#include <iostream>
#include <vector>
#include <utility>
#include <fstream>
#include <iomanip>

#include <dai/alldai.h>
#include <dai/evidence.h>

using namespace dai;
using namespace std;

class Observation{
private:
	vector<vector<pair<size_t, size_t> > > data;

public:
	Observation(const char* filename){

		ifstream file(filename);
		size_t var, value, numObs, dataLen;

		file >> dataLen;
		data.reserve(dataLen);

		for(size_t i=0; i<dataLen; i++){
			file >> numObs;
			data.push_back(vector<pair<size_t, size_t> >(numObs, make_pair(0, 0)));

			for(size_t j=0; j<numObs; j++){
				file >> var;
				data.back().at(j).first = var;
			}

			for(size_t j=0; j<numObs; j++){
				file >> value;
				data.back().at(j).second = value;
			}

		}
	}

	vector<vector<pair<size_t, size_t> > >& getData(){
		return data;
	}
};


class HMMParam{
private:
	vector<Var> trans, obs;
	Permute permTrans, permObs;

public:

	vector<Factor> init, dist;

	HMMParam(const char* filename){
		ifstream is;
		string line;
		size_t nr_Factors;

		is.open( filename );
		if( !is.is_open() ) {
			DAI_THROWE(CANNOT_READ_FILE,"Cannot read from file " + std::string(filename));
		}


		while( (is.peek()) == '#' )
			getline(is,line);
		is >> nr_Factors;

		DAI_ASSERT(nr_Factors == 3);

		getline (is,line);

		if( is.fail() || line.size() > 0 )
			DAI_THROWE(INVALID_FACTORGRAPH_FILE,"Expecting empty line");

		map<long,size_t> vardims;

		for( size_t I = 0; I < nr_Factors; I++ ) {

			size_t nr_members;
			while( (is.peek()) == '#' ) getline(is,line);
			is >> nr_members;

			vector<long> labels;
			for( size_t mi = 0; mi < nr_members; mi++ ) {
				long mi_label;
				while( (is.peek()) == '#' ) getline(is,line);
				is >> mi_label;
				labels.push_back(mi_label);
			}

			vector<size_t> dims;
			for( size_t mi = 0; mi < nr_members; mi++ ) {
				size_t mi_dim;
				while( (is.peek()) == '#' ) getline(is,line);
				is >> mi_dim;
				dims.push_back(mi_dim);
			}

			// add the Factor
			vector<Var> Ivars;
			Ivars.reserve( nr_members );

			for( size_t mi = 0; mi < nr_members; mi++ ) {
				map<long,size_t>::iterator vdi = vardims.find( labels[mi] );
				if( vdi != vardims.end() ) {
					// check whether dimensions are consistent
					if( vdi->second != dims[mi] )
						DAI_THROWE(INVALID_FACTORGRAPH_FILE,"Variable with label " + boost::lexical_cast<string>(labels[mi]) + " has inconsistent dimensions.");
				} else
					vardims[labels[mi]] = dims[mi];
				Ivars.push_back( Var(labels[mi], dims[mi]) );
			}

			if(I < 1){
				init.push_back( Factor( VarSet( Ivars.begin(), Ivars.end(), Ivars.size() ), (Real)0 ) );
			}
			else{
				dist.push_back( Factor( VarSet( Ivars.begin(), Ivars.end(), Ivars.size() ), (Real)0 ) );
			}

			// calculate permutation object
			Permute permindex( Ivars );

			// read values
			size_t nr_nonzeros;
			while( (is.peek()) == '#' ) getline(is,line);
			is >> nr_nonzeros;

			for( size_t k = 0; k < nr_nonzeros; k++ ) {
				size_t li;
				Real val;
				while( (is.peek()) == '#' ) getline(is,line);
				is >> li;

				while( (is.peek()) == '#' ) getline(is,line);
				is >> val;

				// store value, but permute indices first according to internal representation
				if(I < 1){
					init.back().set( permindex.convertLinearIndex( li ), val );
				}
				else{
					dist.back().set( permindex.convertLinearIndex( li ), val );
				}
			}
		}

		is.close();



		//variables for displaying computed distributions
		//internaly the order of variables is in increasing order
		//but here we want in such order so as to print (at|at-1)
		//print row by row, so the fastest varying index corresponds to at-1
		trans.push_back(dist[0].vars().elements().at(1));
		trans.push_back(dist[0].vars().elements().at(0));

		obs.push_back(dist[1].vars().elements().at(0));
		obs.push_back(dist[1].vars().elements().at(1));

		Permute permTransition(trans), permObservation(obs);
		permTrans = permTransition;
		permObs = permObservation;
	}

	void saveHMMparam(const char* filename){
		ofstream os;

		os.open( filename, ios::trunc );

		//create fake Factor graph and write it to file
		FactorGraph graph = FactorGraph();
		graph.createHMMFactorGraph(init, dist, 1);

		os << graph;

		os.close();
	}

	void printHMMparam(const char* filename){

		ofstream myfile;
		myfile.open (filename, ios::trunc);

		myfile << "==============================\n";
		myfile << "P(A0): \n" << init[0] << "\n";
		myfile << "==============================\n";


		size_t li = 0;
		myfile << "P(At|At-1 = " << "): \n";
		for(size_t jj=0; jj<trans[1].states(); jj++){
			for(size_t kk=0; kk<trans[0].states(); kk++){
				myfile << setw(6) << dist[0].p().get(permTrans.convertLinearIndex(li)) << " ";
				li++;
			}
			myfile << "\n";
		}
		myfile << "================================\n";


		li = 0;
		myfile << "P(Ot|At = " << "): \n";
		for(size_t ii=0; ii<obs[1].states(); ii++){
			for(size_t jj=0; jj<obs[0].states(); jj++){
				myfile << setw(8) << dist[1].p().get(permObs.convertLinearIndex(li)) << " ";
				li++;
			}
			myfile << "\n";
		}
		myfile << "================================\n";

		myfile.close();
	}

};

int main(){

	// Read FactorGraph from the file specified by the first command line argument
	HMMParam param("hmm_factor_graph_init.fg");

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
	Observation obs( "/Users/igor/Documents/Projects/anomaly/code/HSMM/libdai/examples/HMMtraining.txt" );
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

	//cout << "train likelihood: "<<"\n";
	//for(size_t i=0; i<likelihood.size(); i++){
	//	cout << likelihood.at(i)<<"\n";
	//}

	param.printHMMparam("hmm_param_learnt.txt");
	param.saveHMMparam("hmm_factor_graph_learnt.fg");

	//===========================================================================================================
	//now evaluate test sequences for anomaly

	//HMMParam param("hmm_factor_graph_learnt.fg");

	// Read test data
	Observation test( "/Users/igor/Documents/Projects/anomaly/code/HSMM/libdai/examples/HMMtesting.txt" );
	vector<vector<pair<size_t, size_t> > > test_data = test.getData();

	vector<Real> likelihood_test;

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

	//cout << "test likelihood: "<<"\n";
	//for(size_t i=0; i<likelihood_test.size(); i++){
	//	cout << likelihood_test.at(i)<<"\n";
	//}

	ofstream os;
	os.open("HMMlikelihood_test.txt", ios::trunc);

	for(size_t i=0; i<likelihood_test.size(); i++){
		os << likelihood_test.at(i)<<"\n";
	}

};












