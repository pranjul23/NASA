#include <dai/HMMparam.h>

#include <iostream>
#include <vector>
#include <utility>
#include <fstream>
#include <iomanip>
#include <limits>


namespace dai{

using namespace std;

HMMparam::HMMparam(const char* filename){
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


		// calculate permutation object
		Permute permindex( Ivars );

		// read # of nonzeros that are currently present in the factor
		size_t nr_nonzeros_curr;
		while( (is.peek()) == '#' ) getline(is,line);
		is >> nr_nonzeros_curr;


		// read # of nonzeros in the factor that are possible
		size_t nr_nonzeros_possib;
		while( (is.peek()) == '#' ) getline(is,line);
		is >> nr_nonzeros_possib;

		// read the type of factor: 0-dense, 1-sparse
		size_t factor_type;
		while( (is.peek()) == '#' ) getline(is,line);
		is >> factor_type;

		if(I < 1){
			init.push_back( Factor( VarSet( Ivars.begin(), Ivars.end(), Ivars.size() ), (Real)0, nr_nonzeros_possib, factor_type ) );
		}
		else{
			dist.push_back( Factor( VarSet( Ivars.begin(), Ivars.end(), Ivars.size() ), (Real)0, nr_nonzeros_possib, factor_type ) );
		}


		for( size_t k = 0; k < nr_nonzeros_curr; k++ ) {
			size_t li;
			Real val;
			//while( (is.peek()) == '#' ) getline(is,line);
			is >> li;

			//while( (is.peek()) == '#' ) getline(is,line);
			is >> val;

			// store value, but permute indices first according to internal representation
			if(I < 1){
				init.back().set( permindex.convertLinearIndex( li ), val );
			}
			else{
				dist.back().set( permindex.convertLinearIndex( li ), val );
			}
		}


		for( size_t k = 0; k < nr_nonzeros_possib; k++ ) {
			size_t ind;
			while( (is.peek()) == '#' ) getline(is,line);
			is >> ind;

			// store value, but permute indices first according to internal representation
			if(I < 1){
				init.back().addNonZeroInd(permindex.convertLinearIndex( ind ));
			}
			else{
				dist.back().addNonZeroInd(permindex.convertLinearIndex( ind ));
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

void HMMparam::saveHMMparam(const char* filename){
	ofstream os;

	os.open( filename, ios::trunc );

	//create fake Factor graph and write it to file
	FactorGraph graph = FactorGraph();
	graph.createHMMFactorGraph(init, dist, 1);

	//os << graph;
	graph.saveFactorGraph(os);

	os.close();
}

void HMMparam::printHMMparam(const char* filename){

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

}
