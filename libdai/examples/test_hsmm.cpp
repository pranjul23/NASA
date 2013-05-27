#include <iostream>
#include <map>
#include <dai/alldai.h>  // Include main libDAI header file
#include <dai/jtree.h>
#include <dai/varset.h>

using namespace dai;
using namespace std;

int main(){
	
	// Read FactorGraph from the file specified by the first command line argument
	FactorGraph fg;
	vector<Factor> init;
	vector<Factor> dist;

	size_t T = 2;

	fg.ReadHMMFromFile(init, dist, "hsmm_factor_graph.fg");
	fg.createFG(init, dist, T);

	// Set some constants
	size_t maxstates = 1000000;
	size_t maxiter = 10000;
	Real   tol = 1e-9;
	size_t verb = 3;
	
	// Store the constants in a PropertySet object
	PropertySet opts;
	opts.set("maxiter",maxiter);  // Maximum number of iterations
	opts.set("tol",tol);          // Tolerance for convergence
	opts.set("verbose",verb);     // Verbosity (amount of output generated)
		
	
	// Bound treewidth for junctiontree
	bool do_jt = true;
	try {
		boundTreewidth(fg, &eliminationCost_MinFill, maxstates );
	} catch( Exception &e ) {
		if( e.getCode() == Exception::OUT_OF_MEMORY ) {
			do_jt = false;
			cout << "Skipping junction tree (need more than " << maxstates << " states)." << endl;
		}
		else
			throw;
	}
	
	JTree jt;
	
	if( do_jt ) {
		jt = JTree( fg, opts("updates",string("HUGIN"))("heuristic",string("MINWEIGHT")) );
		jt.init();


		cout << "Factors' variables and their values before clamping:\n";
		for( size_t i = 0; i < jt.nrFactors(); i++ ){
					cout << jt.factor(i).vars() << endl;
					cout << jt.factor(i).p() << endl;
		}
		cout << "\n\n";

//		jt.clamp(4, 2);
//		jt.clamp(7, 2);
//		jt.clamp(10, 0);

		cout << "Factors' variables and their values after clamping:\n";
		for( size_t i = 0; i < jt.nrFactors(); i++ ){
					cout << jt.factor(i).vars() << endl;
					cout << jt.factor(i).p() << endl;
		}
		cout << "\n\n";

		jt.run();

	}
		
	
	/*
	if( do_jt ) {
		// Report variable marginals for fg, calculated by the junction tree algorithm
		cout << "Exact variable marginals:" << endl;
		for( size_t i = 0; i < fg.nrVars(); i++ ) // iterate over all variables in fg
			cout << jt.belief(fg.var(i)) << endl; // display the "belief" of jt for that variable
	}
	 */		

	VarSet X, tmp;
	std::vector<VarSet> vsNumer;

	X |= fg.var(0);
	X |= fg.var(1);
	X |= fg.var(3);
	cout << "X = " << X << "\n";

	vsNumer.push_back(X);
	tmp = X;
	X /= tmp;

	for(size_t i=2; i<=T; i++){
		X |= fg.var(3*i-4);
		X |= fg.var(3*i-3);
		X |= fg.var(3*i);
		vsNumer.push_back(X);

		cout << "X = " << X << "\n";

		tmp = X;
		X /= tmp;
	}


	cout << "Initial A0 distribution: " << jt.calcInitDistrib( VarSet(fg.var(1)) ) << "\n";
	cout << "Initial D0 distribution: " << jt.calcInitDistrib( VarSet(fg.var(0)) ) << "\n";

	cout << "Transition distribution: " << jt.calcDistrib(vsNumer, 2) << "\n";

//	X |= fg.var(3);
//	X |= fg.var(5);
	
//	cout << jt.belief(fg.var(2)) << endl;
	
//	for( size_t I = 0; I < fg.nrFactors(); I++ ) // iterate over all factors in fg
//		cout << jt.calcMarginal(fg.factor(I).vars()) << endl; // display the belief of bp for the variables in that factor

	
/* 
	if( do_jt ){
		for( size_t i=0; i < N; i++ )			
			X |= fg.var(3*i+2);	
		
		cout << "X = " << X << "\n";
		cout << "Joint marginal of output variables: \n";
		cout << jt.calcMarginal(X) << "\n";
	}		 

*/ 
 
}










