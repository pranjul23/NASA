#include <iostream>
#include <vector>
#include <utility>
#include <fstream>
#include <iomanip>

#include <dai/alldai.h>
#include <dai/trainHMM.h>
#include <dai/testHMM.h>

using namespace dai;
using namespace std;


int main(){
	//===========================================================
	//train model using normal sequences

	TrainHMM model("HMMtraining.txt");
	model.train();


	//===========================================================
	//now evaluate test sequences for anomaly

	TestHMM evaluation("HMMtesting.txt");
	evaluation.test();
};












