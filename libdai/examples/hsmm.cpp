#include <iostream>
#include <vector>
#include <utility>
#include <fstream>
#include <iomanip>

#include <dai/alldai.h>
#include <dai/trainHSMM.h>
#include <dai/testHSMM.h>

using namespace dai;
using namespace std;


int main(){

	//===========================================================
	//train model using normal sequences

	TrainHSMM model("HSMMtraining.txt");
    model.train();

	//===========================================================
	//now evaluate test sequences for anomaly

	TestHSMM evaluation("HSMMtesting.txt");
	evaluation.test();
};












