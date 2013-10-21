# pragma once

#include "tensor.hpp"

using namespace std;

class TensorFactor
{
	private:
		vector<int> front_vars;
		vector<int> back_vars;

		int num_front_vars; 
		int num_back_vars; 

		Tensor factor;



	void Set(vector<int> index_array, double val);
	static TensorFactor FactorMultiply(TensorFactor factor1, TensorFactor factor2);

}
