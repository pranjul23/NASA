# pragma once

#include <vector>
#include "tensor.hpp"
#include "CPT.hpp"
#include "BayesianNetwork.hpp"

class TensorVarElim
{
	public:
		static double VE(BayesianNetwork& bNet, vector<int>& evidence_vars, vector<int>& evidence_vals, vector<int>& elim_order);
};
