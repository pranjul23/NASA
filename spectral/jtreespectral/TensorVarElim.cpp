# pragma once

#include <vector>
#include "tensor.hpp"
#include "CPT.hpp"
#include "BayesianNetwork.hpp"
#include "TensorCPT.hpp"
#include "TensorVarElim.hpp"
#include "VectorPlus.hpp"
#include <assert.h>
#include <map>
#include <boost/shared_ptr.hpp>

using namespace std;

// currently assumes all evidence variables have no descendants for simplicity
double TensorVarElim::VE(BayesianNetwork& bNet, vector<int>& evidence_vars, vector<int>& evidence_vals, vector<int>& elim_order)
{

	map< int, vector< boost::shared_ptr<TensorCPT> >* > vars_to_CPTs; // variables and the tensor factors they are associated with
	double joint_val = 0;
	// initial fill
	for (int i = 0; i < bNet.NumNodes(); ++i)
		vars_to_CPTs[i] = new vector< boost::shared_ptr<TensorCPT> >();

	for (int i = 0; i < bNet.NumNodes(); ++i)
	{
		CPT& node_cpt = bNet.GetCPT(i);
		vector<int>& vars = node_cpt.Vars();
		vector<int> mode_counts;
		for (int j = 0; j < vars.size(); ++j)
		{
			if (vars[j] == i)
				mode_counts.push_back(1);
			else
				mode_counts.push_back(2);
		}

		boost::shared_ptr<TensorCPT> t_node_cpt(new TensorCPT(bNet.GetCPT(i), vars, mode_counts));
		int evidence_index = VectorPlus::Find(evidence_vars, i);
		if (evidence_index >= 0)
		{
			boost::shared_ptr<TensorCPT> new_t_cpt(new TensorCPT());
			vector<int>& i_modes = t_node_cpt->Modes(evidence_vars[evidence_index]);
			t_node_cpt->Slice(*new_t_cpt, i_modes[0], evidence_vals[evidence_index]);
			t_node_cpt = new_t_cpt;
		}
		
		int first_var = VectorPlus::FindNextElement(vars, elim_order, 0);
		vars_to_CPTs[first_var]->push_back(t_node_cpt);
	}

	// run elimination
	for (int i = 0; i < elim_order.size(); ++i)
	{
		int e_var = elim_order[i];
		vector< boost::shared_ptr<TensorCPT> >& e_CPTs = *(vars_to_CPTs[e_var]);
		//if (e_CPTs.size() == 0)
		//	continue;
		assert(e_CPTs.size() >= 1);

		boost::shared_ptr<TensorCPT> curr_tensor(new TensorCPT(*e_CPTs[0]));

		
		for (int j = 1; j < e_CPTs.size(); ++j)
		{
			TensorCPT* new_tensor = new TensorCPT();
			TensorCPT::Multiply(*new_tensor, *curr_tensor, *(e_CPTs[j]));
			curr_tensor.reset(new_tensor);
		}

		if (!VectorPlus::Contains(evidence_vars, e_var))
		{
			TensorCPT* new_tensor = new TensorCPT();
			vector<int> e_var_vec = VectorPlus::CreateSingleton(e_var);
			TensorCPT::OneReduce(*new_tensor, *curr_tensor, e_var_vec);
			curr_tensor.reset(new_tensor);
		}

		if (i < elim_order.size() - 1)
		{
			vector<int>& new_vars = curr_tensor->Vars();
			int next_var = VectorPlus::FindNextElement(new_vars, elim_order, i + 1);
			vars_to_CPTs[next_var]->push_back(curr_tensor);
		}
		else
		{
			assert(curr_tensor->Order() == 0);
			joint_val =  curr_tensor->At(0);
			break;
		}
	}

	map< int, vector< boost::shared_ptr<TensorCPT> >* >::iterator iter;
	for (iter = vars_to_CPTs.begin(); iter != vars_to_CPTs.end(); ++iter)
	{
		delete(iter->second);
	}
	return joint_val;
}
