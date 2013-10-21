# pragma once

#include <vector>
#include "Matrix.hpp"
#include "TensorCPT.hpp"
#include <map>
#include "BayesianNetwork.hpp"
#include "JTNode.hpp"
#include "MultiVector.hpp"
#include "VectorPlus.hpp"
#include "JTree.hpp"
#include <assert.h>
#include <queue>
#include "TensorVarElim.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <assert.h>
#include <cmath>

using namespace std;


/*JTree::JTree(JTree& copy)
{
	bNet = new BayesianNetwork(*(copy.bNet));
	this->training_samples = copy.training_samples;
	this->all_obs_vars = new vector<int>(*copy.all_obs_vars);


		Matrix* training_samples;
		vector<int>* all_obs_vars; 

		Matrix* tree_matrix;
		vector<JTNode*>* node_list;
		int root;

		vector<int>* elimination_order;

		vector<TensorCPT*>* template_cpts;

		double thresh;

		bool clique_params_flag;


}*/

// default constructor
JTree::JTree(BayesianNetwork* bNet, Matrix& tree_matrix, vector<vector<int>*>& nodes_to_vars,
	         Matrix* training_samples, double thresh, bool clique_params_flag)
{
	this->clique_params_flag = clique_params_flag;
	this->thresh = thresh;
	elimination_order = NULL;
	this->bNet = bNet;
	this->tree_matrix = new Matrix(tree_matrix);
	this->all_obs_vars = new vector<int>();
	for (int i = 0; i < bNet->NumNodes(); ++i)
	{
		if (bNet->is_observed(i))
		{
			this->all_obs_vars->push_back(i);
		}
	}


	this->node_list = new vector<JTNode*>();
	for (int i = 0; i < tree_matrix.Dim(0); ++i)
	{
		node_list->push_back(new JTNode(this, i, clique_params_flag));
	}

	// store parents and children

	for (int i = 0; i < NumCliques(); ++i)
	{
		vector<int> parents;
		vector<int> children;

		this->tree_matrix->RowFind(parents, i);
		this->tree_matrix->ColFind(children, i);

		if (parents.size() == 0)
			this->root = i;
		else
			node_list->at(i)->SetParent(parents[0]);

		node_list->at(i)->SetChildren(children);

		node_list->at(i)->SetCVars(*(nodes_to_vars[i]));

	}

	// store C, S, and R vars
	for (int i = 0; i < NumCliques(); ++i)
	{
		vector<int> i_Sep;
		vector<int> i_Res;
		vector<int>& i_vars = node_list->at(i)->GetCVars();

		int parent = GetParent(i);

		if (parent >= 0)
		{
			vector<int>& parent_vars = node_list->at(parent)->GetCVars();
			VectorPlus::Intersect(i_Sep, i_vars, parent_vars);
			VectorPlus::SetDiff(i_Res, i_vars, i_Sep);
		}
		else
		{
			VectorPlus::Copy(i_Res, i_vars);
		}

		node_list->at(i)->SetRVars(i_Res);
		node_list->at(i)->SetSVars(i_Sep);
	}

	ComputeEliminationOrder();

	// assign factors to cliques

	for (int i = 0; i < bNet->NumNodes(); ++i)
	{
		CPT& i_cpt = bNet->GetCPT(i);
		vector<int>& cpt_vars = i_cpt.Vars();

		for (int j = 0; j < NumCliques(); ++j)
		{
			vector<int>& j_vars = GetCVars(j);
			bool is_subset = VectorPlus::IsSubset(j_vars, cpt_vars);
			if (is_subset)
			{
				node_list->at(j)->AddFactor(i_cpt);
				break;
			}
			else if (j == NumCliques() - 1)
			assert("template index not found");
		}
	}
	
	for (int i = 0; i < NumCliques(); ++i)
	{
		node_list->at(i)->CreateCliqueFactor();
	}
	

	this->training_samples = training_samples;
}


bool JTree::IsLeaf(int clique_index)
{
	vector<int>& children = GetChildren(clique_index);
	if (children.size() == 0)
		return false;
	else
		return true;
}

void JTree::ComputeEliminationOrder()
{
	assert(elimination_order == NULL);
	elimination_order = new vector<int>();
	vector<int> temp_vector;
	queue<int> desc_queue;
	desc_queue.push(root);

	while (!desc_queue.empty())
	{
		int node_index = desc_queue.front();
		desc_queue.pop();
		temp_vector.push_back(node_index);
		vector<int>& children = GetChildren(node_index);
		for (int i = 0; i < children.size(); ++i)
			desc_queue.push(children[i]);
	}

	VectorPlus::Reverse(*elimination_order, temp_vector);
}


double JTree::LearnOnlineEMParameters(double learning_rate, int oracle_flag, int max_iter)
{
	int iter = 0;
	int weight_iter = 0;
	double prev_like = 0;
	vector<int> orig_seq;
	VectorPlus::Seq(orig_seq, 0, 1, training_samples->Dim(0));

	vector<TensorCPT*> total_marginals(NumCliques(), NULL);
	vector<int> temp1;
	vector<int> temp2;
	double temp_prob = ComputeEmpiricalMarginals(total_marginals, temp1, temp2);
	cout << "temp_prob: ";
	cout << temp_prob;
	cout << "\n";
	for (iter = 0; iter >= 0; ++iter)
	{
		cout << iter;
		cout << "\n";

		vector<int> permuted_examples;
		VectorPlus::Permute(permuted_examples, orig_seq, iter);

		double likelihood =  0;
		for (int n = 0; n < permuted_examples.size(); ++n)
		{
			++weight_iter;

			vector<int> evidence_vals(all_obs_vars->size(), 0);

			for (int i = 0; i < evidence_vals.size(); ++i)
			{
				evidence_vals.at(i) = (int)(training_samples->At(permuted_examples[n],all_obs_vars->at(i)));
			}

			vector<TensorCPT*> marginals;
		//	double oracle_val = TensorVarElim::VE(*bNet, *all_obs_vars, evidence_vals, bNet->ReverseTopologicalOrder());
			double joint_val = ComputeEmpiricalMarginals(marginals, *all_obs_vars, evidence_vals);
		
			for (int i = 0; i < marginals.size(); ++i)
			{
				marginals[i]->Divide(joint_val);
			}

			double weight = pow(weight_iter + 2.0, -1.0 * learning_rate);
			for (int i = 0; i < total_marginals.size(); ++i)
			{
				if (total_marginals[i] == NULL)
				{
					assert(0);
					total_marginals[i] = new TensorCPT(*marginals[i]);
					total_marginals[i]->Add(0.0000000001);
				}
				else
				{
					TensorCPT::WeightedAdd(*total_marginals[i], *marginals[i], 1 - weight, weight);
				}
				delete(marginals[i]);
			}
			likelihood += log(joint_val);

			for (int i = 0; i < total_marginals.size(); ++i)
			{
				total_marginals[i]->Divide(total_marginals[i]->Sum());
				node_list->at(i)->UpdateFactorWrapper(*(total_marginals[i]));
			}
		}
	//	if (z == max_iter - 1)
	//	{
			cout << "loglike: ";
			cout << likelihood;
			cout << "\n\n";
	//	}

	


	/*	if (iter > 0)
		{
			// test for convergence
			bool converged = TestConvergence(prev_like, likelihood, thresh);
			if (converged)
				break;
		} */

		if (iter > 0)
		{
			// test for convergence
			bool converged = TestConvergence(prev_like, likelihood, thresh);
			if (converged)
			{
				for (int i = 0; i < total_marginals.size(); ++i)
				{
					delete(total_marginals.at(i));
				}

				return likelihood;
			}
		}

		prev_like = likelihood;
	}

	for (int i = 0; i < total_marginals.size(); ++i)
	{
		delete(total_marginals.at(i));
	}

	return prev_like;
}


double JTree::LearnEMParameters(int max_iter, int oracle_flag)
{
	// initialize all CPTs randomly
/*	for (int i = 0; i < NumCliques(); ++i)
	{
		node_list->at(i)->InitializeRandomCPT();
	} */


	if (max_iter == 0)
		return -2;

	int iter = 0;
	double prev_like = 0;
	while (true)
	{
		cout << iter;
		cout << "\n";

		vector<TensorCPT*> total_marginals(NumCliques(), NULL);

		double likelihood =  0;

		for (int n = 0; n < training_samples->Dim(0); ++n)
		{

			vector<int> evidence_vals(all_obs_vars->size(), 0);

			for (int i = 0; i < evidence_vals.size(); ++i)
			{
				evidence_vals.at(i) = (int)(training_samples->At(n,all_obs_vars->at(i)));
			}

			vector<TensorCPT*> marginals;
	//		double oracle_val = TensorVarElim::VE(*bNet, *all_obs_vars, evidence_vals, bNet->ReverseTopologicalOrder());
			double joint_val = ComputeEmpiricalMarginals(marginals, *all_obs_vars, evidence_vals);

			for (int i = 0; i < marginals.size(); ++i)
			{
				marginals[i]->Divide(joint_val);
			}

			for (int i = 0; i < total_marginals.size(); ++i)
			{
				if (total_marginals[i] == NULL)
					total_marginals[i] = new TensorCPT(*marginals[i]);
				else
					TensorCPT::Add(*total_marginals[i], *marginals[i]);

				delete(marginals[i]);
			}
			likelihood += log(joint_val);

		}
	//	if (z == max_iter - 1)
	//	{
			cout << "loglike: ";
			cout << likelihood;
			cout << "\n\n";
	//	}

		for (int i = 0; i < total_marginals.size(); ++i)
		{
			total_marginals[i]->Add(0.0000000001);
			total_marginals[i]->Divide(total_marginals[i]->Sum());
			node_list->at(i)->UpdateFactorWrapper(*(total_marginals[i]));
			delete(total_marginals.at(i));
		}

		if (iter > 0)
		{
			// test for convergence
			bool converged = TestConvergence(prev_like, likelihood, thresh);
			if (converged)
			{
				return likelihood;
			}
		}
		++iter;
		prev_like = likelihood;
	}
	return -1;
}

double JTree::ComputeEmpiricalJointProb(vector<TensorCPT*>& marginals, vector<int>& evidence_vars, vector<int>& evidence_vals)
{
	cout << "yee";
	cout << "\n";
	double joint_val = 0;
	vector<vector<TensorCPT*>*> upward_message_map;
	vector<TensorCPT*> downward_message_map;
	vector<TensorCPT*> integrated_message_map;

	marginals.assign(NumCliques(), NULL);
	for (int i = 0; i < NumCliques(); ++i)
	{
		upward_message_map.push_back(new vector<TensorCPT*>());
		downward_message_map.push_back(NULL);
		integrated_message_map.push_back(NULL);
	}

	for (int i = 0; i < NumCliques(); ++i)
	{
		int clique_index = elimination_order->at(i);

		vector<int> evidence_i_vars;
		vector<int> evidence_i_vals;

		VectorPlus::Intersect(evidence_i_vars, GetCVars(clique_index), evidence_vars);
		VectorPlus::MatchSub(evidence_i_vals, evidence_vals, evidence_vars, evidence_i_vars);

		TensorCPT* outgoing_message;

		if (node_list->at(clique_index)->GetNumFactors() == 0)
		{
			vector<TensorCPT*>& message_vec = *(upward_message_map[clique_index]);
			for (int j = 0; j < message_vec.size(); ++j)
			{
				if (j == 0)
				{
					outgoing_message = new TensorCPT(*message_vec[j]);
				}
				else
				{
					TensorCPT* new_outgoing = new TensorCPT();
					TensorCPT::ElementwiseMultiply(*new_outgoing, *outgoing_message, *message_vec[j]);
					delete(outgoing_message);
					outgoing_message = new_outgoing;
				}
			}
		}
		else
		{
			outgoing_message = new TensorCPT();
			TensorCPT& old_cpt = GetAgglomCliqueFactor(clique_index);

			old_cpt.Select(*outgoing_message, evidence_i_vars, evidence_i_vals); 

			vector<TensorCPT*>& message_vec = *(upward_message_map[clique_index]);
			for (int j = 0; j < message_vec.size(); ++j)
			{
				TensorCPT* new_outgoing = new TensorCPT();
				TensorCPT::ElementwiseMultiply(*new_outgoing, *outgoing_message, *message_vec[j]);
				delete(outgoing_message);
				outgoing_message = new_outgoing;
			}
		}
		
		downward_message_map[clique_index] = outgoing_message;
		TensorCPT *integrated_message = new TensorCPT();

		vector<int>& one_reduce_vars = GetRVars(clique_index);
		TensorCPT::OneReduce(*integrated_message, *outgoing_message, one_reduce_vars);
		integrated_message_map[clique_index] = integrated_message;
		if (clique_index != root)
		{
			upward_message_map[GetParent(clique_index)]->push_back(integrated_message); 
		}
		if (clique_index == root)
		{
		//	downward_message_map[clique_index] = integrated_message;
			vector<int>& dims = integrated_message->Dims();
			assert(dims.size() == 0);
			joint_val = integrated_message->At(0);
			break;
		}
	}

	for (int i = 0; i < upward_message_map.size(); ++i)
	{
		//for (int j = 0; j < upward_message_map.at(i)->size(); ++j)
	//	{
	//		delete(upward_message_map.at(i)->at(j));
	//	}
		delete(upward_message_map.at(i));
	}

	for (int i = 0; i < downward_message_map.size(); ++i)
	{
		delete(downward_message_map.at(i));
		delete(integrated_message_map.at(i));
	}
	

	return joint_val;
}
double JTree::ComputeEmpiricalMarginals(vector<TensorCPT*>& marginals, vector<int>& evidence_vars, vector<int>& evidence_vals)
{
	double joint_val = 0;
	vector<vector<TensorCPT*>*> upward_message_map;
	vector<TensorCPT*> downward_message_map;
	vector<TensorCPT*> integrated_message_map;

	marginals.assign(NumCliques(), NULL);
	for (int i = 0; i < NumCliques(); ++i)
	{
		upward_message_map.push_back(new vector<TensorCPT*>());
		downward_message_map.push_back(NULL);
		integrated_message_map.push_back(NULL);
	}

	for (int i = 0; i < NumCliques(); ++i)
	{
		int clique_index = elimination_order->at(i);

		vector<int> evidence_i_vars;
		vector<int> evidence_i_vals;

		VectorPlus::Intersect(evidence_i_vars, GetCVars(clique_index), evidence_vars);
		VectorPlus::MatchSub(evidence_i_vals, evidence_vals, evidence_vars, evidence_i_vars);

		TensorCPT* outgoing_message;

		if (node_list->at(clique_index)->GetNumFactors() == 0)
		{
			vector<TensorCPT*>& message_vec = *(upward_message_map[clique_index]);
			for (int j = 0; j < message_vec.size(); ++j)
			{
				if (j == 0)
				{
					outgoing_message = new TensorCPT(*message_vec[j]);
				}
				else
				{
					TensorCPT* new_outgoing = new TensorCPT();
					TensorCPT::ElementwiseMultiply(*new_outgoing, *outgoing_message, *message_vec[j]);
					delete(outgoing_message);
					outgoing_message = new_outgoing;
				}
			}
		}
		else
		{
			outgoing_message = new TensorCPT();
			TensorCPT& old_cpt = GetAgglomCliqueFactor(clique_index);

			old_cpt.Select(*outgoing_message, evidence_i_vars, evidence_i_vals); 

			vector<TensorCPT*>& message_vec = *(upward_message_map[clique_index]);
			for (int j = 0; j < message_vec.size(); ++j)
			{
				TensorCPT* new_outgoing = new TensorCPT();
				TensorCPT::ElementwiseMultiply(*new_outgoing, *outgoing_message, *message_vec[j]);
				delete(outgoing_message);
				outgoing_message = new_outgoing;
			}
		}
		
		downward_message_map[clique_index] = outgoing_message;
		TensorCPT *integrated_message = new TensorCPT();

		vector<int>& one_reduce_vars = GetRVars(clique_index);
		TensorCPT::OneReduce(*integrated_message, *outgoing_message, one_reduce_vars);
		integrated_message_map[clique_index] = integrated_message;
		if (clique_index != root)
		{
			upward_message_map[GetParent(clique_index)]->push_back(integrated_message); 
		}
		if (clique_index == root)
		{
		//	downward_message_map[clique_index] = integrated_message;
			vector<int>& dims = integrated_message->Dims();
			assert(dims.size() == 0);
			joint_val = integrated_message->At(0);
			break;
		}
	}

	for (int i = 0; i < NumCliques(); ++i)
	{
		int clique_index = elimination_order->at(NumCliques() - i - 1);
		TensorCPT*  marginal_tensor = new TensorCPT(); 
		TensorCPT& i_tensor = *(downward_message_map[clique_index]);
		

		if (clique_index == root)
		{
			marginal_tensor->Initialize(i_tensor);
		}
		else
		{
			TensorCPT& parent_tensor = *(marginals[GetParent(clique_index)]);
			TensorCPT integrated_parent_message;
			vector<int> parent_reduce_vars;
			VectorPlus::SetDiff(parent_reduce_vars, GetCVars(GetParent(clique_index)), GetSVars(clique_index));
			TensorCPT::OneReduce(integrated_parent_message, parent_tensor, parent_reduce_vars);
			TensorCPT& integrated_message = *(integrated_message_map[clique_index]);;
			TensorCPT divided_message;
		//	vector<int>& one_reduce_vars = GetRVars(clique_index);
		//	TensorCPT::OneReduce(integrated_message, i_tensor, one_reduce_vars);
			integrated_message.ElementwiseInvert();
			TensorCPT::ElementwiseMultiply(divided_message, integrated_message, integrated_parent_message);

			TensorCPT::ElementwiseMultiply(*marginal_tensor,  i_tensor, divided_message);
		}

		marginals[clique_index] = marginal_tensor;

	}

	for (int i = 0; i < upward_message_map.size(); ++i)
	{
		//for (int j = 0; j < upward_message_map.at(i)->size(); ++j)
	//	{
	//		delete(upward_message_map.at(i)->at(j));
	//	}
		delete(upward_message_map.at(i));
	}

	for (int i = 0; i < downward_message_map.size(); ++i)
	{
		delete(downward_message_map.at(i));
		delete(integrated_message_map.at(i));
	}
	

	return joint_val;
}

TensorCPT& JTree::GetAgglomCliqueFactor(int clique_index) { return node_list->at(clique_index)->GetAgglomCliqueFactor(); }
vector<int>& JTree::GetChildren(int clique_index) { return node_list->at(clique_index)->GetChildren(); }
int JTree::GetParent(int clique_index) {return node_list->at(clique_index)->GetParent(); }
vector<int>& JTree::GetCVars(int clique_index) {return node_list->at(clique_index)->GetCVars(); }
vector<int>& JTree::GetRVars(int clique_index) {return node_list->at(clique_index)->GetRVars(); }
vector<int>& JTree::GetSVars(int clique_index) {return node_list->at(clique_index)->GetSVars(); }


bool JTree::TestConvergence(double prev_likelihood, double curr_likelihood, double thresh)
{
	double avg = abs(prev_likelihood + curr_likelihood)  / 2.0;
	double diff = abs(curr_likelihood - prev_likelihood);

	if ((diff / avg) < thresh)
	{
		cout << avg;
		cout << "\n";
		cout << diff;
		cout << "\n";
		cout << diff / avg;
		cout << "\n";
		return true;
	}
	else
		return false;
}

JTree::~JTree()
{
	delete(tree_matrix);
	for (int i = 0; i < node_list->size(); ++i)
	{
		delete(node_list->at(i));
	}
	delete(node_list);
	delete(all_obs_vars);
	
	delete(elimination_order);
}
