# pragma once

#include "BayesianNetwork.hpp"
#include "VectorPlus.hpp"
#include "StatsLibrary.hpp"
#include <queue>
#include <assert.h>
#include <iostream>
#include "MultiVector.hpp"
#include "TensorVarElim.hpp"
#include "FastIndexer.hpp"
#include <set>

using namespace std;

// default constructor
BayesianNetwork::BayesianNetwork()
{
	parent_lists = NULL;
	children_lists = NULL;
	descendant_lists = NULL;
	dims = NULL;
	nodes = NULL;
	CPTs = NULL;
	type_vector = NULL;
	topological_order = NULL; // some topological order
	reverse_topological_order = NULL;
	elimination_order = NULL;
	isInitialized = false;

	CPTs = new vector<CPT*>();

	num_hidden_states = -1;
	samples = NULL;

	empirical_prob_map = new map<vector<int>, Tensor*>();
	empirical_multi_prob_map = new map<MultiVector<int>, Tensor*>();
}

BayesianNetwork::BayesianNetwork(vector<vector<int>*>& parent_lists, vector<int>& type_vector, vector<int>& dims)
{
	this->parent_lists = NULL;
	this->children_lists = NULL;
	this->descendant_lists = NULL;
	this->dims = NULL;
	this->nodes = NULL;
	this->CPTs = new vector<CPT*>();
	this->type_vector = NULL;
	this->topological_order = NULL;
	reverse_topological_order = NULL;
	isInitialized = false;
	empirical_multi_prob_map = new map<MultiVector<int>, Tensor*>();
	empirical_prob_map = new  map<vector<int>, Tensor*>();
	samples = NULL;
	Initialize(parent_lists, type_vector, dims);

	elimination_order = NULL;

	GenerateRandomCPTs();
}

BayesianNetwork::~BayesianNetwork()
{
	for (int i = 0; i < parent_lists->size(); ++i)
	{
		delete(parent_lists->at(i));
	}
	delete(parent_lists);

	for (int i = 0; i < children_lists->size(); ++i)
	{
		delete(children_lists->at(i));
	}
	delete(children_lists);

	for (int i = 0; i < descendant_lists->size(); ++i)
	{
		delete(descendant_lists->at(i));
	}
	delete(descendant_lists);
	
	for (int i = 0; i < CPTs->size(); ++i)
	{
		delete(CPTs->at(i));
	}
	delete(CPTs);

	delete(type_vector);
	delete(dims);
	delete(nodes);

	delete(topological_order);
	delete(reverse_topological_order);

	map<vector<int>, Tensor*>::iterator prob_iter;

	for (prob_iter = empirical_prob_map->begin(); prob_iter != empirical_prob_map->end(); ++prob_iter)
	{
		delete(prob_iter->second);
	}
	delete(empirical_prob_map);

	map<MultiVector<int>, Tensor*>::iterator multi_prob_iter;
	for (multi_prob_iter = empirical_multi_prob_map->begin(); multi_prob_iter != empirical_multi_prob_map->end(); ++multi_prob_iter)
	{
		delete(multi_prob_iter->second);
	}
	delete(empirical_multi_prob_map);

	delete(samples);
}

void BayesianNetwork::Initialize(vector<vector<int>*>& parent_lists, vector<int>& type_vector, vector<int>& dims)
{
	int num_nodes = type_vector.size();
	this->nodes = new vector<int>();
	VectorPlus::Seq(*nodes, 0, 1, num_nodes);
	this->dims = new vector<int>(dims);
	this->type_vector = new vector<int>(type_vector);
	this->parent_lists = new vector< vector<int>* >(num_nodes, NULL);
	this->children_lists = new vector<vector<int>*>(num_nodes, NULL);
	this->descendant_lists = new vector< vector<int>* >(num_nodes, NULL);
	topological_order = NULL;

	for (int i = 0; i < num_nodes; ++i)
	{
		this->parent_lists->at(i) = new vector<int>(*parent_lists[i]);
		this->children_lists->at(i) = new vector<int>();
		this->descendant_lists->at(i) = new vector<int>();
	}

	for (int i = 0; i < num_nodes; ++i)
	{
		for (int j = 0; j < parent_lists[i]->size(); ++j)
		{
			int parent = this->parent_lists->at(i)->at(j);
			children_lists->at(parent)->push_back(i);
		}
	}
/*	
	for (int i = 0; i < num_nodes; ++i)
	{
		cout << i;
		cout << "\n";
		queue<int> ascendant_queue;
		ascendant_queue.push(i);
		while(!ascendant_queue.empty())
		{
			int node = ascendant_queue.front();
			ascendant_queue.pop();
			for (int j = 0; j < this->parent_lists->at(node)->size(); ++j)
			{
				ascendant_queue.push(this->parent_lists->at(node)->at(j));
			}
			if (!VectorPlus::Contains(*descendant_lists->at(node), i) && node != i)
			{
				descendant_lists->at(node)->push_back(i);
			}
		}
	} */

	ComputeTopologicalOrder();
	
	int var = VectorPlus::Find(*(this->type_vector), HIDDEN_FLAG);
	num_hidden_states = this->dims->at(var);

	CPTs->assign(num_nodes, NULL);
}

// generate a random CPT
void BayesianNetwork::GenerateRandomCPT(int node_id)
{	
	vector<int> i_nodes;
	vector<int> i_dims;
	vector<int> node_id_vec =  VectorPlus::CreateSingleton(node_id);
	VectorPlus::Concat(i_nodes, node_id_vec, *parent_lists->at(node_id));
	VectorPlus::Subset(i_dims, *dims, i_nodes);
	CPTs->at(node_id) = new CPT(i_nodes, i_dims);
	CPTs->at(node_id)->FillWithRandom();
}


// generate a random CPT
void BayesianNetwork::GenerateRandomCPTs()
{	
	for (int i = 0; i < NumNodes(); ++i)
	{
		GenerateRandomCPT(i);
	}
}

// Brute Force topological order, takes cubic time
void BayesianNetwork::ComputeTopologicalOrder()
{
	assert(topological_order == NULL);
	topological_order = new vector<int>();
	reverse_topological_order = new vector<int>();
	for (int i = 0; i < NumNodes(); ++i)
	{
		int node = -1; 
		for (int j = 0; j < NumNodes(); ++j)
		{
			if (!VectorPlus::Contains(*topological_order, j))
			{
				bool top_flag = true;
				for (int k = 0; k < parent_lists->at(j)->size(); ++k)
				{
					if (!VectorPlus::Contains(*topological_order, parent_lists->at(j)->at(k)))
					{
						top_flag = false;
					}
				}
				if (top_flag == true)
				{
					node = j;
					break;
				}
			}
		}
		assert(node != -1);
		topological_order->push_back(node); 
	}
	VectorPlus::Reverse(*reverse_topological_order, *topological_order);
}

void BayesianNetwork::GenerateSamples(Matrix& samples, int num_samples)
{
	vector<int> dims = VectorPlus::CreatePair(num_samples, NumNodes());

	samples.Initialize(dims);

	for (int n = 0; n < num_samples; ++n)
	{
		for (int i = 0; i < NumNodes(); ++i)
		{
			int node = topological_order->at(i);
			vector<int>& parent_vec = *(parent_lists->at(node));
			vector<double> parent_vals;
			samples.Sub(parent_vals, n, parent_vec);
			vector<int> i_parent_vals;
			for (int j = 0; j < parent_vals.size(); ++j)
				i_parent_vals.push_back((int)parent_vals[j]);

			double val = CPTs->at(node)->Sample(i_parent_vals);
			samples.Set(n,node, val);
		}
	}
}

void BayesianNetwork::GenerateTrainSamples(Matrix& samples, int num_samples)
{
	GenerateSamples(samples, num_samples);
	assert(this->samples == NULL);
	this->samples = new Matrix(samples);
}

int BayesianNetwork::GetNumHiddenStates()
{
	return num_hidden_states;
}

void BayesianNetwork::ClearEmpiricalProbabilityMap()
{
	empirical_prob_map->clear();
}

double BayesianNetwork::ComputeEmpiricalProbVal(vector<int>& vars, vector<int>& vals)
{
	double count = 0;
	for (int n = 0; n < samples->Dim(0); ++n)
	{
		vector<int> sample_vars;
		for (int i = 0; i < vars.size(); ++i)
		{
			sample_vars.push_back((int)(samples->At(n, vars[i])));
		}
		if (VectorPlus::Equals(sample_vars, vals))
			count++;
	}
	return count / samples->NumElements();

}

void BayesianNetwork::ComputeConditionalProbTensor(Tensor& conditional_prob, vector<int>& front_vars, vector<int>& back_vars)
{
	assert(0);
	vector<int> all_vars;
	VectorPlus::Concat(all_vars, front_vars, back_vars);

	vector<int> var_dims;
	VectorPlus::Subset(var_dims, *dims, all_vars);

	Tensor joint_prob;
	Tensor marginal_prob;
	ComputeEmpiricalProbTensor(joint_prob, all_vars);
	ComputeEmpiricalProbTensor(marginal_prob, back_vars);
	conditional_prob.Initialize(var_dims);

	vector<int> back_modes; 
	VectorPlus::Seq(back_modes, front_vars.size(), 1, all_vars.size());
	for (int i = 0; i < conditional_prob.NumElements(); ++i)
	{
		vector<int> indices;
		vector<int> back_indices;
		conditional_prob.ComputeIndexArray(indices, i);
		VectorPlus::Subset(back_indices, indices, back_modes);
		double numer = joint_prob.At(indices);
		double denom = marginal_prob.At(back_indices);
		conditional_prob.Set(i, numer / denom);
	}
}

void BayesianNetwork::ComputeConditionalMultiProbTensor(Tensor& conditional_multi_prob, MultiVector<int>& front_multi_vars, MultiVector<int>& back_multi_vars)
{
	assert(0);
	MultiVector<int> all_multi_vars;
	MultiVector<int>::Concat(all_multi_vars, front_multi_vars, back_multi_vars);

	vector<int> multi_var_dims;
	for (int i = 0; i < all_multi_vars.Size(); ++i)
	{
		vector<int>& stacked_vars = all_multi_vars.GetVec(i);
		vector<int> stacked_dims;
		VectorPlus::Subset(stacked_dims, *dims, stacked_vars);
		multi_var_dims.push_back(VectorPlus::Sum(stacked_dims));
	}
	
	conditional_multi_prob.Initialize(multi_var_dims);

	vector<int>& var_combos = all_multi_vars.GetVecSizes();
	vector<int> combo_offsets;
	Tensor::ComputeOffsets(combo_offsets, var_combos);

	int num_var_combos = VectorPlus::Product(var_combos);

	vector<int> front_modes;
	vector<int> back_modes;
	VectorPlus::Seq(front_modes, 0, 1, front_multi_vars.Size());
	VectorPlus::Seq(back_modes, front_modes.size(), 1, all_multi_vars.Size());
	for (int i = 0; i < num_var_combos; ++i)
	{
		vector<int> combo_indices;
		Tensor::ComputeIndexArray(combo_indices, combo_offsets, i); 

		vector<int> front_indices;
		vector<int> back_indices;
		VectorPlus::Subset(front_indices, combo_indices, front_modes);
		VectorPlus::Subset(back_indices, combo_indices, back_indices);
		

		vector<int> front_vars;
		vector<int> back_vars;
		vector<int> all_vars;
		front_multi_vars.At(front_vars, front_indices);
		back_multi_vars.At(back_vars, back_indices);
		VectorPlus::Concat(all_vars, front_vars, back_vars);

		Tensor prob_tensor;
		ComputeConditionalProbTensor(prob_tensor, front_vars, back_vars);

		for (int j = 0; j < prob_tensor.NumElements(); ++j)
		{
			vector<int> indices;
			vector<int> multi_indices;
			prob_tensor.ComputeIndexArray(indices, j);
			ConvertIndices(multi_indices, indices, all_multi_vars, all_vars);
			conditional_multi_prob.Set(multi_indices, prob_tensor.At(j));
		}
	}
}

void BayesianNetwork::ComputeOracleProbTensor(Tensor& prob_tensor, vector<int>& vars)
{
	for (int i = 0; i < vars.size(); ++i)
		assert(type_vector->at(vars[i]) == OBSERVED_FLAG);

	if (empirical_prob_map->count(vars) == 0)
	{
		vector<int> var_dims;
		VectorPlus::Subset(var_dims, *dims, vars);
		prob_tensor.Initialize(var_dims);
		prob_tensor.FillWithConst(0);

		for (int i = 0; i < prob_tensor.NumElements(); ++i)
		{
			vector<int> var_vals;
			prob_tensor.ComputeIndexArray(var_vals, i);
			double prob = TensorVarElim::VE(*this, vars, var_vals, *reverse_topological_order); 
			prob_tensor.Set(i, prob);
		}

	//	return *prob_tensor;
	}
	else
	{
	//	return *((*empirical_prob_map)[vars]);
	}



}

void BayesianNetwork::ComputeEmpiricalProbTensor(Tensor& prob_tensor, vector<int>& vars)
{
	if (empirical_prob_map->count(vars) == 0)
	{
		vector<int> var_dims;
		VectorPlus::Subset(var_dims, *dims, vars);
		prob_tensor.Initialize(var_dims);
		prob_tensor.FillWithConst(0);
		for (int n = 0; n < samples->Dim(0); ++n)
		{	
			vector<int> var_vals;
			for (int i = 0; i < vars.size(); ++i)
			{
				var_vals.push_back((int)(samples->At(n, vars[i])));
			}
			double curr_val = prob_tensor.At(var_vals) * samples->Dim(0);
			prob_tensor.Set(var_vals, ((curr_val + 1) / samples->Dim(0)));
		}
	//	prob_tensor.Add(0.0000000001);
	//	(*empirical_prob_map)[vars] = prob_tensor;
		//return *prob_tensor;
	}
	else
	{
		//return *((*empirical_prob_map)[vars]);
	}
}

void BayesianNetwork::ComputeEmpiricalMultiProbTensor(Tensor& multi_prob_tensor, MultiVector<int>& multi_vars)
{
	if (empirical_multi_prob_map->count(multi_vars) == 0)
	{
		vector<int> multi_var_dims;
		for (int i = 0; i < multi_vars.Size(); ++i)
		{
			vector<int>& stacked_vars = multi_vars.GetVec(i);
			vector<int> stacked_dims;
			VectorPlus::Subset(stacked_dims, *dims, stacked_vars);
			multi_var_dims.push_back(VectorPlus::Sum(stacked_dims));
		}
		multi_prob_tensor.Initialize(multi_var_dims);

		vector<int>& var_combos = multi_vars.GetVecSizes();
		vector<int> combo_offsets;
		Tensor::ComputeOffsets(combo_offsets, var_combos);

		int num_var_combos = VectorPlus::Product(var_combos);

		FastIndexer combo_indexer(var_combos);
		for (int i = 0; i < num_var_combos; ++i)
		{
			vector<int>& combo_indices = combo_indexer.GetNext();
		//	Tensor::ComputeIndexArray(combo_indices, combo_offsets, i); 

			vector<int> vars;
			multi_vars.At(vars, combo_indices); 

			Tensor prob_tensor;
			ComputeEmpiricalProbTensor(prob_tensor, vars);
		//	ComputeOracleProbTensor(prob_tensor, vars);
			
			FastIndexer indexer(prob_tensor.Dims());
		//	for (int j = 0; j < prob_tensor.NumElements(); ++j)
			int j = 0;
			while (indexer.HasNext())
			{
				vector<int>& indices = indexer.GetNext();
			//	vector<int> multi_indices;
			//	vector<int> indices;
			//	prob_tensor.ComputeIndexArray(indices, j);
			//	ConvertIndices(multi_indices, indices, multi_vars, vars);
			//	assert(VectorPlus::Equals(indices, multi_indices));
			//	multi_prob_tensor.Set(multi_indices, prob_tensor.At(j));
				multi_prob_tensor.Set(indices, prob_tensor.At(j++));
			}
		}

	//	(*empirical_multi_prob_map)[multi_vars] = multi_prob_tensor;
	//	return *multi_prob_tensor;
	}
	else
	{
	//	return *((*empirical_multi_prob_map)[multi_vars]);
	}
}

void BayesianNetwork::ConvertIndices(vector<int>& multi_indices, vector<int>& indices, MultiVector<int>& multi_vars, vector<int>& vars)
{
	vector<int> offset;
	for (int i = 0; i < multi_vars.Size(); ++i)
	{
		vector<int>& i_vars = multi_vars.GetVec(i);

		int off_val = 0;
		for (int j = 0; j < i_vars.size(); ++j)
		{
			if (i_vars[j] != vars[i])
			{
				off_val += dims->at(i_vars[j]);
			}
			else
			{
				break;
			}
		}

		offset.push_back(off_val);
	}

	VectorPlus::Add(multi_indices, offset, indices);
}

bool BayesianNetwork::IsAdj(int node1, int node2)
{
	vector<int>& parents = *(parent_lists->at(node1));
	vector<int>& children = *(children_lists->at(node1));

	if (VectorPlus::Contains(parents, node2) || VectorPlus::Contains(children, node2))
		return true;
	else
		return false;
}

// assumes max clique size is four
void BayesianNetwork::FindCliques(vector<vector<int>*>& nodes_to_vars)
{
	
	set< vector<int> > clique_set;




	// find quintets

/*	for (int i = 0; i < NumNodes(); ++i)
	{
		for (int j = i + 1; j < NumNodes(); ++j)
		{
			for (int k = j + 1; k < NumNodes(); ++k)
			{
				for (int l = k + 1; l < NumNodes(); ++l)
				{
					for (int m = l + 1; m < NumNodes(); ++m)
					{
					// every quartet of nodes, test if there is an edge among all four of them
					if (IsAdj(i,j) && IsAdj(i,k) && IsAdj(i,l) && IsAdj(i,m) && IsAdj(j,k) && IsAdj(j,l) && IsAdj(j,m) && IsAdj(k,l) && IsAdj(k,m) && IsAdj(l,m))
					{
						vector<int>* new_vec = new vector<int>();
						new_vec->push_back(i);
						new_vec->push_back(j);
						new_vec->push_back(k);
						new_vec->push_back(l);
						new_vec->push_back(m);
						nodes_to_vars.push_back(new_vec);
						clique_set.insert(VectorPlus::CreateQuartet(i,j,k,l));
						clique_set.insert(VectorPlus::CreateQuartet(j,k,l,m));
						clique_set.insert(VectorPlus::CreateQuartet(i,k,l,m));
						clique_set.insert(VectorPlus::CreateQuartet(i,j,k,m));
						clique_set.insert(VectorPlus::CreateQuartet(i,j,l,m));


						clique_set.insert(VectorPlus::CreateTriple(i,j,k));
						clique_set.insert(VectorPlus::CreateTriple(i,j,l));
						clique_set.insert(VectorPlus::CreateTriple(i,j,m));
						clique_set.insert(VectorPlus::CreateTriple(i,k,l));
						clique_set.insert(VectorPlus::CreateTriple(i,k,m));
						clique_set.insert(VectorPlus::CreateTriple(i,l,m));

						clique_set.insert(VectorPlus::CreateTriple(j,k,l));
						clique_set.insert(VectorPlus::CreateTriple(j,k,m));
						clique_set.insert(VectorPlus::CreateTriple(j,l,m));
						clique_set.insert(VectorPlus::CreateTriple(k,l,m));

						clique_set.insert(VectorPlus::CreatePair(i,j));
						clique_set.insert(VectorPlus::CreatePair(i,k));
						clique_set.insert(VectorPlus::CreatePair(i,l));
						clique_set.insert(VectorPlus::CreatePair(i,m));

						clique_set.insert(VectorPlus::CreatePair(j,k));
						clique_set.insert(VectorPlus::CreatePair(j,l));
						clique_set.insert(VectorPlus::CreatePair(j,m));

						clique_set.insert(VectorPlus::CreatePair(j,k));
						clique_set.insert(VectorPlus::CreatePair(k,l));
					}
					}
				}
			}
		}
	} */


	// find quartets
	for (int i = 0; i < NumNodes(); ++i)
	{
		for (int j = i + 1; j < NumNodes(); ++j)
		{
			for (int k = j + 1; k < NumNodes(); ++k)
			{
				for (int l = k + 1; l < NumNodes(); ++l)
				{
					
					// every quartet of nodes, test if there is an edge among all four of them
					if (IsAdj(i,j) && IsAdj(i,k) && IsAdj(i,l) && IsAdj(j,k) && IsAdj(j,l) && IsAdj(k,l))
					{
						vector<int>* new_vec = new vector<int>();
						new_vec->push_back(i);
						new_vec->push_back(j);
						new_vec->push_back(k);
						new_vec->push_back(l);
						nodes_to_vars.push_back(new_vec);
						clique_set.insert(VectorPlus::CreateTriple(i,j,k));
						clique_set.insert(VectorPlus::CreateTriple(i,j,l));
						clique_set.insert(VectorPlus::CreateTriple(i,k,l));
						clique_set.insert(VectorPlus::CreateTriple(j,k,l));
						clique_set.insert(VectorPlus::CreatePair(i,j));
						clique_set.insert(VectorPlus::CreatePair(i,k));
						clique_set.insert(VectorPlus::CreatePair(i,l));
						clique_set.insert(VectorPlus::CreatePair(j,k));
						clique_set.insert(VectorPlus::CreatePair(j,l));
						clique_set.insert(VectorPlus::CreatePair(k,l));
					}
				}
			}
		}
	}


	for (int i = 0; i < NumNodes(); ++i)
	{
		for (int j = i + 1; j < NumNodes(); ++j)
		{
			for (int k = j + 1; k < NumNodes(); ++k)
			{

				// every quartet of nodes, test if there is an edge among all four of them
				if (IsAdj(i,j) && IsAdj(j,k) && IsAdj(i,k) && clique_set.count(VectorPlus::CreateTriple(i,j, k)) == 0)
				{
					vector<int>* new_vec = new vector<int>();
					new_vec->push_back(i);
					new_vec->push_back(j);
					new_vec->push_back(k);
					nodes_to_vars.push_back(new_vec);
	
					clique_set.insert(VectorPlus::CreatePair(i,j));
					clique_set.insert(VectorPlus::CreatePair(j,k));
					clique_set.insert(VectorPlus::CreatePair(i,k));
				}
			}
		}
	}

	// finds pairs
	for (int i = 0; i < NumNodes(); ++i)
	{
		for (int j = i + 1; j < NumNodes(); ++j)
		{
			
			if (IsAdj(i,j) && clique_set.count(VectorPlus::CreatePair(i,j)) == 0)
			{
				vector<int>* new_vec = new vector<int>();
				new_vec->push_back(i);
				new_vec->push_back(j);
				nodes_to_vars.push_back(new_vec);
			}
		}
	}

}


void BayesianNetwork::FindJunctionTree(Matrix& jtree_matrix, vector<vector<int>*> cliques)
{
	// create adjacency matrix
	int num_cliques = cliques.size();
	vector<int> tree_matrix_dims(2,num_cliques);
	Matrix adj_matrix(tree_matrix_dims);
	jtree_matrix.Initialize(tree_matrix_dims);
	for (int i = 0; i < num_cliques; ++i)
	{
		for (int j = 0; j < num_cliques; ++j)
		{
			if (i == j)
				adj_matrix.Set(i,j, 0);
			else
			{
				vector<int> intersect_vars;
				VectorPlus::Intersect(intersect_vars, *(cliques[i]), *(cliques[j]));
				adj_matrix.Set(i,j, intersect_vars.size());
			}
		}
	}

	set<int> curr_set;
	set<int> other_set;
	
	curr_set.insert(0);

	for (int i = 1; i < num_cliques; ++i)
	{
		other_set.insert(i);
	}

	while (other_set.size() > 0)
	{
		int max_row = -1;
		int max_col = -1;
		double max_val = -1;
		for (int i = 0; i < num_cliques; ++i)
		{
			for (int j = 0; j < num_cliques; ++j)
			{
				if (curr_set.count(i) == 1 && other_set.count(j) == 1 && i != j)
				{
					if (adj_matrix.At(i,j) > max_val)
					{
						max_val = adj_matrix.At(i,j);
						max_row = i;
						max_col = j;
					}

				}
			}
		}
		assert(max_row >= 0 && max_col >= 0);
		curr_set.insert(max_col);
		other_set.erase(max_col);
		jtree_matrix.Set(max_col, max_row, 1.0);
	}
}
