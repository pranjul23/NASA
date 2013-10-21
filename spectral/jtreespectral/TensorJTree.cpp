# pragma once

#include "MultiVector.hpp"
#include "TensorCPT.hpp"
#include "tensor.hpp"
#include <assert.h>
#include <algorithm>
#include <math.h>
#include <queue>
#include "VectorPlus.hpp"
#include <assert.h>
#include "StatsLibrary.hpp"
#include <map>
#include "TensorJTree.hpp"
#include <iostream>
#include <queue>
#include <set>
#include <Eigen/Dense>
#include <Eigen/SVD>

#include <boost/shared_ptr.hpp>

using namespace std;

// default constructor
TensorJTree::TensorJTree()
{
	bNet = NULL;
	tree_matrix = NULL;
	node_list = NULL;
	elimination_order = NULL;
	root = -1; 
}

TensorJTree::~TensorJTree()
{
	delete(tree_matrix);
	for (int i = 0; i < node_list->size(); ++i)
	{
		delete(node_list->at(i));
	}
	delete(node_list);
	delete(all_obs_vars);
	delete(elimination_order);

	delete(clique_to_template);

	for (int i = 0; i < template_to_cliques->size(); ++i)
	{
		delete(template_to_cliques->at(i));
	}
	
	delete(template_to_cliques);
}
// default constructor
TensorJTree::TensorJTree(BayesianNetwork& bNet)
{

}

// default constructor
TensorJTree::TensorJTree(BayesianNetwork* bNet, Matrix& tree_matrix, vector<vector<int>*>& nodes_to_vars, int max_obs_per_mode, vector<int>& template_vec, int max_runs)
{
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


	this->node_list = new vector<TensorJTNode*>();
	for (int i = 0; i < tree_matrix.Dim(0); ++i)
	{
		node_list->push_back(new TensorJTNode(this, i, max_obs_per_mode, max_runs));
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


	}

	// store C, S, and R vars
	for (int i = 0; i < NumCliques(); ++i)
	{
		vector<int> i_Sep;
		vector<int>& i_vars = *(nodes_to_vars[i]);

		int parent = GetParent(i);

		if (parent >= 0)
		{
			vector<int>& parent_vars = *(nodes_to_vars[parent]);
			VectorPlus::Intersect(i_Sep, i_vars, parent_vars);
		}

		node_list->at(i)->SetSVars(i_Sep);
	}

	for (int i = 0; i < NumCliques(); ++i)
	{
		vector<int> reduced_C_vars;
		VectorPlus::Union(reduced_C_vars, GetSVars(i));
		
		for (int j = 0; j < nodes_to_vars[i]->size(); ++j)
		{
			if (bNet->is_observed(nodes_to_vars[i]->at(j)))
			{
				if (!VectorPlus::Contains(reduced_C_vars, nodes_to_vars[i]->at(j)))
				{
					reduced_C_vars.push_back(nodes_to_vars[i]->at(j));
				}
			}
		}

		vector<int>& children = GetChildren(i);
		for (int j = 0; j < children.size(); ++j)
		{
			VectorPlus::Union(reduced_C_vars, GetSVars(children[j]));
		}
		node_list->at(i)->SetCVars(reduced_C_vars);
		
		int parent = GetParent(i);
		if (parent >= 0)
		{
			vector<int>& parent_vars = *(nodes_to_vars[parent]);
			vector<int> i_Res;
			VectorPlus::SetDiff(i_Res, reduced_C_vars, GetSVars(i));
			node_list->at(i)->SetRVars(i_Res);
		}
		else
		{
			node_list->at(i)->SetRVars(reduced_C_vars);
		}

	}
	
	for (int i = 0; i < NumCliques(); ++i)
	{
		int parent = GetParent(i);
		vector<int>& children = GetChildren(i);

		for (int j = 0; j < children.size(); ++j)
		{
			vector<int>& j_sep = node_list->at(children[j])->GetSVars();
			node_list->at(i)->UniqueChildSepAdd(j_sep, children[j]);
		}
	}


	ComputeEliminationOrder();

	clique_to_template = new vector<int>(template_vec);

	vector<int> unique_templates;
	VectorPlus::Unique(unique_templates, *clique_to_template);
	template_to_cliques = new vector<vector<int>*>();
	for (int i = 0; i <  unique_templates.size(); ++i)
	{
		template_to_cliques->push_back(new vector<int>());
	}
	for (int i = 0; i < clique_to_template->size(); ++i)
	{
		template_to_cliques->at(clique_to_template->at(i))->push_back(i);
	}
	vector<int> dist_dims(2, NumCliques());
	dist_vars = new vector<vector<int>*>();
	for (int i = 0; i < NumCliques(); ++i)
	{
		dist_vars->push_back(new vector<int>());
	}
	ComputeDistances();
}

void TensorJTree::LearnTemplateSpectralParameters()
{
	for (int i = 0; i < node_list->size(); ++i)
	{
		node_list->at(i)->MarkMultModes();
		FindObservationPartitions(i);
	}

	for (int i = 0; i < node_list->size(); ++i)
	{
		node_list->at(i)->MarkObservableRepresentation();
	}

	for (int i = 0; i < node_list->size(); ++i)
	{
		cout << i;
		cout << "\n";
		node_list->at(i)->ComputeObservableRepresentation();
	}
		cout << "yeee";
	cout << "\n";
	for (int i = 0; i < node_list->size(); ++i)
	{
				cout << i;
		cout << "\n";
		node_list->at(i)->FinalizeTemplateObservableRepresentation();
	}

		cout << "yeee";
	cout << "\n";
/*	template_to_cliques = new vector<vector<int>*>();
	for (int i = 0; i < NumCliques(); ++i)
	{
		template_to_cliques->push_back(new vector<int>(1,i));
	} */
	cout << "yeee";
	cout << "\n";
	for (int i = 0; i < template_to_cliques->size(); ++i)
	{
		vector<int>& i_cliques = *(template_to_cliques->at(i));

		vector<double> big_forward_vec;
		Matrix big_back_matrix;
		vector<int> result_dims;

		bool invModesExist = node_list->at(i_cliques[0])->InvModesExists();
		if (invModesExist)
		{

			for (int j = 0; j < i_cliques.size(); ++j)
			{
				cout << j;
				cout << "\n";
				if (node_list->at(i_cliques[0])->is_invalid())
				{
					assert(0);
					continue;
				}


				int index = i_cliques[j];
				vector<Tensor*>& forward_list = node_list->at(index)->GetForwardList();
				vector<Tensor*>& extra_list = node_list->at(index)->GetExtraList();

				Tensor& X = node_list->at(index)->GetTransformTensor().GetTensor();

				if (j == 0)
				{
					result_dims = X.Dims();
				}

				assert(forward_list.size() == extra_list.size());
		
				vector<int> inv_group_modes =  node_list->at(index)->GetInvModeGroup();
				vector<int> inv_modes = node_list->at(index)->GetInvModes();

				for (int k = 0; k < forward_list.size(); ++k)
				{

					Tensor& forward_tensor = *(forward_list[k]);
					Tensor& back_tensor = *(extra_list[k]);
	
					if (j == 0 && k == 0)
					{
						Tensor::CreateLinearSystem(big_forward_vec, big_back_matrix, X, back_tensor, forward_tensor, inv_group_modes, inv_modes); 
					}
					else
					{
						vector<double> forward_vec;
						Matrix back_matrix;

						cout << "yeebeforesystem";
						Tensor::CreateLinearSystem(forward_vec, back_matrix, X, back_tensor, forward_tensor, inv_group_modes, inv_modes); 
						cout << "yeeaftersystem";
						VectorPlus::Concat(big_forward_vec, forward_vec);
						big_back_matrix.MatrixConcat(back_matrix);
						cout << "yeeafterconcat";
					}
				}
				node_list->at(index)->DeleteForwardList();
				node_list->at(index)->DeleteExtraList();
			}
			cout << big_back_matrix.Dim(0);
			cout << "\n";
			cout << big_back_matrix.Dim(1);
			cout << "\n";
			Eigen::MatrixXd A(big_back_matrix.Dim(0), big_back_matrix.Dim(1));

			cout << "allocatedbigone";
			for (int j = 0; j < big_back_matrix.Dim(0); ++j)
			{
				for (int k = 0; k < big_back_matrix.Dim(1); ++k)
				{
					(A)(j,k) = big_back_matrix.At(j,k);
				}
			}

			Eigen::VectorXd b(big_forward_vec.size());
			for (int j = 0; j < big_forward_vec.size(); ++j)
			{
				b(j) = big_forward_vec.at(j);
			}
			cout << "before solve";
			cout << "\n";
			Eigen::MatrixXd newA = A.transpose() * (A);
			Eigen::MatrixXd newb = A.transpose() * b;
			Eigen::VectorXd res = newA.ldlt().solve(newb);
			cout << "after solve";
			cout << "\n";		
			// solve the linear system here
			Tensor result_tensor(result_dims);
			for (int j = 0; j < res.size(); ++j)
			{
				result_tensor.Set(j, res(j));
			}
			
			Tensor& curr_result = node_list->at(i_cliques[0])->GetTransformTensor().GetTensor();

		//	assert(Tensor::Diff(result_tensor, curr_result) < .0001);
			for (int j = 0; j < i_cliques.size(); ++j)
			{
				node_list->at(i_cliques[j])->SetTransformTensor(result_tensor);
			}
		}
	}
}



void TensorJTree::LearnSpectralParameters()
{
	for (int i = 0; i < node_list->size(); ++i)
	{
		node_list->at(i)->MarkMultModes();
		FindObservationPartitions(i);
	}

	for (int i = 0; i < node_list->size(); ++i)
	{
		node_list->at(i)->MarkObservableRepresentation();
	}

	for (int i = 0; i < node_list->size(); ++i)
	{
		cout << i;
		cout << "\n";
		node_list->at(i)->ComputeObservableRepresentation();
	}

	for (int i = 0; i < node_list->size(); ++i)
	{
		node_list->at(i)->FinalizeObservableRepresentation();
	}


	for (int i = 0; i < node_list->size(); ++i)
	{
		cout << i;
		cout << "\n";
		vector<double> big_forward_vec;
	//	big_forward_vec.reserve(10000);
		Matrix big_back_matrix;
	//	big_back_matrix.Reserve(10000);  	

		bool invModesExist = node_list->at(i)->InvModesExists();
		
		if (invModesExist)
		{
			vector<Tensor*>& forward_list = node_list->at(i)->GetForwardList();
			vector<Tensor*>& extra_list = node_list->at(i)->GetExtraList();

			Tensor& X = node_list->at(i)->GetTransformTensor().GetTensor();
			vector<int> result_dims = X.Dims();
			assert(forward_list.size() == extra_list.size());
		
			vector<int> inv_group_modes =  node_list->at(i)->GetInvModeGroup();
			vector<int> inv_modes = node_list->at(i)->GetInvModes();

			for (int j = 0; j < forward_list.size(); ++j)
			{
				Tensor& forward_tensor = *(forward_list[j]);
				Tensor& back_tensor = *(extra_list[j]);
	
				if (j == 0)
				{
					Tensor::CreateLinearSystem(big_forward_vec, big_back_matrix, X, back_tensor, forward_tensor, inv_group_modes, inv_modes); 
				}
				else
				{
					vector<double> forward_vec;
					Matrix back_matrix;

					Tensor::CreateLinearSystem(forward_vec, back_matrix, X, back_tensor, forward_tensor, inv_group_modes, inv_modes); 

					VectorPlus::Concat(big_forward_vec, forward_vec);
					big_back_matrix.MatrixConcat(back_matrix);
				}
			}

			Eigen::MatrixXd A(big_back_matrix.Dim(0), big_back_matrix.Dim(1));
			for (int j = 0; j < big_back_matrix.Dim(0); ++j)
			{
				for (int k = 0; k < big_back_matrix.Dim(1); ++k)
				{
					A(j,k) = big_back_matrix.At(j,k);
				}
			}

			Eigen::VectorXd b(big_forward_vec.size());
			for (int j = 0; j < big_forward_vec.size(); ++j)
			{
				b(j) = big_forward_vec.at(j);
			}
			//cout << "yee";
			Eigen::MatrixXd newA = A.transpose() * A;
			Eigen::MatrixXd newb = A.transpose() * b;
			Eigen::VectorXd res = newA.ldlt().solve(newb);
			//cout << "yeeldlt solved";
			//Eigen::VectorXd res = A.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
			//cout << "yeesvdsolved";
			// solve the linear system here
			Tensor result_tensor(result_dims);
			for (int j = 0; j < res.size(); ++j)
			{
				result_tensor.Set(j, res(j));
			}
			
			Tensor& curr_result = node_list->at(i)->GetTransformTensor().GetTensor();
		//	assert(Tensor::Diff(result_tensor, curr_result) < .0001);
			node_list->at(i)->SetTransformTensor(result_tensor);
		}
	}

	// average by computing template and see what happens
	// hack it so we change the template
/*	template_to_cliques = new vector<vector<int>*>();
	for (int i = 0; i < NumCliques(); ++i)
	{
		template_to_cliques->push_back(new vector<int>(1,i));
	} */

	/*for (int i = 0; i < template_to_cliques->size(); ++i)
	{
		vector<int>& i_cliques = *(template_to_cliques->at(i));

		vector<double> big_forward_vec;
		Matrix big_back_matrix;
		vector<int> result_dims;

		bool invModesExist = node_list->at(i_cliques[0])->InvModesExists();
		if (invModesExist)
		{
			for (int j = 0; j < i_cliques.size(); ++j)
			{
				if (node_list->at(i_cliques[0])->is_invalid())
					continue;

				Tensor& forward_tensor = node_list->at(i_cliques[j])->GetForwardTensor();
				Tensor& back_tensor = node_list->at(i_cliques[j])->GetUExtraTensor();
				Tensor& X = node_list->at(i_cliques[j])->GetTransformTensor().GetTensor();
			
				vector<int>& inv_group_modes =  node_list->at(i_cliques[j])->GetInvModeGroup();
				vector<int>& inv_modes = node_list->at(i_cliques[j])->GetInvModes();
				if (j == 0)
				{
					Tensor::CreateLinearSystem(big_forward_vec, big_back_matrix, X, back_tensor, forward_tensor, inv_group_modes, inv_modes); 
					result_dims = X.Dims();
				}
				else
				{
					vector<double> forward_vec;
					Matrix back_matrix;

					Tensor::CreateLinearSystem(forward_vec, back_matrix, X, back_tensor, forward_tensor, inv_group_modes, inv_modes); 

					VectorPlus::Concat(big_forward_vec, forward_vec);
					big_back_matrix.MatrixConcat(back_matrix);
				}
			}

			Eigen::MatrixXd A(big_back_matrix.Dim(0), big_back_matrix.Dim(1));
			for (int i = 0; i < big_back_matrix.Dim(0); ++i)
			{
				for (int j = 0; j < big_back_matrix.Dim(1); ++j)
				{
					A(i,j) = big_back_matrix.At(i,j);
				}
			}

			Eigen::VectorXd b(big_forward_vec.size());
			for (int i = 0; i < big_forward_vec.size(); ++i)
			{
				b(i) = big_forward_vec.at(i);
			}

			Eigen::VectorXd res = A.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
			// solve the linear system here
			Tensor result_tensor(result_dims);
			for (int i = 0; i < res.size(); ++i)
			{
				result_tensor.Set(i, res(i));
			}
			
			Tensor& curr_result = node_list->at(i_cliques[0])->GetTransformTensor().GetTensor();

		//	assert(Tensor::Diff(result_tensor, curr_result) < .0001);
			for (int j = 0; j < i_cliques.size(); ++j)
			{
				node_list->at(i_cliques[j])->SetTransformTensor(result_tensor);
			}
		}
	} */

	
/*	for (int i = 0; i < template_to_cliques->size(); ++i)
	{
		vector<int>& i_cliques = *(template_to_cliques->at(i));
		TensorCPT* avg_tensor =  new TensorCPT(node_list->at(i_cliques[0])->GetTransformTensor());
		for (int j = 1; j < i_cliques.size(); ++j)
		{
			TensorCPT::Add(*avg_tensor, node_list->at(i_cliques[j])->GetTransformTensor());
		}
		avg_tensor->Divide(i_cliques.size());

		for (int j = 0; j < i_cliques.size(); ++j)
		{
			node_list->at(i_cliques[j])->SetTransformTensor(avg_tensor->GetTensor());
		}
		delete(avg_tensor);
	} */
}

void TensorJTree::FindObservationPartitions(int clique_index)
{
	map<vector<int>, vector<int>*>& unique_child_seps = GetChildSeps(clique_index);
	map<vector<int>, vector<int>*>::iterator iter;

//	map<vector<int>, vector<int>*> local_obs_vars;
	vector<int> total_union;
	for (iter = unique_child_seps.begin(); iter != unique_child_seps.end(); ++iter)
	{
		vector<int> sep = iter->first;
		vector<int> sep_nodes = *iter->second;
	//	vector<int>* union_obs_vars = new vector<int>();
		vector<int> union_obs_vars;
		for (int i = 0; i < sep_nodes.size(); ++i)
		{
			vector<int> descendants;
			FindObsDescendants(descendants, sep_nodes.at(i));
			node_list->at(clique_index)->AddObsPartition(sep, descendants, i);
			VectorPlus::Union(union_obs_vars, descendants);
		}
	//	local_obs_vars[sep] = union_obs_vars;
	
		//if (clique_index != root)
	//	VectorPlus::Union(total_union, *union_obs_vars);
		vector<int> other_nodes;
		vector<int>& sorted_vars = *(dist_vars->at(clique_index));
		VectorPlus::SetDiff(other_nodes, sorted_vars, union_obs_vars);

		//vector<int> reverse_other_nodes;
		//VectorPlus::Reverse(reverse_other_nodes, other_nodes);
		node_list->at(clique_index)->AddObsPartition(sep, other_nodes, sep_nodes.size());
	}

/*	vector<int> extra(2, -1);
	vector<int>* extra_junk = new vector<int>();
	VectorPlus::SetDiff(*extra_junk, *all_obs_vars, total_union);
	local_obs_vars[extra] = extra_junk;

	for (iter = unique_child_seps.begin(); iter != unique_child_seps.end(); ++iter)
	{
		vector<int> sep = iter->first;
		vector<int> sep_nodes = *iter->second;

		vector<int> extra_obs;
		int curr_index = 0;
		int done_flag = false;
		while(done_flag == false)
		{
			done_flag = true;
			map<vector<int>, vector<int>*>::iterator inner_iter;
			for (inner_iter = local_obs_vars.begin(); inner_iter != local_obs_vars.end(); ++inner_iter)
			{
				vector<int> fun = inner_iter->first;
				if (VectorPlus::Equals(sep, fun))
					continue;

				vector<int>& all_other_obs = *inner_iter->second;
				if (curr_index < all_other_obs.size())
				{
					extra_obs.push_back(all_other_obs[curr_index]);
					done_flag = false;
				}
			}
			curr_index++;
		}
		vector<int> reverse_obs;
		VectorPlus::Reverse(reverse_obs, extra_obs);
		node_list->at(clique_index)->AddObsPartition(sep, reverse_obs, sep_nodes.size());
	}

	map<vector<int>, vector<int>*>::iterator inner_iter;
	for (inner_iter = local_obs_vars.begin(); inner_iter != local_obs_vars.end(); ++inner_iter)
	{
		delete(inner_iter->second);
	} */
}

// Finds observed descendant variables of a given clique node
void TensorJTree::FindObsDescendants(vector<int>& desc_vec, int clique_index)
{
	vector<int> non_unique_vec;
	queue<int> active_queue;
	active_queue.push(clique_index);
	while(!active_queue.empty())
	{
		int node = active_queue.front();
		active_queue.pop();
		vector<int>& children = GetChildren(node);
		for (int i = 0; i < children.size(); ++i)
		{
			active_queue.push(children.at(i));
		}


		vector<int>& c_vars = GetCVars(node);
		for (int i = 0; i < c_vars.size(); ++i)
		{
			if (bNet->is_observed(c_vars[i]))
				non_unique_vec.push_back(c_vars[i]);
		}
		
	}

	VectorPlus::Unique(desc_vec, non_unique_vec);
}


bool TensorJTree::IsLeaf(int clique_index)
{
	vector<int>& children = GetChildren(clique_index);
	if (children.size() == 0)
		return false;
	else
		return true;
}

void TensorJTree::ComputeEliminationOrder()
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


double TensorJTree::ComputeMarginalEmpiricalProbability(vector<int> evidence_vars, vector<int> evidence_vals)
{
	map<int, MultiVector<TensorCPT*>*> message_map;

	for (int i = 0; i < NumCliques(); ++i)
	{
		message_map[i] = new MultiVector<TensorCPT*>();
	}

	for (int i = 0; i < NumCliques(); ++i)
	{
		int clique_index = elimination_order->at(i);

		vector<int> evidence_i_vars;
		vector<int> evidence_i_vals;

		VectorPlus::Intersect(evidence_i_vars, GetCVars(clique_index), evidence_vars);
		VectorPlus::MatchSub(evidence_i_vals, evidence_vals, evidence_vars, evidence_i_vars);

		TensorCPT* outgoing_message = new TensorCPT();
		TensorCPT& old_cpt = GetTransformTensor(clique_index);

		old_cpt.Slice(*outgoing_message, evidence_i_vars, evidence_i_vals); 
		int num_seps = GetNumUniqueChildSeps(clique_index);

		MultiVector<TensorCPT*>* message_multi_vec = message_map[clique_index];

		for (int j = num_seps - 1; j >= 0; j--)
		{
	
			vector<TensorCPT*>& message_vec = message_multi_vec->GetVec(j);
			TensorCPT* curr_message = new TensorCPT(*message_vec[0]);
			for (int k = 1; k < message_vec.size(); ++k)
			{
				TensorCPT* new_t_cpt = new TensorCPT();
				TensorCPT::Multiply(*new_t_cpt, *curr_message, *(message_vec[k]));
				delete(curr_message);
				curr_message = new_t_cpt;
			}

			TensorCPT* new_outgoing = new TensorCPT();
			TensorCPT::Multiply(*new_outgoing, *outgoing_message, *curr_message);
			delete(outgoing_message);
			delete(curr_message);
			outgoing_message = new_outgoing;
		}

		if (clique_index == root)
		{
			vector<int>& dims = outgoing_message->Dims();
			assert(dims.size() == 0);
			return outgoing_message->At(0);
		}
		else
		{

			MultiVector<TensorCPT*>& out_message_multi_vec = *(message_map[GetParent(clique_index)]);
			int s_index = GetSIndex(clique_index);
			int same_sibling_index = GetSameSiblingIndex(clique_index);

			out_message_multi_vec.Set(s_index, same_sibling_index, outgoing_message); 
	
		}
	}

	for (int i = 0; i < NumCliques(); ++i)
	{
		delete(message_map[i]);
	}

	assert(0);
	return 0;
}

bool TensorJTree::IsFirstInTemplate(int clique_index)
{
	if (GetFirstInTemplate(clique_index) == clique_index)
		return true;
	else
		return false;
}

int TensorJTree::GetFirstInTemplate(int clique_index)
{
	int template_index = clique_to_template->at(clique_index);
	//return 1;
	if (template_index == 0)
		return template_to_cliques->at(template_index)->at(0);
	else
	{
		int prev_size = template_to_cliques->at(template_index - 1)->size();
		if (prev_size > 1)
			return template_to_cliques->at(template_index - 1)->at(0);
		else
			return template_to_cliques->at(template_index)->at(0);
	}
}

int TensorJTree::GetLastInTemplate(int clique_index)
{
	int template_index = clique_to_template->at(clique_index);
	//return 1;
	int index = template_to_cliques->at(template_index)->size() - 1;
	return template_to_cliques->at(template_index)->at(index);

}


void TensorJTree::ComputeDistances()
{
	for (int i = 0; i < NumCliques(); ++i)
	{
		queue<int> dist_queue;
		dist_queue.push(i);
		int curr_dist = 0;
		set<int> covered_nodes;
		covered_nodes.insert(i);
		while(dist_queue.size() > 0)
		{
			int node = dist_queue.front();

			if (node != i)
			{
				vector<int>& r_vars = GetRVars(node);
				for (int j = 0; j < r_vars.size(); ++j)
				{
					if (bNet->is_observed(r_vars[j]))
						dist_vars->at(i)->push_back(r_vars[j]);
				}
			}
			dist_queue.pop();
			int parent = GetParent(node);
			if (parent >= 0 && covered_nodes.count(parent) == 0)
			{
				dist_queue.push(parent);
				covered_nodes.insert(parent);
			}

			vector<int>& children = GetChildren(node);
			for (int j = 0; j < children.size(); ++j)
			{
				if (covered_nodes.count(children[j]) == 0)
				{
					dist_queue.push(children[j]);
					covered_nodes.insert(children[j]);
				}
			}

			curr_dist++;
		}
	}
	vector<int> extra_obs;
	VectorPlus::SetDiff(extra_obs, *all_obs_vars, *(dist_vars->at(0)));
	for (int i = 0; i < NumCliques(); ++i)
	{
		vector<int>& r_vars = GetRVars(i);
		vector<int> i_extra_obs;
		VectorPlus::SetDiff(i_extra_obs, extra_obs, r_vars);
		VectorPlus::Concat(*(dist_vars->at(i)), i_extra_obs);
	}
}

TensorCPT& TensorJTree::GetTransformTensor(int clique_index)
{ 
	return node_list->at(clique_index)->GetTransformTensor(); 
}

int TensorJTree::GetOtherNode(int node_id, vector<int>& mode_group) { return node_list->at(node_id)->GetOtherNode(mode_group); }

int TensorJTree::GetNumDescendantPartitions(int node_id, vector<int>& svars) {return node_list->at(node_id)->GetNumDescendantPartitions(svars); }

Tensor& TensorJTree::GetU(int clique_index) { return node_list->at(clique_index)->GetU(); }
int TensorJTree::GetSIndex(int clique_index) { return node_list->at(clique_index)->GetSIndex(); }
int TensorJTree::GetSameSiblingIndex(int clique_index) { return node_list->at(clique_index)->GetSameSiblingIndex(); }
int TensorJTree::GetChildSepIndex(int clique_index, vector<int>& S_vars) { return node_list->at(clique_index)->GetChildSepIndex(S_vars); }
int TensorJTree::GetNumUniqueChildSeps(int clique_index) {return node_list->at(clique_index)->GetNumUniqueChildSeps(); }
vector<int>& TensorJTree::GetChildren(int clique_index) { return node_list->at(clique_index)->GetChildren(); }
int TensorJTree::GetParent(int clique_index) {return node_list->at(clique_index)->GetParent(); }
vector<int>& TensorJTree::GetCVars(int clique_index) {return node_list->at(clique_index)->GetCVars(); }
vector<int>& TensorJTree::GetRVars(int clique_index) {return node_list->at(clique_index)->GetRVars(); }
vector<int>& TensorJTree::GetSVars(int clique_index) {return node_list->at(clique_index)->GetSVars(); }
vector<int>& TensorJTree::GetObsVector(int node, vector<int>& modes) {return node_list->at(node)->GetObsVector(modes); }
vector<int>& TensorJTree::GetObsVector(int node, int other_node) {return node_list->at(node)->GetObsVector(other_node); }
map<vector<int>, vector<int>*>& TensorJTree::GetChildSeps(int clique_index) {return node_list->at(clique_index)->GetChildSeps(); }
vector<int>& TensorJTree::GetDescendants(int clique_index, vector<int>& mode_group, int index) { return node_list->at(clique_index)->GetDescendants(mode_group, index); }
