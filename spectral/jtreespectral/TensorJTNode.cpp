# pragma once

#include "TensorJTNode.hpp"
#include "VectorPlus.hpp"
#include <assert.h>
#include "MultiVector.hpp"
#include <boost/shared_ptr.hpp>
#include <limits>

TensorJTNode::TensorJTNode()
{
	jtree = NULL;

	parent_id = -1;
	node_id = -1;
	max_obs_per_mode = -1;
	children = NULL;
	same_siblings = NULL;

	max_runs = -1;

	S_index = -1;
	same_sibling_index = -1;

	S_vars = NULL;
	R_vars = NULL;
	C_vars = NULL;
	modes_to_vars = new vector<int>();
	inv_mode_group = new vector<int>();

	child_seps_to_children = new map<vector<int>, vector<int>*>();
	child_seps_to_obs_vars = new map<vector<int>, vector<vector<int>*>*>();

	mode_groups_to_nodes = new map<vector<int>, int>();
	mode_groups_to_types = new map<vector<int>, int>();
	mode_groups_to_obs = new map<vector<int>, vector<int>*>();

	extra_tensor_list = new vector<Tensor*>();
	forward_tensor_list = new vector<Tensor*>();

	//mode_groups_to_U = new map<vector<int>, Tensor*>();
	U = NULL;
	inv_mode_group = new vector<int>();

	extra_tensor = NULL;
	U_extra_tensor = NULL;
	transform_tensor = NULL;
	B_transform_tensor = NULL;

	invalid_flag = false;
}

TensorJTNode::TensorJTNode(TensorJTree* jtree, int node_id, int max_obs_per_mode, int max_runs)
{
	this->jtree = jtree;

	parent_id = -1;
	this->max_obs_per_mode = max_obs_per_mode;
	this->node_id = node_id;
	children = NULL;
	same_siblings = NULL;

	this->max_runs = max_runs;

	S_index = -1;
	same_sibling_index = -1;

	S_vars = NULL;
	R_vars = NULL;
	C_vars = NULL;
	modes_to_vars = new vector<int>();
	inv_mode_group = new vector<int>();
	child_seps_to_children = new map<vector<int>, vector<int>*>();
	child_seps_to_obs_vars = new map<vector<int>, vector<vector<int>*>*>();

	mode_groups_to_nodes = new map<vector<int>, int>();
	mode_groups_to_types = new map<vector<int>, int>();
	mode_groups_to_obs = new map<vector<int>, vector<int>*>();

	extra_tensor_list = new vector<Tensor*>();
	forward_tensor_list = new vector<Tensor*>();
//	mode_groups_to_U = new map<vector<int>, Tensor*>();

	inv_mode_group = new vector<int>();
	U = NULL;
	extra_tensor = NULL;
	U_extra_tensor = NULL;
	transform_tensor = NULL;
	B_transform_tensor = NULL;

	invalid_flag = false;
}

TensorJTNode::~TensorJTNode()
{
	delete(children);
	delete(same_siblings);

	delete(S_vars);
	delete(R_vars);
	delete(C_vars);
	delete(modes_to_vars);

	map<vector<int>, vector<int>*>::iterator iter;
	for (iter =  child_seps_to_children->begin(); iter != child_seps_to_children->end(); ++iter)
	{
		delete(iter->second);
	}
	delete(child_seps_to_children);

	for (iter =  mode_groups_to_obs->begin(); iter != mode_groups_to_obs->end(); ++iter)
	{
		delete(iter->second);
	}
	delete(mode_groups_to_obs);


	map<vector<int>, vector<vector<int>*>*>::iterator d_iter;
	for (d_iter = child_seps_to_obs_vars->begin(); d_iter != child_seps_to_obs_vars->end(); ++d_iter)
	{
		vector<vector<int>*>* vec = d_iter->second;
		for (int i = 0; i < vec->size(); ++i)
		{
			delete(vec->at(i));
		}
		delete(vec);
	}
	delete(child_seps_to_obs_vars);


	delete(mode_groups_to_nodes);
	delete(mode_groups_to_types);

/*	map<vector<int>, Tensor*>::iterator U_iter;
	for (U_iter =  mode_groups_to_U->begin(); U_iter != mode_groups_to_U->end(); ++U_iter)
	{
		delete(U_iter->second);
	}
	delete(mode_groups_to_U); */
	delete(transform_tensor);

}

void TensorJTNode::UniqueChildSepAdd(vector<int>& child_sep, int child_id)
{
	if (child_seps_to_children->count(child_sep) == 0)
	{
		(*child_seps_to_children)[child_sep] = new vector<int>();
	}
	(*child_seps_to_children)[child_sep]->push_back(child_id);
}

void TensorJTNode::AddObsPartition(vector<int>& sep, vector<int>& descendants, int index)
{
	if (child_seps_to_obs_vars->count(sep) == 0)
	{
		(*child_seps_to_obs_vars)[sep] = new vector<vector<int>*>((*child_seps_to_children)[sep]->size() + 1, 0);
	}
	vector<int>* desc_copy = new vector<int>(descendants);
	((*child_seps_to_obs_vars)[sep])->at(index) = desc_copy;
}

void TensorJTNode::MarkMultModes()
{
	int curr_mode = 0;

	if (parent_id == -1)
		same_siblings = new vector<int>();

	if (parent_id != -1)
	{
		map<vector<int>, vector<int>*>& sibling_seps = jtree->GetChildSeps(parent_id);
		
		same_siblings = new vector<int>(*sibling_seps[*S_vars]);
		
		int node_index = VectorPlus::Find(*same_siblings, node_id);	
		
		if (same_siblings->size() == 1)
		{
			vector<int> mode_group;

			for (int i = 0; i < S_vars->size(); ++i)
			{
				modes_to_vars->push_back((*S_vars)[i]);
				mode_group.push_back(curr_mode++);
			}
			
			(*mode_groups_to_nodes)[mode_group] =  parent_id;
			(*mode_groups_to_types)[mode_group] =  INVERSE;
			VectorPlus::Copy(*inv_mode_group, mode_group);
		}
		else if (node_index == 0)
		{
			vector<int> mode_group1;
			for (int i = 0; i < S_vars->size(); ++i)
			{
				modes_to_vars->push_back(S_vars->at(i));
				mode_group1.push_back(curr_mode++);
			}
			(*mode_groups_to_nodes)[mode_group1] =  parent_id;
			(*mode_groups_to_types)[mode_group1] =  INVERSE;
			VectorPlus::Copy(*inv_mode_group, mode_group1);

			vector<int> mode_group2;
			for (int i = 0; i < S_vars->size(); ++i)
			{
				modes_to_vars->push_back(S_vars->at(i));
				mode_group2.push_back(curr_mode++);
			}

			(*mode_groups_to_nodes)[mode_group2] = same_siblings->at(node_index + 1);
			(*mode_groups_to_types)[mode_group2] = NORMAL;
		}
		else if (node_index == same_siblings->size() - 1)
		{
			vector<int> mode_group;
			for (int i = 0; i < S_vars->size(); ++i)
			{
				modes_to_vars->push_back(S_vars->at(i));
				mode_group.push_back(curr_mode++);
			}

			(*mode_groups_to_nodes)[mode_group] = same_siblings->at(node_index - 1);
			(*mode_groups_to_types)[mode_group] = INVERSE;
			VectorPlus::Copy(*inv_mode_group, mode_group);

		}
		else
		{
			vector<int> mode_group1;
			for (int i = 0; i < S_vars->size(); ++i)
			{
				modes_to_vars->push_back(S_vars->at(i));
				mode_group1.push_back(curr_mode++);
			}
			(*mode_groups_to_nodes)[mode_group1] = same_siblings->at(node_index - 1);
			(*mode_groups_to_types)[mode_group1] = INVERSE;
			VectorPlus::Copy(*inv_mode_group, mode_group1);

			vector<int> mode_group2;
			for (int i = 0; i < S_vars->size(); ++i)
			{
				modes_to_vars->push_back(S_vars->at(i));
				mode_group2.push_back(curr_mode++);
			}

			(*mode_groups_to_nodes)[mode_group2] = same_siblings->at(node_index + 1);
			(*mode_groups_to_types)[mode_group2] = NORMAL;
		}
	}

	// front part, applies to all non leaves in junction tree
	if (children->size() > 0)
	{
		map<vector<int>, vector<int>*>::iterator iter;

		for (iter = child_seps_to_children->begin(); iter != child_seps_to_children->end(); ++iter)
		{
			vector<int> sep = iter->first;
			vector<int> sep_nodes = *iter->second;

			vector<int> mode_group;
			for (int i = 0; i < sep.size(); ++i)
			{
				modes_to_vars->push_back(sep[i]);
				mode_group.push_back(curr_mode++);
			}

			(*mode_groups_to_nodes)[mode_group] = sep_nodes[0];
			(*mode_groups_to_types)[mode_group] = NORMAL;

		}

	}
	//else
	//{
		for (int i = 0; i < R_vars->size(); ++i)
		{  
			if (!VectorPlus::Contains(*modes_to_vars, R_vars->at(i)))
			{
				modes_to_vars->push_back(R_vars->at(i));
			}
		}
	//}
}

void TensorJTNode::MarkObservableRepresentation()
{	
	map<vector<int>, int>::iterator iter;

	for (iter = mode_groups_to_nodes->begin(); iter != mode_groups_to_nodes->end(); ++iter)
	{
		vector<int> mode_group =  iter->first;

		vector<int> sep_vars;
		VectorPlus::Subset(sep_vars, *modes_to_vars, mode_group);
		int mult_node = iter->second;
		int mult_type = (*mode_groups_to_types)[mode_group];

		int child_index = VectorPlus::Find(*children, mult_node);

		if (child_index >= 0) // child case
		{
			vector<int>& same_children = *((*child_seps_to_children)[sep_vars]);
			int same_child_index = VectorPlus::Find(same_children, mult_node);

			int numPartitions =(*child_seps_to_obs_vars)[sep_vars]->size();
		//	if(parent_id == -1)
		//		numPartitions = numPartitions - 1;

			int num_3_groups = (int)(numPartitions / 3);
			if (num_3_groups == 0)
				num_3_groups = 1;

			int num_obs_per_partition = numeric_limits<int>::max();
			for (int i = 0; i < numPartitions; ++i)
			{
				vector<int>& descendant_obs = *((*child_seps_to_obs_vars)[sep_vars]->at(i));
				if (descendant_obs.size() < num_obs_per_partition)
				{
					assert(descendant_obs.size() > 0);
					num_obs_per_partition = descendant_obs.size();
				}
			}


			int num_obs_available = num_3_groups * num_obs_per_partition;
			int num_obs_available_per_mode = (int)(num_obs_available) / mode_group.size();
			int num_obs_needed = min(num_obs_available_per_mode, max_obs_per_mode) * mode_group.size();
			assert(same_child_index == 0);
			
			
			(*mode_groups_to_obs)[mode_group] = new vector<int>();
			int curr_same_child_index = 0;
			int curr_obs_index = 0;
			for (int i = 0; i < num_obs_needed; ++i)
			{
				vector<int>& descendant_obs = *((*child_seps_to_obs_vars)[sep_vars]->at(curr_same_child_index));
				(*mode_groups_to_obs)[mode_group]->push_back(descendant_obs[curr_obs_index]);

				if (numPartitions >= 3)
				{
					curr_same_child_index += 3;
					if (curr_same_child_index >= numPartitions)
					{
						curr_same_child_index = curr_same_child_index % (numPartitions);
						curr_obs_index++;
					}
				}
				else
				{
					curr_obs_index++;
				}
			}
		}
		else // parent,sibling case
		{
			int mod_obs_index = -1;
			int mod_sibling_size = (int) min((double)same_siblings->size(), 3.0);
			int index = VectorPlus::Find(*same_siblings, node_id);
			int mod_index = index % mod_sibling_size;

			if (mult_type == INVERSE)
			{
				if (NumSameSiblings() == 1)
					mod_obs_index = 1;
				else if (mod_index ==  0)
					mod_obs_index = mod_sibling_size - 1;
				else
					mod_obs_index = mod_index - 1;
			}
			else
			{
				mod_obs_index = (mod_index + 1) % mod_sibling_size;
			}
	
			int curr_obs_index = 0;
			int curr_sibling_index = mod_obs_index;
			int numPartitions = jtree->GetNumDescendantPartitions(parent_id, *S_vars);
	//		if(jtree->GetDescendants(parent_id, *S_vars, numPartitions-1).size() == 0)
	//			numPartitions = numPartitions - 1;

			int num_3_groups = (int)((numPartitions - 1) / 3);
			if (num_3_groups == 0)
				num_3_groups = 1;

			int num_obs_per_partition = numeric_limits<int>::max();
			for (int i = 0; i < numPartitions; ++i)
			{
				vector<int>& descendant_obs = jtree->GetDescendants(parent_id, *S_vars, i);
				if (descendant_obs.size() < num_obs_per_partition && descendant_obs.size() != 0)
					num_obs_per_partition = descendant_obs.size();
			}
			
			int num_obs_available = num_3_groups * num_obs_per_partition;
			if (num_obs_available < mode_group.size())
			{
				invalid_flag = true;
				assert(0);
				continue;
			}

			int num_obs_available_per_mode = (int)(num_obs_available) / mode_group.size();

			int num_obs_needed = 0;
			if (mult_type == NORMAL)
			{
				num_obs_needed = min(num_obs_available_per_mode, max_obs_per_mode) * mode_group.size();

				assert(curr_sibling_index != numPartitions - 1);

				(*mode_groups_to_obs)[mode_group] = new vector<int>();
				for (int i = 0; i < num_obs_needed; ++i)
				{
					assert(curr_sibling_index != numPartitions - 1);
					vector<int>& descendant_obs = jtree->GetDescendants(parent_id, *S_vars, curr_sibling_index);

					(*mode_groups_to_obs)[mode_group]->push_back(descendant_obs[curr_obs_index]);
				
					if (numPartitions - 1 >= 3)
					{
						curr_sibling_index += 3;
						if (curr_sibling_index >= (numPartitions - 1))
						{
							curr_obs_index++;
						}
						curr_sibling_index = curr_sibling_index % (numPartitions - 1);
					}
					else
					{
						curr_obs_index++;
					}
				}

			}
			else
			{
				int curr_index = curr_sibling_index;
				(*mode_groups_to_obs)[mode_group] = new vector<int>();

				if (!(curr_index == numPartitions - 1))
				{
					while(true)
					{
						num_obs_needed += jtree->GetDescendants(parent_id, *S_vars, curr_index).size();
						curr_index += 3;
						if (curr_index >= numPartitions - 1)
							break;
					}

					for (int i = 0; i < num_obs_needed; ++i)
					{
						assert(curr_sibling_index != numPartitions - 1);
						vector<int>& descendant_obs = jtree->GetDescendants(parent_id, *S_vars, curr_sibling_index);

						if (descendant_obs.size() > curr_obs_index)
							(*mode_groups_to_obs)[mode_group]->push_back(descendant_obs[curr_obs_index]);
				
						if (numPartitions - 1 >= 3)
						{
							curr_sibling_index += 3;
							if (curr_sibling_index >= (numPartitions - 1))
							{
								curr_obs_index++;
							}
							curr_sibling_index = curr_sibling_index % (numPartitions - 1);
						}
						else
						{
							curr_obs_index++;
						}
					}

				}

				vector<int>& descendant_obs = jtree->GetDescendants(parent_id, *S_vars, numPartitions - 1);
				for (int i = 0; i < descendant_obs.size(); ++i)
				{
					(*mode_groups_to_obs)[mode_group]->push_back(descendant_obs[i]);
				}
			}
		}
	}
}

void TensorJTNode::ComputeObservableRepresentation()
{
	if (invalid_flag)
		return;

	int num_runs = 0;
	int num_inv_obs = 0;
	if (inv_mode_group->size() == 0)
	{
		num_runs = 1;
	}
	else
	{
		int other_node = (*mode_groups_to_nodes)[*inv_mode_group];
		vector<int>& inv_obs = jtree->GetObsVector(other_node, node_id);
		vector<int>& extra_obs = *((*mode_groups_to_obs)[*inv_mode_group]);
		num_inv_obs = inv_obs.size();
		num_runs = (int)(extra_obs.size() / inv_obs.size());
		num_runs = (int)min((double)num_runs, (double)max_runs);
	}


	int r_curr_index = 0;
	vector<int>* all_extra_modes = NULL;
	for (int r = 0; r < num_runs; ++r)
	{
		map<vector<int>, vector<int>*>::iterator iter;

		MultiVector<int> forward_vars; // multi vector

		vector<int> remaining_modes;
		vector<int> new_modes_to_vars(*modes_to_vars);
		VectorPlus::Seq(remaining_modes, 0, 1, modes_to_vars->size());

		for (iter = mode_groups_to_obs->begin(); iter != mode_groups_to_obs->end(); ++iter)
		{
			vector<int> mode_group =  iter->first;

			vector<int> temp_modes;
			VectorPlus::SetDiff(temp_modes, remaining_modes, mode_group);
			remaining_modes = temp_modes;

			vector<int>& obs_vars = *(iter->second);
			if ((*mode_groups_to_types)[mode_group] == INVERSE)
			{
				int num_obs_per_mode = num_inv_obs / mode_group.size();
				int curr_index = r_curr_index;
				for (int i = 0; i < mode_group.size(); ++i)
				{
					for (int j = 0; j < num_obs_per_mode; ++j)
					{
						forward_vars.Set(mode_group[i], j, obs_vars[curr_index++]);
					}
				}
			}
			else
			{
				int num_obs_per_mode = obs_vars.size() / mode_group.size();
				int curr_index = 0;
				for (int i = 0; i < mode_group.size(); ++i)
				{
					for (int j = 0; j < num_obs_per_mode; ++j)
					{
						forward_vars.Set(mode_group[i], j, obs_vars[curr_index++]);
					}
				}
			}
			
		}

		for (int i = 0; i < remaining_modes.size(); ++i)
		{
			forward_vars.Set(remaining_modes[i], 0, modes_to_vars->at(remaining_modes[i]));
		}
	
		Tensor* curr_tensor = new Tensor();
		jtree->ComputeEmpiricalMultiProbTensor(*curr_tensor, forward_vars);
	//	B_transform_tensor = new TensorCPT();
	//	B_transform_tensor->Initialize(*curr_tensor, *modes_to_vars);

		forward_tensor_list->push_back(new Tensor(*curr_tensor));

		map<vector<int>, int>::iterator type_iter;

		for (type_iter = mode_groups_to_types->begin(); type_iter != mode_groups_to_types->end(); ++type_iter)
		{
			if (type_iter->second == INVERSE)
			{
				vector<int> mode_group = type_iter->first;
				int other_node = (*mode_groups_to_nodes)[mode_group];
				vector<int>& inv_obs = jtree->GetObsVector(other_node, node_id);
				vector<int>& extra_obs = *((*mode_groups_to_obs)[mode_group]);

		
				// fill up with some extra obs
				vector<int> r_extra_obs;
				VectorPlus::Subset(r_extra_obs, extra_obs, r_curr_index, r_curr_index + inv_obs.size());
				r_curr_index += inv_obs.size();
				assert(inv_obs.size() == r_extra_obs.size());
				// fill up extra tensor, first modes come from other, last modes from from extra
				vector<int> inv_modes;
				vector<int> extra_modes;

				int total_modes = 2*mode_group.size();
				VectorPlus::Seq(inv_modes, 0, 1, mode_group.size());
				VectorPlus::Seq(extra_modes, mode_group.size(), 1, total_modes);

				if (all_extra_modes == NULL)
				{
					all_extra_modes = new vector<int>(extra_modes);
				}

				MultiVector<int> backward_vars;
				int curr_index1 = 0;
				int curr_index2 = 0;

				int num_obs_per_mode = inv_obs.size() / mode_group.size();

				for (int i = 0; i < total_modes; ++i)
				{
					if (i < inv_modes.size())
					{
						for (int j = 0; j < num_obs_per_mode; ++j)
							backward_vars.Set(i, j, inv_obs[curr_index1++]);
	
					}
					else
					{
						for (int j = 0; j < num_obs_per_mode; ++j)
							backward_vars.Set(i, j, r_extra_obs[curr_index2++]);
					}
				}

				Tensor extra_tensor;
				jtree->ComputeEmpiricalMultiProbTensor(extra_tensor, backward_vars);
				Tensor inv_tensor;
				extra_tensor_list->push_back(new Tensor(extra_tensor));

				if (r == 0)
				{
					assert(U == NULL);
					U = new Tensor();

					bool succ = extra_tensor.Inverse(*U, inv_tensor, extra_modes, jtree->GetNumHiddenStates()); 
					Tensor* new_tensor = new Tensor();
					Tensor::Multiply2(*new_tensor, *curr_tensor, inv_tensor, mode_group, inv_modes);
					delete(curr_tensor);
					curr_tensor = new_tensor;	
					assert(succ == true);
				}
				
			}
		}

		if (r == 0)
		{
			transform_tensor = new TensorCPT();
			transform_tensor->Initialize(*curr_tensor, *modes_to_vars);
		}
		delete(curr_tensor);
	}

	if (U != NULL)
	{
		delete(U);
		U = NULL;
	}

	if (all_extra_modes != NULL)
	{
		assert(U == NULL);
		U = new Tensor();
		bool succ = Tensor::ComputeJointSVD(*U, *extra_tensor_list, *all_extra_modes, jtree->GetNumHiddenStates()); 
	} 
}

void TensorJTNode::FinalizeTemplateObservableRepresentation()
{
	if (invalid_flag)
	{
		assert(0);
		return;
	}

	map<vector<int>, int>::iterator type_iter;
	boost::shared_ptr<Tensor> curr_tensor(new Tensor(transform_tensor->GetTensor()));
	for (type_iter = mode_groups_to_types->begin(); type_iter != mode_groups_to_types->end(); ++type_iter)
	{
		if (type_iter->second == NORMAL)
		{
			vector<int> mode_group = type_iter->first;
			int other_node = (*mode_groups_to_nodes)[mode_group];
			int template_node = jtree->GetFirstInTemplate(other_node);
			//int template_other_node = jtree->GetOtherNode(template_node, mode_group);
			Tensor* new_tensor = new Tensor();
			vector<int> U_modes;
			VectorPlus::Seq(U_modes, 0, 1, mode_group.size());
			Tensor& U = jtree->GetU(template_node);
			Tensor::Multiply2(*new_tensor, *curr_tensor, U, mode_group, U_modes);
			curr_tensor.reset(new_tensor);

			for (int i = 0; i < forward_tensor_list->size(); ++i)
			{
				Tensor* new_i_tensor = new Tensor();
				Tensor::Multiply2(*new_i_tensor, *(forward_tensor_list->at(i)), U, mode_group, U_modes);
				delete(forward_tensor_list->at(i));
				forward_tensor_list->at(i) = new_i_tensor;
			}

		}
	}

	transform_tensor->SetTensor(*curr_tensor);


	if (InvModesExists())
	{
		for (int i = 0; i < extra_tensor_list->size(); ++i)
		{
			Tensor* new_tensor = new Tensor();
			vector<int> i_modes = GetInvModes();

			int template_node = jtree->GetFirstInTemplate(node_id);
			Tensor& U = jtree->GetU(template_node);

			Tensor::Multiply(*new_tensor, U, *(extra_tensor_list->at(i)), i_modes, i_modes);
			delete(extra_tensor_list->at(i));
			extra_tensor_list->at(i) = new_tensor;
		}
	}
}


void TensorJTNode::FinalizeObservableRepresentation()
{
	if (invalid_flag)
	{
		assert(0);
	}

	map<vector<int>, int>::iterator type_iter;
	boost::shared_ptr<Tensor> curr_tensor(new Tensor(transform_tensor->GetTensor()));
	for (type_iter = mode_groups_to_types->begin(); type_iter != mode_groups_to_types->end(); ++type_iter)
	{
		if (type_iter->second == NORMAL)
		{
			vector<int> mode_group = type_iter->first;
			int other_node = (*mode_groups_to_nodes)[mode_group];
		//	int template_node = jtree->GetFirstInTemplate(other_node);
		//	int template_other_node = jtree->GetOtherNode(template_node, mode_group);
			Tensor* new_tensor = new Tensor();
			vector<int> U_modes;
			VectorPlus::Seq(U_modes, 0, 1, mode_group.size());
			Tensor& U = jtree->GetU(other_node);
			Tensor::Multiply2(*new_tensor, *curr_tensor, U, mode_group, U_modes);
			curr_tensor.reset(new_tensor);

			for (int i = 0; i < forward_tensor_list->size(); ++i)
			{
				Tensor* new_i_tensor = new Tensor();
				Tensor::Multiply2(*new_i_tensor, *(forward_tensor_list->at(i)), U, mode_group, U_modes);
				delete(forward_tensor_list->at(i));
				forward_tensor_list->at(i) = new_i_tensor;
			}

		}
	}

	transform_tensor->SetTensor(*curr_tensor);

	// do it again for B_transform_tensor

	/*curr_tensor.reset(new Tensor(B_transform_tensor->GetTensor()));
	for (type_iter = mode_groups_to_types->begin(); type_iter != mode_groups_to_types->end(); ++type_iter)
	{
		if (type_iter->second == NORMAL)
		{
			vector<int> mode_group = type_iter->first;
			int other_node = (*mode_groups_to_nodes)[mode_group];
			int template_node = jtree->GetFirstInTemplate(other_node);
		//	int template_other_node = jtree->GetOtherNode(template_node, mode_group);
			Tensor* new_tensor = new Tensor();
			vector<int> U_modes;
			VectorPlus::Seq(U_modes, 0, 1, mode_group.size());
			Tensor& U = jtree->GetU(template_node);
			Tensor::Multiply2(*new_tensor, *curr_tensor, U, mode_group, U_modes);
			curr_tensor.reset(new_tensor);
		}
	}

	B_transform_tensor->SetTensor(*curr_tensor); */


	if (InvModesExists())
	{
		for (int i = 0; i < extra_tensor_list->size(); ++i)
		{
			Tensor* new_tensor = new Tensor();
			vector<int> i_modes = GetInvModes();
			Tensor::Multiply(*new_tensor, *U, *(extra_tensor_list->at(i)), i_modes, i_modes);
			delete(extra_tensor_list->at(i));
			extra_tensor_list->at(i) = new_tensor;
		}

		/*if (!jtree->IsFirstInTemplate(node_id))
		{
			if (U_extra_tensor != NULL)
				delete(U_extra_tensor);

			int template_node = jtree->GetFirstInTemplate(node_id);
			Tensor& U = jtree->GetU(template_node);
			U_extra_tensor = new Tensor();
			Tensor::Multiply(*U_extra_tensor, U, *extra_tensor, GetInvModes(), GetInvModes());
		}*/
	} 

}

// with few modes it is just faster to do sequential search
vector<int>& TensorJTNode::GetObsVector(int node)
{
	map<vector<int>, int>::iterator iter;
	for (iter = mode_groups_to_nodes->begin(); iter != mode_groups_to_nodes->end(); ++iter)
	{
		if (iter->second == node)
		{
			return *((*mode_groups_to_obs)[iter->first]);
		}
	}
	assert(0);
}

int TensorJTNode::GetSameSiblingIndex()
{
	if (same_sibling_index >= 0)
	{
		return same_sibling_index; 
	}
	else
	{
		same_sibling_index = VectorPlus::Find(*same_siblings, node_id);
		return same_sibling_index;
	}
}

Tensor& TensorJTNode::GetU()
{ 
	return *U;
}

int TensorJTNode::GetChildSepIndex(vector<int> sep_nodes)
{
	map<vector<int>, vector<int>*>::iterator iter;
	int index = 0;
	for (iter = child_seps_to_children->begin(); iter != child_seps_to_children->end(); ++iter)
	{
		vector<int> sep = iter->first; 
		if (VectorPlus::Equals(sep, sep_nodes))
			return index;
		else
			index++;
	}
	return -1;
}

int TensorJTNode::GetSIndex()
{
	if (S_index >= 0)
	{
		return S_index;
	}
	else
	{
		S_index = jtree->GetChildSepIndex(parent_id, *S_vars);
		return S_index;
	}
}
