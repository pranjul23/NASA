# pragma once

#include "CPT.hpp"
#include "tensor.hpp"
#include <assert.h>
#include <algorithm>
#include <math.h>
#include <queue>
#include "VectorPlus.hpp"
#include <assert.h>
#include "StatsLibrary.hpp"

using namespace std;

// default constructor
CPT::CPT()
{
	prob_tensor = NULL;
	
	all_modes = NULL;
	all_vars = NULL;

	vars_to_modes = NULL;
	modes_to_vars = NULL;

	parent_vars = NULL;
	parent_modes = NULL;
	parent_dims = NULL;

	front_vars = NULL;
	front_modes = NULL;
	front_dims = NULL;

	isInitialized = false;
}

// copy constructor
CPT::CPT(CPT& other)
{
	prob_tensor = new Tensor(*other.prob_tensor);
	
	all_modes = new vector<int>(*other.all_modes);
	all_vars = new vector<int>(*other.all_vars);

	vars_to_modes = new map<int,int>(*other.vars_to_modes);
	modes_to_vars = new vector<int>(*other.modes_to_vars);

	parent_vars = new vector<int>(*other.parent_vars);
	parent_modes = new vector<int>(*other.parent_modes);
	parent_dims = new vector<int>(*other.parent_dims);

	front_vars = new vector<int>(*other.front_vars);
	front_modes = new vector<int>(*other.front_modes);
	front_dims = new vector<int>(*other.front_dims);

	isInitialized = false;
}

CPT::CPT(vector<int>& vars, vector<int>& dims)
{
	prob_tensor = NULL;
	
	all_modes = NULL;
	all_vars = NULL;

	vars_to_modes = NULL;
	modes_to_vars = NULL;

	parent_vars = NULL;
	parent_modes = NULL;
	parent_dims = NULL;

	front_vars = NULL;
	front_modes = NULL;
	front_dims = NULL;

	isInitialized = false;

	Initialize(vars, dims);
}

CPT::CPT(vector<int>& front_vars, vector<int>& back_vars, vector<int>& dims)
{
	prob_tensor = NULL;
	
	all_modes = NULL;
	all_vars = NULL;

	vars_to_modes = NULL;
	modes_to_vars = NULL;

	parent_vars = NULL;
	parent_modes = NULL;
	parent_dims = NULL;

	this->front_vars = NULL;
	front_modes = NULL;
	front_dims = NULL;

	isInitialized = false;

	Initialize(front_vars, back_vars, dims);
}

// destructor
CPT::~CPT()
{
	delete(prob_tensor);

	delete(all_modes);
	delete(all_vars);

	delete(vars_to_modes);
	delete(modes_to_vars);

	delete(parent_vars);
	delete(parent_modes);
	delete(parent_dims);

	delete(front_vars);
	delete(front_modes);
	delete(front_dims);
}


// Initialize
// assumes first variable is the only front variable by default
void CPT::Initialize(vector<int>& vars, vector<int>& dims)
{
	assert(!isInitialized);

	all_vars = new vector<int>(vars);
	all_modes = new vector<int>();
	VectorPlus::Seq(*all_modes, 0, 1, vars.size());
	modes_to_vars = new vector<int>(vars);
	prob_tensor = new Tensor(dims);

	vars_to_modes = new map<int, int>();
	for (int i = 0; i < modes_to_vars->size(); ++i)
	{
		(*vars_to_modes)[modes_to_vars->at(i)] = i;
	}

	front_vars = new vector<int>();
	front_modes = new vector<int>();
	front_dims = new vector<int>();

	front_vars->push_back(vars[0]);
	front_modes->push_back(0);
	front_dims->push_back(dims[0]);

	parent_vars = new vector<int>();
	VectorPlus::Subset(*parent_vars, *all_vars, 1, all_vars->size());
	parent_modes = new vector<int>(vars.size() - 1, 0);
	parent_dims = new vector<int>(vars.size() -1, 0);

	for (int i = 0; i < parent_modes->size(); ++i)
	{
		parent_modes->at(i) = i+1;
		parent_dims->at(i) = dims[i+1];
	}

	isInitialized = true;
}

void CPT::Initialize(vector<int>& front_vars, vector<int>& back_vars, vector<int>& dims)
{
	assert(!isInitialized);

	all_vars = new vector<int>();
	VectorPlus::Union(*all_vars, front_vars, back_vars);

	all_modes = new vector<int>();
	VectorPlus::Seq(*all_modes, 0, 1, all_vars->size());
	modes_to_vars = new vector<int>(*all_vars);
	prob_tensor = new Tensor(dims);

	vars_to_modes = new map<int, int>();
	for (int i = 0; i < modes_to_vars->size(); ++i)
	{
		(*vars_to_modes)[modes_to_vars->at(i)] = i;
	}

	this->front_vars = new vector<int>();
	front_modes = new vector<int>();
	front_dims = new vector<int>();

	for (int i = 0; i < front_vars.size(); ++i)
	{
		this->front_vars->push_back(front_vars[i]);
		front_modes->push_back((*vars_to_modes)[front_vars[i]]);
		front_dims->push_back(dims[i]);
	}

	parent_vars = new vector<int>();
	parent_modes = new vector<int>();
	parent_dims = new vector<int>();

	for (int i = 0; i < back_vars.size(); ++i)
	{
		this->parent_vars->push_back(back_vars[i]);
		parent_modes->push_back((*vars_to_modes)[back_vars[i]]);
		parent_dims->push_back(dims[i + front_vars.size()]);
	}

	isInitialized = true;
}

// Fill CPT randomly such that normalization is satisfied
void CPT::FillWithConst()
{
	do
	{
		prob_tensor->FillWithConst(1);

		if (parent_modes->size() == 0)
		{
			double sum = prob_tensor->Sum();
			for (int i = 0; i < prob_tensor->NumElements(); ++i)
			{
				prob_tensor->Set(i, prob_tensor->At(i) / sum);
			}
			return;
		}

		vector<int> parent_offsets;
		vector<int> front_offsets;
		Tensor::ComputeOffsets(parent_offsets, *parent_dims);
		Tensor::ComputeOffsets(front_offsets, *front_dims);
		int num_front_vals = VectorPlus::Product(*front_dims);
		int num_parent_vals = VectorPlus::Product(*parent_dims);
		for (int i = 0; i < num_parent_vals; ++i)
		{
			vector<int> parent_indices;
			Tensor::ComputeIndexArray(parent_indices, parent_offsets, i);
			double partial_sum = prob_tensor->ComputePartialSum(*parent_modes, parent_indices);

			
			for (int j = 0; j < num_front_vals; ++j)
			{
				vector<int> total_indices;
				vector<int> front_indices;
				Tensor::ComputeIndexArray(front_indices, front_offsets, j);
				Tensor::MergeIndices(total_indices,  *front_modes, *parent_modes, front_indices, parent_indices);  
				prob_tensor->Set(total_indices, prob_tensor->At(total_indices) / partial_sum);
			}
		}
	} while (prob_tensor->Min() < 0.01);
}

// Fill CPT randomly such that normalization is satisfied
void CPT::FillWithRandom()
{
	do
	{
		prob_tensor->FillWithRandom();

		if (parent_modes->size() == 0)
		{
			double sum = prob_tensor->Sum();
			for (int i = 0; i < prob_tensor->NumElements(); ++i)
			{
				prob_tensor->Set(i, prob_tensor->At(i) / sum);
			}
			return;
		}

		vector<int> parent_offsets;
		vector<int> front_offsets;
		Tensor::ComputeOffsets(parent_offsets, *parent_dims);
		Tensor::ComputeOffsets(front_offsets, *front_dims);
		int num_front_vals = VectorPlus::Product(*front_dims);
		int num_parent_vals = VectorPlus::Product(*parent_dims);
		for (int i = 0; i < num_parent_vals; ++i)
		{
			vector<int> parent_indices;
			Tensor::ComputeIndexArray(parent_indices, parent_offsets, i);
			double partial_sum = prob_tensor->ComputePartialSum(*parent_modes, parent_indices);

			
			for (int j = 0; j < num_front_vals; ++j)
			{
				vector<int> total_indices;
				vector<int> front_indices;
				Tensor::ComputeIndexArray(front_indices, front_offsets, j);
				Tensor::MergeIndices(total_indices,  *front_modes, *parent_modes, front_indices, parent_indices);  
				prob_tensor->Set(total_indices, prob_tensor->At(total_indices) / partial_sum);
			}
		}
	} while (prob_tensor->Min() < 0.01);
}

int CPT::Sample(vector<int>& parent_vals)
{
	vector<double> probs;
	assert(front_vars->size() == 1);
	int node_mode = front_modes->at(0);
	for (int i = 0; i < prob_tensor->Dim(node_mode); ++i)
	{
		if (parent_vals.size() == 0)
		{
			probs.push_back(prob_tensor->At(i));
		}
		else
		{
			vector<int> total_indices;
			vector<int> node_mode_vec(1, node_mode);
			vector<int> node_index(1,i);
	
			Tensor::MergeIndices(total_indices,  node_mode_vec, *parent_modes, node_index, parent_vals);
		
			probs.push_back(prob_tensor->At(total_indices));
		}
	}

	return StatsLibrary::sample_multinomial(probs);
}
