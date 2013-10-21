# pragma once

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
#include "FastIndexer.hpp"

using namespace std;

// default constructor
TensorCPT::TensorCPT()
{
	prob_tensor = NULL;

	vars_to_modes = NULL;
	modes_to_vars = NULL;
	
	all_vars = NULL;
	all_modes = NULL;

	isInitialized = false;
}

TensorCPT::~TensorCPT()
{
	delete(prob_tensor);

	delete(modes_to_vars);
	
	delete(all_vars);
	delete(all_modes);

    map<int, vector<int>*>::iterator vec_iter;
    for (vec_iter = vars_to_modes->begin(); vec_iter != vars_to_modes->end(); ++vec_iter) {
        delete(vec_iter->second);
    }

	delete(vars_to_modes);
}


TensorCPT::TensorCPT(TensorCPT& other)
{
	prob_tensor = new Tensor(*(other.prob_tensor));
	modes_to_vars = new vector<int>(*(other.modes_to_vars));
	vars_to_modes = new map<int, vector<int>*>();

	all_vars = new vector<int>(*(other.all_vars));
	all_modes = new vector<int>(*(other.all_modes));


	isInitialized = true;

	for (int i = 0; i < all_vars->size(); ++i)
	{
		(*vars_to_modes)[all_vars->at(i)] = new vector<int>(*((*other.vars_to_modes)[all_vars->at(i)]));
	}
}

TensorCPT::TensorCPT(CPT& table, vector<int>& vars, vector<int>& mode_counts)
{
	prob_tensor = NULL;
	
	all_vars = NULL;
	all_modes = NULL;

	vars_to_modes = NULL;
	modes_to_vars = NULL;


	isInitialized = false;

	Initialize(table, vars, mode_counts);
}

void TensorCPT::Initialize(TensorCPT& other)
{
	prob_tensor = new Tensor(*(other.prob_tensor));
	modes_to_vars = new vector<int>(*(other.modes_to_vars));
	vars_to_modes = new map<int, vector<int>*>();

	all_vars = new vector<int>(*(other.all_vars));
	all_modes = new vector<int>(*(other.all_modes));


	isInitialized = true;

	for (int i = 0; i < all_vars->size(); ++i)
	{
		(*vars_to_modes)[all_vars->at(i)] = new vector<int>(*((*other.vars_to_modes)[all_vars->at(i)]));
	}
}

void TensorCPT::Initialize(Tensor* tensor, vector<int>& modes_to_vars)
{
	assert(!isInitialized);

	vars_to_modes = new map<int, vector<int>*>();
	prob_tensor = tensor;
	this->modes_to_vars = new vector<int>(modes_to_vars);

	all_modes = new vector<int>();
	VectorPlus::Seq(*all_modes, 0, 1, modes_to_vars.size());

	all_vars = new vector<int>();
	
	for (int i = 0; i < this->modes_to_vars->size(); ++i)
	{
		if (vars_to_modes->count(this->modes_to_vars->at(i)) == 0)
		{
			(*vars_to_modes)[this->modes_to_vars->at(i)] = new vector<int>();
			all_vars->push_back(this->modes_to_vars->at(i));
		}

		(*vars_to_modes)[this->modes_to_vars->at(i)]->push_back(i);
	}


	isInitialized = true;
}

void TensorCPT::Initialize(Tensor& tensor, vector<int>& modes_to_vars)
{
	assert(!isInitialized);

	vars_to_modes = new map<int, vector<int>*>();
	prob_tensor = new Tensor(tensor);
	this->modes_to_vars = new vector<int>(modes_to_vars);

	all_modes = new vector<int>();
	VectorPlus::Seq(*all_modes, 0, 1, modes_to_vars.size());

	all_vars = new vector<int>();
	
	for (int i = 0; i < this->modes_to_vars->size(); ++i)
	{
		if (vars_to_modes->count(this->modes_to_vars->at(i)) == 0)
		{
			(*vars_to_modes)[this->modes_to_vars->at(i)] = new vector<int>();
			all_vars->push_back(this->modes_to_vars->at(i));
		}

		(*vars_to_modes)[this->modes_to_vars->at(i)]->push_back(i);
	}


	isInitialized = true;
}


void TensorCPT::Initialize(CPT& table, vector<int>& vars, vector<int>& mode_counts)
{
	assert(!isInitialized);
	
	all_vars = new vector<int>(vars);
	vars_to_modes = new map<int, vector<int>*>();
	modes_to_vars = new vector<int>();

	vector<int> dims;

	int total_modes = 0;
	for (int i = 0; i < vars.size(); ++i)
	{
		(*vars_to_modes)[vars[i]] = new vector<int>();
	}
	for (int i = 0; i < vars.size(); ++i)
	{
		for (int j = 0; j < mode_counts[i]; ++j)
		{
			(*vars_to_modes)[vars[i]]->push_back(total_modes++);
			modes_to_vars->push_back(vars[i]);
			dims.push_back(table.VarDim(vars[i]));
		}
	}
	prob_tensor = new Tensor(dims);
	all_modes = new vector<int>();
	VectorPlus::Seq(*all_modes, 0, 1, total_modes);

	vector<int>& table_vars = table.Vars();
	FastIndexer t_indexer(table.Dims());
	for (int i = 0; i < table.NumElements(); ++i)
	{
		vector<int>& indices = t_indexer.GetNext();
	//	table.ComputeIndexArray(indices, i);
		vector<int> diag_indices;
		ComputeDiagIndices(diag_indices, vars, indices);
		prob_tensor->Set(diag_indices, table.At(i));
	}

	isInitialized = true;
}

void TensorCPT::ComputeDiagIndices(vector<int>& diag_indices, vector<int>& vars, vector<int>& indices)
{
	diag_indices.assign(all_modes->size(), 0);
	for (int i = 0; i < vars.size(); ++i) 
	{
		vector<int>& modes = *((*vars_to_modes)[vars.at(i)]);
		for (int j = 0; j < modes.size(); ++j)
		{
			diag_indices[modes[j]] = indices[i];
		}
	}
}

void TensorCPT::ElementwiseMultiply(TensorCPT& result_CPT, TensorCPT& cpt_A, TensorCPT& cpt_B)
{
	vector<int> mult_vars;
	VectorPlus::Intersect(mult_vars, cpt_A.Vars(), cpt_B.Vars());

	vector<int> A_mult_modes;
	vector<int> B_mult_modes;

	for (int i = 0; i < mult_vars.size(); ++i)
	{
		vector<int>& A_modes = cpt_A.Modes(mult_vars[i]);
		assert(A_modes.size() == 1);
	    A_mult_modes.push_back(A_modes[0]);

		vector<int>& B_modes = cpt_B.Modes(mult_vars[i]);
		assert(B_modes.size() == 1);
		B_mult_modes.push_back(B_modes[0]);
	}

	// create new tensor
	Tensor* result_tensor = new Tensor();
	Tensor::ElementwiseMultiply(*result_tensor, *cpt_A.prob_tensor, *cpt_B.prob_tensor, A_mult_modes, B_mult_modes);

	// create new modes maps
	vector<int> new_modes_to_vars;
	vector<int>& A_all_modes = cpt_A.Modes();

	vector<int> B_free_modes;
	VectorPlus::SetDiff(B_free_modes, cpt_B.Modes(), B_mult_modes);

	for (int i = 0; i < A_all_modes.size(); ++i)
	{
		new_modes_to_vars.push_back(cpt_A.Var(A_all_modes[i]));
	}
	for (int i = 0; i < B_free_modes.size(); ++i)
	{
		new_modes_to_vars.push_back(cpt_B.Var(B_free_modes[i]));
	}

	result_CPT.Initialize(result_tensor, new_modes_to_vars);
}

void TensorCPT::Multiply(TensorCPT& result_CPT, TensorCPT& cpt_A, TensorCPT& cpt_B)
{
	vector<int> mult_vars;
	VectorPlus::Intersect(mult_vars, cpt_A.Vars(), cpt_B.Vars());

	vector<int> A_mult_modes;
	vector<int> B_mult_modes;

	for (int i = 0; i < mult_vars.size(); ++i)
	{
		vector<int>& A_modes = cpt_A.Modes(mult_vars[i]);
	    A_mult_modes.push_back(A_modes[A_modes.size() - 1]);

		vector<int>& B_modes = cpt_B.Modes(mult_vars[i]);
		B_mult_modes.push_back(B_modes[0]);
	}

	// create new tensor
	Tensor* result_tensor = new Tensor();
	Tensor::Multiply(*result_tensor, *cpt_A.prob_tensor, *cpt_B.prob_tensor, A_mult_modes, B_mult_modes);

	// create new modes maps
	vector<int> new_modes_to_vars;
	vector<int> A_free_modes;
	vector<int> B_free_modes;
	VectorPlus::SetDiff(A_free_modes, cpt_A.Modes(), A_mult_modes);
	VectorPlus::SetDiff(B_free_modes, cpt_B.Modes(), B_mult_modes);

	for (int i = 0; i < A_free_modes.size(); ++i)
	{
		new_modes_to_vars.push_back(cpt_A.Var(A_free_modes[i]));
	}
	for (int i = 0; i < B_free_modes.size(); ++i)
	{
		new_modes_to_vars.push_back(cpt_B.Var(B_free_modes[i]));
	}

	result_CPT.Initialize(result_tensor, new_modes_to_vars);
}

void TensorCPT::OneReduce(TensorCPT& result_CPT, TensorCPT& cpt_A, vector<int>& one_vars)
{
	if (one_vars.size() == 0)
	{
		result_CPT.Initialize(cpt_A);
		return;
	}

	vector<int> one_modes;
	for (int i = 0; i < one_vars.size(); ++i)
	{
		vector<int> m = cpt_A.Modes(one_vars[i]);
		assert(m.size() == 1);
		one_modes.push_back(m[0]);
	}

	vector<int> ones_dims;
	VectorPlus::Subset(ones_dims, cpt_A.Dims(), one_modes);
	Tensor ones(ones_dims);
	ones.FillWithConst(1);

	Tensor* result_tensor = new Tensor();
	Tensor::Multiply(*result_tensor, *cpt_A.prob_tensor, ones, one_modes, ones.Modes());

		// create new modes maps
	vector<int> new_modes_to_vars;
	vector<int> A_free_modes;
	VectorPlus::SetDiff(A_free_modes, cpt_A.Modes(), one_modes);
	for (int i = 0; i < A_free_modes.size(); ++i)
	{
		new_modes_to_vars.push_back(cpt_A.Var(A_free_modes[i]));
	}
	result_CPT.Initialize(result_tensor, new_modes_to_vars);
}

void TensorCPT::Slice(TensorCPT& result_CPT, int fixed_mode, int fixed_index)
{
	Tensor* result_tensor = new Tensor();
	prob_tensor->Slice(*result_tensor, fixed_mode, fixed_index);
	vector<int> new_modes_to_vars;
	vector<int> var_vec = VectorPlus::CreateSingleton((*modes_to_vars)[fixed_mode]);
	VectorPlus::SetDiff(new_modes_to_vars, *modes_to_vars, var_vec);
	result_CPT.Initialize(result_tensor, new_modes_to_vars);
}

void TensorCPT::Slice(TensorCPT& result_CPT, vector<int> fixed_vars, vector<int> fixed_vals)
{	
	vector<int> temp_modes;
	vector<int> temp_vals;
	for (int i = 0; i < fixed_vars.size(); ++i)
	{
		vector<int>& modes = Modes(fixed_vars[i]);
		VectorPlus::Union(temp_modes, modes);
		for (int j = 0; j < modes.size(); ++j)
		{
			temp_vals.push_back(fixed_vals[i]);
		}
	}

	vector<int> fixed_modes;
	vector<int> fixed_indices;
	VectorPlus::MultiSort(fixed_modes, fixed_indices, temp_modes, temp_vals);

	Tensor* result_tensor = new Tensor();
	prob_tensor->Slice(*result_tensor, fixed_modes, fixed_indices);

	vector<int> new_modes_to_vars;
	VectorPlus::SetDiff(new_modes_to_vars, *modes_to_vars, fixed_vars);
	result_CPT.Initialize(result_tensor, new_modes_to_vars);
}

void TensorCPT::Select(TensorCPT& result_CPT, vector<int> fixed_vars, vector<int> fixed_vals)
{
	vector<int> temp_modes;
	vector<int> temp_vals;
	for (int i = 0; i < fixed_vars.size(); ++i)
	{
		vector<int>& modes = Modes(fixed_vars[i]);
		VectorPlus::Union(temp_modes, modes);
		for (int j = 0; j < modes.size(); ++j)
		{
			temp_vals.push_back(fixed_vals[i]);
		}
	}

	vector<int> fixed_modes;
	vector<int> fixed_indices;
	VectorPlus::MultiSort(fixed_modes, fixed_indices, temp_modes, temp_vals);

	Tensor* result_tensor = new Tensor();
	prob_tensor->Select(*result_tensor, fixed_modes, fixed_indices);
	result_CPT.Initialize(result_tensor, *modes_to_vars);
}

void TensorCPT::Add(TensorCPT& cpt_A, TensorCPT& cpt_B)
{
	Tensor::Add(*(cpt_A.prob_tensor), *(cpt_B.prob_tensor));

}

void TensorCPT::Divide(TensorCPT& result_CPT, TensorCPT& cpt_A, double val)
{
	Tensor* result_tensor = new Tensor();
	Tensor::Divide(*result_tensor, *(cpt_A.prob_tensor), val);
	result_CPT.Initialize(result_tensor, *(cpt_A.modes_to_vars));
}

void TensorCPT::WeightedAdd(TensorCPT& cpt_A, TensorCPT& cpt_B, double weight_A, double weight_B)
{
	Tensor::WeightedAdd(*(cpt_A.prob_tensor), *(cpt_B.prob_tensor), weight_A, weight_B);
}


// assumes one mode per variable ----- used for EM only
void TensorCPT::ComputeConditionalFactor(TensorCPT& conditional_factor, TensorCPT& factor, vector<int>& front_vars, vector<int>& back_vars)
{
	vector<int>& all_vars = factor.Vars();
	vector<int> extra_vars;
	vector<int> relevant_vars;
	VectorPlus::Union(relevant_vars, front_vars, back_vars);
	VectorPlus::SetDiff(extra_vars, all_vars, relevant_vars);
	
	OneReduce(conditional_factor, factor, extra_vars);

	TensorCPT denom_factor;
	OneReduce(denom_factor, conditional_factor, front_vars);

	vector<int> front_modes;
	for (int i = 0; i < front_vars.size(); ++i)
	{
		vector<int>& temp_modes = factor.Modes(front_vars[i]);
		assert(temp_modes.size() == 1);
		front_modes.push_back(temp_modes[0]);
	}
	
	FastIndexer c_indexer(conditional_factor.Dims());

	for (int i = 0; i < conditional_factor.NumElements(); ++i)
	{
		vector<int> indices = c_indexer.GetNext();
	//	conditional_factor.ComputeIndexArray(indices, i);
		
		vector<int> parent_indices;
		VectorPlus::CSubset(parent_indices, indices, front_modes);
		
		if (parent_indices.size() == 0)
			conditional_factor.Set(i, conditional_factor.At(i));
		else
			conditional_factor.Set(i, conditional_factor.At(i) / (denom_factor.At(parent_indices)));
	}
}

// only works for non diagonal tensors
void TensorCPT::Rearrange(TensorCPT& result_CPT, vector<int>& new_modes_to_vars)
{
	Tensor* result_tensor = new Tensor();

	vector<int> old_to_new_modes(modes_to_vars->size(), 0);
	for (int i = 0; i < modes_to_vars->size(); ++i)
	{
		int index = VectorPlus::Find(new_modes_to_vars, modes_to_vars->at(i));
		old_to_new_modes[i] = index;
	}

	Tensor::Rearrange(*result_tensor, GetTensor(), old_to_new_modes);
	result_CPT.Initialize(result_tensor, new_modes_to_vars);
}
