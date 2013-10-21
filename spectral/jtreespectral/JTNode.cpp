# pragma once

#include "JTNode.hpp"
#include "VectorPlus.hpp"
#include <vector>
#include <assert.h>
#include "MultiVector.hpp"
#include <boost/shared_ptr.hpp>
#include "CPT.hpp"

JTNode::JTNode()
{
	jtree = NULL;

	parent_id = -1;
	node_id = -1;
	children = NULL;

	S_vars = NULL;
	R_vars = NULL;
	C_vars = NULL;
	agglomerate_tensor = NULL;
	factor_list = new vector<TensorCPT*>();
	front_vars_list = new vector<vector<int>*>();

	clique_params_flag = false;
}

JTNode::JTNode(JTree* jtree, int node_id, bool clique_params_flag)
{
	this->jtree = jtree;

	parent_id = -1;
	this->node_id = node_id;
	children = NULL;

	S_vars = NULL;
	R_vars = NULL;
	C_vars = NULL;

	factor_list = new vector<TensorCPT*>();
	front_vars_list = new vector<vector<int>*>();
	agglomerate_tensor = NULL;

	this->clique_params_flag = clique_params_flag;
}

JTNode::JTNode(JTree* jtree, int node_id)
{
	this->jtree = jtree;

	parent_id = -1;
	this->node_id = node_id;
	children = NULL;

	S_vars = NULL;
	R_vars = NULL;
	C_vars = NULL;

	factor_list = new vector<TensorCPT*>();
	front_vars_list = new vector<vector<int>*>();
	agglomerate_tensor = NULL;

	clique_params_flag = false;
}

JTNode::~JTNode()
{
	delete(children);

	delete(S_vars);
	delete(R_vars);
	delete(C_vars);

	for (int i = 0; i < factor_list->size(); ++i)
	{
		delete(factor_list->at(i));
		delete(front_vars_list->at(i));
	}

	delete(factor_list);

	delete(front_vars_list);
	delete(agglomerate_tensor);
}

// basically a way of making the agglomerate factor the permanent factor
void JTNode::CreateCliqueFactor()
{
	// delete everything in factor list, only update the agglomerate tensor
	RecomputeAgglomCliqueFactor();

	if (clique_params_flag)
	{
		/*for (int i = 0; i < factor_list->size(); ++i)
		{
			delete(factor_list->at(i));
			delete(front_vars_list->at(i));
		}

		delete(factor_list);
		delete(front_vars_list);

		factor_list = new vector<TensorCPT*>();
		front_vars_list = new vector<vector<int>*>(); */

		clique_params_flag = true;
	}
	else
	{
//		assert(0);
	}
}

void JTNode::AddFactor(CPT& t_cpt)
{
	vector<int>& cpt_vars = t_cpt.Vars();
	vector<int>& front_vars = t_cpt.GetFrontVars();

	vector<int> mode_counts(cpt_vars.size(), 1);
	CPT* rand_cpt = new CPT(t_cpt);
	rand_cpt->FillWithRandom();
	TensorCPT* new_factor = new TensorCPT(*rand_cpt, cpt_vars, mode_counts);
	delete(rand_cpt);
	factor_list->push_back(new_factor);
	front_vars_list->push_back(new vector<int>(front_vars));
}

void JTNode::UpdateFactorWrapper(TensorCPT& clique_marginal)
{
	if (clique_params_flag == true)
		UpdateCliqueFactor(clique_marginal);
	else
		UpdateFactors(clique_marginal);
}

void JTNode::UpdateCliqueFactor(TensorCPT& clique_marginal)
{
	assert(clique_params_flag == true);
	TensorCPT* new_conditional_factor = new TensorCPT();
	TensorCPT::ComputeConditionalFactor(*new_conditional_factor, clique_marginal, *R_vars, *S_vars);

	delete(agglomerate_tensor);
	agglomerate_tensor = new_conditional_factor;
}

void JTNode::UpdateFactors(TensorCPT& clique_marginal)
{
	assert(clique_params_flag == false);
	for (int i = 0; i < factor_list->size(); ++i)
	{
		vector<int>& all_vars = factor_list->at(i)->Vars();
		vector<int>& front_var = *(front_vars_list->at(i));
		vector<int> back_vars;
		VectorPlus::SetDiff(back_vars, all_vars, front_var);

		TensorCPT* new_conditional_factor = new TensorCPT();
		TensorCPT::ComputeConditionalFactor(*new_conditional_factor, clique_marginal, front_var, back_vars);

		delete(factor_list->at(i));
		factor_list->at(i) = new_conditional_factor;
	}

	RecomputeAgglomCliqueFactor();
}

TensorCPT& JTNode::GetAgglomCliqueFactor()
{
	return *agglomerate_tensor;
}


void JTNode::RecomputeAgglomCliqueFactor()
{
	if (agglomerate_tensor != NULL)
		delete(agglomerate_tensor);

	if (factor_list->size() == 0)
		assert(0);

	agglomerate_tensor = new TensorCPT(*(factor_list->at(0)));

	for (int i = 1; i < factor_list->size(); ++i)
	{
		TensorCPT* temp_tensor = new TensorCPT();
		TensorCPT::ElementwiseMultiply(*temp_tensor, *agglomerate_tensor, *(factor_list->at(i)));
		delete(agglomerate_tensor);
		agglomerate_tensor = temp_tensor;
	}
}
