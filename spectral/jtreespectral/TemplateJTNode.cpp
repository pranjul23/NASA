# pragma once

#include "TemplateJTNode.hpp"
#include "VectorPlus.hpp"
#include <vector>
#include <assert.h>
#include "MultiVector.hpp"
#include <boost/shared_ptr.hpp>
#include "CPT.hpp"

TemplateJTNode::TemplateJTNode()
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
	template_list = new vector<int>();
	front_var_list = new vector<int>();
}

TemplateJTNode::TemplateJTNode(TemplateJTree* jtree, int node_id)
{
	this->jtree = jtree;

	parent_id = -1;
	this->node_id = node_id;
	children = NULL;

	S_vars = NULL;
	R_vars = NULL;
	C_vars = NULL;

	factor_list = new vector<TensorCPT*>();
	front_var_list = new vector<int>();
	agglomerate_tensor = NULL;
	template_list = new vector<int>();
}

TemplateJTNode::~TemplateJTNode()
{
	delete(children);

	delete(S_vars);
	delete(R_vars);
	delete(C_vars);

	for (int i = 0; i < factor_list->size(); ++i)
		delete(factor_list->at(i));

	delete(factor_list);
	delete(front_var_list);
	delete(template_list);
	delete(agglomerate_tensor);
}


void TemplateJTNode::AddFactor(CPT& t_cpt, int front_var)
{
	vector<int>& cpt_vars = t_cpt.Vars();
	vector<int> mode_counts(cpt_vars.size(), 1);
	CPT* rand_cpt = new CPT(t_cpt);
	rand_cpt->FillWithRandom();
	TensorCPT* new_factor = new TensorCPT(*rand_cpt, cpt_vars, mode_counts);
	delete(rand_cpt);
	factor_list->push_back(new_factor);
	front_var_list->push_back(front_var);
}

void TemplateJTNode::UpdateFactors(TensorCPT& clique_marginal)
{
	for (int i = 0; i < factor_list->size(); ++i)
	{
		vector<int>& all_vars = factor_list->at(i)->Vars();
		vector<int> front_var = VectorPlus::CreateSingleton(front_var_list->at(i));
		vector<int> back_vars;
		VectorPlus::SetDiff(back_vars, all_vars, front_var);

		TensorCPT* new_conditional_factor = new TensorCPT();
		TensorCPT::ComputeConditionalFactor(*new_conditional_factor, clique_marginal, front_var, back_vars);

		delete(factor_list->at(i));
		factor_list->at(i) = new_conditional_factor;
	}
}

TensorCPT& TemplateJTNode::GetAgglomCliqueFactor()
{
/*	if (agglomerate_tensor != NULL)
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
	return *agglomerate_tensor; */

	if (agglomerate_tensor != NULL)
	{
		return *agglomerate_tensor;

	}

	if (template_list->size() == 0)
		assert(0);

	int first_template_index = template_list->at(0);
	agglomerate_tensor = new TensorCPT();
	jtree->GetFactor(*agglomerate_tensor, first_template_index, front_var_list->at(0));
	for (int i = 1; i < template_list->size(); ++i)
	{
		int template_index = template_list->at(i);
		TensorCPT* temp_tensor = new TensorCPT();
		TensorCPT new_factor;
		jtree->GetFactor(new_factor, template_index, front_var_list->at(i));
		TensorCPT::ElementwiseMultiply(*temp_tensor, *agglomerate_tensor, new_factor);
		delete(agglomerate_tensor);
		agglomerate_tensor = temp_tensor;
	}
	return *agglomerate_tensor;
}
