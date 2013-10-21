# pragma once

# include "TemplateBayesianNetwork.hpp"
#include "VectorPlus.hpp"

TemplateBayesianNetwork::TemplateBayesianNetwork()
  :BayesianNetwork()
{}

TemplateBayesianNetwork::TemplateBayesianNetwork(vector<vector<int>*>& parent_lists, vector<int>& type_vector, vector<int>& dims, vector<int>& cpt_template)
	:BayesianNetwork(parent_lists, type_vector, dims)
{
	InitializeTemplate(cpt_template);
}

void TemplateBayesianNetwork::Initialize(vector<vector<int>*>& parent_lists, vector<int>& type_vector, vector<int>& dims, vector<int>& cpt_template)
{
	BayesianNetwork::Initialize(parent_lists, type_vector, dims);
	InitializeTemplate(cpt_template);
}


void TemplateBayesianNetwork::InitializeTemplate(vector<int>& cpt_template)
{
	var_to_template = new vector<int>(cpt_template);

	vector<int> unique_templates;
	VectorPlus::Unique(unique_templates, *var_to_template);
	template_to_vars = new vector<vector<int>*>();
	for (int i = 0; i <  unique_templates.size(); ++i)
	{
		template_to_vars->push_back(new vector<int>());
	}
	for (int i = 0; i < var_to_template->size(); ++i)
	{
		template_to_vars->at(var_to_template->at(i))->push_back(i);
	}

	for (int i = 0; i < unique_templates.size(); ++i)
	{
		int first_var = template_to_vars->at(i)->at(0);
		for (int j = 1; j < template_to_vars->at(i)->size(); ++j)
		{
			CPTs->at(template_to_vars->at(i)->at(j))->SetTensor(CPTs->at(first_var)->GetTensor());
		}
	}
}

TemplateBayesianNetwork::~TemplateBayesianNetwork()
{
	delete(var_to_template);
	for (int i = 0; i < template_to_vars->size(); ++i)
		delete(template_to_vars->at(i));
	delete(template_to_vars);
}
