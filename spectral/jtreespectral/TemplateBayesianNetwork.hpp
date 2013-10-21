# pragma once

#include <vector>
#include "Matrix.hpp"
#include "CPT.hpp"
#include "TensorCPT.hpp"
#include "MultiVector.hpp"
#include "BayesianNetwork.hpp"

using namespace std;

class TemplateBayesianNetwork : public BayesianNetwork {

private:
	vector<int>* var_to_template;
	vector<vector<int>*>* template_to_vars;

public:
	TemplateBayesianNetwork();
	TemplateBayesianNetwork(vector<vector<int>*>& parent_lists, vector<int>& type_vector, vector<int>& dims, vector<int>& cpt_template);
	~TemplateBayesianNetwork();

	void Initialize(vector<vector<int>*>& parent_lists, vector<int>& type_vector, vector<int>& dims, vector<int>& cpt_template);
	void InitializeTemplate(vector<int>& cpt_template);

	int GetTemplateIndex(int node_id) { return var_to_template->at(node_id); }
	vector<int>& GetNodes(int template_index) { return *(template_to_vars->at(template_index)); }

	int GetNumTemplates() {return template_to_vars->size(); }

	CPT& GetTemplateCPT(int template_index) { return *(CPTs->at(template_to_vars->at(template_index)->at(0))); }
};
