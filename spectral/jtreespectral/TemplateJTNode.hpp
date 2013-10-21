# pragma once

#include <vector>
#include <map>
#include "tensor.hpp"
#include "TemplateJTree.hpp"
#include "TensorCPT.hpp"
#include "VectorPlus.hpp"

using namespace std;


class TemplateJTree;

class TemplateJTNode {

private:

	TemplateJTree* jtree;

	int parent_id;
	int node_id;

	vector<int>* children;

	vector<int>* S_vars;
	vector<int>* R_vars;
	vector<int>* C_vars;

	vector<TensorCPT*>* factor_list;
	TensorCPT* agglomerate_tensor;

	vector<int>* template_list;

	vector<int>* front_var_list;

public:

	TemplateJTNode();
	~TemplateJTNode();
	TemplateJTNode(TemplateJTree* jtree, int node_id);

	void AddFactor(CPT& cpt, int front_var);
	int GetFrontVar(int template_index) { int index = VectorPlus::Find(*template_list, template_index);  return front_var_list->at(index); }
	void AddTemplate(int template_index, int front_var) { template_list->push_back(template_index); front_var_list->push_back(front_var); }
	int GetNumFactors() { return template_list->size(); }
	bool ContainsTemplate(int template_index) { return VectorPlus::Contains(*template_list, template_index); }
	void SetParent(int parent) { this->parent_id = parent; }
	void SetChildren(vector<int>& children) {this->children = new vector<int>(children); }
	void SetCVars(vector<int>& c_vars) {this->C_vars = new vector<int>(c_vars); }
	void SetSVars(vector<int>& s_vars) {this->S_vars = new vector<int>(s_vars); }
	void SetRVars(vector<int>& r_vars) {this->R_vars = new vector<int>(r_vars); }

	void UpdateFactors(TensorCPT& clique_marginal);

	int GetParent() { return parent_id; }
	vector<int>& GetChildren() { return *children; }
	vector<int>& GetCVars() { return *C_vars; }
	vector<int>& GetSVars() {return *S_vars; }
	vector<int>& GetRVars() {return *R_vars; }

	TensorCPT& GetAgglomCliqueFactor();
	void ResetAgglomCliqueFactor() { if (agglomerate_tensor != NULL) { delete(agglomerate_tensor); agglomerate_tensor = NULL; } }
};


