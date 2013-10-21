# pragma once

#include <vector>
#include <map>
#include "tensor.hpp"
#include "JTree.hpp"
#include "TensorCPT.hpp"
#include "VectorPlus.hpp"

using namespace std;


class JTree;

class JTNode {

private:

	JTree* jtree;

	int parent_id;
	int node_id;

	vector<int>* children;

	vector<int>* S_vars;
	vector<int>* R_vars;
	vector<int>* C_vars;

	vector<TensorCPT*>* factor_list;
	TensorCPT* agglomerate_tensor;

	vector<vector<int>*>* front_vars_list;

	bool clique_params_flag;

public:

	JTNode();
	~JTNode();
	JTNode(JTree* jtree, int node_id, bool clique_params_flag);
	JTNode(JTree* jtree, int node_id);

	JTNode(JTNode& other);

	void AddFactor(CPT& cpt);
	void CreateCliqueFactor();

	int GetNumFactors() { return factor_list->size(); }
	void SetParent(int parent) { this->parent_id = parent; }
	void SetChildren(vector<int>& children) {this->children = new vector<int>(children); }
	void SetCVars(vector<int>& c_vars) {this->C_vars = new vector<int>(c_vars); }
	void SetSVars(vector<int>& s_vars) {this->S_vars = new vector<int>(s_vars); }
	void SetRVars(vector<int>& r_vars) {this->R_vars = new vector<int>(r_vars); }

	void UpdateFactorWrapper(TensorCPT& clique_marginal);

	void UpdateFactors(TensorCPT& clique_marginal);

	int GetParent() { return parent_id; }
	vector<int>& GetChildren() { return *children; }
	vector<int>& GetCVars() { return *C_vars; }
	vector<int>& GetSVars() {return *S_vars; }
	vector<int>& GetRVars() {return *R_vars; }

	TensorCPT& GetAgglomCliqueFactor();
	void ResetAgglomCliqueFactor() { if (agglomerate_tensor != NULL) { delete(agglomerate_tensor); agglomerate_tensor = NULL; } }

	void UpdateCliqueFactor(TensorCPT& clique_marginal);

	void RecomputeAgglomCliqueFactor();
};


