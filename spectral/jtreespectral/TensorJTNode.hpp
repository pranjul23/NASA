# pragma once

#include <vector>
#include <map>
#include "tensor.hpp"
#include "TensorJTree.hpp"
#include "TensorCPT.hpp"
#include "VectorPlus.hpp"

using namespace std;

enum { NORMAL, INVERSE};

class TensorJTree;

class TensorJTNode {

private:

	TensorJTree* jtree;

	int parent_id;
	int node_id;
	int max_obs_per_mode;
	int max_runs; 

	int S_index;
	int same_sibling_index;

	vector<int>* children;
	vector<int>* same_siblings;

	vector<int>* S_vars;
	vector<int>* R_vars;
	vector<int>* C_vars;
	vector<int>* modes_to_vars;

	map<vector<int>, vector<int>*>* child_seps_to_children;
	map<vector<int>, vector<vector<int>*>*>* child_seps_to_obs_vars;

	map<vector<int>, int>* mode_groups_to_nodes;
	map<vector<int>, int>*	mode_groups_to_types;
	map<vector<int>, vector<int>*>* mode_groups_to_obs;

	Tensor* U;
	TensorCPT* transform_tensor;

	Tensor* U_extra_tensor;
	Tensor* extra_tensor;
	vector<int>* inv_mode_group;

	vector<Tensor*>* extra_tensor_list;
	vector<Tensor*>* forward_tensor_list;

	TensorCPT* B_transform_tensor;

	TensorCPT* final_transform_tensor;

	bool invalid_flag;

public:

	TensorJTNode();
	~TensorJTNode();
	TensorJTNode(TensorJTree* jtree, int node_id, int num_obs_per_mode, int max_runs);
	void SetParent(int parent) { this->parent_id = parent; }
	void SetChildren(vector<int>& children) {this->children = new vector<int>(children); }
	void SetCVars(vector<int>& c_vars) {this->C_vars = new vector<int>(c_vars); }
	void SetSVars(vector<int>& s_vars) {this->S_vars = new vector<int>(s_vars); }
	void SetRVars(vector<int>& r_vars) {this->R_vars = new vector<int>(r_vars); }
	void SetModesToVars(vector<int>& modes_to_vars) {this->modes_to_vars = new vector<int>(modes_to_vars); }
	vector<int>& GetObsVector(vector<int>& modes) { return *((*mode_groups_to_obs)[modes]); }
	vector<int>& GetObsVector(int node);

	int GetParent() { return parent_id; }
	vector<int>& GetChildren() { return *children; }
	vector<int>& GetCVars() { return *C_vars; }
	vector<int>& GetSVars() {return *S_vars; }
	vector<int>& GetRVars() {return *R_vars; }
	map<vector<int>, vector<int>*>& GetChildSeps() {return *child_seps_to_children; }
	void AddObsPartition(vector<int>& sep, vector<int>& descendants, int index);
	vector<int>& GetDescendants(vector<int>& mode_group, int index) { return *((*child_seps_to_obs_vars)[mode_group]->at(index)); }
	void UniqueChildSepAdd(vector<int>& child_sep, int child_id);
	
	void MarkObservableRepresentation();
	void MarkMultModes();
	void ComputeObservableRepresentation();

	int NumSameSiblings() {return same_siblings->size(); }

	int GetSameSiblingIndex();
	int GetSIndex();
	int GetChildSepIndex(vector<int> sep_nodes);

	int GetNumDescendantPartitions(vector<int>& svars) { return ((*child_seps_to_obs_vars)[svars])->size(); }

	Tensor& GetU();
	TensorCPT& GetTransformTensor() { return *transform_tensor; }
	void SetTransformTensor(Tensor& tensor) { transform_tensor->SetTensor(tensor); invalid_flag = false; }
	TensorCPT& GetTransformTensor(TensorCPT& t_cpt, vector<int> evidence_vars, vector<int> evidence_vals);
	Tensor& GetForwardTensor() { return (B_transform_tensor->GetTensor()); }
	Tensor& GetUExtraTensor() { return *U_extra_tensor; }
	Tensor& GetExtraTensor() { return *extra_tensor; }
	int GetNumUniqueChildSeps() { return child_seps_to_children->size(); }
	void FinalizeObservableRepresentation();
	void FinalizeTemplateObservableRepresentation();

	void GetTransformTensor(TensorCPT& t_cpt, vector<int>& evidence_i_vars, vector<int>& evidence_i_vals);

	bool InvModesExists() { return (inv_mode_group->size() > 0); }
	vector<int> GetInvModeGroup() { return *inv_mode_group; }
	vector<int> GetInvModes() {
	
		vector<int> inv_modes;
		VectorPlus::Seq(inv_modes, 0, 1, inv_mode_group->size());
		return inv_modes;
	}

	bool is_invalid() { return invalid_flag; }

	int GetOtherNode(vector<int>& mode_group) { return (*mode_groups_to_nodes)[mode_group]; }

	vector<Tensor*>& GetForwardList() { return *forward_tensor_list; }
	vector<Tensor*>& GetExtraList() { return *extra_tensor_list; }

	void DeleteForwardList() { for (int i = 0; i < forward_tensor_list->size(); ++i) { delete(forward_tensor_list->at(i)); } }
	void DeleteExtraList() { for (int i = 0; i < extra_tensor_list->size(); ++i) { delete(extra_tensor_list->at(i)); } }
};

