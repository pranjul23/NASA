# pragma once

#include <vector>
#include "tensor.hpp"
#include "CPT.hpp"

using namespace std;

class TensorCPT {

private:

	Tensor* prob_tensor;

	map<int, vector<int>*>* vars_to_modes;
	vector<int>* modes_to_vars;

	vector<int>* all_vars;
	vector<int>* all_modes;

	bool isInitialized;

public:

	TensorCPT();
	TensorCPT(TensorCPT& other);
	void Initialize(TensorCPT& other);
	TensorCPT(CPT& table, vector<int>& vars, vector<int>& mode_counts);
	void Initialize(Tensor* tensor, vector<int>& modes_to_vars);
	void Initialize(Tensor& tensor, vector<int>& modes_to_vars);
	void Initialize(CPT& table, vector<int>& vars, vector<int>& mode_counts);

	void SetTensor(Tensor& tensor) { if (prob_tensor != NULL) {delete(prob_tensor);} prob_tensor = new Tensor(tensor); }
	Tensor& GetTensor() { return *prob_tensor; }
	int Order() { return prob_tensor->Order(); }

	double At(vector<int>& index_array) {return prob_tensor->At(index_array); }
	double At(int index) {return prob_tensor->At(index); }
	void Set(vector<int>& index_array, double val) { prob_tensor->Set(index_array, val); }
	void Set(int index, double val) { prob_tensor->Set(index, val); }
	void ComputeIndexArray(vector<int>& indices, int n) { prob_tensor->ComputeIndexArray(indices, n); }
	double NumElements() {return prob_tensor->NumElements(); }
	vector<int>& Vars() {return *all_vars;}
	vector<int>& Dims() {return prob_tensor->Dims(); }
	int Dim(int index) {return prob_tensor->Dim(index); }
	int Var(int mode) { return (*modes_to_vars)[mode]; }
	vector<int>& Modes(int var) { return *((*vars_to_modes)[var]); }
	vector<int>& Modes() {return *all_modes; } 

	int NumVars() {return all_vars->size(); }

	void ComputeDiagIndices(vector<int>& diag_indices, vector<int>& vars, vector<int>& indices);

	void Slice(TensorCPT& result_tensor, int fixed_mode, int fixed_index);
	
	static void Add(TensorCPT& cpt_A, TensorCPT& cpt_B);
	static void Multiply(TensorCPT& result_CPT, TensorCPT& cpt_A, TensorCPT& cpt_B);
	static void ElementwiseMultiply(TensorCPT& result_CPT, TensorCPT& cpt_A, TensorCPT& cpt_B);
	static void Divide(TensorCPT& result_CPT, TensorCPT& cpt_A, double val);
	static void OneReduce(TensorCPT& result_CPT, TensorCPT& cpt_A, vector<int>& one_vars);

	void Slice(TensorCPT& result_CPT, vector<int> fixed_vars, vector<int> fixed_vals);
	void Select(TensorCPT& result_CPT, vector<int> fixed_vars, vector<int> fixed_vals);

	void Divide(double val) { prob_tensor->Divide(val); }
	void ElementwiseInvert() { prob_tensor->ElementwiseInvert(); }

	static void ComputeConditionalFactor(TensorCPT& conditional_factor, TensorCPT& factor, vector<int>& front_vars, vector<int>& back_vars);

	double Sum() { return prob_tensor->Sum(); }

	void Rearrange(TensorCPT& result_CPT, vector<int>& new_modes_to_vars);

	void Add(double val) { prob_tensor->Add(val); }

	static void WeightedAdd(TensorCPT& cpt_A, TensorCPT& cpt_B, double weight_A, double weight_B);

	~TensorCPT();
};
