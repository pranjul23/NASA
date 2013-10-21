# pragma once

#include <vector>
#include <map>
#include "tensor.hpp"

using namespace std;

class CPT {

private:

	Tensor* prob_tensor;

	vector<int>* all_modes;
	vector<int>* all_vars;

	map<int, int>* vars_to_modes;
	vector<int>* modes_to_vars;

	vector<int>* parent_vars;
	vector<int>* parent_modes;
	vector<int>* parent_dims;

	vector<int>* front_vars;
	vector<int>* front_modes;
	vector<int>* front_dims;

	bool isInitialized;

public:

	CPT();
	CPT(CPT& other);
	CPT(vector<int>& vars, vector<int>& dims);
	CPT(vector<int>& front_vars, vector<int>& parent_vars, vector<int>& dims);
	void Initialize(vector<int>& vars, vector<int>& dims);
	void Initialize(vector<int>& front_vars, vector<int>& parent_vars, vector<int>& dims);
	void FillWithRandom();

	double At(vector<int>& index_array) {return prob_tensor->At(index_array); }
	double At(int index) {return prob_tensor->At(index); }

	void Set(int index, double val) { prob_tensor->Set(index, val); }
	void Set(vector<int>& index_array, double val) { prob_tensor->Set(index_array, val); }

	int NumElements() {return prob_tensor->NumElements(); }
	vector<int>& Vars() {return *modes_to_vars;}
	vector<int>& Dims() { return prob_tensor->Dims(); }
	int Var(int mode) {return (*modes_to_vars)[mode]; }
	int Mode(int var) {return (*vars_to_modes)[var]; }
	int NumVars() {return all_modes->size(); }
	int VarDim(int var) { return prob_tensor->Dim((*vars_to_modes)[var]); }
	int ModeDim(int mode) { return prob_tensor->Dim(mode); }
	void ComputeIndexArray(vector<int>& indices, int index) { return prob_tensor->ComputeIndexArray(indices, index); } 

	Tensor& GetTensor() {return *prob_tensor; }
	void SetTensor(Tensor& other) { if (prob_tensor != NULL) { delete prob_tensor; } prob_tensor = new Tensor(other); }

	vector<int>& GetFrontVars() {return *front_vars; }
	vector<int>& GetParentVars() { return *parent_vars; }

	void FillWithConst();

	int Sample(vector<int>& parent_vals);

	~CPT();
};
