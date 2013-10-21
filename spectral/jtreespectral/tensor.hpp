// This is a limited tensor class for the simple times one operation that is needed.

# pragma once

#include <vector>

class Matrix;

using namespace std;

class Tensor {

protected:

	vector<double>* tensor_array;
	vector<int>* dims;
	vector<int>* offsets;
	vector<int>* all_modes;

	bool isInitialized;
	


public:
	Tensor();
	Tensor(Tensor& other);
	Tensor(vector<int>& dim_array);

	void Initialize(vector<int>& dim_array);
	void FillWithRandom();
	void FillWithConst(int val);

	int Order() { return dims->size(); }
	int NumElements() {return tensor_array->size(); }
	vector<int>& Dims() {return *dims;}
	int Dim(int index) {return (*dims)[index]; }

	vector<int>& Modes() { return *all_modes; }

	double At(vector<int>& index_array);
	double At(int index);
	void Set(vector<int>& indices, double val);
	void Set(int index, double val);
	
	void Slice(vector<double>& result_vec, vector<int>& fixed_modes, vector<int>& fixed_indices);
	void Slice(Tensor& result_tensor, int fixed_mode, int fixed_index);
	void Slice(Tensor& result_tensor, vector<int>& fixed_modes, vector<int>& fixed_indices);
	void Select(Tensor& result_tensor, vector<int>& fixed_modes, vector<int>& fixed_indices);
	int ComputeIndex(vector<int>& indices);
	void ComputeIndexArray(vector<int>& indices, int index);

	double Min();
	double Max();

	~Tensor();

	
	// static operations
	static void Add(Tensor& t1, Tensor& t2);
	void Add(double val);
	static void Add(Tensor& result_tensor, Tensor& t1, Tensor& t2);
	static void Multiply(Tensor& result_tensor, Tensor& t1, Tensor& t2, vector<int>& mult_indices1, vector<int>& mult_indices2);
	static void Divide(Tensor& result_tensor, Tensor& t1, double val);

	void Divide(double val);

	static double InnerProduct(Tensor& t1, Tensor& t2);

	bool Inverse(Tensor& Umat, Tensor& result_tensor, vector<int>& mult_modes, int nonzerovals);
	void ElementwiseInvert();

	// static helpers
	static void ComputeOffsets(vector<int>& offsets, vector<int>& dims);
	static void ComputeIndexArray(vector<int>& indices, vector<int>& dims, int index);
	static int ComputeIndex(vector<int>& indices, vector<int>& offsets);

	static void MergeIndices(vector<int>& indices, vector<int>& modes1, vector<int>& modes2, vector<int>& indices1, vector<int>& indices2);
	double ComputePartialSum(vector<int>& fixed_modes, vector<int>& fixed_indices);
	void ComputePartialSum(Tensor& result_tensor, vector<int>& fixed_modes, vector<int>& fixed_indices);
	static void ElementwiseMultiply(Tensor& result_tensor, Tensor& t1, Tensor& t2, vector<int>& mult_indices1, vector<int>& mult_indices2);
	static void Multiply2(Tensor& result_tensor, Tensor& t1, Tensor& t2, vector<int>& mult_modes1, vector<int>& mult_modes2);
	double Sum();

	static void Rearrange(Tensor& result_tensor, Tensor& orig_tensor, vector<int>& old_to_new_modes);

	static void CreateLinearSystem(vector<double>& B_vec, Matrix& A_matrix, 
									Tensor& X, Tensor& A, Tensor& B,
  									vector<int>& mult_modesX, vector<int>& mult_modesA);

	static bool Equals(Tensor& A, Tensor& B);

	static double Diff(Tensor& A, Tensor& B);

	void Reserve(int size_to_reserve) { tensor_array->reserve(size_to_reserve); }

	static void WeightedAdd(Tensor& t1, Tensor& t2, double weight_t1, double weight_t2);

	static bool ComputeJointSVD(Tensor& Umat, vector<Tensor*>& extra_tensor_list, vector<int>& mult_modes, int nonzerovals); 
};
