# pragma once

#include <vector>
#include "tensor.hpp"

using namespace std;

class Matrix : public  Tensor {

public:
	Matrix();
	Matrix(vector<int>& dim_array);
	Matrix(Matrix& other);
	~Matrix();

	double At(int row, int col);
	void Set(int row, int col, double val);

	void Row(vector<double>& result_vec, int n);
	void Col(vector<double>& result_vec, int n);
	void RowFind(vector<int>& result_vec, int r);
	void ColFind(vector<int>& result_vec, int c);

	void Sub(vector<double>& result_vec, int r, vector<int>& col_indices);
	void Sub(vector<double>& result_vec, vector<int> row_indices, int c);

	void MatrixConcat(Matrix& concat_mat);

	void NormalizeColumns();
};
