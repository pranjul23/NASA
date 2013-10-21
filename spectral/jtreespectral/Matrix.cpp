# pragma once

# include "Matrix.hpp"
#include <assert.h>
#include <cmath>

Matrix::Matrix()
  :Tensor()
{}

Matrix::Matrix(vector<int>& dims)
  :Tensor(dims)
{}

Matrix::Matrix(Matrix& other)
	:Tensor(other)
{}

Matrix::~Matrix() {}

double Matrix::At(int row, int col)
{
	vector<int> indices(2,0);
	indices[0] = row;
	indices[1] = col;
	assert(row < Dim(0));
	assert(col < Dim(1));
	return Tensor::At(indices);
}

void Matrix::Set(int row, int col, double val)
{
	vector<int> indices(2,0);
	indices[0] = row;
	indices[1] = col;
	assert(row < Dim(0));
	assert(col < Dim(1));
	Tensor::Set(indices, val);
}

void Matrix::Row(vector<double>& result_vec, int n)
{
	for (int i = 0; i < Tensor::Dim(1); ++i)
	{
		result_vec.push_back(At(n, i));
	}
}

void Matrix::Col(vector<double>& result_vec, int n)
{	
	for (int i = 0; i < Tensor::Dim(0); ++i)
	{
		result_vec.push_back(At(i, n));
	}
}

void Matrix::RowFind(vector<int>& result_vec, int r)
{
	for (int i = 0; i < Tensor::Dim(1); ++i)
	{
		if (At(r,i) != 0)
			result_vec.push_back(i);
	}
}

void Matrix::ColFind(vector<int>& result_vec, int c)
{
	for (int i = 0; i < Tensor::Dim(0); ++i)
	{
		if (At(i, c) != 0)
			result_vec.push_back(i);
	}
}

void Matrix::Sub(vector<double>& result_vec, int r, vector<int>& col_indices)
{
	for (int i = 0; i < col_indices.size(); ++i)
	{
		result_vec.push_back(At(r, col_indices[i]));
	}
}

void Matrix::Sub(vector<double>& result_vec, vector<int> row_indices, int c)
{
	for (int i = 0; i < row_indices.size(); ++i)
	{
		result_vec.push_back(At(row_indices[i], c));
	}
}

void Matrix::MatrixConcat(Matrix& concat_mat)
{
	if (this->NumElements() == 0)
	{
		assert(0);
	/*	this->Initialize(concat_mat.Dims());
		for (int i = 0; i < concat_mat.NumElements(); ++i)
		{
			this->tensor_array->at(i) = concat_mat.tensor_array->at(i);
		}*/

	}
	else
	{
		assert(concat_mat.Dim(1) == this->Dim(1)); // assert they have same number of columns
		for (int i = 0; i  < concat_mat.NumElements(); ++i)
		{
			this->tensor_array->push_back(concat_mat.tensor_array->at(i));
		}
		this->dims->at(0) += concat_mat.Dim(0);
	}
}

void Matrix::NormalizeColumns()
{
	for (int i = 0; i < this->Dim(1); ++i)
	{
		double norm = 0;
		for (int j = 0; j < this->Dim(0); ++j)
		{
			norm += pow(At(j,i), 2.0);
		}

		norm = sqrt(norm);

		for (int j = 0; j < this->Dim(0); ++j)
		{
			Set(j,i, At(j,i) / norm);
		}
	}
}