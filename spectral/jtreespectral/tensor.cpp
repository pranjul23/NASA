# pragma once

#include "tensor.hpp"
#include <assert.h>
#include <algorithm>
#include <math.h>
#include "VectorPlus.hpp"
#include <iostream>
#include "Matrix.hpp"
#include "FastIndexer.hpp"
#include <Eigen/Dense>
#include <Eigen/SVD>

#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>

using namespace std;

// default constructor
Tensor::Tensor()
{
	tensor_array = NULL;
	dims = NULL;
	offsets = NULL;
	all_modes = NULL;

	isInitialized = false;
}

// copy constructor
Tensor::Tensor(Tensor& other)
{
	tensor_array = new vector<double>(*other.tensor_array);
	dims = new vector<int>(*other.dims);
	offsets = new vector<int>(*other.offsets);
	all_modes = new vector<int>(*other.all_modes);

	isInitialized = true;
}


// constructor that initializes tensor
Tensor::Tensor(vector<int>& dims)
{
	tensor_array = NULL;
	this->dims = NULL;
	offsets = NULL;
	all_modes = NULL;

	isInitialized = false;
	Initialize(dims);
}

// destructor
Tensor::~Tensor()
{
	delete(tensor_array);
	delete(dims);
	delete(offsets);
	delete(all_modes);
}

// initialize
void Tensor::Initialize(vector<int>& dims)
{
	assert(!isInitialized);

	if (dims.size() == 0)
	{
		this->dims = new vector<int>();
		offsets = new vector<int>(1, 0);
		all_modes = new vector<int>();
		int total_elements = 1;
		tensor_array = new vector<double>(total_elements, 0);
		isInitialized = true;
	}
	else
	{
		this->dims = new vector<int>(dims);
		offsets = new vector<int>();
		ComputeOffsets(*offsets, dims);

		all_modes = new vector<int>();
		VectorPlus::Seq(*all_modes, 0, 1, dims.size());

		// count total number of elements
		int total_elements = VectorPlus::Product(dims);
		tensor_array = new vector<double>(total_elements, 0);
		isInitialized = true;
	}

}

void Tensor::FillWithRandom()
{
	assert(isInitialized);
	for (int i = 0; i < tensor_array->size(); ++i)
	{
		tensor_array->at(i) = rand();
	}
}

void Tensor::FillWithConst(int val)
{
	assert(isInitialized);
	for (int i = 0; i < tensor_array->size(); ++i)
	{
		tensor_array->at(i) = val;
	}
}

// access tensor element wrapper
double Tensor::At(vector<int>& index_array)
{
	assert(isInitialized);
	int index = ComputeIndex(index_array);
	return At(index);
}

// access tensor element
double Tensor::At(int index)
{
	assert(isInitialized);
	return tensor_array->at(index);
}

// set tensor element wrapper
void Tensor::Set(vector<int>& indices, double val)
{
	assert(isInitialized);
	int index = ComputeIndex(indices);
	Set(index, val);
}

// set tensor element
void Tensor::Set(int index, double val)
{
	assert(isInitialized);
	(*tensor_array)[index] = val;
}


void Tensor::ComputeOffsets(vector<int>& offsets, vector<int>& dims)
{
	assert(offsets.size() == 0);
	offsets.assign(dims.size(), 0);
	// fill in offsets for fast indexing
	for (int i = dims.size() - 1; i >=0; --i)
	{
		if (i == dims.size() - 1)
		{
			offsets[i] = 1;
		}
		else
		{
			offsets[i] = offsets[i+1] * dims[i + 1];
		}
	}
}

double Tensor::Min()
{
	return VectorPlus::Min(*tensor_array);
}

double Tensor::Max()
{
	return VectorPlus::Max(*tensor_array);
}

int Tensor::ComputeIndex(vector<int>& indices)
{
	return ComputeIndex(indices, *offsets);
}

void Tensor::ComputeIndexArray(vector<int>& indices, int index)
{
	ComputeIndexArray(indices, *offsets, index);
}

// add two tensors, store result in result_tensor which is assumed 
// to be uninitialized
void Tensor::Add(Tensor& t1, Tensor& t2)
{
	assert(t1.NumElements() == t2.NumElements());

	for (int i = 0; i < t1.NumElements(); ++i)
	{
		t1.Set(i, t1.At(i) + t2.At(i));
	}
}

// add two tensors, store result in result_tensor which is assumed 
// to be uninitialized
void Tensor::Add(Tensor& result_tensor, Tensor& t1, Tensor& t2)
{
	assert(t1.NumElements() == t2.NumElements());
	result_tensor.Initialize(t1.Dims());

	for (int i = 0; i < result_tensor.NumElements(); ++i)
	{
		result_tensor.Set(i, t1.At(i) + t2.At(i));
	}
}

// inner product among two tensors of the same size
double Tensor::InnerProduct(Tensor& t1, Tensor& t2)
{
	assert(t1.NumElements() == t2.NumElements());
	double val = 0;
	for (int i = 0; i < t1.NumElements(); ++i)
	{
		val += t1.At(i) * t2.At(i);
	}
	return val;
}

void Tensor::ElementwiseMultiply(Tensor& result_tensor, Tensor& t1, Tensor& t2, vector<int>& mult_modes1, vector<int>& mult_modes2)
{

	assert(mult_modes1.size() == mult_modes2.size());


	int numMultElements = 1;	
	vector<int> mult_dims(mult_modes1.size(), 0);

	for (int i = 0; i < mult_modes1.size(); ++i) 
	{
		assert(t1.Dim(mult_modes1[i]) == t2.Dim(mult_modes2[i]));
		mult_dims[i] = t1.Dim(mult_modes1[i]);
		numMultElements = numMultElements * mult_dims[i];
	}
	vector<int> mult_offsets;
	ComputeOffsets(mult_offsets, mult_dims);
	int result_order = t1.Order() + t2.Order() - mult_modes2.size(); 

	if (result_order == 0)
		assert(0);

	vector<int> result_dims;

	vector<int> free_modes1;
	vector<int> free_modes2;


	// find free indices from t1
	for (int i = 0; i < t1.Order(); ++i)
	{
		if (!VectorPlus::Contains(mult_modes1, i))
		{
			free_modes1.push_back(i);
		}
	}

	// find free indices from t2
	for (int i = 0; i < t2.Order(); ++i)
	{
		if (!VectorPlus::Contains(mult_modes2, i))
		{
			free_modes2.push_back(i);
		}
	}

	int curr_index = 0; 
	for (int i = 0; i < result_order; ++i)
	{
		if (i < t1.Order())
		{
			result_dims.push_back(t1.Dim(i));
		}
		else
		{
			result_dims.push_back(t2.Dim(free_modes2[curr_index++]));
		}
	}



	// initialize result_tensor
	result_tensor.Initialize(result_dims);

	// fill in elements from result tensor
	FastIndexer indexer(result_dims);
	int n = 0;
	while (indexer.HasNext())
	{
		vector<int>& indices = indexer.GetNext();
		vector<int> indices1(t1.Order(), 0);
		vector<int> free_indices2;
	//	result_tensor.ComputeIndexArray(indices, n);

		for (int i = 0; i < result_tensor.Order(); ++i)
		{
			if (i < t1.Order())
				indices1[i] = indices[i];
			else
				free_indices2.push_back(indices[i]);
		}

		
		vector<int> mult_indices;
		VectorPlus::Subset(mult_indices, indices1, mult_modes1);


		vector<int> indices2;
		MergeIndices(indices2, mult_modes2, free_modes2, mult_indices, free_indices2);

		double val = 0;
		double val1 = t1.At(indices1);
		if (val1 != 0)
		{
			val = val1 * t2.At(indices2);
		}
		result_tensor.Set(n++, val);
	}
}

/*void Tensor::ElementwiseMultiply(Tensor& result_tensor, Tensor& t1, Tensor& t2, vector<int>& mult_modes1, vector<int>& mult_modes2)
{

	assert(mult_modes1.size() == mult_modes2.size());


	int numMultElements = 1;	
	vector<int> mult_dims(mult_modes1.size(), 0);

	for (int i = 0; i < mult_modes1.size(); ++i) 
	{
		assert(t1.Dim(mult_modes1[i]) == t2.Dim(mult_modes2[i]));
		mult_dims[i] = t1.Dim(mult_modes1[i]);
		numMultElements = numMultElements * mult_dims[i];
	}
	vector<int> mult_offsets;
	ComputeOffsets(mult_offsets, mult_dims);
	int result_order = t1.Order() + t2.Order() - mult_modes2.size(); 

	if (result_order == 0)
		assert(0);

	vector<int> result_dims;

	vector<int> free_modes1;
	vector<int> free_modes2;


	// find free indices from t1
	for (int i = 0; i < t1.Order(); ++i)
	{
		if (!VectorPlus::Contains(mult_modes1, i))
		{
			free_modes1.push_back(i);
		}
	}

	// find free indices from t2
	for (int i = 0; i < t2.Order(); ++i)
	{
		if (!VectorPlus::Contains(mult_modes2, i))
		{
			free_modes2.push_back(i);
		}
	}

	int curr_index = 0; 
	for (int i = 0; i < result_order; ++i)
	{
		if (i < t1.Order())
		{
			result_dims.push_back(t1.Dim(i));
		}
		else
		{
			result_dims.push_back(t2.Dim(free_modes2[curr_index++]));
		}
	}



	// initialize result_tensor
	result_tensor.Initialize(result_dims);

	// fill in elements from result tensor
	for (int n = 0; n < result_tensor.NumElements(); ++n)
	{
		vector<int> indices;
		vector<int> indices1;
		vector<int> free_indices2;
		result_tensor.ComputeIndexArray(indices, n);

		for (int i = 0; i < result_tensor.Order(); ++i)
		{
			if (i < t1.Order())
				indices1.push_back(indices[i]);
			else
				free_indices2.push_back(indices[i]);
		}

		
		vector<int> mult_indices;
		VectorPlus::Subset(mult_indices, indices1, mult_modes1);


		vector<int> indices2;
		MergeIndices(indices2, mult_modes2, free_modes2, mult_indices, free_indices2);

		double val = t1.At(indices1) * t2.At(indices2);
		result_tensor.Set(n, val);
	}
}*/

// multiply two tensors, store result in result_tensor which is assumed to be uninitialized
void Tensor::Multiply2(Tensor& result_tensor, Tensor& t1, Tensor& t2, vector<int>& mult_modes1, vector<int>& mult_modes2)
{
	assert(mult_modes1.size() == mult_modes2.size());

	if (t1.Order() == mult_modes1.size() && t2.Order() == mult_modes2.size())
	{
		double val = InnerProduct(t1, t2);
		vector<int> fake_dims;
		result_tensor.Initialize(fake_dims);
		result_tensor.Set(0, val);
		return;
	}

	int numMultElements = 1;	
	vector<int> mult_dims(mult_modes1.size(), 0);

	for (int i = 0; i < mult_modes1.size(); ++i) 
	{
		assert(t1.Dim(mult_modes1[i]) == t2.Dim(mult_modes2[i]));
		mult_dims[i] = t1.Dim(mult_modes1[i]);
		numMultElements = numMultElements * mult_dims[i];
	}
	vector<int> mult_offsets;
	ComputeOffsets(mult_offsets, mult_dims);
	int result_order = t1.Order() + t2.Order() - mult_modes1.size() - mult_modes2.size(); 

	if (result_order == 0)
		assert(0);

	vector<int> result_dims;

	vector<int> free_modes1;
	vector<int> free_modes2;


	// find free indices from t1
	for (int i = 0; i < t1.Order(); ++i)
	{
		if (!VectorPlus::Contains(mult_modes1, i))
		{
			free_modes1.push_back(i);
		}
	}

	// find free indices from t2
	for (int i = 0; i < t2.Order(); ++i)
	{
		if (!VectorPlus::Contains(mult_modes2, i))
		{
			free_modes2.push_back(i);
		}
	}

	int curr_index = 0; 
	for (int i = 0; i < result_order; ++i)
	{
		if (VectorPlus::Contains(mult_modes1, i))
		{
			result_dims.push_back(t2.Dim(free_modes2[curr_index++]));
		}
		else
		{
			result_dims.push_back(t1.Dim(i));
		}
	}



	// initialize result_tensor
	result_tensor.Initialize(result_dims);

	FastIndexer result_indexer(result_dims);
	// fill in elements from result tensor
	for (int n = 0; n < result_tensor.NumElements(); ++n)
	{
		
		vector<int> free_indices1;
		vector<int> free_indices2;

//		result_tensor.ComputeIndexArray(indices, n);

		vector<int>& indices = result_indexer.GetNext();

		for (int i = 0; i < result_tensor.Order(); ++i)
		{
			if (!VectorPlus::Contains(mult_modes1, i))
				free_indices1.push_back(indices[i]);
			else
				free_indices2.push_back(indices[i]);
		}

		vector<int> indices1(mult_modes1.size() + free_modes1.size(), 0); 
		vector<int> indices2(mult_modes2.size() + free_modes2.size(), 0);

		for (int i = 0; i < free_modes1.size(); ++i)
		{
			indices1[free_modes1[i]] = free_indices1[i];
		}

		for (int i = 0; i < free_modes2.size(); ++i)
		{
			indices2[free_modes2[i]] = free_indices2[i];
		}

		// sum over elementwise products of mult-mode elements
		double temp_sum = 0;
		FastIndexer indexer(mult_dims);
		while (indexer.HasNext())
		{
			vector<int>& mult_indices = indexer.GetNext();
			//ComputeIndexArray(mult_indices, mult_offsets, k);

			for (int i = 0; i < mult_indices.size(); ++i)
			{
				indices1[mult_modes1[i]] = mult_indices[i];
				indices2[mult_modes2[i]] = mult_indices[i];
			}


		//	MergeIndices(indices1, mult_modes1, free_modes1, mult_indices, free_indices1);
		//	MergeIndices(indices2, mult_modes2, free_modes2, mult_indices, free_indices2);

			double val1 = t1.At(indices1);
			if (val1 == 0)
				temp_sum += 0;
			else
				temp_sum += val1 * t2.At(indices2);
		}

		result_tensor.Set(n, temp_sum);
	}
}




// multiply two tensors, store result in result_tensor which is assumed to be uninitialized
void Tensor::Multiply(Tensor& result_tensor, Tensor& t1, Tensor& t2, vector<int>& mult_modes1, vector<int>& mult_modes2)
{
	assert(mult_modes1.size() == mult_modes2.size());

	if (t1.Order() == mult_modes1.size() && t2.Order() == mult_modes2.size())
	{
		double val = InnerProduct(t1, t2);
		vector<int> fake_dims;
		result_tensor.Initialize(fake_dims);
		result_tensor.Set(0, val);
		return;
	}

	int numMultElements = 1;	
	vector<int> mult_dims(mult_modes1.size(), 0);

	for (int i = 0; i < mult_modes1.size(); ++i) 
	{
		assert(t1.Dim(mult_modes1[i]) == t2.Dim(mult_modes2[i]));
		mult_dims[i] = t1.Dim(mult_modes1[i]);
		numMultElements = numMultElements * mult_dims[i];
	}
	vector<int> mult_offsets;
	ComputeOffsets(mult_offsets, mult_dims);
	int result_order = t1.Order() + t2.Order() - mult_modes1.size() - mult_modes2.size(); 

	if (result_order == 0)
		assert(0);

	vector<int> result_dims;

	vector<int> free_modes1;
	vector<int> free_modes2;


	// find free indices from t1
	for (int i = 0; i < t1.Order(); ++i)
	{
		if (!VectorPlus::Contains(mult_modes1, i))
		{
			result_dims.push_back(t1.Dim(i));
			free_modes1.push_back(i);
		}
	}

	// find free indices from t2
	for (int i = 0; i < t2.Order(); ++i)
	{
		if (!VectorPlus::Contains(mult_modes2, i))
		{
			result_dims.push_back(t2.Dim(i));
			free_modes2.push_back(i);
		}
	}

	// initialize result_tensor
	result_tensor.Initialize(result_dims);

	// fill in elements from result tensor

	FastIndexer result_indexer(result_dims);

	for (int n = 0; n < result_tensor.NumElements(); ++n)
	{
		vector<int>& indices = result_indexer.GetNext();
		vector<int> free_indices1;
		vector<int> free_indices2;
	//	result_tensor.ComputeIndexArray(indices, n);

		for (int i = 0; i < result_tensor.Order(); ++i)
		{
			if (i < free_modes1.size())
				free_indices1.push_back(indices[i]);
			else
				free_indices2.push_back(indices[i]);
		}

		// sum over elementwise products of mult-mode elements
		double temp_sum = 0;
		FastIndexer mult_indexer(mult_dims);
		for (int k = 0; k < numMultElements; ++k)
		{
			vector<int>& mult_indices = mult_indexer.GetNext();
		//	ComputeIndexArray(mult_indices, mult_offsets, k);

			vector<int> indices1; 
			vector<int> indices2;

			MergeIndices(indices1, mult_modes1, free_modes1, mult_indices, free_indices1);
			MergeIndices(indices2, mult_modes2, free_modes2, mult_indices, free_indices2);

			temp_sum += t1.At(indices1) * t2.At(indices2);
		}

		result_tensor.Set(n, temp_sum);
	}
}

// bool Tensor::Inverse(Tensor& Umat, Tensor& result_tensor, vector<int>& mult_modes, int nonzerovals)
//{	
//	int num_sing_vals = (int)pow((double)nonzerovals, (double)mult_modes.size());
//
//	vector<int> free_modes;
//	vector<int> mult_dims;
//	vector<int> free_dims;
//	vector<int> mult_offsets;
//	vector<int> free_offsets;
//
//	VectorPlus::SetDiff(free_modes, *all_modes, mult_modes);  
//	VectorPlus::Subset(mult_dims, *dims, mult_modes);
//	VectorPlus::CSubset(free_dims, *dims, mult_modes);
//	ComputeOffsets(mult_offsets, mult_dims);
//	ComputeOffsets(free_offsets, free_dims);
//
//	vector<int> usmalldims(free_modes.size(), nonzerovals);
//	vector<int> udims;
//	vector<int> usmall_offsets;
//	ComputeOffsets(usmall_offsets, usmalldims);
//
//	int numMultElements = VectorPlus::Product(mult_dims);
//	int numFreeElements = VectorPlus::Product(free_dims);
//
//	assert(numMultElements == numFreeElements);
//
//    Eigen::MatrixXd matricized_tensor(numFreeElements,numMultElements);
//
//	cout << "copy start 1\n";
//	for (int i = 0; i < numFreeElements; ++i)
//	{
//		vector<int> free_indices;
//		ComputeIndexArray(free_indices, free_offsets, i);
//		for (int j = 0; j < numMultElements; ++j)
//		{
//			vector<int> mult_indices;
//			ComputeIndexArray(mult_indices, mult_offsets, j);
//			vector<int> total_indices;
//			VectorPlus::Concat(total_indices, free_indices, mult_indices);
//			matricized_tensor(i,j) = this->At(total_indices);
//		}
//	}
//	cout << "copy end 1\n";	
////	MatrixXd matricized_inverse = matricized_tensor.inverse();
////	cout << matricized_inverse;
////	cout << "\n";
//	//compute pseudoinverse
//	cout << "svd start 1\n";	
//		Eigen::JacobiSVD<Eigen::MatrixXd> svd(matricized_tensor, Eigen::ComputeFullU);
//		Eigen::MatrixXd U = svd.matrixU();
//		cout << "svd end 1\n";	
//		Eigen::MatrixXd thinU = U.leftCols(num_sing_vals);
//		
//		Eigen::MatrixXd reduced_matrix = thinU.transpose() * matricized_tensor;
//
//		Eigen::JacobiSVD<Eigen::MatrixXd> svd2(reduced_matrix, Eigen::ComputeFullU | Eigen::ComputeFullV);
//		U = svd2.matrixU();
//
//
//		Eigen::MatrixXd Uleft = U.leftCols(num_sing_vals);
//		Eigen::MatrixXd V = svd2.matrixV();
//		Eigen::MatrixXd Vleft = V.leftCols(num_sing_vals);
//		Eigen::VectorXd s_vals = svd2.singularValues();
//
//		Eigen::MatrixXd Sigma(num_sing_vals, num_sing_vals);
//		for (int i = 0; i < num_sing_vals; ++i)
//		{
//			for (int j = 0; j < num_sing_vals; ++j)
//			{
//				Sigma(i,j) = 0;
//			}
//		}
//
//		for (int s = 0; s < num_sing_vals;  ++s)
//		{
//			cout << s_vals(s);
//			cout << "\n";
//			if (abs(s_vals(s)) < .00000000001)
//				assert(0);
//			Sigma(s,s) = 1.0 / s_vals(s);
//			
//		}
//
//		Eigen::MatrixXd pinverse_temp = Vleft * Sigma;
//	 Eigen::MatrixXd matricized_inverse = pinverse_temp * Uleft.transpose();
//	 
//	
//
//		// rearrange into tensor
//	 	vector<int> result_dims;
//		VectorPlus::Concat(result_dims, free_dims, usmalldims);
//		result_tensor.Initialize(result_dims);
//
//
//		for (int i = 0; i < matricized_inverse.rows(); ++i)
//		{
//			vector<int> left_indices;
//			ComputeIndexArray(left_indices, free_offsets, i);
//
//			for (int j = 0; j < matricized_inverse.cols(); ++j)
//			{
//				vector<int> right_indices;
//				ComputeIndexArray(right_indices, usmall_offsets, j);
//				vector<int> total_indices;
//
//				VectorPlus::Concat(total_indices, left_indices, right_indices);
//				result_tensor.Set(total_indices, matricized_inverse(i, j));
//			}
//		}
//
//	
//		VectorPlus::Concat(udims, free_dims, usmalldims);
//		Umat.Initialize(udims);
//
//		
//
//		vector<int> semi_dims;
//		semi_dims.push_back(thinU.rows());
//		semi_dims.push_back(thinU.cols());
//
//
//		/*boost::mt19937 rng; // I don't seed it on purpouse (it's not relevant)
//
//		boost::normal_distribution<> nd(0.0, 1.0);
//
//		boost::variate_generator<boost::mt19937&, 
//                           boost::normal_distribution<> > var_nor(rng, nd);
//
//
//		Matrix rand_matrix(semi_dims);
//		for (int i = 0; i < rand_matrix.Dim(0); ++i)
//		{
//			for (int j = 0; j < rand_matrix.Dim(1); ++j)
//			{
//				double x = var_nor();
//				rand_matrix.Set(i, j ,x);
//			}
//		}
//
//		rand_matrix.NormalizeColumns();*/
//		cout << "copy start 2 \n";	
//		for (int i = 0; i < thinU.rows(); ++i)
//		{
//			vector<int> left_indices;
//			ComputeIndexArray(left_indices, free_offsets, i);
//
//			for (int j = 0; j < thinU.cols(); ++j)
//			{
//				vector<int> right_indices;
//				ComputeIndexArray(right_indices, usmall_offsets, j);
//				vector<int> indices;
//				VectorPlus::Concat(indices, left_indices, right_indices);
//			//	Umat.Set(indices, rand_matrix.At(i,j));
//				Umat.Set(indices, thinU(i,j));
//
//				if (thinU.rows() == thinU.cols())
//				{
//					if (VectorPlus::Equals(left_indices, right_indices))
//					{	
//						Umat.Set(indices, 1);
//					}
//					else
//					{
//						Umat.Set(indices, 0);
//					}
//				}
//			}
//		}
//		cout << "copy end 2 \n";	
//	return true;
//}


bool Tensor::Inverse(Tensor& Umat, Tensor& result_tensor, vector<int>& mult_modes, int nonzerovals)
{	
	int num_sing_vals = (int)pow((double)nonzerovals, (double)mult_modes.size());

	vector<int> free_modes;
	vector<int> mult_dims;
	vector<int> free_dims;
	vector<int> mult_offsets;
	vector<int> free_offsets;

	VectorPlus::SetDiff(free_modes, *all_modes, mult_modes);  
	VectorPlus::Subset(mult_dims, *dims, mult_modes);
	VectorPlus::CSubset(free_dims, *dims, mult_modes);
	ComputeOffsets(mult_offsets, mult_dims);
	ComputeOffsets(free_offsets, free_dims);

	vector<int> usmalldims(free_modes.size(), nonzerovals);
	vector<int> udims;
	vector<int> usmall_offsets;
	ComputeOffsets(usmall_offsets, usmalldims);

	VectorPlus::Concat(udims, free_dims, usmalldims);
	Umat.Initialize(udims);

	vector<int> result_dims;
		VectorPlus::Concat(result_dims, free_dims, usmalldims);
		result_tensor.Initialize(result_dims);

/*	int numMultElements = VectorPlus::Product(mult_dims);
	int numFreeElements = VectorPlus::Product(free_dims);

	assert(numMultElements == numFreeElements);

    Eigen::MatrixXd matricized_tensor(numFreeElements,numMultElements);

	cout << "copy start 1\n";
	for (int i = 0; i < numFreeElements; ++i)
	{
		vector<int> free_indices;
		ComputeIndexArray(free_indices, free_offsets, i);
		for (int j = 0; j < numMultElements; ++j)
		{
			vector<int> mult_indices;
			ComputeIndexArray(mult_indices, mult_offsets, j);
			vector<int> total_indices;
			VectorPlus::Concat(total_indices, free_indices, mult_indices);
			matricized_tensor(i,j) = this->At(total_indices);
		}
	}
	cout << "copy end 1\n";	
//	MatrixXd matricized_inverse = matricized_tensor.inverse();
//	cout << matricized_inverse;
//	cout << "\n";
	//compute pseudoinverse
	cout << "svd start 1\n";	
		Eigen::JacobiSVD<Eigen::MatrixXd> svd(matricized_tensor, Eigen::ComputeFullU);
		Eigen::MatrixXd U = svd.matrixU();
		cout << "svd end 1\n";	
		Eigen::MatrixXd thinU = U.leftCols(num_sing_vals);
		
	/*	MatrixXd reduced_matrix = thinU.transpose() * matricized_tensor;

		JacobiSVD<MatrixXd> svd2(reduced_matrix, ComputeFullU | ComputeFullV);
		U = svd2.matrixU();


		MatrixXd Uleft = U.leftCols(num_sing_vals);
		MatrixXd V = svd2.matrixV();
		MatrixXd Vleft = V.leftCols(num_sing_vals);
		VectorXd s_vals = svd2.singularValues();

		MatrixXd Sigma(num_sing_vals, num_sing_vals);
		for (int i = 0; i < num_sing_vals; ++i)
		{
			for (int j = 0; j < num_sing_vals; ++j)
			{
				Sigma(i,j) = 0;
			}
		}

		for (int s = 0; s < num_sing_vals;  ++s)
		{
			cout << s_vals(s);
			cout << "\n";
			if (abs(s_vals(s)) < .00000000001)
				assert(0);
			Sigma(s,s) = 1.0 / s_vals(s);
			
		}

		MatrixXd pinverse_temp = Vleft * Sigma;
	 MatrixXd matricized_inverse = pinverse_temp * Uleft.transpose();
	 */
	

		// rearrange into tensor
/*	 	vector<int> result_dims;
		VectorPlus::Concat(result_dims, free_dims, usmalldims);
		result_tensor.Initialize(result_dims);


	/*	for (int i = 0; i < matricized_inverse.rows(); ++i)
		{
			vector<int> left_indices;
			ComputeIndexArray(left_indices, free_offsets, i);

			for (int j = 0; j < matricized_inverse.cols(); ++j)
			{
				vector<int> right_indices;
				ComputeIndexArray(right_indices, usmall_offsets, j);
				vector<int> total_indices;

				VectorPlus::Concat(total_indices, left_indices, right_indices);
				result_tensor.Set(total_indices, matricized_inverse(i, j));
			}
		} */

	
	/*	VectorPlus::Concat(udims, free_dims, usmalldims);
		Umat.Initialize(udims);

		

		vector<int> semi_dims;
		semi_dims.push_back(thinU.rows());
		semi_dims.push_back(thinU.cols());

		cout << "copy start 2 \n";	
		for (int i = 0; i < thinU.rows(); ++i)
		{
			vector<int> left_indices;
			ComputeIndexArray(left_indices, free_offsets, i);

			for (int j = 0; j < thinU.cols(); ++j)
			{
				vector<int> right_indices;
				ComputeIndexArray(right_indices, usmall_offsets, j);
				vector<int> indices;
				VectorPlus::Concat(indices, left_indices, right_indices);
			//	Umat.Set(indices, rand_matrix.At(i,j));
				Umat.Set(indices, thinU(i,j));

				if (thinU.rows() == thinU.cols())
				{
					if (VectorPlus::Equals(left_indices, right_indices))
					{	
						Umat.Set(indices, 1);
					}
					else
					{
						Umat.Set(indices, 0);
					}
				}
			}
		}
		cout << "copy end 2 \n"; */
	return true;
}



int Tensor::ComputeIndex(vector<int>& indices, vector<int>& offsets)
{
	int index = 0;
	for (int i = 0; i < offsets.size(); ++i)
	{
		index = index + indices[i] * offsets[i];
	}
	return index;
}

void Tensor::ComputeIndexArray(vector<int>& indices, vector<int>& offsets, int index)
{
	int curr_val = index;
	assert(indices.size() == 0);
	indices.reserve(offsets.size());
	for (int i = 0; i < offsets.size(); ++i)
	{
		indices.push_back((int)(floor(((double)curr_val) / offsets[i])));
		curr_val = curr_val % offsets[i];
	}
}

// Merges Indices
void Tensor::MergeIndices(vector<int>& indices, vector<int>& modes1, vector<int>& modes2, vector<int>& indices1, vector<int>& indices2)
{
	indices.reserve(modes1.size() + modes2.size());
	for (int i = 0; i < modes1.size() + modes2.size(); ++i)
	{
		int first_index = VectorPlus::Find(modes1, i);
		if (first_index >= 0)
			indices.push_back(indices1[first_index]);
		else
		{
			int second_index = VectorPlus::Find(modes2, i);
			indices.push_back(indices2[second_index]);
		}
	}
}

// Computes sum over tensor where certain modes are held fixed
double Tensor::ComputePartialSum(vector<int>& fixed_modes, vector<int>& fixed_indices)
{
	vector<int> free_modes;
	VectorPlus::SetDiff(free_modes, *all_modes, fixed_modes);
	
	vector<int> free_offsets;

	vector<int> free_dims;
	VectorPlus::Subset(free_dims, *dims, free_modes);
	ComputeOffsets(free_offsets, free_dims);

	int num_vals = VectorPlus::Product(free_dims);

	double sum = 0;
	FastIndexer indexer(free_dims);
	while (indexer.HasNext())
//	for (int i = 0; i < num_vals; ++i)
	{
		vector<int>& free_indices = indexer.GetNext();
	//	vector<int> free_indices;
	//	ComputeIndexArray(free_indices, free_offsets, i);
		vector<int> indices;
		MergeIndices(indices, fixed_modes, free_modes, fixed_indices, free_indices);
		sum += At(indices);
	}
	return sum;
}

void Tensor::Slice(Tensor& result_tensor, vector<int>& fixed_modes, vector<int>& fixed_indices)
{
	vector<int> free_modes;
	vector<int> free_dims;
	VectorPlus::SetDiff(free_modes, *all_modes, fixed_modes);
	VectorPlus::Subset(free_dims, *dims, free_modes);
	result_tensor.Initialize(free_dims);
	int i = 0;
	FastIndexer indexer(free_dims);
	while (indexer.HasNext())
//	for (int i = 0; i < result_tensor.NumElements(); ++i)
	{
	//	vector<int> result_indices;
	//	result_tensor.ComputeIndexArray(result_indices, i);
		vector<int>& result_indices = indexer.GetNext();
		vector<int> orig_indices;
		MergeIndices(orig_indices, free_modes, fixed_modes, result_indices, fixed_indices);
		result_tensor.Set(i++, this->At(orig_indices));
	}
}

void Tensor::Slice(vector<double>& result_vec, vector<int>& fixed_modes, vector<int>& fixed_indices)
{
	assert(fixed_modes.size() == Order() - 1);
	assert(result_vec.size() == 0);
	vector<int> free_modes;
	VectorPlus::SetDiff(free_modes, *all_modes, fixed_modes);
	assert(free_modes.size() == 1);

	for (int i = 0; i < dims->at(free_modes[0]); ++i)
	{
		vector<int> indices;
		vector<int> free_index = VectorPlus::CreateSingleton(i);
		MergeIndices(indices, fixed_modes, free_modes, fixed_indices, free_index);
		result_vec.push_back(At(indices));
	}
}

void Tensor::Slice(Tensor& result_tensor, int fixed_mode, int fixed_index)
{
	vector<int> result_dims;
	vector<int> free_modes;


	vector<int> fixed_mode_vec  = VectorPlus::CreateSingleton(fixed_mode);
	vector<int> fixed_index_vec = VectorPlus::CreateSingleton(fixed_index);

	VectorPlus::SetDiff(free_modes, *all_modes,  fixed_mode_vec);
	VectorPlus::Subset(result_dims, *dims, free_modes);
	
	result_tensor.Initialize(result_dims);

	if (result_tensor.NumElements() == 1)
	{
		assert(fixed_mode == 1);
		result_tensor.Set(0, this->At(fixed_index));
		return;
	}
	FastIndexer indexer(result_dims);

	int i = 0;
	while (indexer.HasNext())
//	for (int i = 0; i < result_tensor.NumElements(); ++i)
	{
		vector<int>& indices = indexer.GetNext();
	//	vector<int> indices;
		vector<int> total_indices;
	//	result_tensor.ComputeIndexArray(indices, i); 
		MergeIndices(total_indices, free_modes, fixed_mode_vec, indices, fixed_index_vec); 

		result_tensor.Set(i++, this->At(total_indices));
	}
}

void Tensor::Select(Tensor& result_tensor, vector<int>& fixed_modes, vector<int>& fixed_indices)
{
	result_tensor.Initialize(*dims);

	int i = 0;
	FastIndexer indexer(*dims);
	while (indexer.HasNext())
//	for (int i = 0; i < result_tensor.NumElements(); ++i)
	{
		vector<int>& i_indices = indexer.GetNext();
	//	vector<int> i_indices;
	//	result_tensor.ComputeIndexArray(i_indices, i);
		vector<int> fixed_i_indices;
		VectorPlus::Subset(fixed_i_indices, i_indices, fixed_modes);
		if (VectorPlus::Equals(fixed_i_indices, fixed_indices))
			result_tensor.Set(i, this->At(i));
		else
			result_tensor.Set(i, 0);
		i++;
	}

}

double Tensor::Sum()
{
	return VectorPlus::Sum(*tensor_array);
}

void Tensor::Divide(Tensor& result_tensor, Tensor& t1, double val)
{
	result_tensor.Initialize(t1.Dims());

	for (int i = 0; i < t1.NumElements(); ++i)
	{
		result_tensor.Set(i, t1.At(i) / val);
	}
}

void Tensor::Divide(double val)
{
	for (int i = 0; i < NumElements(); ++i)
	{
		Set(i, At(i) / val);
	}
}

void Tensor::ElementwiseInvert()
{
	for (int i = 0; i < NumElements(); ++i)
	{
		Set(i, 1.0 / At(i));
	}
}

void Tensor::Rearrange(Tensor& result_tensor, Tensor& orig_tensor, vector<int>& old_to_new_modes)
{
	vector<int> new_dims;
	VectorPlus::Rearrange(new_dims, orig_tensor.Dims(), old_to_new_modes);

	result_tensor.Initialize(new_dims);

	int i = 0;
	FastIndexer indexer(orig_tensor.Dims());
	while (indexer.HasNext())
	//for (int i = 0; i < orig_tensor.NumElements(); ++i)
	{
		vector<int>& orig_indices = indexer.GetNext();
//		vector<int> orig_indices;
//		orig_tensor.ComputeIndexArray(orig_indices, i);
		vector<int> new_indices;
		VectorPlus::Rearrange(new_indices, orig_indices, old_to_new_modes);

		result_tensor.Set(new_indices, orig_tensor.At(i++));
	}
}

void Tensor::CreateLinearSystem(vector<double>& B_vec, Matrix& A_matrix, 
						Tensor& X, Tensor& A, Tensor& B,
						vector<int>& mult_modesX, vector<int>& mult_modesA)
{
	// fake multiply x and A together to create B, creating the linear system in the process

	assert(mult_modesX.size() == mult_modesA.size());

	if (X.Order() == mult_modesX.size() && A.Order() == mult_modesA.size())
	{
		assert(0);
	}

	int numMultElements = 1;	
	vector<int> mult_dims(mult_modesX.size(), 0);

	for (int i = 0; i < mult_modesX.size(); ++i) 
	{
		assert(X.Dim(mult_modesX[i]) == A.Dim(mult_modesA[i]));
		mult_dims[i] = X.Dim(mult_modesX[i]);
		numMultElements = numMultElements * mult_dims[i];
	}
	vector<int> mult_offsets;
	ComputeOffsets(mult_offsets, mult_dims);
	int result_order = X.Order() + A.Order() - mult_modesX.size() - mult_modesA.size(); 

	if (result_order == 0)
		assert(0);

	vector<int> result_dims;

	vector<int> free_modesX;
	vector<int> free_modesA;


	// find free indices from X
	for (int i = 0; i < X.Order(); ++i)
	{
		if (!VectorPlus::Contains(mult_modesX, i))
		{
			free_modesX.push_back(i);
		}
	}

	// find free indices from A
	for (int i = 0; i < A.Order(); ++i)
	{
		if (!VectorPlus::Contains(mult_modesA, i))
		{
			free_modesA.push_back(i);
		}
	}
	vector<int> a_mat_dims = VectorPlus::CreatePair(B.NumElements(), X.NumElements());
	A_matrix.Initialize(a_mat_dims);
	B_vec.reserve(B.NumElements());
	// fill in elements from result tensor

	FastIndexer B_indexer(B.Dims());

	for (int n = 0; n < B.NumElements(); ++n)
	{
		B_vec.push_back(B.At(n));

		vector<int>& indices = B_indexer.GetNext();
		vector<int> free_indicesX;
		vector<int> free_indicesA;
	//	B.ComputeIndexArray(indices, n);

		for (int i = 0; i < B.Order(); ++i)
		{
			if (!VectorPlus::Contains(mult_modesX, i))
				free_indicesX.push_back(indices[i]);
			else
				free_indicesA.push_back(indices[i]);
		}

		// sum over elementwise products of mult-mode elements
		double temp_sum = 0;
		FastIndexer mult_indexer(mult_dims);
		for (int k = 0; k < numMultElements; ++k)
		{
			vector<int>& mult_indices = mult_indexer.GetNext();
		//	ComputeIndexArray(mult_indices, mult_offsets, k);

			vector<int> indicesX; 
			vector<int> indicesA;

			MergeIndices(indicesX, mult_modesX, free_modesX, mult_indices, free_indicesX);
			MergeIndices(indicesA, mult_modesA, free_modesA, mult_indices, free_indicesA);

			
			A_matrix.Set(n, X.ComputeIndex(indicesX), A.At(indicesA));
		}
	}
}

bool Tensor::Equals(Tensor& A, Tensor& B)
{
	assert(A.NumElements() == B.NumElements());

	for (int i = 0; A.NumElements(); ++i)
	{
		if (A.At(i) != B.At(i))
			return false;
	}
	return true;
}

double Tensor::Diff(Tensor& A, Tensor& B)
{
	assert(A.NumElements() == B.NumElements());

	double diff = 0;
	for (int i = 0; i < A.NumElements(); ++i)
	{
		diff += abs(A.At(i) - B.At(i));
	}
	return diff;
}

void Tensor::Add(double val)
{
	for (int i = 0; i < this->NumElements(); ++i)
	{
		this->tensor_array->at(i) += val;
	}
}

void Tensor::WeightedAdd(Tensor& t1, Tensor& t2, double weight_t1, double weight_t2)
{
	assert(t1.NumElements() == t2.NumElements());

	for (int i = 0; i < t1.NumElements(); ++i)
	{
		t1.Set(i, weight_t1 * t1.At(i) + weight_t2 * t2.At(i));
	}
}

bool Tensor::ComputeJointSVD(Tensor& Umat, vector<Tensor*>& extra_tensor_list, vector<int>& mult_modes, int nonzerovals)
{
	int num_sing_vals = (int)pow((double)nonzerovals, (double)mult_modes.size());

	vector<int> free_modes;
	vector<int> mult_dims;
	vector<int> free_dims;
	vector<int> mult_offsets;
	vector<int> free_offsets;

	Tensor& temp_tensor = *(extra_tensor_list.at(0));
	VectorPlus::SetDiff(free_modes, temp_tensor.Modes(), mult_modes);  
	VectorPlus::Subset(mult_dims, temp_tensor.Dims(), mult_modes);
	VectorPlus::CSubset(free_dims, temp_tensor.Dims(), mult_modes);
	ComputeOffsets(mult_offsets, mult_dims);
	ComputeOffsets(free_offsets, free_dims);

	vector<int> usmalldims(free_modes.size(), nonzerovals);
	vector<int> udims;
	vector<int> usmall_offsets;
	ComputeOffsets(usmall_offsets, usmalldims);

	int numMultElements = VectorPlus::Product(mult_dims);
	int numFreeElements = VectorPlus::Product(free_dims);

	assert(numMultElements == numFreeElements);

    Eigen::MatrixXd matricized_tensor(numFreeElements,extra_tensor_list.size() * numMultElements);

	//cout << "copy start 1\n";
	for (int z = 0; z < extra_tensor_list.size(); ++z)
	{
		int z_offset = z * numMultElements;
		FastIndexer i_indexer(free_dims);
		for (int i = 0; i < numFreeElements; ++i)
		{
			vector<int>& free_indices = i_indexer.GetNext();
		//	ComputeIndexArray(free_indices, free_offsets, i);
			FastIndexer j_indexer(mult_dims);
			for (int j = 0; j < numMultElements; ++j)
			{
				vector<int>& mult_indices = j_indexer.GetNext();
		//		ComputeIndexArray(mult_indices, mult_offsets, j);
				vector<int> total_indices;
				VectorPlus::Concat(total_indices, free_indices, mult_indices);
				matricized_tensor(i,z_offset + j) = extra_tensor_list.at(z)->At(total_indices);
			}
		}
	}
//	cout << "copy end 1\n";	
//	MatrixXd matricized_inverse = matricized_tensor.inverse();
//	cout << matricized_inverse;
//	cout << "\n";
	//compute pseudoinverse
	//cout << "svd start 1\n";	
		Eigen::JacobiSVD<Eigen::MatrixXd> svd(matricized_tensor, Eigen::ComputeFullU);
		Eigen::MatrixXd U = svd.matrixU();
	//	cout << "svd end 1\n";	
		Eigen::MatrixXd thinU = U.leftCols(num_sing_vals);
		
	
		VectorPlus::Concat(udims, free_dims, usmalldims);
		Umat.Initialize(udims);

		

		vector<int> semi_dims;
		semi_dims.push_back(thinU.rows());
		semi_dims.push_back(thinU.cols());

	//	cout << "copy start 2 \n";	
		FastIndexer i_indexer(free_dims);
		for (int i = 0; i < thinU.rows(); ++i)
		{
			vector<int>& left_indices = i_indexer.GetNext();
		//	ComputeIndexArray(left_indices, free_offsets, i);

			FastIndexer j_indexer(usmalldims);
			for (int j = 0; j < thinU.cols(); ++j)
			{
				vector<int>& right_indices = j_indexer.GetNext();
			//	ComputeIndexArray(right_indices, usmall_offsets, j);
				vector<int> indices;
				VectorPlus::Concat(indices, left_indices, right_indices);
			//	Umat.Set(indices, rand_matrix.At(i,j));
				Umat.Set(indices, thinU(i,j));

				if (thinU.rows() == thinU.cols())
				{
					if (VectorPlus::Equals(left_indices, right_indices))
					{	
						Umat.Set(indices, 1);
					}
					else
					{
						Umat.Set(indices, 0);
					}
				}
			}
		}
	//	cout << "copy end 2 \n";	
	return true;
}