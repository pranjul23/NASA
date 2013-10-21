# pragma once

#include <vector>
#include <algorithm>
#include "VectorPlus.hpp"
#include <assert.h>
#include <set>
#include <limits>

using namespace std;

bool VectorPlus::IsSubset(vector<int>& vec_A, vector<int>& vec_B)
{
	for (int i = 0; i < vec_B.size(); ++i)
	{
		if (!Contains(vec_A, vec_B[i]))
			return false;
	}

	return true;

}

bool VectorPlus::Equals(vector<int>& vec_A, vector<int>& vec_B)
{
	if (vec_A.size() != vec_B.size())
		return false;

	for (int i = 0; i < vec_A.size(); ++i)
	{
		if (vec_A[i] != vec_B[i])
			return false;
	}

	return true;

}

void VectorPlus::Add(vector<int>& sum_vec, vector<int>& vec_A, vector<int>& vec_B)
{
	assert(sum_vec.size() == 0);
	assert(vec_A.size() == vec_B.size());
	sum_vec.reserve(vec_A.size());
	for (int i = 0; i < vec_A.size(); ++i)
	{
		sum_vec.push_back(vec_A[i] + vec_B[i]);
	}
}

void VectorPlus::MatchSub(vector<int>& result_vec, vector<int> vec, vector<int> match_vec, vector<int> sub_match_vec)
{
	result_vec.reserve(sub_match_vec.size());
	for (int i = 0; i < sub_match_vec.size(); ++i)
	{
		int index = Find(match_vec, sub_match_vec[i]);
		result_vec.push_back(vec[index]);
	}
}


bool VectorPlus::Contains(vector<int>& vec, int element)
{
	vector<int>::iterator iter = find(vec.begin(), vec.end(), element);
	if (iter == vec.end())
		return false;
	else
		return true;
}

void VectorPlus::Unique(vector<int>& result_vec, vector<int>& vec)
{
	set<int> unique_set;
	for (int i = 0; i < vec.size(); ++i)
		unique_set.insert(vec[i]);

	set<int>::iterator iter;
	result_vec.reserve(unique_set.size());
	for (iter = unique_set.begin(); iter != unique_set.end(); ++iter)
		result_vec.push_back(*iter);
}

int VectorPlus::Find(vector<int>& vec, int element)
{
	for (int i = 0; i < vec.size(); ++i)
	{
		if (vec[i] == element)
			return i;
	}
	return -1;
}

void VectorPlus::Find(vector<int>& result_vec, vector<int>& vec, int element)
{
	for (int i = 0; i < vec.size(); ++i)
	{
		if (vec[i] == element)
			result_vec.push_back(i);
	}
}

void VectorPlus::Union(vector<int>& result_vec, vector<int>& vec_A, vector<int>& vec_B)
{
	Copy(result_vec, vec_A);
	Union(result_vec, vec_B);
}

void VectorPlus::Union(vector<int>& vec_A, vector<int>& vec_B)
{
	for (int i = 0; i < vec_B.size(); ++i)
	{
		if (!Contains(vec_A, vec_B[i]))
			vec_A.push_back(vec_B[i]);
	}
}


void VectorPlus::Concat(vector<int>& result_vec, vector<int>& vec_A, vector<int>& vec_B)
{
	assert(result_vec.size() == 0);
	result_vec.reserve(vec_A.size() + vec_B.size());
	for (int i = 0; i <vec_A.size(); ++i)
	{
		result_vec.push_back(vec_A[i]);
	}

	for (int i = 0; i < vec_B.size(); ++i)
	{
		result_vec.push_back(vec_B[i]);
	}
}

int VectorPlus::Product(vector<int>& vec)
{
	if (vec.size() == 0)
		assert(0);

	int prod = 1;
	for (int i = 0; i < vec.size(); ++i)
	{
		prod = prod * vec[i];
	}
	return prod;
}

int VectorPlus::Sum(vector<int>& vec)
{
	if (vec.size() == 0)
		assert(0);

	int sum = 0;
	for (int i = 0; i < vec.size(); ++i)
	{
		sum += vec[i];
	}
	return sum;
}

double VectorPlus::Sum(vector<double>& vec)
{
	if (vec.size() == 0)
		assert(0);

	double sum = 0;
	for (int i = 0; i < vec.size(); ++i)
	{
		sum += vec[i];
	}
	return sum;
}


int VectorPlus::AbsSum(vector<int>& vec)
{
	if (vec.size() == 0)
		assert(0);

	int abs_sum = 0;
	for (int i = 0; i < vec.size(); ++i)
	{
		abs_sum += abs(vec[i]);
	}
	return abs_sum;
}

void VectorPlus::SetDiff(vector<int>& result_vec, vector<int>& vec_A, vector<int>& vec_B)
{
	assert(result_vec.size() == 0);
	result_vec.reserve(vec_A.size());
	for (int i = 0; i < vec_A.size(); ++i)
	{
		if (!Contains(vec_B, vec_A[i]))
		{
			result_vec.push_back(vec_A[i]);
		}
	}
}

void VectorPlus::Intersect(vector<int>& result_vec, vector<int>& vec_A, vector<int>& vec_B)
{
	assert(result_vec.size() == 0);
	result_vec.reserve(vec_A.size());
	for (int i = 0; i < vec_A.size(); ++i)
	{		
		if (Contains(vec_B, vec_A[i]))
		{
			result_vec.push_back(vec_A[i]);
		}
	}
}

void VectorPlus::Subset(vector<int>& result_vec, vector<int>& vec, vector<int>& indices)
{
	assert(result_vec.size() == 0);
	result_vec.reserve(indices.size());
	for (int i = 0; i < indices.size(); ++i)
	{
		result_vec.push_back(vec[indices[i]]);
	}
}

void VectorPlus::CSubset(vector<int>& result_vec, vector<int>& vec, vector<int>& indices)
{
	assert(result_vec.size() == 0);
	result_vec.reserve(vec.size() - indices.size());
	
	vector<int> c_indices;
	vector<int> all_indices;
	VectorPlus::Seq(all_indices, 0, 1, vec.size());
	VectorPlus::SetDiff(c_indices, all_indices, indices);

	for (int i = 0; i < c_indices.size(); ++i)
	{
		result_vec.push_back(vec[c_indices[i]]);
	}
}

void VectorPlus::Subset(vector<int>& result_vec, vector<int>& vec, int start_index, int end_index)
{
	assert(result_vec.size() == 0);
	result_vec.reserve(end_index - start_index);
	for (int i = start_index; i < end_index; ++i)
	{
		result_vec.push_back(vec[i]);
	}
}

double VectorPlus::Min(vector<double>& vec)
{
	double min = numeric_limits<double>::max();
	for (int i = 0; i < vec.size(); ++i)
	{
		if (vec[i] < min)
		{
			min = vec[i];
		}
	}
	return min;
}

double VectorPlus::Max(vector<double>& vec)
{
	double max = numeric_limits<double>::min();

	for (int i = 0; i < vec.size(); ++i)
	{
		if (vec[i] > max)
		{
			max = vec[i];
		}
	}

	return max;
}

void VectorPlus::Seq(vector<int>& result_vec, int start, int step, int end)
{
	assert(result_vec.size() == 0);
	result_vec.reserve((end - start) / step + 1);
	for (int i = start; i < end; i+=step)
	{
		result_vec.push_back(i);
	}
}

int VectorPlus::FindNextElement(vector<int>& vec_A, vector<int>& vec_B, int index)
{
	for (int i = index; i < vec_B.size(); ++i)
	{
		if (Contains(vec_A, vec_B[i]))
			return vec_B[i];
	}
	return -1;
}

vector<int> VectorPlus::CreateSingleton(int val)
{
	vector<int> vec(1, val);
	return vec;
}

vector<int> VectorPlus::CreatePair(int val1, int val2)
{
	vector<int> vec(2);
	vec[0] = val1;
	vec[1] = val2;

	return vec;
}

void VectorPlus::Reverse(vector<int>& result_vec, vector<int>& vec)
{
	assert(result_vec.size() == 0);
	result_vec.reserve(vec.size());
	for (int i = vec.size() - 1; i >= 0; --i)
	{
		result_vec.push_back(vec[i]);
	}
}

void VectorPlus::Copy(vector<int>& result_vec, vector<int>& vec)
{
	assert(result_vec.size() == 0);
	result_vec.reserve(vec.size());
	for (int i = 0; i < vec.size(); ++i)
	{
		result_vec.push_back(vec[i]);
	}
}

int VectorPlus::GetCyclicIndex(vector<int>& vec, int index)
{
	if (index >= 0 && index < vec.size())
		return vec[index];
	else if (index == -1)
		return vec[vec.size() - 1];
	else
		return vec[0];
}

// very slow (quadratic), only for small vectors
void VectorPlus::MultiSort(vector<int>& result_vec1, vector<int>& result_vec2, vector<int>& vec_to_sort, vector<int>& other_vec)
{
	vector<int> vec_copy;
	VectorPlus::Copy(vec_copy, vec_to_sort);
	for (int i = 0; i < vec_to_sort.size(); ++i)
	{
		int minIndex = MinIndex(vec_copy);
		result_vec1.push_back(vec_to_sort[minIndex]);
		result_vec2.push_back(other_vec[minIndex]);

		vec_copy[minIndex] = numeric_limits<int>::max();
	}
}
int VectorPlus::MinIndex(vector<int>& vec)
{
	int min = numeric_limits<int>::max();
	int minIndex = -1;
	for (int i = 0; i < vec.size(); ++i)
	{
		if (vec[i] < min)
		{
			min = vec[i];
			minIndex = i;
		}
	}
	return minIndex;
}

int VectorPlus::MaxIndex(vector<int>& vec)
{
	int max = numeric_limits<int>::min();
	int maxIndex = -1;
	for (int i = 0; i < vec.size(); ++i)
	{
		if (vec[i] > max)
		{
			max = vec[i];
			maxIndex = i;
		}
	}

	return maxIndex;

}

void VectorPlus::Rearrange(vector<int>& result_vec, vector<int>& vec_A, vector<int>& old_to_new_indices)
{
	result_vec.assign(vec_A.size(), 0);
	for (int i = 0; i < vec_A.size(); ++i)
	{
		result_vec[old_to_new_indices[i]] = vec_A[i];
	}
}

void VectorPlus::Concat(vector<int>& vec_A, vector<int>& vec_B)
{
	for (int i = 0; i <vec_B.size(); ++i)
	{
		vec_A.push_back(vec_B[i]);
	}
}


void VectorPlus::Concat(vector<double>& vec_A, vector<double>& vec_B)
{
	for (int i = 0; i <vec_B.size(); ++i)
	{
		vec_A.push_back(vec_B[i]);
	}
}

vector<int> VectorPlus::CreateTriple(int val1, int val2, int val3)
{
	vector<int> result_vec(3);
	result_vec[0] = val1;
	result_vec[1] = val2;
	result_vec[2] = val3;

	return result_vec;
}

vector<int> VectorPlus::CreateQuartet(int val1, int val2, int val3, int val4)
{
	vector<int> result_vec(4);
	result_vec[0] = val1;
	result_vec[1] = val2;
	result_vec[2] = val3;
	result_vec[3] = val4;

	return result_vec;
}

void VectorPlus::Permute(vector<int>& permuted_seq, vector<int>& orig_seq, int s)
{
	VectorPlus::Copy(permuted_seq, orig_seq);
	
	int curr_min_index = 0;
	for (int i = 0; i < orig_seq.size(); ++i)
	{
		int rand_index = curr_min_index + rand() % (orig_seq.size() - curr_min_index);
		
		// swap elements at rand and curr_min_index and update curr_min_index
		int temp = permuted_seq[curr_min_index];
		permuted_seq[curr_min_index] = permuted_seq[rand_index];
		permuted_seq[rand_index] = temp;
		curr_min_index++;
	}
	assert(curr_min_index == orig_seq.size());
}