# pragma once

#include <vector>

using namespace std;

template <class T>

class MultiVector
{
	private:

		vector<int>* vec_sizes;
		void DeleteIfPointer(T& item);
		void DeleteIfPointer(T* item);

	public:

		vector<vector<T>*>* multi_vector;

		MultiVector();
		~MultiVector();

		int Size() {return multi_vector->size(); }

		T At(int index1, int index2);
		void At(vector<T>& result, vector<int>& indices);
		vector<int>& GetVecSizes();
		void Set(int index1, int index2, T val);

		vector<T>& GetVec(int index1);

		static void Concat(MultiVector<T>& result_multi_vec, MultiVector<T>& multi_vec1, MultiVector<T>& multi_vec2);
};

template <class T>
MultiVector<T>::MultiVector()
{
	multi_vector = new vector<vector<T>*>();
	vec_sizes = new vector<int>();
}

template <class T>
void MultiVector<T>::DeleteIfPointer(T& item) {}

template <class T>
void MultiVector<T>::DeleteIfPointer(T* item) { delete item; }

template <class T>
MultiVector<T>::~MultiVector()
{
	for (int i = 0; i < multi_vector->size(); ++i)
	{
		vector<T>* vec = multi_vector->at(i);
		for (int j = 0; j < vec->size(); ++j)
		{
			DeleteIfPointer(vec->at(j));
		}
		delete(vec);
	}
	delete(multi_vector);
	delete(vec_sizes);
}

template <class T>
T MultiVector<T>::At(int index1, int index2)
{
	return multi_vector->at(index1)->at(index2);
}

template <class T>
void MultiVector<T>::Set(int index1, int index2, T val)
{
	if (multi_vector->size() <= index1)
	{
		while (multi_vector->size() <= index1)
		{
			multi_vector->push_back(new vector<T>());
		}
	}

	vector<T>* vec = multi_vector->at(index1);
	if (vec->size() <= index2)
	{
		while (vec->size() <= index2)
		{
			vec->push_back(0);
		}
	}

	multi_vector->at(index1)->at(index2) = val;
}

template <class T>
vector<T>& MultiVector<T>::GetVec(int index1)
{
	return *(multi_vector->at(index1));
}

template <class T>
void MultiVector<T>::At(vector<T>& result, vector<int>& indices)
{
	for (int i = 0; i < Size(); ++i)
	{
		result.push_back(multi_vector->at(i)->at(indices[i]));
	}
}

template <class T>
vector<int>& MultiVector<T>::GetVecSizes()
{
	vec_sizes->clear();
	for (int i = 0; i < Size(); ++i)
	{
		vec_sizes->push_back(multi_vector->at(i)->size());
	}
	return *vec_sizes;
}

template <class T>
void MultiVector<T>::Concat(MultiVector<T>& result_multi_vec, MultiVector<T>& multi_vec1, MultiVector<T>& multi_vec2)
{
	int curr_index = 0;
	for (int i = 0; i < multi_vec1.Size(); ++i)
	{
		vector<T>& vec = multi_vec1.GetVec(i);
		for (int j = 0; j < vec.size(); ++j)
		{
			result_multi_vec.Set(curr_index, j, vec[j]);
		}
		curr_index++;
	}

	for (int i = 0; i < multi_vec2.Size(); ++i)
	{
		vector<T>& vec = multi_vec2.GetVec(i);
		for (int j = 0; j < vec.size(); ++j)
		{
			result_multi_vec.Set(curr_index, j, vec[j]);
		}
		curr_index++;
	}
}

template <class T>
const bool operator < ( const MultiVector<T>& multi_vec1, const  MultiVector<T>& multi_vec2 )
{
	if (multi_vec1.multi_vector->size() < multi_vec2.multi_vector->size())
	{
		return true;
	}
	if (multi_vec1.multi_vector->size() > multi_vec2.multi_vector->size())
	{
		return false;
	}

	for (int i = 0; i < multi_vec1.multi_vector->size(); ++i)
	{
		if (multi_vec1.multi_vector->at(i)->size() < multi_vec2.multi_vector->at(i)->size())
			return true;
		else if (multi_vec1.multi_vector->at(i)->size() > multi_vec2.multi_vector->at(i)->size())
			return false;

		for (int j = 0; j < multi_vec1.multi_vector->at(i)->size(); ++j)
		{
			if (multi_vec1.multi_vector->at(i)->at(j) < multi_vec2.multi_vector->at(i)->at(j))
				return true;
			else if (multi_vec1.multi_vector->at(i)->at(j) < multi_vec2.multi_vector->at(i)->at(j))
				return false;
		}

	}

	return false;
}
