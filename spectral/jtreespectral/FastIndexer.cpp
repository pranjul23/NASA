#include "FastIndexer.hpp"
#include "VectorPlus.hpp"

using namespace std;

// default constructor
FastIndexer::FastIndexer()
{
	indices = NULL;
	dims = NULL;

}

// default constructor
FastIndexer::FastIndexer(vector<int>& dims)
{
	indices = new vector<int>(dims.size(), 0);
	this->dims = new vector<int>(dims);

	curr_index = dims.size() - 1;

	indexedCount = 0;
	totalCount = VectorPlus::Product(dims);
}

bool FastIndexer::HasNext() { return (indexedCount < totalCount); }

vector<int>& FastIndexer::GetNext()
{
	if (indexedCount == 0)
	{
		indexedCount++;
		return *indices;
	}

	if (indices->at(curr_index) < dims->at(curr_index) - 1)
	{
		indices->at(curr_index)++;
	}
	else
	{
		indices->at(curr_index) = 0;
		curr_index--;
		while(true)
		{
			if (indices->at(curr_index) < dims->at(curr_index) - 1)
			{
				indices->at(curr_index)++;
				break;
			}
			else
			{
				indices->at(curr_index) = 0;
				curr_index--;
			}
		}
		curr_index = dims->size() - 1;
	}
	indexedCount++;
	return *indices;
}

FastIndexer::~FastIndexer()
{
	delete(dims);
	delete(indices);
}