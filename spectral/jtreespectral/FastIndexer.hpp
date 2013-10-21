// This is a limited tensor class for the simple times one operation that is needed.

# pragma once

#include <vector>

using namespace std;

class FastIndexer {

protected:

	vector<int>* indices;
	vector<int>* dims;
	
	int curr_index;
	int indexedCount;
	int totalCount;

public:
	FastIndexer();
	FastIndexer(vector<int>& dims);

	bool HasNext();
	vector<int>& GetNext();
	
	~FastIndexer();
};
