# pragma once

#include <vector>

using namespace std;

class VectorPlus
{
public:
		static bool Contains(vector<int>& vec, int element);
		static bool Equals(vector<int>& vec_A, vector<int>& vec_B);
		static bool IsSubset(vector<int>& vec_A, vector<int>& vec_B);
		static int Find(vector<int>& vec, int element);
		static void Find(vector<int>& result_vec, vector<int>& vec, int element);
		static int Product(vector<int>& vec);
		static int Sum(vector<int>& vec);
		static double Sum(vector<double>& vec);
		static void Copy(vector<int>& result_vec, vector<int>& vec);

		static int GetCyclicIndex(vector<int>& vec, int index);
		static void Unique(vector<int>& result_vec, vector<int>& vec);
		static void Union(vector<int>& result_vec, vector<int>& vec_A, vector<int>& vec_B);
		static void Union(vector<int>& vec_A, vector<int>& vec_B);
		static int AbsSum(vector<int>& vec);
		static void SetDiff(vector<int>& result_vec, vector<int>& vec_A, vector<int>& vec_B);
		static void Intersect(vector<int>& result_vec, vector<int>& vec_A, vector<int>& vec_B);
		static void Reverse(vector<int>& result_vec, vector<int>& vec);
		static void Concat(vector<int>& result_vec, vector<int>& vec_A, vector<int>& vec_B);

		static void Concat(vector<int>& vec_A, vector<int>& vec_B);
		static void Concat(vector<double>& vec_A, vector<double>& vec_B);

		static void Subset(vector<int>& result_vec, vector<int>& vec, vector<int>& indices);
		static void CSubset(vector<int>& result_vec, vector<int>& vec, vector<int>& indices);
		static void Subset(vector<int>& result_vec, vector<int>& vec, int start_index, int end_index);

		static void Seq(vector<int>& result_vec, int start, int step, int end);
		static int FindNextElement(vector<int>& vec_A, vector<int>& vec_B, int index);
		static double Min(vector<double>& vec);
		static double Max(vector<double>& vec);
		static int MinIndex(vector<int>& vec);
		static int MaxIndex(vector<int>& vec);

		static void MatchSub(vector<int>& result_vec, vector<int> vec, vector<int> match_vec, vector<int> sub_match_vec);

		static void MultiSort(vector<int>& result_vec1, vector<int>& result_vec2, vector<int>& vec_to_sort, vector<int>& other_vec);

		static vector<int> CreateSingleton(int val);
		static vector<int> CreatePair(int val1, int val2);
		static vector<int> CreateTriple(int val1, int val2, int val3);
		static vector<int> CreateQuartet(int val1, int val2, int val3, int val4);

		static void Add(vector<int>& sum_vec, vector<int>& vec_A, vector<int>& vec_B);

		static void Permute(vector<int>& permuted_seq, vector<int>& orig_seq, int s);

		static void Rearrange(vector<int>& result_vec, vector<int>& vec_A, vector<int>& old_to_new_indices);
};
