# pragma once

#include <vector>
#include "Matrix.hpp"
#include "CPT.hpp"
#include "TensorCPT.hpp"
#include "MultiVector.hpp"
#include <assert.h>

using namespace std;

enum { HIDDEN_FLAG, OBSERVED_FLAG };

class BayesianNetwork {

private:

	int num_hidden_states;

	vector<vector<int>*>* parent_lists;
	vector<vector<int>*>* children_lists;
	vector<vector<int>*>* descendant_lists;
	
	vector<int>* type_vector;
	vector<int>* dims;
	vector<int>* nodes;

	vector<int>* topological_order;
	vector<int>* reverse_topological_order;

	vector<int>* elimination_order;

	map<vector<int>, Tensor*>* empirical_prob_map; 
	map<MultiVector<int>, Tensor*>* empirical_multi_prob_map;

	Matrix* samples;

	bool isInitialized;

	void ComputeTopologicalOrder();
	void GenerateRandomCPTs();

public:
	
	vector<CPT*>* CPTs;

	BayesianNetwork();
	BayesianNetwork(vector<vector<int>*>& parent_lists, vector<int>& type_vector, vector<int>& dims);
	

	void Initialize(vector<vector<int>*>& parent_lists, vector<int>& type_vector, vector<int>& dims);
	
	void GenerateSamples(Matrix& samples, int num_samples);
	void GenerateTrainSamples(Matrix& samples, int num_samples);
	void GenerateTestSamples(Matrix& samples, int num_samples) { GenerateSamples(samples, num_samples); }

	int NumNodes() { return nodes->size(); }
	CPT& GetCPT(int node) { return *(CPTs->at(node)); }
	vector<int>& parents(int node) { return *(parent_lists->at(node)); }
	vector<int>& children(int node) { return *(children_lists->at(node)); }
	vector<int>& descendants(int node) {return *(descendant_lists->at(node)); }
	vector<int>& TopologicalOrder() {return *topological_order; }
	vector<int>& ReverseTopologicalOrder() {return *reverse_topological_order; }
	bool is_hidden(int node)  { return (type_vector->at(node) == HIDDEN_FLAG); }
	bool is_observed(int node) {return (type_vector->at(node) == OBSERVED_FLAG); }	

	void SetCPT(int node_id, int index, double val) { CPTs->at(node_id)->Set(index, val); }
	void SetCPT(int node_id, vector<int>& index_array, double val) { CPTs->at(node_id)->Set(index_array, val); }
	vector<int>& GetDims() { return *dims; }
	void ClearEmpiricalProbabilityMap();
	void ComputeOracleProbTensor(Tensor& prob_tensor, vector<int>& vars);
	double ComputeEmpiricalProbVal(vector<int>& vars, vector<int>& vals);
	void ComputeEmpiricalProbTensor(Tensor& prob_tensor, vector<int>& vars);

	void ConvertIndices(vector<int>& multi_indices, vector<int>& indices, MultiVector<int>& multi_vars, vector<int>& vars);
	void ComputeEmpiricalMultiProbTensor(Tensor& multi_prob_tensor, MultiVector<int>& multi_vars);

	void ComputeConditionalProbTensor(Tensor& conditional_prob, vector<int>& front_vars, vector<int>& back_vars);
	void ComputeConditionalMultiProbTensor(Tensor& conditional_multi_prob, MultiVector<int>& multi_front_vars, MultiVector<int>& multi_back_vars);

	int GetNumHiddenStates();

	void GenerateRandomCPT(int node_id);

	void FindCliques(vector<vector<int>*>& nodes_to_vars);

	bool IsAdj(int node1, int node2);

	void FindJunctionTree(Matrix& jtree_matrix, vector<vector<int>*> cliques);

	void SetEliminationOrder(vector<int>& order) { elimination_order = new vector<int>(order); }

	void SetTrainingSamples(Matrix& training_samples) { assert(samples == NULL); samples = new Matrix(training_samples); } 
	~BayesianNetwork();
};
