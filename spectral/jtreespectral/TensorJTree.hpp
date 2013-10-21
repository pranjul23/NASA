# pragma once

#include <vector>
#include "Matrix.hpp"
#include "TensorCPT.hpp"
#include <map>
#include "BayesianNetwork.hpp"
#include "TensorJTNode.hpp"
#include "MultiVector.hpp"

using namespace std;

class TensorJTNode;

class TensorJTree
{
	private:

		BayesianNetwork* bNet;
		
		vector<int>* all_obs_vars; 

		Matrix* tree_matrix;
		vector<TensorJTNode*>* node_list;
		int root;

		vector<int>* elimination_order;

		void FindObsDescendants(vector<int>& desc_vec, int clique_index);

		vector<int>* clique_to_template;
		vector<vector<int>*>* template_to_cliques;

		vector<vector<int>*>* dist_vars;

	public:

		TensorJTree();
		~TensorJTree();
		TensorJTree(BayesianNetwork& bNet);
		TensorJTree(BayesianNetwork* bNet, Matrix& tree_matrix, vector<vector<int>*>& nodes_to_vars, int max_obs_per_mode, vector<int>& template_vec, int max_runs);

		void ComputeDistances();

		int NumCliques() { return node_list->size(); }
		vector<int>& GetChildren(int clique_index);
		int GetParent(int clique_index);
		vector<int>& GetCVars(int clique_index);
		vector<int>& GetRVars(int clique_index);
		vector<int>& GetSVars(int clique_index);

		map<vector<int>, vector<int>*>& GetChildSeps(int clique_index);
		vector<int>& GetDescendants(int clique_index, vector<int>& mode_group, int index);
		bool IsLeaf(int clique_index);

		void LearnSpectralParameters();
		vector<int>& GetObsVector(int node, vector<int>& modes);
		void ComputeEmpiricalMultiProbTensor(Tensor& multi_prob_tensor, MultiVector<int>& multi_vars) { return bNet->ComputeEmpiricalMultiProbTensor(multi_prob_tensor, multi_vars); }
		void FindObservationPartitions(int clique_index);

		void ComputeEliminationOrder();
		double ComputeMarginalEmpiricalProbability(vector<int> evidence_vars, vector<int> evidence_vals);

		bool IsFirstInTemplate(int clique_index);
		int GetFirstInTemplate(int clique_index);
		int GetLastInTemplate(int clique_index);

		Tensor& GetU(int clique_index);
		int GetNumUniqueChildSeps(int clique_index);

		int GetSIndex(int clique_index);
		int GetSameSiblingIndex(int clique_index);
		int GetChildSepIndex(int clique_index, vector<int>& S_vars);

		int GetNumHiddenStates() { return bNet->GetNumHiddenStates(); }

		TensorCPT& GetTransformTensor(int clique_index);

		vector<int>& GetObsVector(int node, int other_node);

		int GetNumDescendantPartitions(int node_id, vector<int>& svars);

		int GetOtherNode(int node, vector<int>& mode_group);

		void LearnTemplateSpectralParameters();
};
