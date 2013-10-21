# pragma once

#include <vector>
#include "Matrix.hpp"
#include "TensorCPT.hpp"
#include <map>
#include "BayesianNetwork.hpp"
#include "JTNode.hpp"
#include "MultiVector.hpp"
#include "TemplateBayesianNetwork.hpp"

using namespace std;

class JTNode;

class JTree
{
	private:

		BayesianNetwork* bNet;
		Matrix* training_samples;
		vector<int>* all_obs_vars; 

		Matrix* tree_matrix;
		vector<JTNode*>* node_list;
		int root;

		vector<int>* elimination_order;

		double thresh;

		bool clique_params_flag;

	public:

		JTree();
		~JTree();
		JTree(BayesianNetwork& bNet);
//		JTree(JTree& other);
		JTree(BayesianNetwork* bNet, Matrix& tree_matrix, vector<vector<int>*>& nodes_to_vars, Matrix* samples, double thresh, bool clique_params_flag);

		int NumCliques() { return node_list->size(); }
		vector<int>& GetChildren(int clique_index);
		int GetParent(int clique_index);
		vector<int>& GetCVars(int clique_index);
		vector<int>& GetRVars(int clique_index);
		vector<int>& GetSVars(int clique_index);

		bool IsLeaf(int clique_index);
		TensorCPT& GetAgglomCliqueFactor(int clique_index);

		double LearnEMParameters(int max_iter, int oracle_flag);
		double ComputeEmpiricalMarginals(vector<TensorCPT*>& marginals, vector<int>& evidence_vars, vector<int>& evidence_vals);

		void ComputeEliminationOrder();
		int GetNumHiddenStates() { return bNet->GetNumHiddenStates(); }

		double LearnOnlineEMParameters(double learning_rate, int oracle_flag, int max_iter);

		bool TestConvergence(double prev_likelihood, double curr_likelihood, double thresh);

		double ComputeEmpiricalJointProb(vector<TensorCPT*>& marginals, vector<int>& evidence_vars, vector<int>& evidence_vals);
};
