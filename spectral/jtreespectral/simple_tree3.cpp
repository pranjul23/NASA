#include <iostream>
#include "tensor.hpp"
#include "BayesianNetwork.hpp"
#include "TensorVarElim.hpp"
#include "TensorJTree.hpp"
#include "VectorPlus.hpp"

using namespace std;

int main ()
{
	
  // node names
  int A = 0;
  int O1 = 1;
  int O2 = 2;
  int O3 = 3;
  int O4 = 4;
  int O5 = 5;
  int O6 = 6;

  int num_nodes = 7;
  vector<vector<int>*> parent_list;
  for (int i = 0; i < num_nodes; ++i)
	parent_list.push_back(new vector<int>());

  // hidden variables
  parent_list[O1]->push_back(A);
  parent_list[O2]->push_back(A);
  parent_list[O3]->push_back(A);
  parent_list[O4]->push_back(A);
  parent_list[O5]->push_back(A);
  parent_list[O6]->push_back(A);

  vector<int> type_vector(num_nodes, 0);
  type_vector[0] = HIDDEN_FLAG;
  for (int i = 1; i < num_nodes; ++i)
  {
	type_vector[i] = OBSERVED_FLAG;
  }


  vector<int> dims(num_nodes, 2);
  dims[0] = 2;
  BayesianNetwork* bayesNet = new BayesianNetwork(parent_list, type_vector, dims);
  Matrix samples;
  bayesNet->GenerateSamples(samples, 100);

  // test the junction tree
  int jtree_nodes = 7;
  vector<int> tree_matrix_dims(2,jtree_nodes);
  
  //int A = 0;
  int O1A = 1;
  int O2A = 2;
  int O3A = 3;
  int O4A = 4;
  int O5A = 5;
  int O6A = 6;

  Matrix tree_matrix(tree_matrix_dims);
  tree_matrix.Set(O1A, A, 1);
  tree_matrix.Set(O2A, A, 1);
  tree_matrix.Set(O3A, A, 1);
  tree_matrix.Set(O4A, A, 1);
  tree_matrix.Set(O5A, A, 1);
  tree_matrix.Set(O6A, A, 1);

  map<int, vector<int>*> nodes_to_vars;
  for (int i = 0; i < jtree_nodes; ++i)
	nodes_to_vars[i] = new vector<int>();
  
  nodes_to_vars[A]->push_back(A);
  nodes_to_vars[O1A]->push_back(O1);
  nodes_to_vars[O2A]->push_back(O2);
  nodes_to_vars[O3A]->push_back(O3);
  nodes_to_vars[O1A]->push_back(A);
  nodes_to_vars[O2A]->push_back(A);
  nodes_to_vars[O3A]->push_back(A);
  nodes_to_vars[O4A]->push_back(O4);
  nodes_to_vars[O5A]->push_back(O5);
  nodes_to_vars[O6A]->push_back(O6);
  nodes_to_vars[O4A]->push_back(A);
  nodes_to_vars[O5A]->push_back(A);
  nodes_to_vars[O6A]->push_back(A);

  TensorJTree jtree(*bayesNet, tree_matrix, nodes_to_vars);
  jtree.LearnSpectralParameters();
  vector<int> evidence_vars;
  VectorPlus::Seq(evidence_vars, O1, 1, O6 + 1);
  vector<int> evidence_vals(evidence_vars.size(), 1);

  double spectral_prob = jtree.ComputeMarginalEmpiricalProbability(evidence_vars, evidence_vals);
  double oracle_prob = TensorVarElim::VE(*bayesNet, evidence_vars, evidence_vals, bayesNet->ReverseTopologicalOrder());

  return 0;
}