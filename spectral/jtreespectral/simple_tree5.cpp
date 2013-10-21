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
  int B = 1;
  int C = 2;
  int D = 3;

  int O1 = 4;
  int O2 = 5;
  int O3 = 6;
  int O4 = 7;
  int O5 = 8;
  int O6 = 9;

  int num_nodes = 10;
  vector<vector<int>*> parent_list;
  for (int i = 0; i < num_nodes; ++i)
	parent_list.push_back(new vector<int>());

  // hidden variables
  parent_list[B]->push_back(A);
  parent_list[C]->push_back(A);
  parent_list[D]->push_back(A);

  parent_list[O1]->push_back(B);
  parent_list[O2]->push_back(B);
  parent_list[O3]->push_back(C);
  parent_list[O4]->push_back(C);
  parent_list[O5]->push_back(D);
  parent_list[O6]->push_back(D);
 
  vector<int> type_vector(num_nodes, 0);
  for (int i = 0; i < num_nodes; ++i)
  {
	  if (i < 4)
		type_vector[i] = HIDDEN_FLAG;
	  else
		type_vector[i] = OBSERVED_FLAG;
  }


  vector<int> dims(num_nodes, 2);
  BayesianNetwork* bayesNet = new BayesianNetwork(parent_list, type_vector, dims);
  Matrix samples;
  bayesNet->GenerateSamples(samples, 1000);

  int jtree_nodes = 10;
  vector<int> tree_matrix_dims(2,jtree_nodes);
  
  int AB = 1;
  int AC = 2;
  int AD = 3;

  int O1B = 4;
  int O2B = 5;
  int O3C = 6;
  int O4C = 7;
  int O5D = 8;
  int O6D = 9;
  
  Matrix tree_matrix(tree_matrix_dims);
  tree_matrix.Set(B, A, 1);
  tree_matrix.Set(C, A, 1);
  tree_matrix.Set(D, A, 1);

  tree_matrix.Set(O1, B, 1);
  tree_matrix.Set(O2, B, 1);
  tree_matrix.Set(O3, C, 1);
  tree_matrix.Set(O4, C, 1);
  tree_matrix.Set(O5, D, 1);
  tree_matrix.Set(O6, D, 1);

  map<int, vector<int>*> nodes_to_vars;
  for (int i = 0; i < jtree_nodes; ++i)
	nodes_to_vars[i] = new vector<int>();
  
  nodes_to_vars[A]->push_back(A);
  nodes_to_vars[AB]->push_back(A);
  nodes_to_vars[AB]->push_back(B);
  nodes_to_vars[AC]->push_back(A);
  nodes_to_vars[AC]->push_back(C);
  nodes_to_vars[AD]->push_back(A);
  nodes_to_vars[AD]->push_back(D);

  nodes_to_vars[O1B]->push_back(O1);
  nodes_to_vars[O1B]->push_back(B);
  nodes_to_vars[O2B]->push_back(O2);
  nodes_to_vars[O2B]->push_back(B);
  nodes_to_vars[O3C]->push_back(O3);
  nodes_to_vars[O3C]->push_back(C);
  nodes_to_vars[O4C]->push_back(O4);
  nodes_to_vars[O4C]->push_back(C);
  nodes_to_vars[O5D]->push_back(O5);
  nodes_to_vars[O5D]->push_back(D);
  nodes_to_vars[O6D]->push_back(O6);
  nodes_to_vars[O6D]->push_back(D);

 
  vector<int> evidence_vars;
  VectorPlus::Seq(evidence_vars, O1, 1, O6 + 1);
  vector<int> evidence_vals(evidence_vars.size(), 0);
  
  double prob = TensorVarElim::VE(*bayesNet, evidence_vars, evidence_vals, bayesNet->ReverseTopologicalOrder());
  TensorJTree jtree(*bayesNet, tree_matrix, nodes_to_vars);
  jtree.LearnSpectralParameters();
  double val = jtree.ComputeMarginalEmpiricalProbability(evidence_vars, evidence_vals);
  return 0;
}