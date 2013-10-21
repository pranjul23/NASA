#include <iostream>
#include "tensor.hpp"
#include "BayesianNetwork.hpp"
#include "TensorVarElim.hpp"
#include "TensorJTree.hpp"
#include "JTree.hpp"
#include "VectorPlus.hpp"
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

int main ()
{
	
  // node names
  int A = 0;
  int O1 = 1;
  int O2 = 2;
  int O3 = 3;

  int num_nodes = 4;
  vector<vector<int>*> parent_list;
  for (int i = 0; i < num_nodes; ++i)
	parent_list.push_back(new vector<int>());

  // hidden variables
  parent_list[O1]->push_back(A);
  parent_list[O2]->push_back(A);
  parent_list[O3]->push_back(A);

  vector<int> type_vector(num_nodes, 0);
  type_vector[0] = HIDDEN_FLAG;
  for (int i = 1; i < num_nodes; ++i)
  {
	type_vector[i] = OBSERVED_FLAG;
  }


  vector<int> dims(num_nodes, 2);
  dims[0] = 2;
  BayesianNetwork* bayesNet = new BayesianNetwork(parent_list, type_vector, dims);

    string line;
  ifstream param_file("cptfile.txt");
  for (int i = 0; i < num_nodes; ++i)
  {
	  if (i == A)
	  {
		for (int j = 0; j < 2; ++j)
		{
			getline(param_file, line);
			double val = atof(line.c_str());
			bayesNet->SetCPT(i, j, val); 
		}
	  }
	  else
	  {
		  for (int j = 0; j < 2; ++j)
		  {
			getline(param_file, line);
			double val1 = atof(line.c_str());
			bayesNet->SetCPT(i,j, val1);

			getline(param_file, line);
			double val2 = atof(line.c_str());
			bayesNet->SetCPT(i, j + 2, val2);
		  }
	  }
  }

  Matrix samples;
  bayesNet->GenerateSamples(samples, 100);

  // test the junction tree
  int jtree_nodes = 4;
  vector<int> tree_matrix_dims(2,jtree_nodes);
  
  //int A = 0;
  int O1A = 1;
  int O2A = 2;
  int O3A = 3;

  Matrix tree_matrix(tree_matrix_dims);
  tree_matrix.Set(O1A, A, 1);
  tree_matrix.Set(O2A, A, 1);
  tree_matrix.Set(O3A, A, 1);


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

 // JTree EMtree(*bayesNet, tree_matrix, nodes_to_vars, samples);
 // EMtree.LearnEMParameters();
  vector<int> jtree_vec;
  VectorPlus::Seq(jtree_vec, 0, 1, jtree_nodes);
  TensorJTree jtree(bayesNet, tree_matrix, nodes_to_vars, 1, jtree_vec);
  jtree.LearnSpectralParameters();
  vector<int> evidence_vars;
  VectorPlus::Seq(evidence_vars, O1, 1, O3 + 1);
  vector<int> evidence_vals;
  evidence_vals.push_back(0);
  evidence_vals.push_back(0);
  evidence_vals.push_back(0);
  double spectral_prob = jtree.ComputeMarginalEmpiricalProbability(evidence_vars, evidence_vals);
  double oracle_prob = TensorVarElim::VE(*bayesNet, evidence_vars, evidence_vals, bayesNet->ReverseTopologicalOrder());

  return 0;
}