#include "tensor.hpp"
#include "BayesianNetwork.hpp"
#include "TensorVarElim.hpp"
#include "TensorJTree.hpp"
#include "VectorPlus.hpp"
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

int main ()
{
	
  // node names
  int A = 0;
  int B = 1;
  int C = 2;
  int D = 3;
  int E = 4;

  int O1 = 5;
  int O2 = 6;
  int O3 = 7;
  int O4 = 8;
  int O5 = 9;
  int O6 = 10;
  int O7 = 11;
  int O8 = 12;
  int O9 = 13;

  int num_nodes = 14;
  vector<vector<int>*> parent_list;
  for (int i = 0; i < num_nodes; ++i)
	parent_list.push_back(new vector<int>());

  parent_list[B]->push_back(A);
  parent_list[C]->push_back(A);
  parent_list[D]->push_back(B);
  parent_list[D]->push_back(A);
  parent_list[E]->push_back(B);

  // hidden variables
  parent_list[O1]->push_back(C);
  parent_list[O2]->push_back(C);
  parent_list[O3]->push_back(C);
  parent_list[O4]->push_back(D);
  parent_list[O5]->push_back(D);
  parent_list[O6]->push_back(D);
  parent_list[O7]->push_back(E);
  parent_list[O8]->push_back(E);
  parent_list[O9]->push_back(E);


  vector<int> type_vector(num_nodes, 0);

  for (int i = 1; i < num_nodes; ++i)
  {
	  if (i <= 4)
		type_vector[i] = HIDDEN_FLAG;
	  else
		type_vector[i] = OBSERVED_FLAG;
  }


  vector<int> dims(num_nodes, 2);
  dims[0] = 2;
  BayesianNetwork* bayesNet = new BayesianNetwork(parent_list, type_vector, dims);

 /*   string line;
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
	  else if (i == D)
	  {
		  for (int j = 0; j < 4; ++j)
		  {
			 getline(param_file, line);
			double val1 = atof(line.c_str());
			bayesNet->SetCPT(i,j, val1);

			getline(param_file, line);
			double val2 = atof(line.c_str());
			bayesNet->SetCPT(i, j+4, val2);
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
  param_file.close();*/

  Matrix samples;
  bayesNet->GenerateSamples(samples, 100);

  // test the junction tree
  int jtree_nodes = 12;
  vector<int> tree_matrix_dims(2,jtree_nodes);
  
  //int A = 0;
  int ABD = 0;
  int AC = 1;
  int BE = 2;

  int O1C = 3;
  int O2C = 4;
  int O3C = 5;
  int O4D = 6;
  int O5D = 7;
  int O6D = 8;
  int O7E = 9;
  int O8E = 10;
  int O9E = 11;

  Matrix tree_matrix(tree_matrix_dims);
  tree_matrix.Set(BE, ABD, 1);
  tree_matrix.Set(AC, ABD, 1);
  tree_matrix.Set(O1C, AC, 1);
  tree_matrix.Set(O2C, AC, 1);
  tree_matrix.Set(O3C, AC, 1);
  tree_matrix.Set(O4D, ABD, 1);
  tree_matrix.Set(O5D, ABD, 1);
  tree_matrix.Set(O6D, ABD, 1);
  tree_matrix.Set(O7E, BE, 1);
  tree_matrix.Set(O8E, BE, 1);
  tree_matrix.Set(O9E, BE, 1);
  map<int, vector<int>*> nodes_to_vars;
  for (int i = 0; i < jtree_nodes; ++i)
	nodes_to_vars[i] = new vector<int>();
  

  nodes_to_vars[ABD]->push_back(A);
  nodes_to_vars[ABD]->push_back(B);
  nodes_to_vars[ABD]->push_back(D);

  nodes_to_vars[AC]->push_back(A);
  nodes_to_vars[AC]->push_back(C);

  nodes_to_vars[BE]->push_back(B);
  nodes_to_vars[BE]->push_back(E);

  nodes_to_vars[O1C]->push_back(O1);
  nodes_to_vars[O1C]->push_back(C);
  nodes_to_vars[O2C]->push_back(O2);
  nodes_to_vars[O2C]->push_back(C);
  nodes_to_vars[O3C]->push_back(O3);
  nodes_to_vars[O3C]->push_back(C);

  nodes_to_vars[O4D]->push_back(O4);
  nodes_to_vars[O4D]->push_back(D);
  nodes_to_vars[O5D]->push_back(O5);
  nodes_to_vars[O5D]->push_back(D);
  nodes_to_vars[O6D]->push_back(O6);
  nodes_to_vars[O6D]->push_back(D);

  nodes_to_vars[O7E]->push_back(O7);
  nodes_to_vars[O7E]->push_back(E);
  nodes_to_vars[O8E]->push_back(O8);
  nodes_to_vars[O8E]->push_back(E);
  nodes_to_vars[O9E]->push_back(O9);
  nodes_to_vars[O9E]->push_back(E);

  vector<int> evidence_vars;
  VectorPlus::Seq(evidence_vars, O1, 1, O9 + 1);
  vector<int> evidence_vals(evidence_vars.size(), 0);

  TensorJTree jtree(*bayesNet, tree_matrix, nodes_to_vars);
  jtree.LearnSpectralParameters();
  double spectral_prob = jtree.ComputeMarginalEmpiricalProbability(evidence_vars, evidence_vals);

  double oracle_prob = TensorVarElim::VE(*bayesNet, evidence_vars, evidence_vals, bayesNet->ReverseTopologicalOrder());

   cout << spectral_prob;
  cout << "\n";
  cout << oracle_prob;

  return 0;
}