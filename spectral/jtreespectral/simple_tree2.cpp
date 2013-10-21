#include <iostream>
#include "tensor.hpp"
#include "BayesianNetwork.hpp"
#include "TensorVarElim.hpp"
#include "TensorJTree.hpp"
#include "VectorPlus.hpp"
#include "JTree.hpp"

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
  int O7 = 10;
  int O8 = 11;
  int O9 = 12;
  int O10 = 13;
  int O11 = 14;
  int O12 = 15;
  int O13 = 16;
  int O14 = 17;
  int O15 = 18;
  int O16 = 19;
  int O17 = 20;
  int O18 = 21;

  int num_nodes = 22;
  vector<vector<int>*> parent_list;
  for (int i = 0; i < num_nodes; ++i)
	parent_list.push_back(new vector<int>());

  // hidden variables
  parent_list[B]->push_back(A);
  parent_list[C]->push_back(A);
  parent_list[D]->push_back(A);

  parent_list[O1]->push_back(B);
  parent_list[O2]->push_back(B);
  parent_list[O3]->push_back(B);
  parent_list[O4]->push_back(B);
  parent_list[O5]->push_back(B);
  parent_list[O6]->push_back(B);
  parent_list[O7]->push_back(C);
  parent_list[O8]->push_back(C);
  parent_list[O9]->push_back(C);
  parent_list[O10]->push_back(C);
  parent_list[O11]->push_back(C);
  parent_list[O12]->push_back(C);
  parent_list[O13]->push_back(D);
  parent_list[O14]->push_back(D);
  parent_list[O15]->push_back(D);
  parent_list[O16]->push_back(D);
  parent_list[O17]->push_back(D);
  parent_list[O18]->push_back(D);

  vector<int> type_vector(num_nodes, 0);
    vector<int> dims(num_nodes, 0);
  for (int i = 0; i < num_nodes; ++i)
  {
	  if (i < 4)
	  {
		  dims[i] = 2;
		type_vector[i] = HIDDEN_FLAG;
	  }
	  else
	  {
		  dims[i] = 4;
		type_vector[i] = OBSERVED_FLAG;
	  }
  }

  vector<int> template_vec;
  VectorPlus::Seq(template_vec, 0, 1, num_nodes);
  BayesianNetwork* bayesNet = new BayesianNetwork(parent_list, type_vector, dims);
  Matrix samples;
  Matrix test_samples;
  bayesNet->GenerateTrainSamples(samples, 200);
  bayesNet->GenerateTestSamples(test_samples, 1000); 
  int jtree_nodes = 22;
  vector<int> tree_matrix_dims(2,jtree_nodes);
  
  int AB = 1;
  int AC = 2;
  int AD = 3;

  int O1B = 4;
  int O2B = 5;
  int O3B = 6;
  int O4B = 7;
  int O5B = 8;
  int O6B = 9;
  int O7C = 10;
  int O8C = 11;
  int O9C = 12;
  int O10C = 13;
  int O11C = 14;
  int O12C = 15;
  int O13D = 16;
  int O14D = 17;
  int O15D = 18;
  int O16D = 19;
  int O17D = 20;
  int O18D = 21;

  Matrix tree_matrix(tree_matrix_dims);
  tree_matrix.Set(B, A, 1);
  tree_matrix.Set(C, A, 1);
  tree_matrix.Set(D, A, 1);

  tree_matrix.Set(O1, B, 1);
  tree_matrix.Set(O2, B, 1);
  tree_matrix.Set(O3, B, 1);
  tree_matrix.Set(O4, B, 1);
  tree_matrix.Set(O5, B, 1);
  tree_matrix.Set(O6, B, 1);
  tree_matrix.Set(O7, C, 1);
  tree_matrix.Set(O8, C, 1);
  tree_matrix.Set(O9, C, 1);
  tree_matrix.Set(O10, C, 1);
  tree_matrix.Set(O11, C, 1);
  tree_matrix.Set(O12, C, 1);
  tree_matrix.Set(O13, D, 1);
  tree_matrix.Set(O14, D, 1);
  tree_matrix.Set(O15, D, 1);
  tree_matrix.Set(O16, D, 1);
  tree_matrix.Set(O17, D, 1);
  tree_matrix.Set(O18, D, 1);

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
  nodes_to_vars[O3B]->push_back(O3);
  nodes_to_vars[O3B]->push_back(B);
  nodes_to_vars[O4B]->push_back(O4);
  nodes_to_vars[O4B]->push_back(B);
  nodes_to_vars[O5B]->push_back(O5);
  nodes_to_vars[O5B]->push_back(B);
  nodes_to_vars[O6B]->push_back(O6);
  nodes_to_vars[O6B]->push_back(B);
  nodes_to_vars[O7C]->push_back(O7);
  nodes_to_vars[O7C]->push_back(C);
  nodes_to_vars[O8C]->push_back(O8);
  nodes_to_vars[O8C]->push_back(C);
  nodes_to_vars[O9C]->push_back(O9);
  nodes_to_vars[O9C]->push_back(C);
  nodes_to_vars[O10C]->push_back(O10);
  nodes_to_vars[O10C]->push_back(C);
  nodes_to_vars[O11C]->push_back(O11);
  nodes_to_vars[O11C]->push_back(C);
  nodes_to_vars[O12C]->push_back(O12);
  nodes_to_vars[O12C]->push_back(C);
  nodes_to_vars[O13D]->push_back(O13);
  nodes_to_vars[O13D]->push_back(D);
  nodes_to_vars[O14D]->push_back(O14);
  nodes_to_vars[O14D]->push_back(D);
  nodes_to_vars[O15D]->push_back(O15);
  nodes_to_vars[O15D]->push_back(D);
  nodes_to_vars[O16D]->push_back(O16);
  nodes_to_vars[O16D]->push_back(D);
  nodes_to_vars[O17D]->push_back(O17);
  nodes_to_vars[O17D]->push_back(D);
  nodes_to_vars[O18D]->push_back(O18);
  nodes_to_vars[O18D]->push_back(D);
 
  vector<int> evidence_vars;
  VectorPlus::Seq(evidence_vars, O1, 1, O18 + 1);
  vector<int> evidence_vals(evidence_vars.size(), 0);

double total_error = 0;
double EM_error1 = 0;
double EM_error2 = 0;
double EM_error3 = 0;

  JTree EMtree(bayesNet, EM_tree_matrix, EM_nodes_to_vars, train_samples);
  EMtree.LearnEMParameters(1, 1);

  JTree EMtree2(bayesNet, EM_tree_matrix, EM_nodes_to_vars, train_samples);
  EMtree2.LearnEMParameters(10, 0);

  JTree EMtree3(bayesNet, EM_tree_matrix, EM_nodes_to_vars, train_samples);
  EMtree3.LearnEMParameters(1, 0);

  TensorJTree jtree(*bayesNet, tree_matrix, nodes_to_vars, 6);
  jtree.LearnSpectralParameters();
  vector<TensorCPT*> marginals;
  for (int n = 0; n < test_samples.Dim(0); ++n)
  {
	cout << n;
	cout << "\n";
	for (int i = 0; i < evidence_vals.size(); ++i)
	{
		evidence_vals.at(i) = test_samples.At(n,evidence_vars[i]);
	}
  
	double prob = TensorVarElim::VE(*bayesNet, evidence_vars, evidence_vals, bayesNet->ReverseTopologicalOrder());
	
	double EM_val = EMtree.ComputeEmpiricalMarginals(marginals, evidence_vars, evidence_vals);
	double spectral_val = jtree.ComputeMarginalEmpiricalProbability(evidence_vars, evidence_vals);
	if (spectral_val < 0)
		spectral_val = 0;
	if (spectral_val > 1)
		spectral_val = 1;
	//cout << val;
	//cout << "\n";
	total_error += abs(prob - spectral_val) / prob; 
	double EM_val = EMtree.ComputeEmpiricalMarginals(marginals, evidence_vars, evidence_vals);
	double EM_val2 = EMtree2.ComputeEmpiricalMarginals(marginals, evidence_vars, evidence_vals);
	double EM_val3 = EMtree3.ComputeEmpiricalMarginals(marginals, evidence_vars, evidence_vals);
  }
  cout << "total error:\n";
  cout << total_error / test_samples.Dim(0);
  cout << "\n";
  cout << EM_error / test_samples.Dim(0);

  return 0;

}