#include "tensor.hpp"
#include "BayesianNetwork.hpp"
#include "TensorVarElim.hpp"
#include "TensorJTree.hpp"
#include "VectorPlus.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include "JTree.hpp"
#include <time.h>
#include <assert.h>

using namespace std;

int main ()
{
  srand(time(NULL));
  // node names
  int A = 0;
  int B = 1;
 

  int O1 = 2;
  int O2 = 3;
  int O3 = 4;
  int O4 = 5;
  int O5 = 6;
  int O6 = 7;
  int O7 = 8;
  int O8 = 9;
  int O9 = 10;
  int O10 = 11;
  int O11 = 12;
  int O12 = 13;
  int O13 = 14;
  int O14 = 15;
  int O15 = 16;
  int O16 = 17;
  int O17 = 18;
  int O18 = 19;

  int num_nodes = 20;
  vector<vector<int>*> parent_list;
  for (int i = 0; i < num_nodes; ++i)
	parent_list.push_back(new vector<int>());

 // parent_list[B]->push_back(A);
 // parent_list[C]->push_back(A);
 // parent_list[D]->push_back(B);
 // parent_list[D]->push_back(A);
 // parent_list[E]->push_back(B);

  // hidden variables
  parent_list[O1]->push_back(A);
  parent_list[O2]->push_back(A);
  parent_list[O3]->push_back(A);
  parent_list[O4]->push_back(A);
  parent_list[O5]->push_back(A);
  parent_list[O6]->push_back(A);
  parent_list[O7]->push_back(A);
  parent_list[O8]->push_back(A);
  parent_list[O9]->push_back(A);
  parent_list[O10]->push_back(A);
  parent_list[O11]->push_back(A);
  parent_list[O12]->push_back(A);
  parent_list[O13]->push_back(A);
  parent_list[O14]->push_back(A);
  parent_list[O15]->push_back(A);
  parent_list[O16]->push_back(A);
  parent_list[O17]->push_back(A);
  parent_list[O18]->push_back(A);

  parent_list[O1]->push_back(B);
  parent_list[O2]->push_back(B);
  parent_list[O3]->push_back(B);
  parent_list[O4]->push_back(B);
  parent_list[O5]->push_back(B);
  parent_list[O6]->push_back(B);
  parent_list[O7]->push_back(B);
  parent_list[O8]->push_back(B);
  parent_list[O9]->push_back(B);
  parent_list[O10]->push_back(B);
  parent_list[O11]->push_back(B);
  parent_list[O12]->push_back(B);
  parent_list[O13]->push_back(B);
  parent_list[O14]->push_back(B);
  parent_list[O15]->push_back(B);
  parent_list[O16]->push_back(B);
  parent_list[O17]->push_back(B);
  parent_list[O18]->push_back(B);

  vector<int> type_vector(num_nodes, 0);

  for (int i = 1; i < num_nodes; ++i)
  {
	  if (i <= 1)
		type_vector[i] = HIDDEN_FLAG;
	  else
		type_vector[i] = OBSERVED_FLAG;
  }


  vector<int> dims(num_nodes, 4);
  dims[0] = 2;
  dims[1] = 2;
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

   for (int i = 0; i < num_nodes; ++i)
  {
	  double noise = ((double)rand()) / (8 *RAND_MAX);
//	  assert(noise < .15);
	  if (i < 2)
	  {
		  bayesNet->SetCPT(i, 0, .5 + noise);
		  bayesNet->SetCPT(i, 1, .5 - noise);

	  }
	  else
	  {
		  bayesNet->SetCPT(i, 0, (1 - noise) / 2);
		  bayesNet->SetCPT(i, 1, noise / 2);
		  bayesNet->SetCPT(i, 2, noise / 2);
		  bayesNet->SetCPT(i, 3, (1 - noise) / 2);
		  bayesNet->SetCPT(i, 4, (1 - noise) / 2);
		  bayesNet->SetCPT(i, 5, noise / 2);
		  bayesNet->SetCPT(i, 6, noise / 2);
		  bayesNet->SetCPT(i, 7, (1 - noise) / 2);


		  bayesNet->SetCPT(i, 8, noise / 2);
		  bayesNet->SetCPT(i, 9, (1 - noise) / 2);
		  bayesNet->SetCPT(i, 10, (1 - noise) / 2);
		  bayesNet->SetCPT(i, 11, noise / 2);
		  bayesNet->SetCPT(i, 12, noise / 2);
		  bayesNet->SetCPT(i, 13, (1 - noise) / 2);
		  bayesNet->SetCPT(i, 14, (1 - noise) / 2);
		  bayesNet->SetCPT(i, 15, noise / 2);
	  }
   }

  Matrix train_samples;
  Matrix test_samples;
  bayesNet->GenerateTrainSamples(train_samples, 500);
  bayesNet->GenerateTestSamples(test_samples, 500);

  // test the junction tree
  int jtree_nodes = 19;
  vector<int> tree_matrix_dims(2,jtree_nodes);
  
  //int A = 0;
  int AB = 0;
  
  int O1AB = 1;
  int O2AB = 2;
  int O3AB = 3;
  int O4AB = 4;
  int O5AB = 5;
  int O6AB = 6;
  int O7AB = 7;
  int O8AB = 8;
  int O9AB = 9;
  int O10AB = 10;
  int O11AB = 11;
  int O12AB = 12;
  int O13AB = 13;
  int O14AB = 14;
  int O15AB = 15;
  int O16AB = 16;
  int O17AB = 17;
  int O18AB = 18;

  Matrix tree_matrix(tree_matrix_dims);
  tree_matrix.Set(O1AB , AB, 1);
  tree_matrix.Set(O2AB , AB, 1);
  tree_matrix.Set(O3AB, AB, 1);
  tree_matrix.Set(O4AB, AB, 1);
  tree_matrix.Set(O5AB, AB, 1);
  tree_matrix.Set(O6AB, AB, 1);
  tree_matrix.Set(O7AB , AB, 1);
  tree_matrix.Set(O8AB , AB, 1);
  tree_matrix.Set(O9AB, AB, 1);
  tree_matrix.Set(O10AB, AB, 1);
  tree_matrix.Set(O11AB, AB, 1);
  tree_matrix.Set(O12AB, AB, 1);
  tree_matrix.Set(O13AB , AB, 1);
  tree_matrix.Set(O14AB , AB, 1);
  tree_matrix.Set(O15AB, AB, 1);
  tree_matrix.Set(O16AB, AB, 1);
  tree_matrix.Set(O17AB, AB, 1);
  tree_matrix.Set(O18AB, AB, 1);

  map<int, vector<int>*> nodes_to_vars;
  for (int i = 0; i < jtree_nodes; ++i)
	nodes_to_vars[i] = new vector<int>();
  

  nodes_to_vars[AB]->push_back(A);
  nodes_to_vars[AB]->push_back(B);

  nodes_to_vars[O1AB]->push_back(O1);
  nodes_to_vars[O1AB]->push_back(A);
  nodes_to_vars[O1AB]->push_back(B);

  nodes_to_vars[O2AB]->push_back(O2);
  nodes_to_vars[O2AB]->push_back(A);
  nodes_to_vars[O2AB]->push_back(B);

  nodes_to_vars[O3AB]->push_back(O3);
  nodes_to_vars[O3AB]->push_back(A);
  nodes_to_vars[O3AB]->push_back(B);

  nodes_to_vars[O4AB]->push_back(O4);
  nodes_to_vars[O4AB]->push_back(A);
  nodes_to_vars[O4AB]->push_back(B);

  nodes_to_vars[O5AB]->push_back(O5);
  nodes_to_vars[O5AB]->push_back(A);
  nodes_to_vars[O5AB]->push_back(B);

  nodes_to_vars[O6AB]->push_back(O6);
  nodes_to_vars[O6AB]->push_back(A);
  nodes_to_vars[O6AB]->push_back(B);

  nodes_to_vars[O7AB]->push_back(O7);
  nodes_to_vars[O7AB]->push_back(A);
  nodes_to_vars[O7AB]->push_back(B);

  nodes_to_vars[O8AB]->push_back(O8);
  nodes_to_vars[O8AB]->push_back(A);
  nodes_to_vars[O8AB]->push_back(B);

  nodes_to_vars[O9AB]->push_back(O9);
  nodes_to_vars[O9AB]->push_back(A);
  nodes_to_vars[O9AB]->push_back(B);

  nodes_to_vars[O10AB]->push_back(O10);
  nodes_to_vars[O10AB]->push_back(A);
  nodes_to_vars[O10AB]->push_back(B);

  nodes_to_vars[O11AB]->push_back(O11);
  nodes_to_vars[O11AB]->push_back(A);
  nodes_to_vars[O11AB]->push_back(B);

  nodes_to_vars[O12AB]->push_back(O12);
  nodes_to_vars[O12AB]->push_back(A);
  nodes_to_vars[O12AB]->push_back(B);

    nodes_to_vars[O13AB]->push_back(O13);
  nodes_to_vars[O13AB]->push_back(A);
  nodes_to_vars[O13AB]->push_back(B);

  nodes_to_vars[O14AB]->push_back(O14);
  nodes_to_vars[O14AB]->push_back(A);
  nodes_to_vars[O14AB]->push_back(B);

  nodes_to_vars[O15AB]->push_back(O15);
  nodes_to_vars[O15AB]->push_back(A);
  nodes_to_vars[O15AB]->push_back(B);

  nodes_to_vars[O16AB]->push_back(O16);
  nodes_to_vars[O16AB]->push_back(A);
  nodes_to_vars[O16AB]->push_back(B);

  nodes_to_vars[O17AB]->push_back(O17);
  nodes_to_vars[O17AB]->push_back(A);
  nodes_to_vars[O17AB]->push_back(B);

  nodes_to_vars[O18AB]->push_back(O18);
  nodes_to_vars[O18AB]->push_back(A);
  nodes_to_vars[O18AB]->push_back(B);

  vector<int> evidence_vars;
  VectorPlus::Seq(evidence_vars, O1, 1, O18 + 1);
  vector<int> evidence_vals(evidence_vars.size(), 0);

  JTree EMtree(bayesNet, tree_matrix, nodes_to_vars, train_samples);
  EMtree.LearnEMParameters();

  TensorJTree jtree(bayesNet, tree_matrix, nodes_to_vars, 3);
  jtree.LearnSpectralParameters();
vector<TensorCPT*> marginals;
double total_error = 0;
double EM_error = 0;
  for (int n = 0; n < test_samples.Dim(0); ++n)
  {

	for (int i = 0; i < evidence_vals.size(); ++i)
	{
		evidence_vals.at(i) = (int)test_samples.At(n,evidence_vars[i]);
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
	cout << EM_val;
	cout << "\t";
	cout << spectral_val;
	cout << "\t";
	cout << prob; 
	cout << "\n";

	total_error += min(abs(prob - spectral_val) / prob, 1.0); 
	EM_error += min(abs(prob - EM_val) / prob, 1.0);
  }
    cout << "total error:\n";
  cout << total_error / test_samples.Dim(0);
  cout << "\n";
  cout << EM_error / test_samples.Dim(0);
  return 0;
}