#include "tensor.hpp"
#include "BayesianNetwork.hpp"
#include "TemplateBayesianNetwork.hpp"
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
  srand(30);

  int num_hidden_nodes = 15;
  int num_hidden_states = 2;

  int num_observed_nodes = 17;
  int num_observed_states = 6;
  // node names

  vector<int> hidden_nodes;
  vector<int> observed_nodes;
  VectorPlus::Seq(hidden_nodes, 0,1, num_hidden_nodes); 
  VectorPlus::Seq(observed_nodes, num_hidden_nodes, 1, num_hidden_nodes + num_observed_nodes);
  
  
  int num_nodes = num_hidden_nodes + num_observed_nodes;
  vector<vector<int>*> parent_list;
  vector<int> type_vector(num_nodes, 0);
    vector<int> dims(num_nodes, 0);
  for (int i = 0; i < num_nodes; ++i)
	parent_list.push_back(new vector<int>());

  vector<int> template_vec(num_nodes, 0);
  for (int i = 0; i < num_hidden_nodes; ++i)
  {
	  type_vector[hidden_nodes[i]] = HIDDEN_FLAG;
	  dims[hidden_nodes[i]] = num_hidden_states;
	  if (i >= 1)
	  {
		parent_list[hidden_nodes[i]]->push_back(hidden_nodes[i-1]);
		template_vec[hidden_nodes[i]] = 1;
	  }
	  if (i >= 2)
	  {
		  parent_list[hidden_nodes[i]]->push_back(hidden_nodes[i-2]);
		  template_vec[hidden_nodes[i]] = 2;
	  }
  }

  for (int i = 0; i < num_observed_nodes; ++i)
  {
	  type_vector[observed_nodes[i]] = OBSERVED_FLAG; 
	  dims[observed_nodes[i]] = num_observed_states; 
	 
	  

	  if (i == num_observed_nodes - 2)
	  {
		  parent_list[observed_nodes[i]]->push_back(hidden_nodes[0]);
		  parent_list[observed_nodes[i]]->push_back(hidden_nodes[1]);
		  template_vec[observed_nodes[i]] = 3;
	  }
	  else if ( i == num_observed_nodes - 1)
	  {
		  parent_list[observed_nodes[i]]->push_back(hidden_nodes[hidden_nodes.size() - 1]);
		  template_vec[observed_nodes[i]] = 4;
	  }
	  else
	  {
		parent_list[observed_nodes[i]]->push_back(hidden_nodes[i]);
		template_vec[observed_nodes[i]] = 4;
	  }
  }


  TemplateBayesianNetwork* bayesNet = new TemplateBayesianNetwork(parent_list, type_vector, dims, template_vec);

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

 /*  for (int i = 0; i < num_nodes; ++i)
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
   } */

  Matrix train_samples;
  Matrix test_samples;
  bayesNet->GenerateTrainSamples(train_samples, 100);
  bayesNet->GenerateTestSamples(test_samples, 500);

  // test the junction tree
  int jtree_nodes = num_observed_nodes - 2;
  vector<int> tree_matrix_dims(2,jtree_nodes);
  Matrix tree_matrix(tree_matrix_dims);

  for (int i = 1; i < jtree_nodes; ++i)
  {
	  tree_matrix.Set(i, i - 1, 1);
  }
 

    map<int, vector<int>*> nodes_to_vars;
  for (int i = 0; i < jtree_nodes; ++i)
	nodes_to_vars[i] = new vector<int>();

  for (int i = 0; i < jtree_nodes; ++i)
  {
	 // if (i != 0)
		nodes_to_vars[i]->push_back(hidden_nodes[i]);
	  if (i + 1 < hidden_nodes.size())
		  nodes_to_vars[i]->push_back(hidden_nodes[i+1]);
	  if (i + 2 < hidden_nodes.size())
		  nodes_to_vars[i]->push_back(hidden_nodes[i+2]);
	  nodes_to_vars[i]->push_back(observed_nodes[i]);

	  if ( i == 0)
		  nodes_to_vars[i]->push_back(observed_nodes[hidden_nodes.size()]);
	  if (i == jtree_nodes - 1)
		  nodes_to_vars[i]->push_back(observed_nodes[hidden_nodes.size() + 1]);
  }
  
  vector<int> jtemplate_vec(jtree_nodes, 0);
  jtemplate_vec[0] = 0;
  jtemplate_vec[1] = 1;
  for (int i = 2; i < 6; ++i)
  {
	jtemplate_vec[i] = 2;
  }
  jtemplate_vec[6] = 3;
  jtemplate_vec[7] = 4; 
//  vector<int> jtemplate_vec;
 // VectorPlus::Seq(jtemplate_vec, 0, 1, jtree_nodes);

  vector<int> evidence_vars;
  VectorPlus::Copy(evidence_vars, observed_nodes);
  vector<int> evidence_vals(evidence_vars.size(), 0);

  JTree EMtree2(bayesNet, tree_matrix, nodes_to_vars, &train_samples);
  EMtree2.LearnEMParameters(1000, 0);

  JTree EMtree3(bayesNet, tree_matrix, nodes_to_vars, &train_samples);
  EMtree3.LearnEMParameters(1000, 0);

  JTree EMtree(bayesNet, tree_matrix, nodes_to_vars, &train_samples);
  EMtree.LearnEMParameters(200, 1);

vector<TensorCPT*> marginals;
double total_error = 0;
double EM_error1 = 0;
double EM_error2 = 0;
double EM_error3 = 0;

 //TensorJTree jtree(bayesNet, tree_matrix, nodes_to_vars, 1, jtemplate_vec);
 //jtree.LearnSpectralParameters();

for (int n = 0; n < test_samples.Dim(0); ++n)
{
	for (int i = 0; i < evidence_vals.size(); ++i)
	{
		evidence_vals.at(i) = (int)test_samples.At(n,evidence_vars[i]);
	}
  
	double prob = TensorVarElim::VE(*bayesNet, evidence_vars, evidence_vals, bayesNet->ReverseTopologicalOrder());
//	double spectral_val = jtree.ComputeMarginalEmpiricalProbability(evidence_vars, evidence_vals);
//	if (spectral_val < 0)
	//	spectral_val = 0;
	//if (spectral_val > 1)
	//	spectral_val = 1;

	double EM_val = EMtree.ComputeEmpiricalMarginals(marginals, evidence_vars, evidence_vals);
	double EM_val2 = EMtree2.ComputeEmpiricalMarginals(marginals, evidence_vars, evidence_vals);
	double EM_val3 = EMtree3.ComputeEmpiricalMarginals(marginals, evidence_vars, evidence_vals);
	
//	cout << spectral_val;
//	cout << "\t";
//	cout << prob;
//	cout << "\n";
	//cout << val;
	//cout << "\n";
	/*cout << EM_val;
	cout << "\t";
	cout << EM_val2;
	cout << "\t";
	cout << EM_val3; 
	cout << "\n";*/

	//total_error += min(abs(prob - spectral_val) / prob, 1.0);
	EM_error1 += min(abs(prob - EM_val) / prob, 1.0);
	EM_error2 += min(abs(prob - EM_val2) / prob, 1.0);
	EM_error3 += min(abs(prob - EM_val3) / prob, 1.0);

  }
 //   cout << "total error:\n";
//	cout << total_error / test_samples.Dim(0);
//	cout << "\n";
  cout << EM_error1 / test_samples.Dim(0);
  cout << "\n";
  cout << EM_error2 / test_samples.Dim(0);
  cout << "\n";
  cout << EM_error3 / test_samples.Dim(0);
  cout << "\n"; 
  return 0;
}