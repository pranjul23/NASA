#include "TensorJTree.hpp"
#include "VectorPlus.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include "JTree.hpp"
#include <time.h>
#include <assert.h>
#include "TensorVarElim.hpp"
#include "tensor.hpp"
#include "BayesianNetwork.hpp"
#include "TemplateBayesianNetwork.hpp"

using namespace std;

int main ()
{
  srand(30);

  // model parameters
  int num_hidden_level1 = 10;
  int num_hidden_level2 = num_hidden_level1 + 1;
  int num_obs_per_node = 2;

  int num_hidden_states = 2;
  int num_observed_states = 4;

  // EM iterations
  int num_EM_iterations = 50;

  // spectral observations
  int stack_count = 1;

  vector<int> hidden_nodes;
  vector<int> observed_nodes;
  VectorPlus::Seq(hidden_nodes, 0,1, num_hidden_level1 + num_hidden_level2);
  int num_hidden_nodes = hidden_nodes.size();
  VectorPlus::Seq(observed_nodes, num_hidden_nodes, 1, num_hidden_nodes + num_hidden_level2 * num_obs_per_node);
  int num_observed_nodes = observed_nodes.size();
  
  int num_nodes = num_hidden_nodes + num_observed_nodes;
  vector<vector<int>*> parent_list;
  vector<int> type_vector(num_nodes, 0);
  vector<int> dims(num_nodes, 0);
  for (int i = 0; i < num_nodes; ++i)
  {
	parent_list.push_back(new vector<int>());
	if (i < num_hidden_nodes)
	{
		type_vector[i] = HIDDEN_FLAG;
		dims[i] = num_hidden_states;
	}
	else
	{
		type_vector[i] = OBSERVED_FLAG;
		dims[i] = num_observed_states;
	}
  }

  vector<int> template_vec;
  VectorPlus::Seq(template_vec, 0, 1, num_nodes);

  for (int i = 1; i < num_hidden_level1; ++i)
  {
	  parent_list.at(i)->push_back(i-1);
  }

  for (int i = 0; i < num_hidden_level2; ++i)
  {
	  if (i == 0)
	  {
		parent_list.at(num_hidden_level1 + i)->push_back(i);
	  }
	  else if (i == num_hidden_level2 - 1)
	  {
		  parent_list.at(num_hidden_level1 + i)->push_back(i-1);
	  }
		
	  else
	  {
		  parent_list.at(num_hidden_level1 + i)->push_back(i - 1);
		  parent_list.at(num_hidden_level1 + i)->push_back(i);
	  }
  }


  for (int i = 0; i < num_observed_nodes; ++i)
  {
	  int hidden_parent = (int)(i / num_obs_per_node) + num_hidden_level1;
	  parent_list.at(observed_nodes[i])->push_back(hidden_parent);
  }



  TemplateBayesianNetwork* bayesNet = new TemplateBayesianNetwork(parent_list, type_vector, dims, template_vec);

 


  Matrix* train_samples = new Matrix();
  Matrix* test_samples = new Matrix();
  bayesNet->GenerateTrainSamples(*train_samples, 500);
  bayesNet->GenerateTestSamples(*test_samples, 500);

  // test the junction tree
//  int jtree_nodes = num_nodes;
//  vector<int> tree_matrix_dims(2,jtree_nodes);
  Matrix tree_matrix;
  vector<vector<int>*> nodes_to_vars;

  bayesNet->FindCliques(nodes_to_vars);
  bayesNet->FindJunctionTree(tree_matrix, nodes_to_vars);


 /* for (int i = 0; i < jtree_nodes; ++i)
	nodes_to_vars[i] = new vector<int>();

  for (int i = 0; i < num_hidden_level1; ++i)
  {
	  if (i == 0)
	  {
		  nodes_to_vars[i]->push_back(i);
		  nodes_to_vars[i]->push_back(i+1);		
		  tree_matrix.Set(i, i+1, 1);
	  }
	  else if (i == num_hidden_level1 -1)
	  {
		  nodes_to_vars[i]->push_back(i);
	  }
	  else
	  {
		  nodes_to_vars[i]->push_back(i);
		  nodes_to_vars[i]->push_back(i+1);
		  tree_matrix.Set(i, i+1, 1);
	  }
  }

  for (int i = 0; i < num_hidden_level2; ++i)
  {
	  int node_id = i + num_hidden_level1;
	  nodes_to_vars[node_id]->push_back(node_id);
	  for (int j = 0; j < parent_list.at(node_id)->size(); ++j)
	  {
		  nodes_to_vars[node_id]->push_back(parent_list.at(node_id)->at(j));
	  }

	  if (i == 0)
	  {
		  tree_matrix.Set(node_id, i, 1);
	  }
	  else if (i == num_hidden_level2 - 1)
	  {
		  tree_matrix.Set(node_id, i - 1, 1);
	  }
	  else
	  {
		  tree_matrix.Set(node_id, i - 1, 1);
	  }
  }

  for (int i = 0; i < num_observed_nodes;  ++i)
  {
	  int parent = parent_list.at(observed_nodes[i])->at(0);
	  tree_matrix.Set(observed_nodes[i], parent, 1);

	  nodes_to_vars[observed_nodes[i]]->push_back(observed_nodes[i]);
	  nodes_to_vars[observed_nodes[i]]->push_back(parent);
  }

 cout << tree_matrix.At(0, 1);
 cout << tree_matrix.At(2, 0);
 cout << tree_matrix.At(3, 0);
 cout << tree_matrix.At(4, 1);
  */

 vector<int> jtemplate_vec;
 VectorPlus::Seq(jtemplate_vec, 0, 1, nodes_to_vars.size());

 vector<int> evidence_vars;
 VectorPlus::Copy(evidence_vars, observed_nodes);

 vector<int> evidence_vals(evidence_vars.size(), 0);

TensorJTree jtree(bayesNet, tree_matrix, nodes_to_vars, 1, jtemplate_vec);
//jtree.LearnSpectralParameters();

 JTree EMtree(bayesNet, tree_matrix, nodes_to_vars, train_samples);
 EMtree.LearnEMParameters(num_EM_iterations, 1);

 JTree EMtree2(bayesNet, tree_matrix, nodes_to_vars, train_samples);
 EMtree2.LearnEMParameters(num_EM_iterations, 0);

 JTree EMtree3(bayesNet, tree_matrix, nodes_to_vars, train_samples);
 EMtree3.LearnEMParameters(num_EM_iterations, 0);

vector<TensorCPT*> marginals;
double total_error = 0;
double EM_error1 = 0;
double EM_error2 = 0;
double EM_error3 = 0;



for (int n = 0; n < test_samples->Dim(0); ++n)
{
	for (int i = 0; i < evidence_vals.size(); ++i)
	{
		evidence_vals.at(i) = (int)test_samples->At(n,evidence_vars[i]);
	}
  
	double prob = TensorVarElim::VE(*bayesNet, evidence_vars, evidence_vals, bayesNet->ReverseTopologicalOrder());
	// double spectral_val = jtree.ComputeMarginalEmpiricalProbability(evidence_vars, evidence_vals);
//	if (spectral_val < 0)
//		spectral_val = 0;
//	if (spectral_val > 1)
//		spectral_val = 1;

	double EM_val = EMtree.ComputeEmpiricalMarginals(marginals, evidence_vars, evidence_vals);
	double EM_val2 = EMtree2.ComputeEmpiricalMarginals(marginals, evidence_vars, evidence_vals);
	double EM_val3 = EMtree3.ComputeEmpiricalMarginals(marginals, evidence_vars, evidence_vals);
	
	/*cout << spectral_val;
	cout << "\t";
	cout << prob;
	cout << "\n"; */
	//cout << val;
	//cout << "\n";
//	cout << EM_val;
//	cout << "\t";
//	cout << EM_val2;
//	cout << "\t";
//	cout << EM_val3; 
//	cout << "\n";

//	total_error += min(abs(prob - spectral_val) / prob, 1.0);
	EM_error1 += min(abs(prob - EM_val) / prob, 1.0);
	EM_error2 += min(abs(prob - EM_val2) / prob, 1.0);
	EM_error3 += min(abs(prob - EM_val3) / prob, 1.0);

  }
    cout << "total error:\n";
	cout << total_error / test_samples->Dim(0);
	cout << "\n";
  cout << EM_error1 / test_samples->Dim(0);
  cout << "\n";
  cout << EM_error2 / test_samples->Dim(0);
  cout << "\n";
  cout << EM_error3 / test_samples->Dim(0);
  cout << "\n";

  return 0;
}