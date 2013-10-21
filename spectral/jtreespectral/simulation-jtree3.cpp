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
#include <cmath>
#include "StatsLibrary.hpp"

using namespace std;

int main (int argc, char *argv[])
{

  assert(argc == 10);
  int hidden_depth = atoi(argv[1]);
  int num_observed_nodes_per_leaf = atoi(argv[2]);
  int num_hidden_states = atoi(argv[3]);
  int num_observed_states = atoi(argv[4]);
  int num_training_examples = atoi(argv[5]);
  double EM_thresh = atof(argv[6]);
  int stack_count = atoi(argv[7]);
  int max_runs = atoi(argv[8]);
  int trial = atoi(argv[9]);

 /* int hidden_depth = 2;
  int num_observed_nodes_per_leaf = 6;
  int num_hidden_states = 2;
  int num_observed_states = 2;
  int num_training_examples = 100;
  double EM_thresh = 0.001;
  int stack_count = 1;
  int max_runs = 2;
  int trial = 1; */

  int rand_seed = trial;
  srand(rand_seed);

  int hidden_R_nodes_per_clique = 3;

  int num_hidden_cliques = (int)pow(2.0, (double)(hidden_depth)) - 1;
  int num_hidden_leaves = (int)(num_hidden_cliques + 1)/2;
  int num_hidden_interior = num_hidden_cliques - num_hidden_leaves;
  int num_hidden_nodes = 3 * num_hidden_interior + 2 * num_hidden_leaves;
  int num_observed_nodes = num_hidden_leaves * num_observed_nodes_per_leaf;

  vector<int> hidden_nodes;
  vector<int> observed_nodes;
  VectorPlus::Seq(hidden_nodes, 0,1, num_hidden_nodes); 
  VectorPlus::Seq(observed_nodes, num_hidden_nodes, 1, num_hidden_nodes + num_observed_nodes);
  int num_nodes = num_hidden_nodes + num_observed_nodes;

  // make junction tree first
  int jtree_nodes = num_hidden_cliques + num_observed_nodes;
  vector<int> tree_matrix_dims(2,jtree_nodes);
  Matrix tree_matrix(tree_matrix_dims);

  for (int i = 1; i < jtree_nodes; ++i)
  {
	  if (i < num_hidden_cliques)
	  {
		tree_matrix.Set(i, (i-1)/2, 1);
	  }
	  else
	  {
		 int obs_index = (i - num_hidden_cliques) / num_observed_nodes_per_leaf;
		 tree_matrix.Set(i, num_hidden_interior + obs_index, 1);
	  }
  }
 
  vector<vector<int>*> parent_list;
  vector<int> type_vector(num_nodes, 0);
  vector<int> dims(num_nodes, 0);
  
  for (int i = 0; i < num_nodes; ++i)
  {
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
  
  for (int i = 0; i < num_nodes; ++i)
	parent_list.push_back(new vector<int>());

  vector<vector<int>*> nodes_to_vars;
  map<int, int> nodes_to_cliques;
  for (int i = 0; i < jtree_nodes; ++i)
	  nodes_to_vars.push_back(new vector<int>());

  int curr_index = 0;
  for (int i = 0; i < jtree_nodes;  ++i)
  {
	  if (i == 0)
	  {
		  for (int j = 0; j < 3; ++j)
		  {
			nodes_to_vars[i]->push_back(curr_index);
			curr_index++;
		  }
	  }

	  else if (i < num_hidden_cliques)
	  {
		  int parent = (i-1)/2;
		  vector<int>& parent_nodes = *(nodes_to_vars[parent]);

		  if (i % 2 == 1) // left child
		  {
			  for (int j = parent_nodes.size() - 3; j < parent_nodes.size() - 1; ++j)
			  {
				  nodes_to_vars[i]->push_back(parent_nodes[j]);
			  }
		  }
		  else
		  {
			  for (int j = parent_nodes.size() - 2; j < parent_nodes.size(); ++j)
			  {
				  nodes_to_vars[i]->push_back(parent_nodes[j]);
			  }
		  }

		  if (i < (num_hidden_cliques - num_hidden_leaves))
		  {
			for (int j = 0; j < 3; ++j)
			{
				nodes_to_vars[i]->push_back(curr_index);
				curr_index++;
			}
		  }
		  else
		  {
			for (int j = 0; j < 2; ++j)
			{
				nodes_to_vars[i]->push_back(curr_index);
				curr_index++;
			}

		  }
	  }
	  else
	  {
		  vector<int> parent_vec;
		  tree_matrix.RowFind(parent_vec, i);
		  int parent = parent_vec[0];
		  vector<int>& parent_nodes = *(nodes_to_vars[parent]);

		  for (int j = parent_nodes.size() - 2; j < parent_nodes.size(); ++j)
		  {
			  nodes_to_vars[i]->push_back(parent_nodes[j]);
		  }

		  nodes_to_vars[i]->push_back(curr_index++);
	  }
  }


  for (int i =0; i < jtree_nodes; ++i)
  {
	  vector<int>& clique_nodes = *(nodes_to_vars)[i];
	  for (int j = 0; j < clique_nodes.size(); ++j)
	  {
		  for (int k = j + 1; k < clique_nodes.size(); ++k)
		  {
			  if (!VectorPlus::Contains(*(parent_list[clique_nodes[k]]), clique_nodes[j]))
				parent_list[clique_nodes[k]]->push_back(clique_nodes[j]);
		  }
	  }
  }
 
  vector<int> template_vec;
  VectorPlus::Seq(template_vec, 0, 1, num_nodes);

  TemplateBayesianNetwork* bayesNet = new TemplateBayesianNetwork(parent_list, type_vector, dims, template_vec);

  Matrix* train_samples = new Matrix();
  Matrix* test_samples = new Matrix();
  bayesNet->GenerateTrainSamples(*train_samples, num_training_examples);
  bayesNet->GenerateTestSamples(*test_samples, 1000);

  vector<int> jtemplate_vec;
  VectorPlus::Seq(jtemplate_vec, 0, 1, jtree_nodes);

  
  vector<int> evidence_vars;
  VectorPlus::Copy(evidence_vars, observed_nodes);
  vector<int> evidence_vals(evidence_vars.size(), 0);

  TensorJTree jtree(bayesNet, tree_matrix, nodes_to_vars, stack_count, jtemplate_vec, max_runs);
  clock_t spectral_init = clock();
  jtree.LearnSpectralParameters();
  clock_t spectral_final = clock() - spectral_init;
  double spectral_time = ((double)spectral_final) / ((double) CLOCKS_PER_SEC);

  JTree EMtree(bayesNet, tree_matrix, nodes_to_vars, train_samples, EM_thresh);
  clock_t EM1_init = clock();
  EMtree.LearnEMParameters(1, 0);
  clock_t EM1_final = clock() - EM1_init;
  double EM1_time = ((double)EM1_final) / ((double) CLOCKS_PER_SEC);

  JTree EMtree2(bayesNet, tree_matrix, nodes_to_vars, train_samples, EM_thresh);
  clock_t EM2_init = clock();
  EMtree2.LearnEMParameters(50, 0);
  clock_t EM2_final = clock() - EM2_init;
  double EM2_time = ((double)EM2_final) / ((double) CLOCKS_PER_SEC);

  JTree EMtree3(bayesNet, tree_matrix, nodes_to_vars, train_samples, EM_thresh);
  clock_t EM3_init = clock();
  EMtree3.LearnEMParameters(10, 0);
  clock_t EM3_final = clock() - EM3_init;
  double EM3_time = ((double)EM3_final) / ((double) CLOCKS_PER_SEC);

  JTree EMtree4(bayesNet, tree_matrix, nodes_to_vars, train_samples, EM_thresh);
  clock_t EM4_init = clock();
  EMtree4.LearnEMParameters(50, 0);
  clock_t EM4_final = clock() - EM4_init;
  double EM4_time = ((double)EM4_final) / ((double) CLOCKS_PER_SEC);

  JTree EMtree5(bayesNet, tree_matrix, nodes_to_vars, train_samples, EM_thresh);
  clock_t EM5_init = clock();
  EMtree5.LearnEMParameters(10, 0);
  clock_t EM5_final = clock() - EM5_init;
  double EM5_time = ((double)EM5_final) / ((double) CLOCKS_PER_SEC);

vector<TensorCPT*> marginals;
double spectral_error = 0;
double EM_error1 = 0;
double EM_error2 = 0;
double EM_error3 = 0;
double EM_error4 = 0;
double EM_error5 = 0;

int win_count = 0;

for (int n = 0; n < test_samples->Dim(0); ++n)
{
	for (int i = 0; i < evidence_vals.size(); ++i)
	{
		evidence_vals.at(i) = (int)test_samples->At(n,evidence_vars[i]);
	}
  
	double prob = TensorVarElim::VE(*bayesNet, evidence_vars, evidence_vals, bayesNet->ReverseTopologicalOrder());
	 double spectral_val = jtree.ComputeMarginalEmpiricalProbability(evidence_vars, evidence_vals);
	if (spectral_val < 0)
		spectral_val = 0;
	if (spectral_val > 1)
		spectral_val = 1;

	double EM_val = EMtree.ComputeEmpiricalMarginals(marginals, evidence_vars, evidence_vals);
	double EM_val2 = EMtree2.ComputeEmpiricalMarginals(marginals, evidence_vars, evidence_vals);
	double EM_val3 = EMtree3.ComputeEmpiricalMarginals(marginals, evidence_vars, evidence_vals);
	double EM_val4 = EMtree4.ComputeEmpiricalMarginals(marginals, evidence_vars, evidence_vals);
	double EM_val5 = EMtree5.ComputeEmpiricalMarginals(marginals, evidence_vars, evidence_vals);

	//cout << spectral_val;
	//cout << "\t";
	//cout << prob;
	//cout << "\n";
	//cout << val;
	//cout << "\n";
//	cout << EM_val;
//	cout << "\t";
//	cout << EM_val2;
//	cout << "\t";
//	cout << EM_val3; 
//	cout << "\n";

	spectral_error += min(abs(prob - spectral_val) / prob, 1.0);
	EM_error1 += min(abs(prob - EM_val) / prob, 1.0);
	EM_error2 += min(abs(prob - EM_val2) / prob, 1.0);
	EM_error3 += min(abs(prob - EM_val3) / prob, 1.0);
	EM_error4 += min(abs(prob - EM_val4) / prob, 1.0);
	EM_error5 += min(abs(prob - EM_val5) / prob, 1.0);

	if (abs(prob - EM_val2) > abs(prob - spectral_val) )
		win_count++;

  }

  spectral_error = spectral_error /  test_samples->Dim(0);
  EM_error1 = EM_error1 / test_samples->Dim(0);
  EM_error2 = EM_error2 / test_samples->Dim(0);
  EM_error3 = EM_error3 / test_samples->Dim(0);
  EM_error4 = EM_error4 / test_samples->Dim(0);
  EM_error5 = EM_error5 / test_samples->Dim(0);


  cout << "total error:\n";
  cout << spectral_error;
  cout << "\n";
  cout << EM_error1;
  cout << "\n";
  cout << EM_error2;
  cout << "\n";
  cout << EM_error3;
  cout << "\n";
  cout << EM_error4;
  cout << "\n";
  cout << EM_error5;
  cout << "\n";
  cout << win_count;
  cout << "\n";

  string outfile_name = "overlapping-jtree-" + to_string(hidden_depth) + "-" + to_string(num_observed_nodes_per_leaf ) 
						+ "-" + to_string(num_hidden_states) + "-" 
						+  to_string(num_observed_states) + "-" + to_string(num_training_examples)
						+ "-" + to_string(EM_thresh) + "-" + to_string(stack_count) + "-" + to_string(max_runs)
						+ "-" + to_string(trial) + ".txt";
//  string outfile_name = "fun.txt";
  ofstream myfile;
  myfile.open(outfile_name.c_str());
  myfile << ("spectral\t" + to_string(spectral_error) + "\t" + to_string(spectral_time) + "\n");
  myfile << ("EM1\t" + to_string(EM_error1) + "\t" + to_string(EM1_time) + "\n");
  myfile << ("EM2\t" + to_string(EM_error2) + "\t" + to_string(EM2_time) + "\n");
  myfile << ("EM3\t" + to_string(EM_error3) + "\t" + to_string(EM3_time) + "\n");
  myfile << ("EM4\t" + to_string(EM_error4) + "\t" + to_string(EM4_time) + "\n");
  myfile << ("EM5\t" + to_string(EM_error5) + "\t" + to_string(EM5_time) + "\n");
  myfile << ("randseed\t" + to_string(rand_seed) + "\t" + to_string(rand_seed) + "\n");

  return 0;
}