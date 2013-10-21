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
#include "StatsLibrary.hpp"
#include <cmath>

using namespace std;

int main (int argc, char *argv[])
{
  //srand(30);
 /* assert(argc == 10);
  int num_hidden_level1 = atoi(argv[1]);
  int num_hidden_level2 = num_hidden_level1;
  int num_obs_per_pair = atoi(argv[2]);

  int num_hidden_states = atoi(argv[3]);
  int num_observed_states = atoi(argv[4]);

  int num_training_examples = atoi(argv[5]);

  // EM iterations
  double EM_thresh = atof(argv[6]);

  // spectral observations
  int stack_count = atoi(argv[7]);
  int max_runs = atoi(argv[8]);

  int trial = atoi(argv[9]); */



  // model parameters
  int num_hidden_level1 = 8;
  int num_hidden_level2 = num_hidden_level1;
  int num_obs_per_pair = 2;

  int num_hidden_states = 2;
  int num_observed_states = 2;

  int num_training_examples = 500;

  // EM iterations
  double EM_thresh = 0.001;

  // spectral observations
  int stack_count = 1;
  int max_runs = 1;

  int trial = 1;

  int rand_seed = 30 * trial;
//  int rand_seed = time(NULL) * trial;
  srand(rand_seed);

  int skip = 0;

  vector<int> hidden_nodes;
  vector<int> observed_nodes;
  VectorPlus::Seq(hidden_nodes, 0,1, num_hidden_level1 + num_hidden_level2);
  int num_hidden_nodes = hidden_nodes.size();
  VectorPlus::Seq(observed_nodes, num_hidden_nodes, 1, num_hidden_nodes + num_hidden_level1 * num_obs_per_pair);
  int num_observed_nodes = observed_nodes.size();
  
  int num_nodes = num_hidden_nodes + num_observed_nodes;
  vector<vector<int>*> j_parent_list;
  vector<vector<int>*> orig_parent_list;
  vector<int> type_vector(num_nodes, 0);
  vector<int> dims(num_nodes, 0);
  for (int i = 0; i < num_nodes; ++i)
  {
	j_parent_list.push_back(new vector<int>());
	orig_parent_list.push_back(new vector<int>());

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
	  j_parent_list.at(i)->push_back(i-1);
	  orig_parent_list.at(i)->push_back(i-1);
	  j_parent_list.at(i)->push_back(num_hidden_level1 + i - 1);
  }

  for (int i = 0; i < num_hidden_level2; ++i)
  {
	  if (i == 0)
	  {
		j_parent_list.at(num_hidden_level1 + i)->push_back(i);
	  }
	  else
	  {
	//	  j_parent_list.at(num_hidden_level1 + i)->push_back(i-1);
		  j_parent_list.at(num_hidden_level1 + i)->push_back(i);
		  j_parent_list.at(num_hidden_level1 + i)->push_back(num_hidden_level1 + i - 1);
		  orig_parent_list.at(num_hidden_level1 + i)->push_back(num_hidden_level1 + i - 1);
	  }
  }


  for (int i = 0; i < num_observed_nodes; ++i)
  {
	  int hidden_parent1 = (int)(i / num_obs_per_pair);
	  int hidden_parent2 = (int)(i / num_obs_per_pair) + num_hidden_level1;
	  j_parent_list.at(observed_nodes[i])->push_back(hidden_parent1);
	  j_parent_list.at(observed_nodes[i])->push_back(hidden_parent2);
	  orig_parent_list.at(observed_nodes[i])->push_back(hidden_parent1);
	  orig_parent_list.at(observed_nodes[i])->push_back(hidden_parent2);

	  if (i % num_obs_per_pair == 1)
	  {
	//	  j_parent_list.at(observed_nodes[i])->push_back(observed_nodes[i-1]);
	//	  orig_parent_list.at(observed_nodes[i])->push_back(observed_nodes[i-1]);
	  }
  }



  TemplateBayesianNetwork* j_bayesNet = new TemplateBayesianNetwork(j_parent_list, type_vector, dims, template_vec);
  TemplateBayesianNetwork* orig_bayesNet = new TemplateBayesianNetwork(orig_parent_list, type_vector, dims, template_vec);
  


  Matrix* train_samples = new Matrix();
  Matrix* test_samples = new Matrix();
  orig_bayesNet->GenerateTrainSamples(*train_samples, num_training_examples);
  orig_bayesNet->GenerateTestSamples(*test_samples, 1000);
  cout << "yee\n";
  // test the junction tree
//  int jtree_nodes = num_nodes;


  
  int first_spectral_nodes = (num_hidden_nodes - 1  - 2*skip);
  int second_spectral_nodes = (num_hidden_level1 - skip);
  int spectral_jtree_nodes =  first_spectral_nodes + second_spectral_nodes;
  vector<int> spectral_tree_matrix_dims(2, spectral_jtree_nodes);
  Matrix spectral_tree_matrix(spectral_tree_matrix_dims);
  vector<vector<int>*> spectral_nodes_to_vars;

  for (int i = 1; i < first_spectral_nodes; ++i)
  {
	  spectral_tree_matrix.Set(i, i - 1, 1);
  }

  for (int i = first_spectral_nodes; i < spectral_jtree_nodes; ++i)
  {
	  int parent = 2 * (i - first_spectral_nodes);
	  spectral_tree_matrix.Set(i, parent , 1);
  }

  for (int i = 0; i < spectral_jtree_nodes; ++i)
	  spectral_nodes_to_vars.push_back(new vector<int>());

  for (int i = 0; i < first_spectral_nodes;  ++i)
  {
	  int index = (int) (i / 2);
	  if (i == first_spectral_nodes - 1)
	  {
		  spectral_nodes_to_vars[i]->push_back(hidden_nodes[index + skip]);
		  spectral_nodes_to_vars[i]->push_back(hidden_nodes[index + skip + num_hidden_level1]);
	  }

	  else if (i % 2 == 0)
	  {
		spectral_nodes_to_vars[i]->push_back(hidden_nodes[index + skip]);
		spectral_nodes_to_vars[i]->push_back(hidden_nodes[index + skip + 1]);
		spectral_nodes_to_vars[i]->push_back(hidden_nodes[index + skip + num_hidden_level1]);
	  }
	  else
	  {
		spectral_nodes_to_vars[i]->push_back(hidden_nodes[index + skip + num_hidden_level1]);
		spectral_nodes_to_vars[i]->push_back(hidden_nodes[index + skip + 1]);
		spectral_nodes_to_vars[i]->push_back(hidden_nodes[index + skip + num_hidden_level1 + 1]);
	  }
  }

  for (int i = 0; i < second_spectral_nodes; ++i)
  {
    spectral_nodes_to_vars[i + first_spectral_nodes]->push_back(observed_nodes[num_obs_per_pair * (i + skip)]);
	spectral_nodes_to_vars[i + first_spectral_nodes]->push_back(observed_nodes[num_obs_per_pair * (i + skip) + 1]);
    spectral_nodes_to_vars[i + first_spectral_nodes]->push_back(hidden_nodes[i + skip]);
	spectral_nodes_to_vars[i + first_spectral_nodes]->push_back(hidden_nodes[i + skip + num_hidden_level1]);
  }

  Matrix EM_tree_matrix;
  vector<vector<int>*> EM_nodes_to_vars;

  j_bayesNet->FindCliques(EM_nodes_to_vars);
  j_bayesNet->FindJunctionTree(EM_tree_matrix, EM_nodes_to_vars);
  cout << "yee2\n";
 vector<int> jtemplate_vec;
 VectorPlus::Seq(jtemplate_vec, 0, 1, spectral_nodes_to_vars.size());

  vector<int> evidence_vars;
  vector<int> omitted_obs;
  for (int i = 0; i < skip; ++i)
  {
	omitted_obs.push_back(observed_nodes[num_obs_per_pair * i]);
	omitted_obs.push_back(observed_nodes[num_obs_per_pair * i  + 1]);
  }

  vector<int> elimination_order;
  for (int i = 0; i < num_hidden_level1; ++i)
  {
	  elimination_order.push_back(observed_nodes[i * num_obs_per_pair]);
	  elimination_order.push_back(observed_nodes[i * num_obs_per_pair + 1]);
	  elimination_order.push_back(hidden_nodes[num_hidden_level1 + i]);
	  elimination_order.push_back(hidden_nodes[i]);

  }
  orig_bayesNet->SetEliminationOrder(elimination_order);
  VectorPlus::SetDiff(evidence_vars, observed_nodes, omitted_obs);
  vector<int> evidence_vals(evidence_vars.size(), 0);

  TensorJTree jtree(orig_bayesNet, spectral_tree_matrix, spectral_nodes_to_vars, stack_count, jtemplate_vec, max_runs);
  clock_t spectral_init = clock();
  jtree.LearnSpectralParameters();
  clock_t spectral_final = clock() - spectral_init;
  double spectral_time = ((double)spectral_final) / ((double) CLOCKS_PER_SEC);

  JTree EMtree(orig_bayesNet, EM_tree_matrix, EM_nodes_to_vars, train_samples, EM_thresh, 0);
  clock_t EM1_init = clock();
  //EMtree.LearnEMParameters(1, 0);
  clock_t EM1_final = clock() - EM1_init;
  double EM1_time = ((double)EM1_final) / ((double) CLOCKS_PER_SEC);

  JTree EMtree2(orig_bayesNet, EM_tree_matrix, EM_nodes_to_vars, train_samples, EM_thresh, 0);
  clock_t EM2_init = clock();
  //EMtree2.LearnEMParameters(50, 0);
  clock_t EM2_final = clock() - EM2_init;
  double EM2_time = ((double)EM2_final) / ((double) CLOCKS_PER_SEC);

  JTree EMtree3(orig_bayesNet, EM_tree_matrix, EM_nodes_to_vars, train_samples, EM_thresh, 1);
  clock_t EM3_init = clock();
 // EMtree3.LearnEMParameters(10, 0);
  clock_t EM3_final = clock() - EM3_init;
  double EM3_time = ((double)EM3_final) / ((double) CLOCKS_PER_SEC);

  JTree EMtree4(orig_bayesNet, EM_tree_matrix, EM_nodes_to_vars, train_samples, EM_thresh, 1);
  clock_t EM4_init = clock();
 // EMtree4.LearnEMParameters(50, 0);
  clock_t EM4_final = clock() - EM4_init;
  double EM4_time = ((double)EM4_final) / ((double) CLOCKS_PER_SEC);

  JTree EMtree5(orig_bayesNet, EM_tree_matrix, EM_nodes_to_vars, train_samples, EM_thresh, 1);
  clock_t EM5_init = clock();
 //EMtree5.LearnEMParameters(10, 0);
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
  
	double prob = TensorVarElim::VE(*orig_bayesNet, evidence_vars, evidence_vals, elimination_order);
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

  string outfile_name = "factorial-hmm-" + to_string(num_hidden_level1) + "-" + to_string(num_hidden_states) + "-" 
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
