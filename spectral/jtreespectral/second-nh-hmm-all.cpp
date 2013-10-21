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
#include <stdlib.h>

using namespace std;

int main (int argc, char *argv[])
{

  assert(argc == 9);
  int num_hidden_nodes = atoi(argv[1]);
  int num_hidden_states = atoi(argv[2]);

  int num_observed_nodes = atoi(argv[1]);
  int num_observed_states = atoi(argv[3]);
  // node names

   int num_training_examples = atoi(argv[4]);

    // EM iterations
  double EM_thresh = atof(argv[5]);

  // spectral observations
  int stack_count = atoi(argv[6]);
  int max_runs = atoi(argv[7]);

  int trial = atoi(argv[8]); 

  int num_EM_trials = 5; 

 /* int num_hidden_nodes = 8;
  int num_hidden_states = 2;
  int num_observed_nodes = 8;
  int num_observed_states = 4;
  int num_training_examples = 10;
  double EM_thresh = 0.01;
  int stack_count = 1;
  int max_runs = 1;

  int trial = 1; */

  double EM_early_thresh = 0.01;
  vector<double> online_rates;
  online_rates.push_back(0.6);
  online_rates.push_back(0.7);
  online_rates.push_back(0.8);
  online_rates.push_back(0.9);
  online_rates.push_back(1);

  int rand_seed = 596 * trial;
  srand(rand_seed);

   int skip = 1;

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
		template_vec[hidden_nodes[i]] = i;
	  }
	  if (i >= 2)
	  {
		  parent_list[hidden_nodes[i]]->push_back(hidden_nodes[i-2]);
		  template_vec[hidden_nodes[i]] = i;
	  }
  }

  for (int i = 0; i < num_observed_nodes; ++i)
  {
	  type_vector[observed_nodes[i]] = OBSERVED_FLAG; 
	  dims[observed_nodes[i]] = num_observed_states; 
	  parent_list[observed_nodes[i]]->push_back(hidden_nodes[i]);
	  template_vec[observed_nodes[i]] = observed_nodes[i];
  }


  TemplateBayesianNetwork* bayesNet = new TemplateBayesianNetwork(parent_list, type_vector, dims, template_vec);

  Matrix* train_samples = new Matrix();
  Matrix* test_samples = new Matrix();
  bayesNet->GenerateTestSamples(*test_samples, 1000);
  bayesNet->GenerateTrainSamples(*train_samples, num_training_examples);

  // test the junction tree

/*     int jtree_nodes = 2 * (num_observed_nodes - skip);
  vector<int> tree_matrix_dims(2,jtree_nodes);
  Matrix tree_matrix(tree_matrix_dims);
 
  vector<vector<int>*> nodes_to_vars;
  for (int i = 0; i < jtree_nodes; ++i)
	  nodes_to_vars.push_back(new vector<int>());

  for (int i = 1; i < jtree_nodes / 2; ++i)
  {
	  tree_matrix.Set(i, i - 1, 1);
  }


  for (int i = jtree_nodes / 2; i < jtree_nodes; ++i)
  {
	  tree_matrix.Set(i, i - (jtree_nodes / 2), 1);
  }


  for (int i = 0; i < (jtree_nodes) / 2;  ++i)
  {
	  nodes_to_vars[i]->push_back(hidden_nodes[i + skip]);
	  if (i + 1 < (jtree_nodes) / 2)
		  nodes_to_vars[i]->push_back(hidden_nodes[i+1 + skip]);
	  if (i + 2 < (jtree_nodes) / 2)
		  nodes_to_vars[i]->push_back(hidden_nodes[i+2 + skip]);
  }

  for (int i = (jtree_nodes) / 2; i < jtree_nodes; ++i)
  {
	  nodes_to_vars[i]->push_back(hidden_nodes[i - (jtree_nodes) / 2 + skip]);
	  nodes_to_vars[i]->push_back(observed_nodes[i - (jtree_nodes) / 2 + skip]);
  } */


      int jtree_nodes = 2 * (num_observed_nodes - skip);
  vector<int> tree_matrix_dims(2,jtree_nodes);
  Matrix tree_matrix(tree_matrix_dims);
 
  vector<vector<int>*> nodes_to_vars;
  for (int i = 0; i < jtree_nodes; ++i)
	  nodes_to_vars.push_back(new vector<int>());

  for (int i = 1; i < jtree_nodes / 2; ++i)
  {
	  tree_matrix.Set(i, i - 1, 1);
  }


  for (int i = jtree_nodes / 2; i < jtree_nodes; ++i)
  {
	  tree_matrix.Set(i, i - (jtree_nodes / 2), 1);
  }


  for (int i = 0; i < (jtree_nodes) / 2;  ++i)
  {
	  nodes_to_vars[i]->push_back(hidden_nodes[i + skip]);
	  if (i + 1 < (jtree_nodes) / 2)
		  nodes_to_vars[i]->push_back(hidden_nodes[i+1 + skip]);
	  if (i + 2 < (jtree_nodes) / 2)
		  nodes_to_vars[i]->push_back(hidden_nodes[i+2 + skip]);
  }

  for (int i = (jtree_nodes) / 2; i < jtree_nodes; ++i)
  {
	  nodes_to_vars[i]->push_back(hidden_nodes[i - (jtree_nodes) / 2 + skip]);
	  nodes_to_vars[i]->push_back(observed_nodes[i - (jtree_nodes) / 2 + skip]);

	  if (i == (jtree_nodes) / 2)
	  {
		  nodes_to_vars[i]->push_back(hidden_nodes[i - (jtree_nodes) / 2 + skip + 1]);
		  for (int j = 0; j < skip; ++j)
		  {
			  nodes_to_vars[i]->push_back(observed_nodes[j]);
		  }
	  }
  } 


/*  int jtree_nodes = num_observed_nodes - skip;
  vector<int> tree_matrix_dims(2,jtree_nodes);
  Matrix tree_matrix(tree_matrix_dims);

  for (int i = 1; i < jtree_nodes; ++i)
  {
	  tree_matrix.Set(i, i - 1, 1);
  }
 

   vector<vector<int>*> nodes_to_vars;
  for (int i = 0; i < jtree_nodes; ++i)
	  nodes_to_vars.push_back(new vector<int>());

  for (int i = 0; i < jtree_nodes;  ++i)
  {
	  nodes_to_vars[i]->push_back(hidden_nodes[i + skip]);
	  nodes_to_vars[i]->push_back(observed_nodes[i + skip]);
	   if (i + 1 < jtree_nodes)
		  nodes_to_vars[i]->push_back(hidden_nodes[i+1 + skip]);
	  if (i + 2 < jtree_nodes)
		  nodes_to_vars[i]->push_back(hidden_nodes[i+2 + skip]);
  } */
  
/* vector<int> jtemplate_vec(jtree_nodes, 0);
  
  for (int i = 1; i < jtree_nodes - 2; ++i)
  {
	jtemplate_vec[i] = 1;
  }
  jtemplate_vec[jtree_nodes - 2] = 2;
  jtemplate_vec[jtree_nodes - 1] = 3; */

 vector<int> jtemplate_vec;
  VectorPlus::Seq(jtemplate_vec, 0, 1, jtree_nodes);

  Matrix EM_tree_matrix;
  vector<vector<int>*> EM_nodes_to_vars;
  bayesNet->FindCliques(EM_nodes_to_vars);
  bayesNet->FindJunctionTree(EM_tree_matrix, EM_nodes_to_vars);

 /* int EM_jtree_nodes = num_observed_nodes;
  vector<int> EM_tree_matrix_dims(2,EM_jtree_nodes);
  Matrix EM_tree_matrix(EM_tree_matrix_dims);

  for (int i = 1; i < EM_jtree_nodes; ++i)
  {
	  EM_tree_matrix.Set(i, i - 1, 1);
  }
 

   vector<vector<int>*> EM_nodes_to_vars;
  for (int i = 0; i < EM_jtree_nodes; ++i)
	  EM_nodes_to_vars.push_back(new vector<int>());

  for (int i = 0; i < EM_jtree_nodes;  ++i)
  {
	  EM_nodes_to_vars[i]->push_back(hidden_nodes[i]);
	  EM_nodes_to_vars[i]->push_back(observed_nodes[i]);
	   if (i + 1 < EM_jtree_nodes)
		  EM_nodes_to_vars[i]->push_back(hidden_nodes[i+1]);
	  if (i + 2 < EM_jtree_nodes)
		  EM_nodes_to_vars[i]->push_back(hidden_nodes[i+2]);
  } */

  vector<int> evidence_vars;
  VectorPlus::Copy(evidence_vars, observed_nodes);
  vector<int> evidence_vals(evidence_vars.size(), 0);
  cout << "yee1";
   clock_t spectral_init = clock();
  TensorJTree jtree(bayesNet, tree_matrix, nodes_to_vars, stack_count, jtemplate_vec, max_runs);
  jtree.LearnSpectralParameters();
  clock_t spectral_final = clock() - spectral_init;
  double spectral_time = ((double)spectral_final) / ((double) CLOCKS_PER_SEC);


  double EM_time = 0;
  JTree* EMtree = NULL;
  double curr_top_like = 0;
  for (int i = 0; i < num_EM_trials; ++i)
  {
		clock_t i_EM_init = clock();
		JTree* i_EMtree = new JTree(bayesNet, EM_tree_matrix, EM_nodes_to_vars, train_samples, EM_thresh, 1);
		double like = i_EMtree->LearnEMParameters(10, 0);
		clock_t i_EM_final = clock() - i_EM_init;
		EM_time += ((double)i_EM_final) / ((double) CLOCKS_PER_SEC);

		if (i == 0 || like > curr_top_like)
		{
			if (EMtree != NULL)
				delete(EMtree);

			EMtree = i_EMtree;
			curr_top_like = like;
		}
		else
		{
			delete(i_EMtree);
		}
  }

  
  

    double EM_low_time = 0;
  JTree* EMlowtree = NULL;
  double curr_top_low_like = 0;
  for (int i = 0; i < num_EM_trials; ++i)
  {
		clock_t i_EM_init = clock();
		JTree* i_EMtree = new JTree(bayesNet, EM_tree_matrix, EM_nodes_to_vars, train_samples, EM_early_thresh, 1);
		double like = i_EMtree->LearnEMParameters(10, 0);
		clock_t i_EM_final = clock() - i_EM_init;
		EM_low_time += ((double)i_EM_final) / ((double) CLOCKS_PER_SEC);

		if (i == 0 || like > curr_top_low_like)
		{
			if (EMlowtree != NULL)
				delete(EMlowtree);

			EMlowtree = i_EMtree;
			curr_top_low_like = like;
		}
		else
		{
			delete(i_EMtree);
		}
  }

  double onlineEM_high_time = 0;
  JTree* EMonlinehightree = NULL;
  double curr_top_online_high_like = 0;
  for (int i = 0; i < num_EM_trials; ++i)
  {
		clock_t i_EM_init = clock();
		JTree* i_EMtree = new JTree(bayesNet, EM_tree_matrix, EM_nodes_to_vars, train_samples, EM_thresh, 1);
		double like = i_EMtree->LearnOnlineEMParameters(online_rates[i], 0, 1000000);
		clock_t i_EM_final = clock() - i_EM_init;
		onlineEM_high_time += ((double)i_EM_final) / ((double) CLOCKS_PER_SEC);

		if (i == 0 || like > curr_top_online_high_like)
		{
			if (EMonlinehightree != NULL)
				delete(EMonlinehightree);

			EMonlinehightree = i_EMtree;
			curr_top_online_high_like = like;
		}
		else
		{
			delete(i_EMtree);
		}
  }

   double onlineEM_low_time = 0;
  JTree* EMonlinelowtree = NULL;
  double curr_top_online_low_like = 0;
  for (int i = 0; i < num_EM_trials; ++i)
  {
		clock_t i_EM_init = clock();
		JTree* i_EMtree = new JTree(bayesNet, EM_tree_matrix, EM_nodes_to_vars, train_samples, EM_early_thresh, 1);
		double like = i_EMtree->LearnOnlineEMParameters(online_rates[i], 0, 1000000);
		clock_t i_EM_final = clock() - i_EM_init;
		onlineEM_low_time += ((double)i_EM_final) / ((double) CLOCKS_PER_SEC);

		if (i == 0 || like > curr_top_online_low_like)
		{
			if (EMonlinelowtree != NULL)
				delete(EMonlinelowtree);

			EMonlinelowtree = i_EMtree;
			curr_top_online_low_like = like;
		}
		else
		{
			delete(i_EMtree);
		}
  }

	vector<TensorCPT*> marginals;
	double spectral_error = 0;
	double EM_error = 0;
	double EM_low_error = 0;
	double onlineEM_high_error = 0;
	double onlineEM_low_error = 0;

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
			spectral_val = 0;

		double EM_val = EMtree->ComputeEmpiricalMarginals(marginals, evidence_vars, evidence_vals);
		double EM_low_val = EMlowtree->ComputeEmpiricalMarginals(marginals, evidence_vars, evidence_vals);
		double EM_online_high_val = EMonlinehightree->ComputeEmpiricalMarginals(marginals, evidence_vars, evidence_vals);
		double EM_online_low_val = EMonlinelowtree->ComputeEmpiricalMarginals(marginals, evidence_vars, evidence_vals);

		//cout << spectral_val;
		//cout << "\t";
		//cout << prob;
		//cout << "\n";
		//cout << val;
		//cout << "\n";
		//cout << EM_val;
		//cout << "\t";
		//cout << EM_val2;
		//	cout << "\t";
		//	cout << EM_val3; 
		//	cout << "\n";

	spectral_error += abs(prob - spectral_val) / prob;
	EM_error += abs(prob - EM_val) / prob;
	EM_low_error += abs(prob - EM_low_val) / prob;
	onlineEM_high_error += abs(prob - EM_online_high_val) / prob;
	onlineEM_low_error += abs(prob - EM_online_low_val) / prob;

  }

  spectral_error = spectral_error /  test_samples->Dim(0);
  EM_error = EM_error / test_samples->Dim(0);
  EM_low_error = EM_low_error / test_samples->Dim(0);
  onlineEM_high_error = onlineEM_high_error / test_samples->Dim(0);
  onlineEM_low_error = onlineEM_low_error / test_samples->Dim(0);

  cout << "total error:\n";
  cout << spectral_error;
  cout << "\n";
  cout << EM_error;
  cout << "\n";
  cout << EM_low_error;
  cout << "\n";
  cout << onlineEM_high_error;
  cout << "\n";
  cout << onlineEM_low_error;
  cout << "\n";

  string outfile_name = "second-nh-hmm-all-" + to_string(num_hidden_nodes) + "-" + to_string(num_hidden_states) + "-" 
						+  to_string(num_observed_states) + "-" + to_string(num_training_examples)
						+ "-" + to_string(EM_thresh) + "-" + to_string(stack_count) + "-" + to_string(max_runs)
						+ "-" + to_string(trial) + ".txt";
//  string outfile_name = "fun.txt";
  ofstream myfile;
  myfile.open(outfile_name.c_str());
  myfile << ("spectral\t" + to_string(spectral_error) + "\t" + to_string(spectral_time) + "\t" + to_string(spectral_time) + "\n");
  myfile << ("EM\t" + to_string(EM_error) + "\t" + to_string(EM_time) + "\t" + to_string(curr_top_like) + "\n");
  myfile << ("EM-low\t" + to_string(EM_low_error) + "\t" + to_string(EM_low_time) + "\t" + to_string(curr_top_low_like) + "\n");
  myfile << ("online-EM-high\t" + to_string(onlineEM_high_error) + "\t" + to_string(onlineEM_high_time) + "\t" + to_string(curr_top_online_high_like) + "\n");
  myfile << ("online-EM-low\t" + to_string(onlineEM_low_error) + "\t" + to_string(onlineEM_low_time) + "\t" + to_string(curr_top_online_low_like) + "\n");
  myfile << ("randseed\t" + to_string(rand_seed) + "\t" + to_string(rand_seed) + to_string(rand_seed) + "\n");

  return 0;
}
