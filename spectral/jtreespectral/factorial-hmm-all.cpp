#include "TensorJTree.hpp"
#include "VectorPlus.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include "JTree.hpp"
#include "TemplateJTree.hpp"
#include <time.h>
#include <assert.h>
#include "TensorVarElim.hpp"
#include "tensor.hpp"
#include "BayesianNetwork.hpp"
#include "TemplateBayesianNetwork.hpp"
#include "StatsLibrary.hpp"
#include <cmath>
#include <stdlib.h>
#include <vector>

using namespace std;

int main (int argc, char *argv[])
{
  //srand(30);
  assert(argc == 10);
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

  int trial = atoi(argv[9]); 

    int num_EM_trials = 5; 
  double EM_early_thresh = 0.01;
  vector<double> online_rates;
  online_rates.push_back(0.6);
  online_rates.push_back(0.7);
  online_rates.push_back(0.8);
  online_rates.push_back(0.9);
  online_rates.push_back(1);

  // model parameters
/*  int num_hidden_level1 = 6;
  int num_hidden_level2 = num_hidden_level1;
  int num_obs_per_pair = 2;

  int num_hidden_states = 2;
  int num_observed_states = 2;

  int num_training_examples = 10;

  // EM iterations
  double EM_thresh = 0.001;

  // spectral observations
  int stack_count = 1;
  int max_runs = 1;

  int trial = 1; */

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
		  j_parent_list.at(observed_nodes[i])->push_back(observed_nodes[i-1]);
		  orig_parent_list.at(observed_nodes[i])->push_back(observed_nodes[i-1]);
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


  /*

 double EM_time = 0;
  JTree* EMtree = NULL;
  double curr_top_like = 0;
  for (int i = 0; i < num_EM_trials; ++i)
  {
		clock_t i_EM_init = clock();
		JTree* i_EMtree = new JTree(orig_bayesNet, EM_tree_matrix, EM_nodes_to_vars, train_samples, EM_thresh, 0);
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
		JTree* i_EMtree = new JTree(orig_bayesNet, EM_tree_matrix, EM_nodes_to_vars, train_samples, EM_early_thresh, 0);
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
		JTree* i_EMtree = new JTree(orig_bayesNet, EM_tree_matrix, EM_nodes_to_vars, train_samples, EM_thresh, 0);
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
		JTree* i_EMtree = new JTree(orig_bayesNet, EM_tree_matrix, EM_nodes_to_vars, train_samples, EM_early_thresh, 0);
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

*/
	vector<TensorCPT*> marginals;
	double spectral_error = 0;
	double EM_error = 0;
	double EM_low_error = 0;
	double onlineEM_high_error = 0;
	double onlineEM_low_error = 0;

	int win_count = 0;

	TemplateJTree oracle_tree(orig_bayesNet, EM_tree_matrix, EM_nodes_to_vars, train_samples, EM_early_thresh);
	oracle_tree.LearnEMParameters(0, 1);

string outfileval_name = "factorial-hmm-vals-" + to_string(num_hidden_nodes) + "-" + to_string(num_hidden_states) + "-" 
						+  to_string(num_observed_states) + "-" + to_string(num_training_examples)
						+ "-" + to_string(EM_thresh) + "-" + to_string(stack_count) + "-" + to_string(max_runs)
						+ "-" + to_string(trial) + ".txt";
//  string outfile_name = "fun.txt";
  ofstream myfileval;
  myfileval.open(outfileval_name.c_str());

	for (int n = 0; n < test_samples->Dim(0); ++n)
	{
		for (int i = 0; i < evidence_vals.size(); ++i)
		{
			evidence_vals.at(i) = (int)test_samples->At(n,evidence_vars[i]);
		}
  
	//	double prob = TensorVarElim::VE(*orig_bayesNet, evidence_vars, evidence_vals, elimination_order);
		double prob = oracle_tree.ComputeEmpiricalMarginals(marginals, evidence_vars, evidence_vals);
	//	cout << "\n";
	//	cout << abs(prob - o_prob) / prob;
	//	cout << "\n";
		double spectral_val = jtree.ComputeMarginalEmpiricalProbability(evidence_vars, evidence_vals);
		if (spectral_val < 0)
			spectral_val = 0;
		if (spectral_val > 1)
			spectral_val = 0;


/*

		double EM_val = EMtree->ComputeEmpiricalMarginals(marginals, evidence_vars, evidence_vals);
		double EM_low_val = EMlowtree->ComputeEmpiricalMarginals(marginals, evidence_vars, evidence_vals);
		double EM_online_high_val = EMonlinehightree->ComputeEmpiricalMarginals(marginals, evidence_vars, evidence_vals);
		double EM_online_low_val = EMonlinelowtree->ComputeEmpiricalMarginals(marginals, evidence_vars, evidence_vals);

				myfileval << (to_string(prob) + "\t" + to_string(spectral_val) + "\t" + to_string(EM_val) + "\t" + to_string(EM_low_val) + "\t" + to_string(EM_online_high_val) + "\t" + to_string(EM_online_low_val) +  "\n");

*/

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

	/*
	EM_error += abs(prob - EM_val) / prob;
	EM_low_error += abs(prob - EM_low_val) / prob;
	onlineEM_high_error += abs(prob - EM_online_high_val) / prob;
	onlineEM_low_error += abs(prob - EM_online_low_val) / prob;

	 */
  }
	myfileval.close();
  spectral_error = spectral_error /  test_samples->Dim(0);

  /*
  EM_error = EM_error / test_samples->Dim(0);
  EM_low_error = EM_low_error / test_samples->Dim(0);
  onlineEM_high_error = onlineEM_high_error / test_samples->Dim(0);
  onlineEM_low_error = onlineEM_low_error / test_samples->Dim(0);

   */


  cout << "total error:\n";
  cout << spectral_error;
  cout << "\n";

  /*
  cout << EM_error;
  cout << "\n";
  cout << EM_low_error;
  cout << "\n";
  cout << onlineEM_high_error;
  cout << "\n";
  cout << onlineEM_low_error;
  cout << "\n";
*/

  string outfile_name = "factorial-hmm-all-" + to_string(num_hidden_level1) + "-" + to_string(num_hidden_states) + "-" 
						+  to_string(num_observed_states) + "-" + to_string(num_training_examples)
						+ "-" + to_string(EM_thresh) + "-" + to_string(stack_count) + "-" + to_string(max_runs)
						+ "-" + to_string(trial) + ".txt";
//  string outfile_name = "fun.txt";
  ofstream myfile;
  myfile.open(outfile_name.c_str());
  myfile << ("spectral\t" + to_string(spectral_error) + "\t" + to_string(spectral_time) + "\t" + to_string(spectral_time) + "\n");

  /*
  myfile << ("EM\t" + to_string(EM_error) + "\t" + to_string(EM_time) + "\t" + to_string(curr_top_like) + "\n");
  myfile << ("EM-low\t" + to_string(EM_low_error) + "\t" + to_string(EM_low_time) + "\t" + to_string(curr_top_low_like) + "\n");
  myfile << ("online-EM-high\t" + to_string(onlineEM_high_error) + "\t" + to_string(onlineEM_high_time) + "\t" + to_string(curr_top_online_high_like) + "\n");
  myfile << ("online-EM-low\t" + to_string(onlineEM_low_error) + "\t" + to_string(onlineEM_low_time) + "\t" + to_string(curr_top_online_low_like) + "\n");
  */
  myfile << ("randseed\t" + to_string(rand_seed) + "\t" + to_string(rand_seed) + to_string(rand_seed) + "\n");

  return 0;
}

