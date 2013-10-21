#include "TensorJTree.hpp"
#include "VectorPlus.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include "TemplateJTree.hpp"
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

  int num_EM_trials = 1; 
  double EM_early_thresh = 0.01;
  vector<double> online_rates;
  online_rates.push_back(0.6);
  online_rates.push_back(0.7);
  online_rates.push_back(0.8);
  online_rates.push_back(0.9);
  online_rates.push_back(1);

 /* int hidden_depth = 2;
  int num_observed_nodes_per_leaf = 2;
  int num_hidden_states = 2;
  int num_observed_states = 2;
  int num_training_examples = 100;
  double EM_thresh = 0.001;
  int stack_count = 1;
  int max_runs = 1;
  int trial = 1;*/
  assert(num_observed_nodes_per_leaf = 2); 
  int rand_seed = 60 * trial;
  srand(rand_seed);

  int hidden_R_nodes_per_clique = 3;

  int num_hidden_cliques = (int)pow(2.0, (double)(hidden_depth)) - 1;
  int num_hidden_leaves = (int)(num_hidden_cliques + 1)/2;
  int num_hidden_interior = num_hidden_cliques - num_hidden_leaves;
  int num_hidden_nodes = 3 * num_hidden_cliques;
  int num_observed_nodes = 2 * num_hidden_leaves * num_observed_nodes_per_leaf;

  vector<int> hidden_nodes;
  vector<int> observed_nodes;
  VectorPlus::Seq(hidden_nodes, 0,1, num_hidden_nodes); 
  VectorPlus::Seq(observed_nodes, num_hidden_nodes, 1, num_hidden_nodes + num_observed_nodes);
  int num_nodes = num_hidden_nodes + num_observed_nodes;

  // make junction tree first
  int jtree_nodes = num_hidden_cliques + num_observed_nodes / 2;
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
		 int obs_index = (i - num_hidden_cliques) / 2;
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

		  for (int j = 0; j < 3; ++j)
		  {
			  nodes_to_vars[i]->push_back(curr_index);
			  curr_index++;
		  }
	  }
	  else
	  {
		  vector<int> parent_vec;
		  tree_matrix.RowFind(parent_vec, i);
		  int parent = parent_vec[0];
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

		  // put the two observations in one clique
		  nodes_to_vars[i]->push_back(curr_index++);
		  nodes_to_vars[i]->push_back(curr_index++);
	  }
  }


  for (int i =0; i < jtree_nodes; ++i)
  {
	  vector<int>& clique_nodes = *(nodes_to_vars)[i];
	  for (int j = 0; j < clique_nodes.size(); ++j)
	  {
		//  if (type_vector[clique_nodes[j]] == OBSERVED_FLAG)
		//	  continue;
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
  bayesNet->GenerateTestSamples(*test_samples, 1000);
  bayesNet->GenerateTrainSamples(*train_samples, num_training_examples);

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

 double EM_time = 0;
  JTree* EMtree = NULL;
  double curr_top_like = 0;
  for (int i = 0; i < num_EM_trials; ++i)
  {
		clock_t i_EM_init = clock();
		JTree* i_EMtree = new JTree(bayesNet, tree_matrix, nodes_to_vars, train_samples, EM_thresh, 1);
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
  for (int i = 0; i < 1; ++i)
  {
		clock_t i_EM_init = clock();
		JTree* i_EMtree = new JTree(bayesNet, tree_matrix, nodes_to_vars, train_samples, EM_early_thresh, 1);
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
		JTree* i_EMtree = new JTree(bayesNet, tree_matrix, nodes_to_vars, train_samples, EM_thresh, 1);
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
  for (int i = 0; i < 1; ++i)
  {
		clock_t i_EM_init = clock();
		JTree* i_EMtree = new JTree(bayesNet, tree_matrix, nodes_to_vars, train_samples, EM_early_thresh, 1);
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

	TemplateJTree oracle_tree(bayesNet, tree_matrix, nodes_to_vars, train_samples, EM_early_thresh);
	oracle_tree.LearnEMParameters(0, 1);

	string outfileval_name = "overlapping-jtree-vals-" + to_string(hidden_depth) + "-" + to_string(num_observed_nodes_per_leaf ) 
						+ "-" + to_string(num_hidden_states) + "-" 
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
  
		double prob = oracle_tree.ComputeEmpiricalMarginals(marginals, evidence_vars, evidence_vals);
		double spectral_val = jtree.ComputeMarginalEmpiricalProbability(evidence_vars, evidence_vals);
		if (spectral_val < 0)
			spectral_val = 0;
		if (spectral_val > 1)
			spectral_val = 0;

		double EM_val = EMtree->ComputeEmpiricalMarginals(marginals, evidence_vars, evidence_vals);
		double EM_low_val = EMlowtree->ComputeEmpiricalMarginals(marginals, evidence_vars, evidence_vals);
		double EM_online_high_val = EMonlinehightree->ComputeEmpiricalMarginals(marginals, evidence_vars, evidence_vals);
		double EM_online_low_val = EMonlinelowtree->ComputeEmpiricalMarginals(marginals, evidence_vars, evidence_vals);

						myfileval << (to_string(prob) + "\t" + to_string(spectral_val) + "\t" + to_string(EM_val) + "\t" + to_string(EM_low_val) + "\t" + to_string(EM_online_high_val) + "\t" + to_string(EM_online_low_val) +  "\n");
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
	myfileval.close();
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

  string outfile_name = "overlapping-jtree-" + to_string(hidden_depth) + "-" + to_string(num_observed_nodes_per_leaf ) 
						+ "-" + to_string(num_hidden_states) + "-" 
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