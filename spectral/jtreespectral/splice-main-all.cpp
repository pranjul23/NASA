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
#include <sstream>

using namespace std;

int main (int argc, char *argv[])
{

  assert(argc == 7);
  int num_hidden_states = atoi(argv[1]);
  double EM_thresh = atof(argv[2]);
  int stack_count = atoi(argv[3]);
  int max_runs = atoi(argv[4]);
  int trial = atoi(argv[5]);
  int num_training_examples = atoi(argv[6]);

 /* int num_hidden_states = 2;
  double EM_thresh = 0.0001;
  int stack_count = 1;
  int max_runs = 1;
  int trial = 1;
  int num_training_examples = 10; */

  int num_EM_trials = 5; 
  double EM_early_thresh = 0.01;
  vector<double> online_rates;
  online_rates.push_back(0.6);
  online_rates.push_back(0.7);
  online_rates.push_back(0.8);
  online_rates.push_back(0.9);
  online_rates.push_back(1);

  int EM_FLAG = 1;
  int num_hidden_nodes = 60;
  int num_observed_nodes = 60;
  int num_observed_states = 4;
  // node names
  vector<int> hidden_nodes;
  vector<int> observed_nodes;
  VectorPlus::Seq(hidden_nodes, 0,1, num_hidden_nodes); 
  VectorPlus::Seq(observed_nodes, num_hidden_nodes, 1, num_hidden_nodes + num_observed_nodes);
  
  
  int num_nodes = num_hidden_nodes + num_observed_nodes;
  // spectral observations

  
//  int num_training_examples = 2675;
  int num_testing_examples = 500;
  int num_columns = 61;
  int num_obs_vars = 60;

  string train_file_name = "splice_train" + to_string(trial) + ".txt";
  string test_file_name = "splice_test" + to_string(trial) + ".txt";

  vector<int> test_dims = VectorPlus::CreatePair(num_testing_examples, num_nodes);
  Matrix* test_samples = new Matrix(test_dims);

  vector<int> test_labels(num_testing_examples, 0);

  ifstream test_file;
  test_file.open(test_file_name.c_str());
  string test_line;
  int curr_test_index = 0;
  while (getline(test_file,test_line))
  {
	  stringstream ss(test_line);
	  string delimited_line;
	  getline(ss, delimited_line, '\t');
	  int class_var = atoi(delimited_line.c_str());
	  test_labels[curr_test_index] = class_var;

	  for (int i = 0; i < num_obs_vars; ++i)
	  {
		getline(ss, delimited_line, '\t');
		int val = atoi(delimited_line.c_str());
		test_samples->Set(curr_test_index, num_hidden_nodes + i, val);
	  }
	  curr_test_index++;
  }
  test_file.close();



  ifstream train_file;
  train_file.open(train_file_name.c_str());

  string line;


  int zero_count = 0;
  int one_count = 0;
  int two_count = 0;

  int curr_line = 0;
  while (getline(train_file,line) && curr_line < num_training_examples)
  {
	  stringstream ss(line);
	  string delimited_line;
	  getline(ss, delimited_line, '\t');
	  int class_var = atoi(delimited_line.c_str());
	  if (class_var == 0)
		  zero_count++;
	  else if (class_var == 1)
		  one_count++;
	  else if (class_var == 2)
		  two_count++;
	  else
		  assert(0);

	  curr_line++;

  }
  train_file.close();
  if (zero_count == 0)
	  assert(0);
  vector<int> zero_dims = VectorPlus::CreatePair(zero_count, num_nodes);
  vector<int> one_dims  = VectorPlus::CreatePair(one_count, num_nodes);
  vector<int> two_dims = VectorPlus::CreatePair(two_count, num_nodes);

  Matrix* train_samples0 = new Matrix(zero_dims);
  Matrix* train_samples1 = new Matrix(one_dims);
  Matrix* train_samples2 = new Matrix(two_dims);

  train_file.open(train_file_name.c_str());

  int curr_row0 = 0;
  int curr_row1 = 0;
  int curr_row2 = 0;
  curr_line = 0;
  while (getline(train_file,line) && curr_line < num_training_examples)
  {
	  stringstream ss(line);
	  string delimited_line;
	  getline(ss, delimited_line, '\t');
	  int class_var = atoi(delimited_line.c_str());
	  
	  for (int i = 0; i < num_obs_vars; ++i)
	  {
		getline(ss, delimited_line, '\t');
		int val = atoi(delimited_line.c_str());
		if (class_var == 0)
		{
			train_samples0->Set(curr_row0, num_hidden_nodes + i, val);
		}
		if (class_var == 1)
		{
			train_samples1->Set(curr_row1, num_hidden_nodes + i, val);
		}
		if (class_var == 2)
		{
			train_samples2->Set(curr_row2, num_hidden_nodes + i, val);
		}
	  }
	  if (class_var == 0)
		  curr_row0++;
	  else if (class_var == 1)
		  curr_row1++;
	  else if (class_var == 2)
		  curr_row2++;
	  else
		  assert(0);

	  curr_line++;
  }
  train_file.close();



  int rand_seed = 10 * trial;
  srand(rand_seed);

   int skip = 1;


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


  TemplateBayesianNetwork* bayesNet0 = new TemplateBayesianNetwork(parent_list, type_vector, dims, template_vec);
  TemplateBayesianNetwork* bayesNet1 = new TemplateBayesianNetwork(parent_list, type_vector, dims, template_vec);
  TemplateBayesianNetwork* bayesNet2 = new TemplateBayesianNetwork(parent_list, type_vector, dims, template_vec);

  bayesNet0->SetTrainingSamples(*train_samples0);
  bayesNet1->SetTrainingSamples(*train_samples1);
  bayesNet2->SetTrainingSamples(*train_samples2);

  // test the junction tree
/* int jtree_nodes = num_observed_nodes - skip;
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
  } 

   vector<int> jtemplate_vec;
  VectorPlus::Seq(jtemplate_vec, 0, 1, jtree_nodes); */

 /* int jtree_nodes = 2 * (num_observed_nodes - skip);
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

  for (int i = 0; i < jtree_nodes / 2;  ++i)
  {
	  nodes_to_vars[i]->push_back(hidden_nodes[i + skip]);
	  if (i + 1 < jtree_nodes / 2)
		  nodes_to_vars[i]->push_back(hidden_nodes[i+1 + skip]);
	  if (i + 2 < jtree_nodes / 2)
		  nodes_to_vars[i]->push_back(hidden_nodes[i+2 + skip]);
  }

  for (int i = jtree_nodes / 2; i < jtree_nodes; ++i)
  {
	  nodes_to_vars[i]->push_back(hidden_nodes[i - jtree_nodes / 2 + skip]);
	  nodes_to_vars[i]->push_back(observed_nodes[i - jtree_nodes / 2 + skip]);
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


  vector<int> jtemplate_vec;
 VectorPlus::Seq(jtemplate_vec, 0, 1, jtree_nodes);

 /* vector<int> jtemplate_vec(jtree_nodes, 0);
  for (int i = 1; i < jtree_nodes/2 - 2; ++i)
  {
	jtemplate_vec[i] = i;
  }
  jtemplate_vec[jtree_nodes/2 - 2] = jtree_nodes/2 - 2;
  jtemplate_vec[jtree_nodes/2 - 1] = jtree_nodes/2 - 1;

  for (int i = jtree_nodes/2; i < jtree_nodes; ++i)
  {
	  int index = (i - jtree_nodes / 2) / 1;
	  jtemplate_vec[i] = jtree_nodes / 2 + index;
  } */

  Matrix EM_tree_matrix;
  vector<vector<int>*> EM_nodes_to_vars;
  bayesNet0->FindCliques(EM_nodes_to_vars);
  bayesNet0->FindJunctionTree(EM_tree_matrix, EM_nodes_to_vars);

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

  vector<int> EM_evidence_vars;
  vector<int> spectral_evidence_vars;
 // vector<int> omitted_obs;
 // for (int i = 0; i < skip; ++i)
 // {
//	omitted_obs.push_back(observed_nodes[i]);
 // }

  VectorPlus::Copy(spectral_evidence_vars, observed_nodes);
  VectorPlus::Copy(EM_evidence_vars, observed_nodes);



  /* learn spectral parameters *****************************************************************/
  clock_t spectral_init0 = clock();
  TensorJTree jtree0(bayesNet0, tree_matrix, nodes_to_vars, stack_count, jtemplate_vec, max_runs);
  jtree0.LearnSpectralParameters();
  clock_t spectral_final0 = clock() - spectral_init0;
  double spectral_time0 = ((double)spectral_final0) / ((double) CLOCKS_PER_SEC);

  clock_t spectral_init1 = clock();
  TensorJTree jtree1(bayesNet1, tree_matrix, nodes_to_vars, stack_count, jtemplate_vec, max_runs);
  jtree1.LearnSpectralParameters();
  clock_t spectral_final1 = clock() - spectral_init1;
  double spectral_time1 = ((double)spectral_final1) / ((double) CLOCKS_PER_SEC);

  clock_t spectral_init2 = clock();
  TensorJTree jtree2(bayesNet2, tree_matrix, nodes_to_vars, stack_count, jtemplate_vec, max_runs);
  jtree2.LearnSpectralParameters();
  clock_t spectral_final2 = clock() - spectral_init2;
  double spectral_time2 = ((double)spectral_final2) / ((double) CLOCKS_PER_SEC);
  /*******************************************************************************************/

  /* learn EM zero parameters *****************************************************************/

   double EM_time0 = 0;
  JTree* EMtree0 = NULL;
  double curr_top_like0 = 0;
  for (int i = 0; i < num_EM_trials; ++i)
  {
		clock_t i_EM_init = clock();
		JTree* i_EMtree = new JTree(bayesNet0, EM_tree_matrix, EM_nodes_to_vars, train_samples0, EM_thresh, 0);
		double like = i_EMtree->LearnEMParameters(10, 0);
		clock_t i_EM_final = clock() - i_EM_init;
		EM_time0 += ((double)i_EM_final) / ((double) CLOCKS_PER_SEC);

		if (i == 0 || like > curr_top_like0)
		{
			if (EMtree0 != NULL)
				delete(EMtree0);

			EMtree0 = i_EMtree;
			curr_top_like0 = like;
		}
		else
		{
			delete(i_EMtree);
		}
  }

  
  

/*    double EM_low_time0 = 0;
  JTree* EMlowtree0 = NULL;
  double curr_top_low_like0 = 0;
  for (int i = 0; i < num_EM_trials; ++i)
  {
		clock_t i_EM_init = clock();
		JTree* i_EMtree = new JTree(bayesNet0, EM_tree_matrix, EM_nodes_to_vars, train_samples0, EM_early_thresh, 1);
		double like = i_EMtree->LearnEMParameters(10, 0);
		clock_t i_EM_final = clock() - i_EM_init;
		EM_low_time0 += ((double)i_EM_final) / ((double) CLOCKS_PER_SEC);

		if (i == 0 || like > curr_top_low_like0)
		{
			if (EMlowtree0 != NULL)
				delete(EMlowtree0);

			EMlowtree0 = i_EMtree;
			curr_top_low_like0 = like;
		}
		else
		{
			delete(i_EMtree);
		}
  } */

  double onlineEM_high_time0 = 0;
  JTree* EMonlinehightree0 = NULL;
  double curr_top_online_high_like0 = 0;
  for (int i = 0; i < num_EM_trials; ++i)
  {
		clock_t i_EM_init = clock();
		JTree* i_EMtree = new JTree(bayesNet0, EM_tree_matrix, EM_nodes_to_vars, train_samples0, EM_thresh, 1);
		double like = i_EMtree->LearnOnlineEMParameters(online_rates[i], 0, 1000000);
		clock_t i_EM_final = clock() - i_EM_init;
		onlineEM_high_time0 += ((double)i_EM_final) / ((double) CLOCKS_PER_SEC);

		if (i == 0 || like > curr_top_online_high_like0)
		{
			if (EMonlinehightree0 != NULL)
				delete(EMonlinehightree0);

			EMonlinehightree0 = i_EMtree;
			curr_top_online_high_like0 = like;
		}
		else
		{
			delete(i_EMtree);
		}
  }

 /*  double onlineEM_low_time0 = 0;
  JTree* EMonlinelowtree0 = NULL;
  double curr_top_online_low_like0 = 0;
  for (int i = 0; i < num_EM_trials; ++i)
  {
		clock_t i_EM_init = clock();
		JTree* i_EMtree = new JTree(bayesNet0, EM_tree_matrix, EM_nodes_to_vars, train_samples0, EM_early_thresh, 1);
		double like = i_EMtree->LearnOnlineEMParameters(online_rates[i], 0, 1000000);
		clock_t i_EM_final = clock() - i_EM_init;
		onlineEM_low_time0 += ((double)i_EM_final) / ((double) CLOCKS_PER_SEC);

		if (i == 0 || like > curr_top_online_low_like0)
		{
			if (EMonlinelowtree0 != NULL)
				delete(EMonlinelowtree0);

			EMonlinelowtree0 = i_EMtree;
			curr_top_online_low_like0 = like;
		}
		else
		{
			delete(i_EMtree);
		}
  } */

  

    /* learn EM One parameters *****************************************************************/
   double EM_time1 = 0;
  JTree* EMtree1 = NULL;
  double curr_top_like1 = 0;
  for (int i = 0; i < num_EM_trials; ++i)
  {
		clock_t i_EM_init = clock();
		JTree* i_EMtree = new JTree(bayesNet1, EM_tree_matrix, EM_nodes_to_vars, train_samples1, EM_thresh, 0);
		double like = i_EMtree->LearnEMParameters(10, 0);
		clock_t i_EM_final = clock() - i_EM_init;
		EM_time1 += ((double)i_EM_final) / ((double) CLOCKS_PER_SEC);

		if (i == 0 || like > curr_top_like0)
		{
			if (EMtree1 != NULL)
				delete(EMtree1);

			EMtree1 = i_EMtree;
			curr_top_like1 = like;
		}
		else
		{
			delete(i_EMtree);
		}
  }

  
  

 /*   double EM_low_time1 = 0;
  JTree* EMlowtree1 = NULL;
  double curr_top_low_like1 = 0;
  for (int i = 0; i < num_EM_trials; ++i)
  {
		clock_t i_EM_init = clock();
		JTree* i_EMtree = new JTree(bayesNet1, EM_tree_matrix, EM_nodes_to_vars, train_samples0, EM_early_thresh, 1);
		double like = i_EMtree->LearnEMParameters(10, 0);
		clock_t i_EM_final = clock() - i_EM_init;
		EM_low_time1 += ((double)i_EM_final) / ((double) CLOCKS_PER_SEC);

		if (i == 0 || like > curr_top_low_like1)
		{
			if (EMlowtree1 != NULL)
				delete(EMlowtree1);

			EMlowtree1 = i_EMtree;
			curr_top_low_like1 = like;
		}
		else
		{
			delete(i_EMtree);
		}
  }*/

  double onlineEM_high_time1 = 0;
  JTree* EMonlinehightree1 = NULL;
  double curr_top_online_high_like1 = 0;
  for (int i = 0; i < num_EM_trials; ++i)
  {
		clock_t i_EM_init = clock();
		JTree* i_EMtree = new JTree(bayesNet1, EM_tree_matrix, EM_nodes_to_vars, train_samples1, EM_thresh, 1);
		double like = i_EMtree->LearnOnlineEMParameters(online_rates[i], 0, 1000000);
		clock_t i_EM_final = clock() - i_EM_init;
		onlineEM_high_time1 += ((double)i_EM_final) / ((double) CLOCKS_PER_SEC);

		if (i == 0 || like > curr_top_online_high_like1)
		{
			if (EMonlinehightree1 != NULL)
				delete(EMonlinehightree1);

			EMonlinehightree1 = i_EMtree;
			curr_top_online_high_like1 = like;
		}
		else
		{
			delete(i_EMtree);
		}
  }

/*   double onlineEM_low_time1 = 0;
  JTree* EMonlinelowtree1 = NULL;
  double curr_top_online_low_like1 = 0;
  for (int i = 0; i < num_EM_trials; ++i)
  {
		clock_t i_EM_init = clock();
		JTree* i_EMtree = new JTree(bayesNet1, EM_tree_matrix, EM_nodes_to_vars, train_samples0, EM_early_thresh, 1);
		double like = i_EMtree->LearnOnlineEMParameters(online_rates[i], 0, 1000000);
		clock_t i_EM_final = clock() - i_EM_init;
		onlineEM_low_time1 += ((double)i_EM_final) / ((double) CLOCKS_PER_SEC);

		if (i == 0 || like > curr_top_online_low_like1)
		{
			if (EMonlinelowtree1 != NULL)
				delete(EMonlinelowtree1);

			EMonlinelowtree1 = i_EMtree;
			curr_top_online_low_like1 = like;
		}
		else
		{
			delete(i_EMtree);
		}
  }*/

  


  /* learn EM two parameters *****************************************************************/
  double EM_time2 = 0;
  JTree* EMtree2 = NULL;
  double curr_top_like2 = 0;
  for (int i = 0; i < num_EM_trials; ++i)
  {
		clock_t i_EM_init = clock();
		JTree* i_EMtree = new JTree(bayesNet2, EM_tree_matrix, EM_nodes_to_vars, train_samples2, EM_thresh, 0);
		double like = i_EMtree->LearnEMParameters(10, 0);
		clock_t i_EM_final = clock() - i_EM_init;
		EM_time2 += ((double)i_EM_final) / ((double) CLOCKS_PER_SEC);

		if (i == 0 || like > curr_top_like0)
		{
			if (EMtree2 != NULL)
				delete(EMtree2);

			EMtree2 = i_EMtree;
			curr_top_like2 = like;
		}
		else
		{
			delete(i_EMtree);
		}
  }

  
  

/*    double EM_low_time2 = 0;
  JTree* EMlowtree2 = NULL;
  double curr_top_low_like2 = 0;
  for (int i = 0; i < num_EM_trials; ++i)
  {
		clock_t i_EM_init = clock();
		JTree* i_EMtree = new JTree(bayesNet2, EM_tree_matrix, EM_nodes_to_vars, train_samples2, EM_early_thresh, 1);
		double like = i_EMtree->LearnEMParameters(10, 0);
		clock_t i_EM_final = clock() - i_EM_init;
		EM_low_time2 += ((double)i_EM_final) / ((double) CLOCKS_PER_SEC);

		if (i == 0 || like > curr_top_low_like2)
		{
			if (EMlowtree2 != NULL)
				delete(EMlowtree2);

			EMlowtree2 = i_EMtree;
			curr_top_low_like2 = like;
		}
		else
		{
			delete(i_EMtree);
		}
  } */

  double onlineEM_high_time2 = 0;
  JTree* EMonlinehightree2 = NULL;
  double curr_top_online_high_like2 = 0;
  for (int i = 0; i < num_EM_trials; ++i)
  {
		clock_t i_EM_init = clock();
		JTree* i_EMtree = new JTree(bayesNet2, EM_tree_matrix, EM_nodes_to_vars, train_samples2, EM_thresh, 1);
		double like = i_EMtree->LearnOnlineEMParameters(online_rates[i], 0, 1000000);
		clock_t i_EM_final = clock() - i_EM_init;
		onlineEM_high_time2 += ((double)i_EM_final) / ((double) CLOCKS_PER_SEC);

		if (i == 0 || like > curr_top_online_high_like2)
		{
			if (EMonlinehightree2 != NULL)
				delete(EMonlinehightree2);

			EMonlinehightree2 = i_EMtree;
			curr_top_online_high_like2 = like;
		}
		else
		{
			delete(i_EMtree);
		}
  }

/*   double onlineEM_low_time2 = 0;
  JTree* EMonlinelowtree2 = NULL;
  double curr_top_online_low_like2 = 0;
  for (int i = 0; i < num_EM_trials; ++i)
  {
		clock_t i_EM_init = clock();
		JTree* i_EMtree = new JTree(bayesNet2, EM_tree_matrix, EM_nodes_to_vars, train_samples2, EM_early_thresh, 1);
		double like = i_EMtree->LearnOnlineEMParameters(online_rates[i], 0, 1000000);
		clock_t i_EM_final = clock() - i_EM_init;
		onlineEM_low_time2 += ((double)i_EM_final) / ((double) CLOCKS_PER_SEC);

		if (i == 0 || like > curr_top_online_low_like2)
		{
			if (EMonlinelowtree2 != NULL)
				delete(EMonlinelowtree2);

			EMonlinelowtree2 = i_EMtree;
			curr_top_online_low_like2 = like;
		}
		else
		{
			delete(i_EMtree);
		}
  } */


  

 

vector<TensorCPT*> marginals;
double spectral_error = 0;
double EM_error = 0;
double online_EM_error = 0;
int win_count = 0;

for (int n = 0; n < test_samples->Dim(0); ++n)
{
	cout << n;
	cout << "\n";
	vector<int> EM_evidence_vals;
	vector<int> spectral_evidence_vals;
	for (int i = 0; i < EM_evidence_vars.size(); ++i)
	{
		EM_evidence_vals.push_back((int)test_samples->At(n,EM_evidence_vars[i]));
	}

	for (int i = 0; i < spectral_evidence_vars.size(); ++i)
	{
		spectral_evidence_vals.push_back((int)test_samples->At(n,spectral_evidence_vars[i]));
	}
  
//	double prob = TensorVarElim::VE(*bayesNet, evidence_vars, evidence_vals, bayesNet->ReverseTopologicalOrder());
	double spectral_val0 = jtree0.ComputeMarginalEmpiricalProbability(spectral_evidence_vars, spectral_evidence_vals);
	double spectral_val1 = jtree1.ComputeMarginalEmpiricalProbability(spectral_evidence_vars, spectral_evidence_vals);
	double spectral_val2 = jtree2.ComputeMarginalEmpiricalProbability(spectral_evidence_vars, spectral_evidence_vals);

	/*cout << spectral_val0;
	cout << "\t";
    cout << spectral_val1;
	cout << "\t";
	cout << spectral_val2;
	cout << "\t";
	cout << test_labels[n];
	cout << "\n"; */
	int spectral_map = -1;

	if (spectral_val0 > spectral_val1 && spectral_val0 > spectral_val2)
	{
		spectral_map = 0;
	}
	else if (spectral_val1 > spectral_val0 && spectral_val1 > spectral_val2)
	{
		spectral_map = 1;
	}
	else
	{
		assert(spectral_val2 >= spectral_val1 && spectral_val2 >= spectral_val0);
		spectral_map = 2;
	}

	double EM_val0 = EMtree0->ComputeEmpiricalMarginals(marginals, EM_evidence_vars, EM_evidence_vals);
	double EM_val1 = EMtree1->ComputeEmpiricalMarginals(marginals, EM_evidence_vars, EM_evidence_vals);
	double EM_val2 = EMtree2->ComputeEmpiricalMarginals(marginals, EM_evidence_vars, EM_evidence_vals);

	double EM_online_val0 = EMonlinehightree0->ComputeEmpiricalMarginals(marginals, EM_evidence_vars, EM_evidence_vals);
	double EM_online_val1 = EMonlinehightree1->ComputeEmpiricalMarginals(marginals, EM_evidence_vars, EM_evidence_vals);
	double EM_online_val2 = EMonlinehightree2->ComputeEmpiricalMarginals(marginals, EM_evidence_vars, EM_evidence_vals);

	int EM_map = -1;
	int online_EM_map = -1;

	/********************** compute likelihood of 0 ************************************************/


	 if (EM_val0 > EM_val1 && EM_val0 > EM_val2)
	{
		EM_map = 0;
	}
	else if (EM_val1 > EM_val0 && EM_val1 > EM_val2)
	{
		EM_map = 1;
	}
	else
	{
		assert(EM_val2 >= EM_val1 && EM_val2 >= EM_val0);
		EM_map = 2;
	}

	if (EM_online_val0 > EM_online_val1 && EM_online_val0 > EM_online_val2)
	{
		online_EM_map = 0;
	}
	else if (EM_online_val1 > EM_online_val0 && EM_online_val1 > EM_online_val2)
	{
		online_EM_map = 1;
	}
	else
	{
		assert(EM_online_val2 >= EM_online_val1 && EM_online_val2 >= EM_online_val0);
		online_EM_map = 2;
	} 



/*	cout << EM_val0;
	cout << "\t";
    cout << EM_val1;
	cout << "\t";
	cout << EM_val2;
	cout << "\t";
	cout << test_labels[n];
	cout << "\n"; */

	if (test_labels[n] != spectral_map)
	{
		spectral_error++;
	}
	if (test_labels[n] != EM_map)
	{
		EM_error++;
	}
	if (test_labels[n] != online_EM_map)
	{
		online_EM_error++;
	}
  }
  cout << "yeeee\n";
  spectral_error = spectral_error / test_samples->Dim(0);
  EM_error = EM_error / test_samples->Dim(0);
  online_EM_error = online_EM_error / test_samples->Dim(0);
  cout << "yeeee2\n";
  double spectral_time = spectral_time0 + spectral_time1 + spectral_time2;
  cout << "yeeee3\n";
  double EM_time =  EM_time0 + EM_time1+ EM_time2;
  double onlineEM_high_time = onlineEM_high_time0 + onlineEM_high_time1+ onlineEM_high_time2;
  cout << "yeeee4\n";
  cout << "total error:\n";
  cout << spectral_error;
  cout << "\n";
  cout << EM_error;
  cout << "\n";
  cout << online_EM_error;
  cout << "\n";

  string outfile_name = "splice2-" + to_string(num_hidden_states) + "-"
						+ to_string(num_training_examples)
						+ "-" + to_string(EM_thresh) 
						+ "-" + to_string(stack_count)
						+ "-" + to_string(max_runs)
						+ "-" + to_string(trial) + ".txt";
//  string outfile_name = "fun.txt";
  ofstream myfile;
  myfile.open(outfile_name.c_str());
  myfile << ("spectral\t" + to_string(spectral_error) + "\t" + to_string(spectral_time) + "\n");
  myfile << ("EM\t" + to_string(EM_error) + "\t" + to_string(EM_time) + "\n");
  myfile << ("online EM\t" + to_string(online_EM_error) + "\t" + to_string(onlineEM_high_time) + "\n");
  myfile << ("randseed\t" + to_string(rand_seed) + "\t" + to_string(rand_seed) + "\n");

  return 0;
}
