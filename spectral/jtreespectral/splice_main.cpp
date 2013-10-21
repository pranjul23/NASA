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

/*  assert(argc == 7);
  int num_hidden_states = atoi(argv[1]);
  double EM_thresh = atof(argv[2]);
  int stack_count = atoi(argv[3]);
  int max_runs = atoi(argv[4]);
  int trial = atoi(argv[5]);
  int num_training_examples = atoi(argv[6]); */

  int num_hidden_states = 2;
  double EM_thresh = 0.5;
  int stack_count = 1;
  int max_runs = 1;
  int trial = 1;
  int num_training_examples = 10;

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

   int skip = 2;


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
  }


  vector<int> jtemplate_vec;
 VectorPlus::Seq(jtemplate_vec, 0, 1, jtree_nodes);

/*  vector<int> jtemplate_vec(jtree_nodes, 0);
  for (int i = 1; i < jtree_nodes/2 - 2; ++i)
  {
	jtemplate_vec[i] = 1;
  }
  jtemplate_vec[jtree_nodes/2 - 2] = 2;
  jtemplate_vec[jtree_nodes/2 - 1] = 3;

  for (int i = jtree_nodes/2; i < jtree_nodes; ++i)
  {
	  jtemplate_vec[i] = 4;
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
  vector<int> omitted_obs;
  for (int i = 0; i < skip; ++i)
  {
	omitted_obs.push_back(observed_nodes[i]);
  }

  VectorPlus::SetDiff(spectral_evidence_vars, observed_nodes, omitted_obs);
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


  clock_t EMZero1_init = clock();
  JTree EMtreeZero1(bayesNet0, EM_tree_matrix, EM_nodes_to_vars, train_samples0, EM_thresh);
  double likeZero1 = EMtreeZero1.LearnEMParameters(EM_FLAG, 0);
  clock_t EMZero1_final = clock() - EMZero1_init;
  double EMZero1_time = ((double)EMZero1_final) / ((double) CLOCKS_PER_SEC);

  clock_t EMZero2_init = clock();
  JTree EMtreeZero2(bayesNet0, EM_tree_matrix, EM_nodes_to_vars, train_samples0, EM_thresh);
  double likeZero2 = EMtreeZero2.LearnEMParameters(EM_FLAG, 0);
  clock_t EMZero2_final = clock() - EMZero2_init;
  double EMZero2_time = ((double)EMZero2_final) / ((double) CLOCKS_PER_SEC);

  clock_t EMZero3_init = clock();
  JTree EMtreeZero3(bayesNet0, EM_tree_matrix, EM_nodes_to_vars, train_samples0, EM_thresh);
   double likeZero3 = EMtreeZero3.LearnEMParameters(EM_FLAG, 0);
  clock_t EMZero3_final = clock() - EMZero3_init;
  double EMZero3_time = ((double)EMZero3_final) / ((double) CLOCKS_PER_SEC);

  clock_t EMZero4_init = clock();
  JTree EMtreeZero4(bayesNet0, EM_tree_matrix, EM_nodes_to_vars, train_samples0, EM_thresh);
   double likeZero4 = EMtreeZero4.LearnEMParameters(EM_FLAG, 0);
  clock_t EMZero4_final = clock() - EMZero4_init;
  double EMZero4_time = ((double)EMZero4_final) / ((double) CLOCKS_PER_SEC);

  clock_t EMZero5_init = clock();
  JTree EMtreeZero5(bayesNet0, EM_tree_matrix, EM_nodes_to_vars, train_samples0, EM_thresh);
  double likeZero5 = EMtreeZero5.LearnEMParameters(EM_FLAG, 0);
  clock_t EMZero5_final = clock() - EMZero5_init;
  double EMZero5_time = ((double)EMZero5_final) / ((double) CLOCKS_PER_SEC);

  /*****************************************************************************************/

    /* learn EM One parameters *****************************************************************/
  clock_t EMOne1_init = clock();
  JTree EMtreeOne1(bayesNet1, EM_tree_matrix, EM_nodes_to_vars, train_samples1, EM_thresh);
  double likeOne1 =  EMtreeOne1.LearnEMParameters(EM_FLAG, 0);
  clock_t EMOne1_final = clock() - EMOne1_init;
  double EMOne1_time = ((double)EMOne1_final) / ((double) CLOCKS_PER_SEC);

  clock_t EMOne2_init = clock();
  JTree EMtreeOne2(bayesNet1, EM_tree_matrix, EM_nodes_to_vars, train_samples1, EM_thresh);
  double likeOne2 = EMtreeOne2.LearnEMParameters(EM_FLAG, 0);
  clock_t EMOne2_final = clock() - EMOne2_init;
  double EMOne2_time = ((double)EMOne2_final) / ((double) CLOCKS_PER_SEC);

  clock_t EMOne3_init = clock();
  JTree EMtreeOne3(bayesNet1, EM_tree_matrix, EM_nodes_to_vars, train_samples1, EM_thresh);
  double likeOne3 = EMtreeOne3.LearnEMParameters(EM_FLAG, 0);
  clock_t EMOne3_final = clock() - EMOne3_init;
  double EMOne3_time = ((double)EMOne3_final) / ((double) CLOCKS_PER_SEC);

  clock_t EMOne4_init = clock();
  JTree EMtreeOne4(bayesNet1, EM_tree_matrix, EM_nodes_to_vars, train_samples1, EM_thresh);
  double likeOne4 = EMtreeOne4.LearnEMParameters(EM_FLAG, 0);
  clock_t EMOne4_final = clock() - EMOne4_init;
  double EMOne4_time = ((double)EMOne4_final) / ((double) CLOCKS_PER_SEC);

  clock_t EMOne5_init = clock();
  JTree EMtreeOne5(bayesNet1, EM_tree_matrix, EM_nodes_to_vars, train_samples1, EM_thresh);
  double likeOne5 = EMtreeOne5.LearnEMParameters(EM_FLAG, 0);
  clock_t EMOne5_final = clock() - EMOne5_init;
  double EMOne5_time = ((double)EMOne5_final) / ((double) CLOCKS_PER_SEC);

  /*****************************************************************************************/


  /* learn EM two parameters *****************************************************************/
  clock_t EMTwo1_init = clock();
  JTree EMtreeTwo1(bayesNet2, EM_tree_matrix, EM_nodes_to_vars, train_samples2, EM_thresh);
  double likeTwo1 = EMtreeTwo1.LearnEMParameters(EM_FLAG, 0);
  clock_t EMTwo1_final = clock() - EMTwo1_init;
  double EMTwo1_time = ((double)EMTwo1_final) / ((double) CLOCKS_PER_SEC);

  clock_t EMTwo2_init = clock();
  JTree EMtreeTwo2(bayesNet2, EM_tree_matrix, EM_nodes_to_vars, train_samples2, EM_thresh);
  double likeTwo2 = EMtreeTwo2.LearnEMParameters(EM_FLAG, 0);
  clock_t EMTwo2_final = clock() - EMTwo2_init;
  double EMTwo2_time = ((double)EMTwo2_final) / ((double) CLOCKS_PER_SEC);

  clock_t EMTwo3_init = clock();
  JTree EMtreeTwo3(bayesNet2, EM_tree_matrix, EM_nodes_to_vars, train_samples2, EM_thresh);
  double likeTwo3 = EMtreeTwo3.LearnEMParameters(EM_FLAG, 0);
  clock_t EMTwo3_final = clock() - EMTwo3_init;
  double EMTwo3_time = ((double)EMTwo3_final) / ((double) CLOCKS_PER_SEC);

  clock_t EMTwo4_init = clock();
  JTree EMtreeTwo4(bayesNet2, EM_tree_matrix, EM_nodes_to_vars, train_samples2, EM_thresh);
  double likeTwo4 = EMtreeTwo4.LearnEMParameters(EM_FLAG, 0);
  clock_t EMTwo4_final = clock() - EMTwo4_init;
  double EMTwo4_time = ((double)EMTwo4_final) / ((double) CLOCKS_PER_SEC);

  clock_t EMTwo5_init = clock();
  JTree EMtreeTwo5(bayesNet2, EM_tree_matrix, EM_nodes_to_vars, train_samples2, EM_thresh);
  double likeTwo5 = EMtreeTwo5.LearnEMParameters(EM_FLAG, 0);
  clock_t EMTwo5_final = clock() - EMTwo5_init;
  double EMTwo5_time = ((double)EMTwo5_final) / ((double) CLOCKS_PER_SEC);

  /*****************************************************************************************/

 

vector<TensorCPT*> marginals;
double spectral_error = 0;
double EM_error = 0;
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
		assert(spectral_val2 >= spectral_val1 && spectral_val2 >= spectral_val1);
		spectral_map = 2;
	}

	double EM_val0 = -1;
	double EM_val1 = -1;
	double EM_val2 = -1;

	int EM_map = -1;

	/********************** compute likelihood of 0 ************************************************/
	if (EM_FLAG)
	{
	if (likeZero1 > likeZero2 && likeZero1 > likeZero3 && likeZero1 > likeZero4 && likeZero1 > likeZero5)
//	if (true)
	{
		EM_val0 = EMtreeZero1.ComputeEmpiricalJointProb(marginals, EM_evidence_vars, EM_evidence_vals);
	}
	else if (likeZero2 > likeZero1 && likeZero2 > likeZero3 && likeZero2 > likeZero4 && likeZero2 > likeZero5)
	{
		EM_val0 = EMtreeZero2.ComputeEmpiricalJointProb(marginals, EM_evidence_vars, EM_evidence_vals);
	}
	else if (likeZero3 > likeZero1 && likeZero3 > likeZero2 && likeZero3 > likeZero4 && likeZero3 > likeZero5)
	{
		EM_val0 = EMtreeZero3.ComputeEmpiricalJointProb(marginals, EM_evidence_vars, EM_evidence_vals);
	}
	else if (likeZero4 > likeZero1 && likeZero4 > likeZero2 && likeZero4 > likeZero3 && likeZero4 > likeZero5)
	{
		EM_val0 = EMtreeZero4.ComputeEmpiricalJointProb(marginals, EM_evidence_vars, EM_evidence_vals);
	}
	else
	{
		assert(likeZero5 > likeZero1 && likeZero5 > likeZero2 && likeZero5 > likeZero3 && likeZero5 > likeZero4);
		EM_val0 = EMtreeZero5.ComputeEmpiricalJointProb(marginals, EM_evidence_vars, EM_evidence_vals);
	}

	/********************** compute likelihood of 1 ************************************************/
	if (likeOne1 > likeOne2 && likeOne1 > likeOne3 && likeOne1 > likeOne4 && likeOne1 > likeOne5)
//	if (true)
	{
		EM_val1 = EMtreeOne1.ComputeEmpiricalJointProb(marginals, EM_evidence_vars, EM_evidence_vals);
	}
	else if (likeOne2 > likeOne1 && likeOne2 > likeOne3 && likeOne2 > likeOne4 && likeOne2 > likeOne5)
	{
		EM_val1 = EMtreeOne2.ComputeEmpiricalJointProb(marginals, EM_evidence_vars, EM_evidence_vals);
	}
	else if (likeOne3 > likeOne1 && likeOne3 > likeOne2 && likeOne3 > likeOne4 && likeOne3 > likeOne5)
	{
		EM_val1 = EMtreeOne3.ComputeEmpiricalMarginals(marginals, EM_evidence_vars, EM_evidence_vals);
	}
	else if (likeOne4 > likeOne1 && likeOne4 > likeOne2 && likeOne4 > likeOne3 && likeOne4 > likeOne5)
	{
		EM_val1 = EMtreeOne4.ComputeEmpiricalMarginals(marginals, EM_evidence_vars, EM_evidence_vals);
	}
	else
	{
		assert(likeOne5 > likeOne1 && likeOne5 > likeOne2 && likeOne5 > likeOne3 && likeOne5 > likeOne4);
		EM_val1 = EMtreeOne5.ComputeEmpiricalMarginals(marginals, EM_evidence_vars, EM_evidence_vals);
	}

	/********************** compute likelihood of 2 ************************************************/
	if (likeTwo1 > likeTwo2 && likeTwo1 > likeTwo3 && likeTwo1 > likeTwo4 && likeTwo1 > likeTwo5)
//	if (true)
	{
		EM_val2 = EMtreeTwo1.ComputeEmpiricalJointProb(marginals, EM_evidence_vars, EM_evidence_vals);
	}
	else if (likeTwo2 > likeTwo1 && likeTwo2 > likeTwo3 && likeTwo2 > likeTwo4 && likeTwo2 > likeTwo5)
	{
		EM_val2 = EMtreeTwo2.ComputeEmpiricalMarginals(marginals, EM_evidence_vars, EM_evidence_vals);
	}
	else if (likeTwo3 > likeTwo1 && likeTwo3 > likeTwo2 && likeTwo3 > likeTwo4 && likeTwo3 > likeTwo5)
	{
		EM_val2 = EMtreeTwo3.ComputeEmpiricalMarginals(marginals, EM_evidence_vars, EM_evidence_vals);
	}
	else if (likeTwo4 > likeTwo1 && likeTwo4 > likeTwo2 && likeTwo4 > likeTwo3 && likeTwo4 > likeTwo5)
	{
		EM_val2 = EMtreeTwo4.ComputeEmpiricalMarginals(marginals, EM_evidence_vars, EM_evidence_vals);
	}
	else
	{
		assert(likeTwo5 > likeTwo1 && likeTwo5 > likeTwo2 && likeTwo5 > likeTwo3 && likeTwo5 > likeTwo4);
		EM_val2 = EMtreeTwo5.ComputeEmpiricalMarginals(marginals, EM_evidence_vars, EM_evidence_vals);
	}
	}

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

  }
  cout << "yeeee\n";
  spectral_error = spectral_error / test_samples->Dim(0);
  EM_error = EM_error / test_samples->Dim(0);
  cout << "yeeee2\n";
  double spectral_time = spectral_time0 + spectral_time1 + spectral_time2;
  cout << "yeeee3\n";
  double EM_time = EMZero1_time + EMZero2_time + EMZero3_time + EMZero4_time + EMZero5_time
	       + EMOne1_time + EMOne2_time + EMOne3_time + EMOne4_time + EMOne5_time 
		   + EMTwo1_time + EMTwo2_time + EMTwo3_time + EMTwo4_time + EMTwo5_time;
  cout << "yeeee4\n";
  cout << "total error:\n";
  cout << spectral_error;
  cout << "\n";
  cout << EM_error;
  cout << "\n";

  string outfile_name = "splice-" + to_string(num_hidden_states) + "-"
						+ to_string(num_training_examples)
						+ "-" + to_string(EM_thresh) 
						+ "-" + to_string(stack_count)
						+ "-" + to_string(max_runs)
						+ "-" + to_string(trial) + ".txt";
//  string outfile_name = "fun.txt";
  ofstream myfile;
  myfile.open(outfile_name.c_str());
  myfile << ("spectral\t" + to_string(spectral_error) + "\t" + to_string(spectral_time) + "\n");
  myfile << ("EMt" + to_string(EM_error) + "\t" + to_string(EM_time) + "\n");
  myfile << ("randseed\t" + to_string(rand_seed) + "\t" + to_string(rand_seed) + "\n");

  return 0;
}
