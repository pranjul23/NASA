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

 /* int num_hidden_nodes = 8;
  int num_hidden_states = 2;
  int num_observed_nodes = 8;
  int num_observed_states = 4;
  int num_training_examples = 10;
  double EM_thresh = 0.001;
  int stack_count = 1;
  int max_runs = 1;

  int trial = 1; */

  int rand_seed = 596 * trial;
  srand(rand_seed);

   int skip = 2;

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
  vector<int> omitted_obs;
  for (int i = 0; i < skip; ++i)
  {
	omitted_obs.push_back(observed_nodes[i]);
  }

  VectorPlus::SetDiff(evidence_vars, observed_nodes, omitted_obs);
  vector<int> evidence_vals(evidence_vars.size(), 0);
  cout << "yee1";
   clock_t spectral_init = clock();
  TensorJTree jtree(bayesNet, tree_matrix, nodes_to_vars, stack_count, jtemplate_vec, max_runs);
  jtree.LearnSpectralParameters();
  clock_t spectral_final = clock() - spectral_init;
  double spectral_time = ((double)spectral_final) / ((double) CLOCKS_PER_SEC);

/*  clock_t onlineEM1_init = clock();
  TemplateJTree EMtreeOnline(bayesNet, EM_tree_matrix, EM_nodes_to_vars, train_samples, EM_thresh);
  double onlinelike1 = EMtreeOnline.LearnOnlineEMParameters(0.9, 1, 1);
  clock_t onlineEM1_final = clock() - onlineEM1_init;
  double onlineEM1_time = ((double)onlineEM1_final) / ((double) CLOCKS_PER_SEC); */

    cout << "yee2";
  clock_t EM1_init = clock();
  JTree EMtree(bayesNet, EM_tree_matrix, EM_nodes_to_vars, train_samples, EM_thresh, 1);
//  double like1 = EMtree.LearnOnlineEMParameters(0.8, 1 , 1000);
  clock_t EM1_final = clock() - EM1_init;
  double EM1_time = ((double)EM1_final) / ((double) CLOCKS_PER_SEC);
    cout << "yee3";
	  cout << "yee4";

  clock_t EM2_init = clock();
  JTree EMtree2(bayesNet, EM_tree_matrix, EM_nodes_to_vars, train_samples, EM_thresh, 1);
// double like2 = EMtree2.LearnOnlineEMParameters(0.8, 1, 1000);
  clock_t EM2_final = clock() - EM2_init;
  double EM2_time = ((double)EM2_final) / ((double) CLOCKS_PER_SEC);
    cout << "yee5";


  clock_t EM3_init = clock();
  JTree EMtree3(bayesNet, EM_tree_matrix, EM_nodes_to_vars, train_samples, EM_thresh, 1);
 // double like3 = EMtree3.LearnOnlineEMParameters(0.7, 0, 1000);
  clock_t EM3_final = clock() - EM3_init;
  double EM3_time = ((double)EM3_final) / ((double) CLOCKS_PER_SEC);

  clock_t EM4_init = clock();
  JTree EMtree4(bayesNet, EM_tree_matrix, EM_nodes_to_vars, train_samples, EM_thresh, 1);
//  double like4 = EMtree4.LearnOnlineEMParameters(0.6,1,100000);
  clock_t EM4_final = clock() - EM4_init;
  double EM4_time = ((double)EM4_final) / ((double) CLOCKS_PER_SEC);

  clock_t EM5_init = clock();
  JTree EMtree5(bayesNet, EM_tree_matrix, EM_nodes_to_vars, train_samples, EM_thresh, 1);
//  double like5 = EMtree5.LearnOnlineEMParameters(0.6, 1, 10000);
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
		spectral_val = 0;

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

	spectral_error += abs(prob - spectral_val) / prob;
	EM_error1 += abs(prob - EM_val) / prob;
	EM_error2 += abs(prob - EM_val2) / prob;
	EM_error3 += abs(prob - EM_val3) / prob;
	EM_error4 += abs(prob - EM_val4) / prob;
	EM_error5 += abs(prob - EM_val5) / prob;

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

  string outfile_name = "second-nh-hmm-" + to_string(num_hidden_nodes) + "-" + to_string(num_hidden_states) + "-" 
						+  to_string(num_observed_states) + "-" + to_string(num_training_examples)
						+ "-" + to_string(EM_thresh) + "-" + to_string(stack_count) + "-" + to_string(max_runs)
						+ "-" + to_string(trial) + ".txt";
//  string outfile_name = "fun.txt";
  ofstream myfile;
  myfile.open(outfile_name.c_str());
 /* myfile << ("spectral\t" + to_string(spectral_error) + "\t" + to_string(spectral_time) + "\t" + to_string(spectral_time) + "\n");
  myfile << ("EM1\t" + to_string(EM_error1) + "\t" + to_string(EM1_time) + "\t" + to_string(like1) + "\n");
  myfile << ("EM2\t" + to_string(EM_error2) + "\t" + to_string(EM2_time) + "\t" + to_string(like2) + "\n");
  myfile << ("EM3\t" + to_string(EM_error3) + "\t" + to_string(EM3_time) + "\t" + to_string(like3) + "\n");
  myfile << ("EM4\t" + to_string(EM_error4) + "\t" + to_string(EM4_time) + "\t" + to_string(like4) + "\n");
  myfile << ("EM5\t" + to_string(EM_error5) + "\t" + to_string(EM5_time) + "\t" + to_string(like5) + "\n");
  myfile << ("randseed\t" + to_string(rand_seed) + "\t" + to_string(rand_seed) + to_string(rand_seed) + "\n"); 
  cout << "\n\n\nat end\n"; */
  return 0;
}
