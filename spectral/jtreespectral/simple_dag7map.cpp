
#include "tensor.hpp"
#include "BayesianNetwork.hpp"
#include "TensorVarElim.hpp"
#include "TensorJTree.hpp"
#include "VectorPlus.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include "JTree.hpp"
#include "time.h"

using namespace std;

int main ()
{
  srand(time(NULL));
	
  // node names
  int A = 0;
  int B = 1;
  int C = 2;
  int D = 3;
  int E = 4;

  int O1 = 5;
  int O2 = 6;
  int O3 = 7;
  int O4 = 8;
  int O5 = 9;
  int O6 = 10;
  int O7 = 11;
  int O8 = 12;
  int O9 = 13;
  int O10 = 14;
  int O11 = 15;
  int O12 = 16;
  int O13 = 17;
  int O14 = 18;
  int O15 = 19;
  int O16 = 20;
  int O17 = 21;
  int O18 = 22;
  int O19 = 23;
  int O20 = 24;
  int O21 = 25;
  int O22 = 26;
  int O23 = 27;
  int O24 = 28;

  int num_nodes = 29;
  vector<vector<int>*> parent_list;
  for (int i = 0; i < num_nodes; ++i)
	parent_list.push_back(new vector<int>());

  parent_list[B]->push_back(A);
  parent_list[C]->push_back(A);
  parent_list[D]->push_back(A);
  parent_list[D]->push_back(B);
  parent_list[E]->push_back(B);

  // hidden variables
  parent_list[O1]->push_back(A);
  parent_list[O2]->push_back(A);
  parent_list[O3]->push_back(A);
  parent_list[O4]->push_back(A);
  parent_list[O5]->push_back(A);
  parent_list[O6]->push_back(A);
  parent_list[O7]->push_back(A);
  parent_list[O8]->push_back(A);
  parent_list[O9]->push_back(A);
  parent_list[O10]->push_back(A);
  parent_list[O11]->push_back(A);
  parent_list[O12]->push_back(A);

  parent_list[O1]->push_back(C);
  parent_list[O2]->push_back(C);
  parent_list[O3]->push_back(C);
  parent_list[O4]->push_back(C);
  parent_list[O5]->push_back(C);
  parent_list[O6]->push_back(C);
  parent_list[O7]->push_back(C);
  parent_list[O8]->push_back(C);
  parent_list[O9]->push_back(C);
  parent_list[O10]->push_back(C);
  parent_list[O11]->push_back(C);
  parent_list[O12]->push_back(C);


  parent_list[O13]->push_back(D);
  parent_list[O14]->push_back(D);
  parent_list[O15]->push_back(D);
  parent_list[O16]->push_back(D);
  parent_list[O17]->push_back(D);
  parent_list[O18]->push_back(D);

  parent_list[O19]->push_back(E);
  parent_list[O20]->push_back(E);
  parent_list[O21]->push_back(E);
  parent_list[O22]->push_back(E);
  parent_list[O23]->push_back(E);
  parent_list[O24]->push_back(E);


  vector<int> type_vector(num_nodes, 0);
   vector<int> dims(num_nodes, 2);
  for (int i = 1; i < num_nodes; ++i)
  {
	  if (i < 5)
	  {
		type_vector[i] = HIDDEN_FLAG;
		dims[i] = 2;
	  }
	  else
	  {
		type_vector[i] = OBSERVED_FLAG;
		dims[i] = 4;
	  }
  }
  vector<int> template_vec;
  VectorPlus::Seq(template_vec, 0, 1, num_nodes);
  TemplateBayesianNetwork* bayesNet = new TemplateBayesianNetwork(parent_list, type_vector, dims, template_vec);
  TemplateBayesianNetwork* bayesNet2 = new TemplateBayesianNetwork(parent_list, type_vector, dims, template_vec);
 /* string line;
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
  Matrix samples;
  Matrix test_samples;
  bayesNet->GenerateTrainSamples(samples, 1000);
  bayesNet2->GenerateTestSamples(test_samples, 500);
  // test the junction tree
  int jtree_nodes = 27;
  vector<int> tree_matrix_dims(2,jtree_nodes);
  
  //int A = 0;
  int BE = 0;
  int ABD = 1;
  int AC = 2;

  int O1AC = 3;
  int O2AC = 4;
  int O3AC = 5;
  int O4AC = 6;
  int O5AC = 7;
  int O6AC = 8;
  int O7AC = 9;
  int O8AC = 10;
  int O9AC = 11;
  int O10AC = 12;
  int O11AC = 13;
  int O12AC = 14;

  int O13D = 15;
  int O14D = 16;
  int O15D = 17;
  int O16D = 18;
  int O17D = 19;
  int O18D = 20;

  int O19E = 21;
  int O20E = 22;
  int O21E = 23;
  int O22E = 24;
  int O23E = 25;
  int O24E = 26;

  Matrix tree_matrix(tree_matrix_dims);
  tree_matrix.Set(ABD, BE, 1);
  tree_matrix.Set(AC, ABD, 1);
  tree_matrix.Set(O1AC, AC, 1);
  tree_matrix.Set(O2AC, AC, 1);
  tree_matrix.Set(O3AC, AC, 1);
  tree_matrix.Set(O4AC, AC, 1);
  tree_matrix.Set(O5AC, AC, 1);
  tree_matrix.Set(O6AC, AC, 1);
  tree_matrix.Set(O7AC, AC, 1);
  tree_matrix.Set(O8AC, AC, 1);
  tree_matrix.Set(O9AC, AC, 1);
  tree_matrix.Set(O10AC, AC, 1);
  tree_matrix.Set(O11AC, AC, 1);
  tree_matrix.Set(O12AC, AC, 1);
  tree_matrix.Set(O13D, ABD, 1);
  tree_matrix.Set(O14D, ABD, 1);
  tree_matrix.Set(O15D, ABD, 1);
  tree_matrix.Set(O16D, ABD, 1);
  tree_matrix.Set(O17D, ABD, 1);
  tree_matrix.Set(O18D, ABD, 1);
  tree_matrix.Set(O19E, BE, 1);
  tree_matrix.Set(O20E, BE, 1);
  tree_matrix.Set(O21E, BE, 1);
  tree_matrix.Set(O22E, BE, 1);
  tree_matrix.Set(O23E, BE, 1);
  tree_matrix.Set(O24E, BE, 1);

  map<int, vector<int>*> nodes_to_vars;
  for (int i = 0; i < jtree_nodes; ++i)
	nodes_to_vars[i] = new vector<int>();
  
  nodes_to_vars[BE]->push_back(B);
  nodes_to_vars[BE]->push_back(E);

  nodes_to_vars[ABD]->push_back(A);
  nodes_to_vars[ABD]->push_back(B);
  nodes_to_vars[ABD]->push_back(D);

  nodes_to_vars[AC]->push_back(A);
  nodes_to_vars[AC]->push_back(C);

  nodes_to_vars[O1AC]->push_back(O1);
  nodes_to_vars[O2AC]->push_back(O2);
  nodes_to_vars[O3AC]->push_back(O3);
  nodes_to_vars[O4AC]->push_back(O4);
  nodes_to_vars[O5AC]->push_back(O5);
  nodes_to_vars[O6AC]->push_back(O6);
  nodes_to_vars[O7AC]->push_back(O7);
  nodes_to_vars[O8AC]->push_back(O8);
  nodes_to_vars[O9AC]->push_back(O9);
  nodes_to_vars[O10AC]->push_back(O10);
  nodes_to_vars[O11AC]->push_back(O11);
  nodes_to_vars[O12AC]->push_back(O12);

  nodes_to_vars[O1AC]->push_back(A);
  nodes_to_vars[O2AC]->push_back(A);
  nodes_to_vars[O3AC]->push_back(A);
  nodes_to_vars[O4AC]->push_back(A);
  nodes_to_vars[O5AC]->push_back(A);
  nodes_to_vars[O6AC]->push_back(A);
  nodes_to_vars[O7AC]->push_back(A);
  nodes_to_vars[O8AC]->push_back(A);
  nodes_to_vars[O9AC]->push_back(A);
  nodes_to_vars[O10AC]->push_back(A);
  nodes_to_vars[O11AC]->push_back(A);
  nodes_to_vars[O12AC]->push_back(A);

  nodes_to_vars[O1AC]->push_back(C);
  nodes_to_vars[O2AC]->push_back(C);
  nodes_to_vars[O3AC]->push_back(C);
  nodes_to_vars[O4AC]->push_back(C);
  nodes_to_vars[O5AC]->push_back(C);
  nodes_to_vars[O6AC]->push_back(C);
  nodes_to_vars[O7AC]->push_back(C);
  nodes_to_vars[O8AC]->push_back(C);
  nodes_to_vars[O9AC]->push_back(C);
  nodes_to_vars[O10AC]->push_back(C);
  nodes_to_vars[O11AC]->push_back(C);
  nodes_to_vars[O12AC]->push_back(C);

  nodes_to_vars[O13D]->push_back(O13);
  nodes_to_vars[O14D]->push_back(O14);
  nodes_to_vars[O15D]->push_back(O15);
  nodes_to_vars[O16D]->push_back(O16);
  nodes_to_vars[O17D]->push_back(O17);
  nodes_to_vars[O18D]->push_back(O18);

  nodes_to_vars[O13D]->push_back(D);
  nodes_to_vars[O14D]->push_back(D);
  nodes_to_vars[O15D]->push_back(D);
  nodes_to_vars[O16D]->push_back(D);
  nodes_to_vars[O17D]->push_back(D);
  nodes_to_vars[O18D]->push_back(D);

  nodes_to_vars[O19E]->push_back(O19);
  nodes_to_vars[O20E]->push_back(O20);
  nodes_to_vars[O21E]->push_back(O21);
  nodes_to_vars[O22E]->push_back(O22);
  nodes_to_vars[O23E]->push_back(O23);
  nodes_to_vars[O24E]->push_back(O24);

  nodes_to_vars[O19E]->push_back(E);
  nodes_to_vars[O20E]->push_back(E);
  nodes_to_vars[O21E]->push_back(E);
  nodes_to_vars[O22E]->push_back(E);
  nodes_to_vars[O23E]->push_back(E);
  nodes_to_vars[O24E]->push_back(E);

  vector<int> jtemplate_vec;
  VectorPlus::Seq(jtemplate_vec, 0, 1, jtree_nodes);

  vector<int> evidence_vars;
  VectorPlus::Seq(evidence_vars, O1, 1, O24 + 1);
  vector<int> evidence_vals(evidence_vars.size(), 0);

  double total_error = 0;
  double EM_error = 0;
  JTree EMtree(bayesNet, tree_matrix, nodes_to_vars, &samples);
  EMtree.LearnEMParameters(20, 0);
  TensorJTree jtree(bayesNet, tree_matrix, nodes_to_vars, 6, jtemplate_vec);
  jtree.LearnSpectralParameters();
  vector<TensorCPT*> marginals;
  for (int n = 0; n < test_samples.Dim(0); ++n)
  {
	int map_var =  (rand() % 24) + O1;

	int oracle_map = -1;
	int EM_map = -1;
	int spectral_map = -1;

	double oracle_map_prob = -1;
	double EM_map_prob = -1;
	double spectral_map_prob = -1;

	for (int k = 0; k < dims[map_var]; ++k)
	{
		for (int i = 0; i < evidence_vals.size(); ++i)
		{
			if (i != map_var)
				evidence_vals.at(i) = test_samples.At(n,evidence_vars[i]);
			else
				evidence_vals.at(i) = k;
		}
		double prob = TensorVarElim::VE(*bayesNet, evidence_vars, evidence_vals, bayesNet->ReverseTopologicalOrder());

		double EM_val = EMtree.ComputeEmpiricalMarginals(marginals, evidence_vars, evidence_vals);

		double spectral_val = jtree.ComputeMarginalEmpiricalProbability(evidence_vars, evidence_vals);


		if (spectral_val < 0)
			spectral_val = 0;
		if (spectral_val > 1)
			spectral_val = 1;

		if (prob > oracle_map_prob)
		{
			oracle_map_prob = prob;
			oracle_map = k;
		}
		if (EM_val > EM_map_prob)
		{
			EM_map_prob = EM_val;
			EM_map = k;
		}
		if (spectral_val > spectral_map_prob)
		{
			spectral_map_prob = spectral_val;
			spectral_map = k;
		}
	}
  


	cout << EM_map;
	cout << "\t";
	cout << spectral_map;
	cout << "\t";
	cout << oracle_map; 
	cout << "\n";

	if (oracle_map != spectral_map)
		total_error += 1; 
	if (oracle_map != EM_map)
		EM_error += 1;
  }
  cout << "total error:\n";
  cout << total_error / test_samples.Dim(0);
  cout << "\n";
  cout << EM_error / test_samples.Dim(0);
  delete(bayesNet);
  return 0;
}