#include <iostream>
#include "tensor.hpp"
#include "BayesianNetwork.hpp"
#include "TensorVarElim.hpp"
#include "TensorJTree.hpp"
#include "VectorPlus.hpp"

using namespace std;

int main ()
{
	
  // node names
  int A = 0;
  int B = 1;
  int C = 2;
  int D = 3;
  int E = 4;
  int F = 5;
  int G = 6;
  int H = 7;
  int I = 8;
  int J = 9;
  int K = 10;
  int L = 11;
  int M = 12;

  int O1 = 13;
  int O2 = 14;
  int O3 = 15;
  int O4 = 16;
  int O5 = 17;
  int O6 = 18;
  int O7 = 19;
  int O8 = 20;
  int O9 = 21;
  int O10 = 22;
  int O11 = 23;
  int O12 = 24;
  int O13 = 25;
  int O14 = 26;
  int O15 = 27;
  int O16 = 28;
  int O17 = 29;
  int O18 = 30;
  int O19 = 31;
  int O20 = 32;
  int O21 = 33;
  int O22 = 34;
  int O23 = 35;
  int O24 = 36;
  int O25 = 37;
  int O26 = 38;
  int O27 = 39;
  int O28 = 40;
  int O29 = 41;
  int O30 = 42;
  int O31 = 43;
  int O32 = 44;
  int O33 = 45;
  int O34 = 46;
  int O35 = 47;
  int O36 = 48;

  int num_nodes = 49;
  vector<vector<int>*> parent_list;
  for (int i = 0; i < num_nodes; ++i)
	parent_list.push_back(new vector<int>());

  // hidden variables
  parent_list[B]->push_back(A);
  parent_list[B]->push_back(D);
  parent_list[C]->push_back(A);
  parent_list[C]->push_back(D);

  parent_list[E]->push_back(A);
  parent_list[E]->push_back(B);
  parent_list[F]->push_back(B);
  parent_list[F]->push_back(C);
  parent_list[G]->push_back(C);
  parent_list[G]->push_back(D);

  parent_list[H]->push_back(A);
  parent_list[H]->push_back(B);
  parent_list[H]->push_back(E);
  parent_list[I]->push_back(A);
  parent_list[I]->push_back(B);
  parent_list[I]->push_back(E);

  parent_list[J]->push_back(B);
  parent_list[J]->push_back(C);
  parent_list[J]->push_back(F);
  parent_list[K]->push_back(B);
  parent_list[K]->push_back(C);
  parent_list[K]->push_back(F);

  parent_list[L]->push_back(C);
  parent_list[L]->push_back(D);
  parent_list[L]->push_back(G);
  parent_list[M]->push_back(C);
  parent_list[M]->push_back(D);
  parent_list[M]->push_back(G);

  parent_list[O1]->push_back(H);
  parent_list[O2]->push_back(H);
  parent_list[O3]->push_back(H);
  parent_list[O4]->push_back(H);
  parent_list[O5]->push_back(H);
  parent_list[O6]->push_back(H);
  
  parent_list[O7]->push_back(I);
  parent_list[O8]->push_back(I);
  parent_list[O9]->push_back(I);
  parent_list[O10]->push_back(I);
  parent_list[O11]->push_back(I);
  parent_list[O12]->push_back(I);

  parent_list[O13]->push_back(J);
  parent_list[O14]->push_back(J);
  parent_list[O15]->push_back(J);
  parent_list[O16]->push_back(J);
  parent_list[O17]->push_back(J);
  parent_list[O18]->push_back(J);

  parent_list[O19]->push_back(K);
  parent_list[O20]->push_back(K);
  parent_list[O21]->push_back(K);
  parent_list[O22]->push_back(K);
  parent_list[O23]->push_back(K);
  parent_list[O24]->push_back(K);

  parent_list[O25]->push_back(L);
  parent_list[O26]->push_back(L);
  parent_list[O27]->push_back(L);
  parent_list[O28]->push_back(L);
  parent_list[O29]->push_back(L);
  parent_list[O30]->push_back(L);

  parent_list[O31]->push_back(M);
  parent_list[O32]->push_back(M);
  parent_list[O33]->push_back(M);
  parent_list[O34]->push_back(M);
  parent_list[O35]->push_back(M);
  parent_list[O36]->push_back(M);

  vector<int> type_vector(num_nodes, 0);
  for (int i = 0; i < num_nodes; ++i)
  {
	  if (i < 13)
		type_vector[i] = HIDDEN_FLAG;
	  else
		type_vector[i] = OBSERVED_FLAG;
  }


  vector<int> dims(num_nodes, 2);
  BayesianNetwork* bayesNet = new BayesianNetwork(parent_list, type_vector, dims);
  Matrix samples;
  bayesNet->GenerateSamples(samples, 100);
//  vector<int> evidence_vars(1,3);
 // vector<int> evidence_vals(1,1);
  //double prob = TensorVarElim::VE(*bayesNet, evidence_vars, evidence_vals, bayesNet->ReverseTopologicalOrder());


  // test the junction tree
  int jtree_nodes = 46;
  vector<int> tree_matrix_dims(2,jtree_nodes);
  
  int ABCD = 0;
  int ABE = 1;
  int BCF = 2;
  int CDG = 3;
  int ABEH = 4;
  int ABEI = 5;
  int BCFJ = 6;
  int BCFK = 7;
  int CDGL = 8;
  int CDGM = 9;
  
  int O1H = 10;
  int O2H = 11;
  int O3H = 12;
  int O4H = 13;
  int O5H = 14;
  int O6H = 15;
 
  int O7I = 16;
  int O8I = 17;
  int O9I = 18;
  int O10I = 19;
  int O11I = 20;
  int O12I = 21;

  int O13J = 22;
  int O14J = 23;
  int O15J = 24;
  int O16J = 25;
  int O17J = 26;
  int O18J = 27;

  int O19K = 28;
  int O20K = 29;
  int O21K = 30;
  int O22K = 31;
  int O23K = 32;
  int O24K = 33;

  int O25L = 34;
  int O26L = 35;
  int O27L = 36;
  int O28L = 37;
  int O29L = 38;
  int O30L = 39;

  int O31M = 40;
  int O32M = 41;
  int O33M = 42;
  int O34M = 43;
  int O35M = 44;
  int O36M = 45;

  Matrix tree_matrix(tree_matrix_dims);
  tree_matrix.Set(ABE, ABCD, 1);
  tree_matrix.Set(BCF, ABCD, 1);
  tree_matrix.Set(CDG, ABCD, 1);

  tree_matrix.Set(ABEH, ABE, 1);
  tree_matrix.Set(ABEI, ABE, 1);

  tree_matrix.Set(BCFJ, BCF, 1);
  tree_matrix.Set(BCFK, BCF, 1);

  tree_matrix.Set(CDGL, CDG, 1);
  tree_matrix.Set(CDGM, CDG, 1);

  tree_matrix.Set(O1H, ABEH, 1);
  tree_matrix.Set(O2H, ABEH, 1);
  tree_matrix.Set(O3H, ABEH, 1);
  tree_matrix.Set(O4H, ABEH, 1);
  tree_matrix.Set(O5H, ABEH, 1);
  tree_matrix.Set(O6H, ABEH, 1);

  tree_matrix.Set(O1H, ABEH, 1);
  tree_matrix.Set(O2H, ABEH, 1);
  tree_matrix.Set(O3H, ABEH, 1);
  tree_matrix.Set(O4H, ABEH, 1);
  tree_matrix.Set(O5H, ABEH, 1);
  tree_matrix.Set(O6H, ABEH, 1);

  tree_matrix.Set(O1H, ABEH, 1);
  tree_matrix.Set(O2H, ABEH, 1);
  tree_matrix.Set(O3H, ABEH, 1);
  tree_matrix.Set(O4H, ABEH, 1);
  tree_matrix.Set(O5H, ABEH, 1);
  tree_matrix.Set(O6H, ABEH, 1);

  tree_matrix.Set(O7I, ABEI, 1);
  tree_matrix.Set(O8I, ABEI, 1);
  tree_matrix.Set(O9I, ABEI, 1);
  tree_matrix.Set(O10I, ABEI, 1);
  tree_matrix.Set(O11I, ABEI, 1);
  tree_matrix.Set(O12I, ABEI, 1);

  tree_matrix.Set(O13J, BCFJ, 1);
  tree_matrix.Set(O14J, BCFJ, 1);
  tree_matrix.Set(O15J, BCFJ, 1);
  tree_matrix.Set(O16J, BCFJ, 1);
  tree_matrix.Set(O17J, BCFJ, 1);
  tree_matrix.Set(O18J, BCFJ, 1);

  tree_matrix.Set(O19K, BCFK, 1);
  tree_matrix.Set(O20K, BCFK, 1);
  tree_matrix.Set(O21K, BCFK, 1);
  tree_matrix.Set(O22K, BCFK, 1);
  tree_matrix.Set(O23K, BCFK, 1);
  tree_matrix.Set(O24K, BCFK, 1);

  tree_matrix.Set(O25L, CDGL, 1);
  tree_matrix.Set(O26L, CDGL, 1);
  tree_matrix.Set(O27L, CDGL, 1);
  tree_matrix.Set(O28L, CDGL, 1);
  tree_matrix.Set(O29L, CDGL, 1);
  tree_matrix.Set(O30L, CDGL, 1);

  tree_matrix.Set(O31M, CDGM, 1);
  tree_matrix.Set(O32M, CDGM, 1);
  tree_matrix.Set(O33M, CDGM, 1);
  tree_matrix.Set(O34M, CDGM, 1);
  tree_matrix.Set(O35M, CDGM, 1);
  tree_matrix.Set(O36M, CDGM, 1);

  map<int, vector<int>*> nodes_to_vars;
  for (int i = 0; i < jtree_nodes; ++i)
	nodes_to_vars[i] = new vector<int>();
  
  nodes_to_vars[ABCD]->push_back(A);
  nodes_to_vars[ABCD]->push_back(B);
  nodes_to_vars[ABCD]->push_back(C);
  nodes_to_vars[ABCD]->push_back(D);

  nodes_to_vars[ABE]->push_back(A);
  nodes_to_vars[ABE]->push_back(B);
  nodes_to_vars[ABE]->push_back(E);

  nodes_to_vars[BCF]->push_back(B);
  nodes_to_vars[BCF]->push_back(C);
  nodes_to_vars[BCF]->push_back(F);

  nodes_to_vars[CDG]->push_back(C);
  nodes_to_vars[CDG]->push_back(D);
  nodes_to_vars[CDG]->push_back(G);

  nodes_to_vars[ABEH]->push_back(A);
  nodes_to_vars[ABEH]->push_back(B);
  nodes_to_vars[ABEH]->push_back(E);
  nodes_to_vars[ABEH]->push_back(H);

  nodes_to_vars[ABEI]->push_back(A);
  nodes_to_vars[ABEI]->push_back(B);
  nodes_to_vars[ABEI]->push_back(E);
  nodes_to_vars[ABEI]->push_back(I);

  nodes_to_vars[BCFJ]->push_back(B);
  nodes_to_vars[BCFJ]->push_back(C);
  nodes_to_vars[BCFJ]->push_back(F);
  nodes_to_vars[BCFJ]->push_back(J);

  nodes_to_vars[BCFK]->push_back(B);
  nodes_to_vars[BCFK]->push_back(C);
  nodes_to_vars[BCFK]->push_back(F);
  nodes_to_vars[BCFK]->push_back(K);

  nodes_to_vars[CDGL]->push_back(C);
  nodes_to_vars[CDGL]->push_back(D);
  nodes_to_vars[CDGL]->push_back(G);
  nodes_to_vars[CDGL]->push_back(L);

  nodes_to_vars[CDGM]->push_back(C);
  nodes_to_vars[CDGM]->push_back(D);
  nodes_to_vars[CDGM]->push_back(G);
  nodes_to_vars[CDGM]->push_back(M);

  nodes_to_vars[O1H]->push_back(O1);
  nodes_to_vars[O1H]->push_back(H);
  nodes_to_vars[O2H]->push_back(O2);
  nodes_to_vars[O2H]->push_back(H);
  nodes_to_vars[O3H]->push_back(O3);
  nodes_to_vars[O3H]->push_back(H);
  nodes_to_vars[O4H]->push_back(O4);
  nodes_to_vars[O4H]->push_back(H);
  nodes_to_vars[O5H]->push_back(O5);
  nodes_to_vars[O5H]->push_back(H);
  nodes_to_vars[O6H]->push_back(O6);
  nodes_to_vars[O6H]->push_back(H);

  nodes_to_vars[O7I]->push_back(O7);
  nodes_to_vars[O7I]->push_back(I);
  nodes_to_vars[O8I]->push_back(O8);
  nodes_to_vars[O8I]->push_back(I);
  nodes_to_vars[O9I]->push_back(O9);
  nodes_to_vars[O9I]->push_back(I);
  nodes_to_vars[O10I]->push_back(O10);
  nodes_to_vars[O10I]->push_back(I);
  nodes_to_vars[O11I]->push_back(O11);
  nodes_to_vars[O11I]->push_back(I);
  nodes_to_vars[O12I]->push_back(O12);
  nodes_to_vars[O12I]->push_back(I);

  nodes_to_vars[O13J]->push_back(O13);
  nodes_to_vars[O13J]->push_back(J);
  nodes_to_vars[O14J]->push_back(O14);
  nodes_to_vars[O14J]->push_back(J);
  nodes_to_vars[O15J]->push_back(O15);
  nodes_to_vars[O15J]->push_back(J);
  nodes_to_vars[O16J]->push_back(O16);
  nodes_to_vars[O16J]->push_back(J);
  nodes_to_vars[O17J]->push_back(O17);
  nodes_to_vars[O17J]->push_back(J);
  nodes_to_vars[O18J]->push_back(O18);
  nodes_to_vars[O18J]->push_back(J);

  nodes_to_vars[O19K]->push_back(O19);
  nodes_to_vars[O19K]->push_back(K);
  nodes_to_vars[O20K]->push_back(O20);
  nodes_to_vars[O20K]->push_back(K);
  nodes_to_vars[O21K]->push_back(O21);
  nodes_to_vars[O21K]->push_back(K);
  nodes_to_vars[O22K]->push_back(O22);
  nodes_to_vars[O22K]->push_back(K);
  nodes_to_vars[O23K]->push_back(O23);
  nodes_to_vars[O23K]->push_back(K);
  nodes_to_vars[O24K]->push_back(O24);
  nodes_to_vars[O24K]->push_back(K);

  nodes_to_vars[O25L]->push_back(O25);
  nodes_to_vars[O25L]->push_back(L);
  nodes_to_vars[O26L]->push_back(O26);
  nodes_to_vars[O26L]->push_back(L);
  nodes_to_vars[O27L]->push_back(O27);
  nodes_to_vars[O27L]->push_back(L);
  nodes_to_vars[O28L]->push_back(O28);
  nodes_to_vars[O28L]->push_back(L);
  nodes_to_vars[O29L]->push_back(O29);
  nodes_to_vars[O29L]->push_back(L);
  nodes_to_vars[O30L]->push_back(O30);
  nodes_to_vars[O30L]->push_back(L);

  nodes_to_vars[O31M]->push_back(O31);
  nodes_to_vars[O31M]->push_back(M);
  nodes_to_vars[O32M]->push_back(O32);
  nodes_to_vars[O32M]->push_back(M);
  nodes_to_vars[O33M]->push_back(O33);
  nodes_to_vars[O33M]->push_back(M);
  nodes_to_vars[O34M]->push_back(O34);
  nodes_to_vars[O34M]->push_back(M);
  nodes_to_vars[O35M]->push_back(O35);
  nodes_to_vars[O35M]->push_back(M);
  nodes_to_vars[O36M]->push_back(O36);
  nodes_to_vars[O36M]->push_back(M);


  TensorJTree jtree(*bayesNet, tree_matrix, nodes_to_vars);
  jtree.LearnSpectralParameters();
  vector<int> evidence_vars;
  VectorPlus::Seq(evidence_vars, O1, 1, O36 + 1);
  vector<int> evidence_vals(evidence_vars.size(), 0);
  double val = jtree.ComputeMarginalEmpiricalProbability(evidence_vars, evidence_vals);
  double prob = TensorVarElim::VE(*bayesNet, evidence_vars, evidence_vals, bayesNet->ReverseTopologicalOrder());

  /* Eigen::MatrixXd m(2,2);
  m(0,0) = 3;
  m(1,0) = 2.5;
  m(0,1) = -1;
  m(1,1) = m(1,0) + m(0,1); */

  vector<int> fun_dims;
  fun_dims.push_back(2);
  fun_dims.push_back(4);
  fun_dims.push_back(6);
  Tensor* t = new Tensor(fun_dims);
  t->FillWithRandom();

  Tensor* u = new Tensor(fun_dims);
  u->FillWithRandom();

  Tensor* v = new Tensor();

 // Tensor::Add(*v, *t, *u);
  vector<int> mult_dims(2,0);
  mult_dims[0] = 1;
  mult_dims[1] = 2;

  Tensor::Multiply(*v, *t, *u, mult_dims, mult_dims); 
  double inner_product = Tensor::InnerProduct(*t, *u);
  vector<int> index_array(8,0);

  for (int i = 0; i < 2; ++i)
  {
	  for (int j = 0; j < 4; ++j)
	  {
		  index_array[0] = i;
		  index_array[1] = j;

		  t->Set(index_array, i + j);
	  }
  }

    for (int i = 0; i < 2; ++i)
	{
	  for (int j = 0; j < 4; ++j)
	  {
		  index_array[0] = i;
		  index_array[1] = j;

		  cout << t->At(index_array);
		  cout << "\n";
	  }
  }

  return 0;
}