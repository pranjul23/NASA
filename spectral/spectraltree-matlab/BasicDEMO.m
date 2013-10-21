% Ankur Parikh
% SAILING Lab
% Carnegie Mellon University

% This is a demo file to demonstrate how to use the spectral and EM
% functions

function [] = BasicDEMO()
 
    % model parameters
    Ko = 6; % number of observed states
    Kh = 2; % number of hidden states
    TREE_DEPTH = 4; % depth of tree
    
    % number of samples
    NUM_TRAINING_SAMPLES = 10000;
    NUM_TEST_SAMPLES = 500;
    
    % EM parameters
    EM_THRESH = 0.0001; % based on likelihood
    NUM_RESTARTS = 5;

    % generate a tree randomly
    input_tree_struct = GenerateTreeMatrixTypeDeep(Ko, Kh, TREE_DEPTH);

    % directed tree matrix, (i,j) = 1 means that i's parent is j
    tree_matrix = input_tree_struct.tree_matrix; 
    
    
    type_indicator = input_tree_struct.type_indicator; % vector indicating which variables are hidden and which observed (0 = hidden
    
    % training and test samples, note that these are N by V matrices, where
    % N is number of samples, and V is TOTAL number of variables (hidden +
    % observed)
    % However, all algorithms will only use the observed variables.
    training_samples = GenerateMatrixTreeSamples(input_tree_struct, NUM_TRAINING_SAMPLES);
    test_samples = GenerateMatrixTreeSamples(input_tree_struct, NUM_TEST_SAMPLES);

    % EM
    %EM_CPT = TreeEM(input_tree_struct.tree_matrix, training_samples, type_indicator, Ko, Kh, EM_THRESH, NUM_RESTARTS);

     
     % Calling spectral to "train"
     % transform_array is cell array of parameters (analogous to CPTs, but
     % the observable representation)
     % The last argument controls the size of the linear system, larger
     % number increases accuracy (but with more computational cost)
     % For description of Linear System trick see Parikh et al. A Spectral
     % Algorithm for Latent Junction Trees, UAI 2012.
     transform_array = SpectralTrainSystem(tree_matrix,  type_indicator, Ko, Kh, training_samples, 10);
     
     EM_error_array =  zeros(NUM_TEST_SAMPLES, 1);
     spectral_error_array = zeros(NUM_TEST_SAMPLES, 1); 
     
     % test
     for n=1:1:NUM_TEST_SAMPLES
            
          test_vars = find(type_indicator); % all leaves of tree
          test_vals = test_samples(n,test_vars);

          % compute true answer
          [oracle_joint_prob] = RunUpwardBeliefPropagation(input_tree_struct.tree_matrix, input_tree_struct.CPT, test_vars, test_vals);         
          
          % compute EM estimate
          [EM_joint_prob] = RunUpwardBeliefPropagation(input_tree_struct.tree_matrix, EM_CPT, test_vars, test_vals);

          % compute spectral estimate, NOTE THAT TEST_VARS IS ASSUMED TO BE
          % ALL LEAVES OF TREE (this is a limitation of the code/implementation, not the
          % actual algorithm)
          spectral_joint_prob = SpectralTest(tree_matrix,  type_indicator, Ko, Kh, transform_array, test_vars, test_vals);
          
          % measure error
          EM_error_array(n) = abs(EM_joint_prob - oracle_joint_prob) / oracle_joint_prob;
          spectral_error_array(n) = abs(spectral_joint_prob - oracle_joint_prob) / oracle_joint_prob;
                   
     end
     
     mean(spectral_error_array)
     mean(EM_error_array)
end
