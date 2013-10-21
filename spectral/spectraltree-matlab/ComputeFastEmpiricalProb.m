% Ankur Parikh
% SAILING Lab
% Carnegie Mellon University

% compute empirical joint distribution of nodes in node_array
function [fast_prob_matrix, EmpiricalProbMap] = ComputeFastEmpiricalProb(node_array, K, samples, type_indicator, tree_matrix, CPT, EmpiricalProbMap)
        
    if (hashmap_contains(EmpiricalProbMap, node_array))
        fast_prob_matrix = hashmap_get(EmpiricalProbMap, node_array);
        return;
    end
    epsilon = 0.000000001; 
    num_vars = length(node_array);
    num_samples = size(samples, 1);
    data_matrix = samples(:,node_array)';
    v = K.^(0:num_vars-1)';
    idx = v' * (data_matrix - 1) + 1; 
    S = sparse(idx, (1:num_samples), ones(1,num_samples), K^num_vars, num_samples);
    meanS = mean(S, 2);
    if (num_vars == 1)
        fast_prob_matrix = full(meanS);
    else
        fast_prob_matrix = reshape(full(meanS), K*ones(1,num_vars));
    end
    fast_prob_matrix = fast_prob_matrix + epsilon;
    fast_prob_matrix = fast_prob_matrix / sum(sum(sum(fast_prob_matrix)));
%     assert(sum(type_indicator(node_array)) == length(node_array));
%     N = size(samples, 1);
%     epsilon = 0.00000001;
%     num_vars = length(node_array);
%     if (num_vars == 1)
%         prob_matrix = zeros(K, 1) + epsilon;
%     else
%        prob_matrix = repmat(epsilon, K*ones(length(node_array),1)'); 
%     end 
%     
%     num_combos = K^num_vars; % assume all variables take on K values
%    
%     if (length(node_array) == 1)
%         for i=1:1:K
%              relevant_samples = (samples(:,node_array(1)) == i);
%              prob =  sum(relevant_samples) / N; 
%              prob_matrix(i) = prob_matrix(i) + prob;
%         end
%     elseif (length(node_array) == 2)
%         for i=1:1:K
%             relevant_samples_i = (samples(:,node_array(1)) == i); % filter out first group
%             for j=1:1:K
%                 relevant_samples = relevant_samples_i & (samples(:,node_array(2)) == j);
%                 prob =  sum(relevant_samples) / N; 
%                 prob_matrix(i,j) = prob_matrix(i,j) + prob;
%             end
%         end
%     elseif (length(node_array) == 3)
%       for i=1:1:K
%             relevant_samples_i = (samples(:,node_array(1)) == i); % filter out first group
%             for j=1:1:K
%                 relevant_samples_j = relevant_samples_i & (samples(:,node_array(2)) == j);
%                 for k=1:1:K
%                     relevant_samples = relevant_samples_j & (samples(:,node_array(3)) == k);
%                     prob =  sum(relevant_samples) / N; 
%                     prob_matrix(i,j,k) = prob_matrix(i,j,k) + prob;
%                 end
%             end
%        end
%     else
%         assert(0);
%     end
%    
% 
%     prob_matrix = prob_matrix / sum(sum(sum(sum(prob_matrix))));
%     assert(sum(sum(sum(abs(prob_matrix - fast_prob_matrix)))) < .0000001);
    EmpiricalProbMap = hashmap_put(EmpiricalProbMap, node_array, fast_prob_matrix);
end

function val = hashmap_contains(hash_map, key)

    str_key = int2str(key(1));
    for i=2:1:length(key)
        str_key = strcat(str_key, '+', int2str(key(i)));
    end
    
    val = hash_map.containsKey(str_key);
end

function val = hashmap_get(hash_map, key)
    str_key = int2str(key(1));
    for i=2:1:length(key)
        str_key = strcat(str_key, '+', int2str(key(i)));
    end
    
    val = hash_map.get(str_key);
end

function hash_map = hashmap_put(hash_map, key, value)
    str_key = int2str(key(1));
    for i=2:1:length(key)
        str_key = strcat(str_key, '+', int2str(key(i)));
    end
    
    hash_map.put(str_key, value);
end

