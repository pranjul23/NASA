% Ankur Parikh
% SAILING Lab
% Carnegie Mellon University

% use direct sampling to generate samples from joint distribution
% assume root is node 1

function samples = GenerateMatrixTreeSamples(input_tree_struct, NUM_SAMPLES)

    
    tree_matrix = input_tree_struct.tree_matrix;
    CPT = input_tree_struct.CPT;
    type_indicator = logical(input_tree_struct.type_indicator);
    
    num_obs = sum(type_indicator);
    total_vars = length(type_indicator);
    
    samples = zeros(NUM_SAMPLES, total_vars);
        
    search_queue = java.util.LinkedList;
    sample_vector = zeros(1, total_vars);
    % add root (assumed to be 1)
      root = find_root(tree_matrix);
    assert(isempty(find(tree_matrix(root,:), 1)));
    search_queue.add(root);
       
       while(search_queue.size() >  0)
            node_id = search_queue.poll();
            parent_id = find(tree_matrix(node_id,:));
            if (isempty(parent_id))
                mult_matrix = mnrnd(1, CPT{node_id}, NUM_SAMPLES);
                [r c] = find(mult_matrix');
                samples(:,node_id) = r;
            else
                cpt_matrix = repmat(CPT{node_id}, 1, NUM_SAMPLES);
                parent_samples = samples(:,parent_id)'; 
                num_parent_vals = size(CPT{node_id}, 2);
                shift_vector = 0:num_parent_vals:(num_parent_vals*(NUM_SAMPLES-1));
                parent_samples_index = parent_samples + shift_vector;
                cond_distr_matrix = cpt_matrix(:, parent_samples_index);
                mult_matrix = mnrnd(ones(NUM_SAMPLES, 1), cond_distr_matrix');
                [r c] = find(mult_matrix');
                samples(:, node_id) = r;
            end
            children = find(tree_matrix(:,node_id));
            for i=1:1:length(children) 
                search_queue.add(children(i));
            end
       end
end
