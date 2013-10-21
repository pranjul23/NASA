% Ankur Parikh
% SAILING Lab
% Carnegie Mellon University

% computes joint of query vars in a tree using CPT as the estimated CPTs
% assumes 1 is the root

function [evidence_marginals, evidence_joint_prob] = RunMatrixBeliefPropagation(tree_matrix, CPT, evidence_vars, evidence_vals_matrix, Ko)

    orig_evidence_vals_matrix = evidence_vals_matrix;
    N = size(evidence_vals_matrix, 1);
    index_scale_vector = 0:Ko:(Ko*(N-1));
    evidence_vals_matrix = evidence_vals_matrix + (index_scale_vector' * ones(1, length(evidence_vars))); 

    total_vars = size(tree_matrix, 1);
    
    % make upward pass first
    upward_CPT = cell(total_vars, 1);
    for v=1:1:total_vars
       upward_CPT{v} = repmat(CPT{v}, N, 1); 
    end
    upward_tree_matrix = tree_matrix;
    outgoing_message_array = cell(total_vars, 1);
    
    global_root = find_root(tree_matrix);

    % upward pass of belief propagation, which is all we need
    while(1)
        leaf_node = find_leaf(upward_tree_matrix);
        upward_tree_matrix(leaf_node,:) = 0; % allows the parent to become a leaf
        upward_tree_matrix(leaf_node,leaf_node) = 1; % so that this won't be picked as a leaf again
        K_leaf = size(CPT{leaf_node}, 1); 
        N_factor = upward_CPT{leaf_node};
                  
        if (isempty(find(evidence_vars == leaf_node, 1)))
           % integrate it out
           new_factor = zeros(N, size(N_factor,2));
           for k=1:1:K_leaf
               new_factor = new_factor + N_factor(k:K_leaf:end,:);
           end
        else
           evidence_index = find(evidence_vars == leaf_node);
           new_factor = N_factor(evidence_vals_matrix(:,evidence_index),:);
        end    
        outgoing_message_array{leaf_node} = new_factor;
        parent = find(tree_matrix(leaf_node,:));
        if (isempty(parent))
            break;
        else
            new_factor = new_factor';
            new_factor = new_factor(:);
                
            if (parent == global_root)
                upward_CPT{parent} = upward_CPT{parent} .* new_factor;
            else
                Kparent = size(upward_CPT{parent}, 2); 
                upward_CPT{parent} = upward_CPT{parent} .* (new_factor * ones(1,Kparent));
            end
        end    
    end
    
    if (isempty(find(evidence_vars == global_root, 1)))
        K_root = length(CPT{global_root});
        evidence_joint_prob = zeros(N,1);
        for k=1:1:K_root
            evidence_joint_prob = evidence_joint_prob + upward_CPT{global_root}(k:K_root:end);
        end
    else
        evidence_index = find(evidence_vars == global_root);
        evidence_joint_prob = upward_CPT{global_root}(evidence_vals_matrix(:,evidence_index));
    end
    
    evidence_marginals = cell(total_vars, 1);
    downward_CPT = cell(total_vars, 1);
    
    % compute downward pass
    downward_CPT{global_root} = upward_CPT{global_root};
    downward_tree_matrix = tree_matrix;
    while(1)
        root_node = find_root(downward_tree_matrix);
        if (root_node == -1)
            break
        end
        evidence_marginals{root_node} = downward_CPT{root_node}; % note that these are pairwise marginals, except the root 
        children = find(tree_matrix(:,root_node));
        downward_tree_matrix(:,root_node) = 0;
        downward_tree_matrix(root_node, root_node) = 1;
        for i=1:1:length(children)
            child_node = children(i);
            factor = downward_CPT{root_node};
            
            if (isempty(find(evidence_vars == root_node, 1)))
                % integrate it out
                if (size(factor, 2) > 1)
                    new_factor = sum(factor, 2); 
                else
                    new_factor = factor;
                end
                outgoing_divider = outgoing_message_array{child_node}';
                outgoing_divider = outgoing_divider(:);
                new_factor = new_factor ./ outgoing_divider;
                 Kchild = size(CPT{child_node}, 1);
                Kparent = size(CPT{child_node}, 2);
                
                new_factor_matrix = zeros(N*Kchild, Kparent);
                for k=1:1:Kparent
                    parent_indices = (k:Kparent:(Kparent*(N)))';
                    for m=1:1:Kchild
                        child_indices = (m:Kchild:(Kchild*(N)))';
                        new_factor_matrix(child_indices,k) = new_factor(parent_indices);
                    end
                end
                
                downward_CPT{child_node} = upward_CPT{child_node} .* new_factor_matrix;
                
            else
                evidence_index = find(evidence_vars == root_node);
                if (size(factor,2) > 1)
                    % select
                    new_factor = sum(factor(evidence_vals_matrix(:,evidence_index),:), 2);
                else
                    new_factor = factor(evidence_vals_matrix(:,evidence_index));
                end
                outgoing_divider = zeros(N,1);
                Kchild = size(CPT{child_node},1);
                upward_CPT_vector = zeros(Kchild*N,1);
                Kparent = size(outgoing_message_array{child_node},2);
                for k=1:1:Kparent
                    temp_array = outgoing_message_array{child_node}(:,k);
                    temp_array(~(orig_evidence_vals_matrix(:,evidence_index) == k)) = 0;
                    outgoing_divider = outgoing_divider + temp_array;
                    
                    temp_array_2 = upward_CPT{child_node}(:,k);
                    for m=1:1:Kchild
                        indices = (m:Kchild:(Kchild*(N)))';
                        sub_indices = indices(~(orig_evidence_vals_matrix(:,evidence_index) == k));
                        temp_array_2((sub_indices)) = 0;
                    end
                    upward_CPT_vector = upward_CPT_vector + temp_array_2;
                end
                
                new_factor = new_factor ./ outgoing_divider;
                
                new_factor_matrix = zeros(Kchild*N, 1);
                for k=1:1:Kchild
                    child_indices = (k:Kchild:(Kchild*(N)))';
                    new_factor_matrix(child_indices) = new_factor;
                end
                
                downward_CPT{child_node} = upward_CPT_vector  .* new_factor_matrix;
            end    
        end
    end
end

function leaf_node = find_leaf(tree_matrix)
    total_vars = size(tree_matrix, 1);
    leaf_node = -1;
    for v=1:1:total_vars
        if (sum(tree_matrix(:,v)) == 0)
           leaf_node = v;
           break;
        end
    end
end

function root_node = find_root(tree_matrix)
    root_node = -1;
    total_vars = size(tree_matrix, 1);
    for v=1:1:total_vars
        if (sum(tree_matrix(v,:)) == 0)
            root_node = v;
            break;
        end
    end
end