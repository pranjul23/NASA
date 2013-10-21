% Ankur Parikh
% SAILING Lab
% Carnegie Mellon University

% computes joint of query vars in a tree using CPT as the estimated CPTs

function [evidence_joint_prob] = RunBeliefPropagation(tree_matrix, CPT, evidence_vars, evidence_vals)

    total_vars = size(tree_matrix, 1);
    
    % make upward pass first
    upward_CPT = CPT;
    upward_tree_matrix = tree_matrix;
    outgoing_message_array = cell(total_vars, 1);
    
    global_root = find_root(tree_matrix);

    % upward pass of belief propagation, which is all we need
    while(1)
        leaf_node = find_leaf(upward_tree_matrix);
        upward_tree_matrix(leaf_node,:) = 0; % allows the parent to become a leaf
        upward_tree_matrix(leaf_node,leaf_node) = 1; % so that this won't be picked as a leaf again
        
         factor = upward_CPT{leaf_node};
                  
        if (isempty(find(evidence_vars == leaf_node, 1)))
           % integrate it out 
           new_factor = sum(factor, 1);   
        else
           evidence_index = find(evidence_vars == leaf_node);
           new_factor = factor(evidence_vals(evidence_index),:);
        end    
        outgoing_message_array{leaf_node} = new_factor';
        parent = find(tree_matrix(leaf_node,:));
        if (isempty(parent))
            break;
        else
            if (parent == global_root)
                upward_CPT{parent} = upward_CPT{parent} .* new_factor';
            else
                Kparent = size(upward_CPT{parent}, 2); 
                upward_CPT{parent} = upward_CPT{parent} .* (new_factor' *ones(1,Kparent));
            end
        end    
    end
    
    if (isempty(find(evidence_vars == global_root, 1)))
        evidence_joint_prob = sum(upward_CPT{global_root});
    else
        evidence_index = find(evidence_vars == global_root);
        evidence_joint_prob = upward_CPT{global_root}(evidence_vals(evidence_index));
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