% Ankur Parikh
% SAILING Lab
% Carnegie Mellon University

function full_tree_matrix = trace_ancestors(tree_matrix)

    [V V] = size(tree_matrix);
    full_tree_matrix = zeros(V,V);
    % not very efficient but assumes tree depth is short
    for v=1:1:V
        curr_node = v;
        count = 1;
        while(1)
            parent_index = find(tree_matrix(curr_node,:) == 1);
            if (size(parent_index,2) == 0)
                break;
            end
            full_tree_matrix(v,parent_index) = count;
            curr_node = parent_index;
            count = count + 1;    
        end
    end
end