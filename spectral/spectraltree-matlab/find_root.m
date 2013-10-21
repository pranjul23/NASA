% Ankur Parikh
% SAILING Lab
% Carnegie Mellon University

function root = find_root(tree_matrix)

    total_vars = size(tree_matrix, 1);
    
    for v=1:1:total_vars
        if (sum(tree_matrix(v,:)) == 0)
            if (tree_matrix(v,v) == -1)
                assert(0);
            end
            root = v;
            break;
        end 
    end
    
end