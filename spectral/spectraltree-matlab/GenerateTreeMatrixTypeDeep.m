% Ankur Parikh
% SAILING Lab
% Carnegie Mellon University

% Generates a tree matrix of type 1

function input_tree_struct = GenerateTreeMatrixTypeDeep(Ko, Kh, num_levels)

    total_subtree_vars = 2^(num_levels) - 1;
    num_subtree_hidden = 2^(num_levels - 1) - 1;
 %   num_subtree_obs = 3 * 2^(num_levels - 2);
  %  total_subtree_vars = num_subtree_hidden + num_subtree_obs;
    total_vars = 2*total_subtree_vars; % two binary trees placed together
    
    tree_matrix = zeros(total_vars, total_vars);
    full_tree_matrix = zeros(total_vars, total_vars); % complete version after running breadth first search 
    
    % make some binary trees
    binary_tree_matrix = zeros(total_subtree_vars, total_subtree_vars);
    for v=2:1:total_subtree_vars
        binary_tree_matrix(v,floor(v/2)) = 1;
    end

    
    % put the binary trees together
    tree_matrix = zeros(total_vars, total_vars);
    tree_matrix(1:total_subtree_vars, 1:total_subtree_vars) = binary_tree_matrix;
    tree_matrix((total_subtree_vars)+1:total_vars, (total_subtree_vars + 1):total_vars) = binary_tree_matrix;
    tree_matrix(total_subtree_vars + 1, 1) = 1; 

    % indicates if variable is observed or hidden
    type_indicator = zeros(total_vars,1);  
    observation_set = (num_subtree_hidden + 1) :total_subtree_vars;
    observation_set = [observation_set (total_subtree_vars + observation_set)];
    type_indicator(observation_set) = 1;
   
    % forward map and backward map of observations
    
    
    full_tree_matrix = trace_ancestors(tree_matrix);

    CPT = cell(total_vars,1);
    for v=1:1:total_vars
        if (sum(full_tree_matrix(v,:)) == 0)
            assert(type_indicator(v) == 0);
           CPT{v} = GenerateRandomCPT(Kh);
        else
            parent_id = find(tree_matrix(v,:));
            if (type_indicator(v) == 0 && type_indicator(parent_id) == 0)
                CPT{v} = GenerateRandomCPT([Kh, Kh]);
            elseif (type_indicator(v) == 0 && type_indicator(parent_id) == 1)
                CPT{v} = GenerateRandomCPT([Kh, Ko]);
            elseif (type_indicator(v) == 1 && type_indicator(parent_id) == 0)
                CPT{v} = GenerateRandomCPT([Ko, Kh]);
            elseif (type_indicator(v) == 1 && type_indicator(parent_id) == 1)
                CPT{v} = GenerateRandomCPT([Ko, Ko]);
            end    
                
        end
    end
    
    input_tree_struct.tree_matrix = tree_matrix;
    input_tree_struct.full_tree_matrix = full_tree_matrix;
    input_tree_struct.CPT = CPT;
    input_tree_struct.type_indicator = type_indicator; 
end



function cpt_matrix = GenerateRandomCPT(Karray)

    while (1)
    if (length(Karray) == 1)
        cpt_matrix = rand(Karray(1),1);
        cpt_matrix = cpt_matrix / sum(cpt_matrix);
    elseif (length(Karray) == 2)
       cpt_matrix = rand(Karray(1),Karray(2));
       for k=1:1:Karray(2)
        cpt_matrix(:,k) = cpt_matrix(:,k) / sum(cpt_matrix(:,k));
       end
    end
      if (min(min(cpt_matrix)) > .01)
         break; 
      end
    end
end




