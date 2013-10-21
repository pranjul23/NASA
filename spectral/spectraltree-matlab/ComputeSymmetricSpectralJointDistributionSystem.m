% Computes joint distribution of tree

function [spectral_joint_val, EmpiricalProbMap, transform_array] = ComputeSymmetricSpectralJointDistributionSystem(tree_matrix, full_tree_matrix, type_indicator,...
                                                                        Ko, Kh, training_samples, test_sample, true_CPT, EmpiricalProbMap, transform_array, ...
                                                                        best_condition_matrix, train_flag, system_size)

total_vars = size(tree_matrix, 1);
global_root = find_root(tree_matrix);

if (train_flag == 1)
    % compute transformed version of each variable
    transform_array = cell(total_vars, 1);



    % mark the correct observations so we know how to compute the transforms
    [transform_array, EmpiricalProbMap] = mark_transform_hidden(global_root, -1, transform_array, ...
                                                        tree_matrix, full_tree_matrix, type_indicator, best_condition_matrix, Ko, Kh,...
                                                        training_samples, true_CPT, EmpiricalProbMap, system_size);


    for v=1:1:total_vars
        if (tree_matrix(v,v) == -1)
            continue
        end
        if (type_indicator(v) == 0)
            % compute hidden node transforms
            [transform_array, EmpiricalProbMap] = compute_transform_hidden(v, transform_array, training_samples, type_indicator, Ko, Kh,...
                                                                       true_CPT, tree_matrix, EmpiricalProbMap);
        else
            % compute observed node transforms
            [transform_array, EmpiricalProbMap] = compute_transform_obs(v, transform_array, training_samples, test_sample,...
                                                                   type_indicator, Ko, Kh, true_CPT, tree_matrix, EmpiricalProbMap);
                                                                   
        end
    end
end


for v=1:1:total_vars

    if (type_indicator(v) == 1)
       transform_array{v}.transform = transform_array{v}.pre_transform(test_sample(v),:); 
    end
end

% compute the joint distribution
spectral_joint_val = compute_spectral_joint(global_root, transform_array, tree_matrix, full_tree_matrix, type_indicator, global_root);

end

% computes the joint recursively
function curr_joint = compute_spectral_joint(node_id, transform_array, tree_matrix, full_tree_matrix, type_indicator, global_root)

    curr_joint = transform_array{node_id}.transform;
    if (type_indicator(node_id) == 1)
        return; % its a leaf
    else
    
        children = find(tree_matrix(:,node_id));

        mult_index = 1;
        for i=1:1:length(children)
            child_node = children(i);
            curr_vector = compute_spectral_joint(child_node, transform_array, tree_matrix, full_tree_matrix, type_indicator, global_root);
            curr_joint = ttv(curr_joint, curr_vector', 1); % tensor vector multiplication
            mult_index = mult_index + 1;
        end
        if (node_id ~= global_root)
            curr_joint = curr_joint.data';
        end
    end
end



function [transform_array, EmpiricalProbMap] = mark_transform_hidden(node_id, backward_obs, transform_array, ...
                                                tree_matrix, full_tree_matrix, type_indicator, best_condition_matrix, Ko, Kh,...
                                                training_samples, CPT, EmpiricalProbMap, system_size)

    transform_array{node_id}.backward = [];
    transform_array{node_id}.backward_extra = [];
    transform_array{node_id}.forward = [];

    % backward marking
    if (backward_obs > 0)
        transform_array{node_id}.backward = backward_obs;
        transform_array{node_id}.backward_extra = find_extra_obs(backward_obs, node_id,  best_condition_matrix, ...
                                                                 full_tree_matrix, type_indicator, tree_matrix, system_size);
 
       extra_obs = transform_array{node_id}.backward_extra;
    
        big_joint = [];
        for z=1:1:length(extra_obs)
            [temp_joint, EmpiricalProbMap] = ComputeFastEmpiricalProb([backward_obs transform_array{node_id}.backward_extra(1)], Ko, training_samples, type_indicator, tree_matrix, CPT, EmpiricalProbMap);
            big_joint = [big_joint temp_joint];
        end
    
        transform_array{node_id}.U = ComputeU(big_joint, Kh);
    end

    % forward marking    
    children = find(tree_matrix(:,node_id));
    for i=1:1:length(children)
		child_node = children(i);
        if (type_indicator(child_node) == 0)
                    
            % pick random observation with this parent
            child_obs_vector = find(full_tree_matrix(:,child_node) & type_indicator);
            transform_array{node_id}.forward = [transform_array{node_id}.forward child_obs_vector(1)];
            [transform_array, EmpiricalProbMap] = mark_transform_hidden(child_node, child_obs_vector(1), transform_array, tree_matrix,...
                                                        full_tree_matrix, type_indicator, best_condition_matrix, Ko, Kh, training_samples, ...
                                                        CPT, EmpiricalProbMap, system_size);
        else
            transform_array{node_id}.forward = [transform_array{node_id}.forward child_node];
            [transform_array, EmpiricalProbMap] = mark_transform_obs(child_node, transform_array, tree_matrix, full_tree_matrix,...
                                                   type_indicator,  best_condition_matrix, Ko, Kh, training_samples, CPT, EmpiricalProbMap, system_size);
        end
    end
end

function [transform_array, EmpiricalProbMap] = mark_transform_obs(node_id, transform_array, tree_matrix, full_tree_matrix, ...
                                                                  type_indicator, best_condition_matrix, Ko, Kh, training_samples, ...
                                                                  CPT, EmpiricalProbMap, system_size)
         
    transform_array{node_id}.backward = node_id;
    transform_array{node_id}.backward_extra = find_extra_obs(node_id, node_id, best_condition_matrix, ...
                                                             full_tree_matrix, type_indicator, tree_matrix, system_size);
    transform_array{node_id}.forward = [];

    extra_obs = transform_array{node_id}.backward_extra;
    
    big_joint = [];
    for z=1:1:length(extra_obs)
        [temp_joint, EmpiricalProbMap] = ComputeFastEmpiricalProb([transform_array{node_id}.backward transform_array{node_id}.backward_extra(z)], Ko, training_samples, type_indicator, tree_matrix, CPT, EmpiricalProbMap);
        big_joint = [big_joint temp_joint];
    end
    
    transform_array{node_id}.U = ComputeU(big_joint, Kh);
end



function  [transform_array, EmpiricalProbMap] = compute_transform_hidden(node_id, transform_array, training_samples, type_indicator, Ko, Kh, true_CPT, tree_matrix, EmpiricalProbMap)
        
        % we need to compute joint distribution of all things in left
        backward_obs = transform_array{node_id}.backward;
        forward_list = transform_array{node_id}.forward;
        extra_obs = transform_array{node_id}.backward_extra;
        
        % hidden node can't be a leaf
        assert(~isempty(forward_list));
        
        
        if (length(extra_obs) == 0)
            [forward_joint, EmpiricalProbMap] = ComputeFastEmpiricalProb([forward_list extra_obs], Ko, training_samples, type_indicator, tree_matrix, true_CPT, EmpiricalProbMap);
            forward_joint = tensor(forward_joint);
            
            children = find(tree_matrix(:,node_id));
            for i=1:1:length(forward_list)
                U = transform_array{children(i)}.U;
                forward_joint = ttm(forward_joint, U', i);
            end
                
            transform_array{node_id}.transform = forward_joint;
        else     
            %forward transform part
            big_Amat = [];
            big_bvec = [];
            fakeX = zeros(Kh * ones(length(forward_list) + 1, 1)');
            for z=1:1:length(extra_obs)
                [forward_joint, EmpiricalProbMap] = ComputeFastEmpiricalProb([forward_list extra_obs(z)], Ko, training_samples, type_indicator, tree_matrix, true_CPT, EmpiricalProbMap);
                forward_joint = tensor(forward_joint);

                children = find(tree_matrix(:,node_id));
                for i=1:1:length(forward_list)
                    U = transform_array{children(i)}.U;
                    forward_joint = ttm(forward_joint, U', i);
                end


                assert(length(backward_obs) == 1);

                assert(backward_obs ~= extra_obs(z));
                U = transform_array{node_id}.U;
                pinv_tensor = ComputeFastEmpiricalProb([backward_obs extra_obs(z)], Ko, training_samples, type_indicator, ...
                                                       tree_matrix, true_CPT, EmpiricalProbMap);
                
                A_tensor = U' * pinv_tensor;
                
                [Amat bvec] = CreateSpectralLinearSystem(A_tensor, forward_joint, fakeX, 1 , length(forward_list) + 1);
                big_Amat = [big_Amat; Amat];
                big_bvec = [big_bvec; bvec];
            end

            total_transform_vec = big_Amat \ big_bvec;
            total_transform = zeros(size(fakeX));
            val_matrix = createAllCombos(size(total_transform));
            
            for n=1:1:size(val_matrix,1)
               cell_index = num2cell(val_matrix(n,:)); 
               total_transform(cell_index{:}) = total_transform_vec(n); 
            end
            
            transform_array{node_id}.transform = tensor(total_transform);
        end
end

% compute transform for whole sequence at once and store it

function [transform_array, EmpiricalProbMap] = compute_transform_obs(node_id, transform_array, training_samples, test_sample, type_indicator, Ko, Kh, tree_matrix, CPT, EmpiricalProbMap)

     backward_obs = transform_array{node_id}.backward;
     extra_obs = transform_array{node_id}.backward_extra;
         
     big_Amat = [];
     big_bvec = [];
     fakeX = zeros([Ko Kh]);
     for z=1:1:length(extra_obs)
         
         [pinv_matrix, EmpiricalProbMap] = ComputeFastEmpiricalProb([backward_obs extra_obs(z)],... 
                                                                   Ko, training_samples, type_indicator, tree_matrix, CPT, EmpiricalProbMap);

          U = transform_array{node_id}.U;
          Atensor = U' * pinv_matrix;
          
          [Amat, bvec] = CreateSpectralLinearSystem(Atensor, pinv_matrix, fakeX, 1, 2);
          big_Amat = [big_Amat; Amat];
          big_bvec = [big_bvec; bvec];
     end
       
     total_transform_vec = big_Amat \ big_bvec;
     total_transform = zeros(size(fakeX));
     val_matrix = createAllCombos(size(total_transform));
            
     for n=1:1:size(val_matrix,1)
        cell_index = num2cell(val_matrix(n,:)); 
        total_transform(cell_index{:}) = total_transform_vec(n); 
     end
     transform_array{node_id}.pre_transform = total_transform;     
end

% find the extra observation, can be any observation that is not a
% descendant of invalid_ancestors set, we would prefer to pick the best
% conditioned one
function extra_obs = find_extra_obs(first_obs, invalid_ancestors, best_condition_matrix, full_tree_matrix, type_indicator, tree_matrix, system_size)

    extra_obs = [];
    best_condition_vector = best_condition_matrix(first_obs,:);
    for v=1:1:length(best_condition_vector)
        v_node = best_condition_vector(v);
        if(type_indicator(v_node) == 0);
            continue;
        end
        if (tree_matrix(v_node, v_node) == -1)
            continue
        end
        if (~isempty(find(invalid_ancestors == v_node, 1)))
            continue
        end
        
        ascendants = find(full_tree_matrix(v_node,:));
        if (~isempty(intersect(ascendants, invalid_ancestors)))
            continue
        end
        
        extra_obs = [extra_obs v_node];
        if (length(extra_obs) >= system_size)
            break
        end
%         if (~isempty(find(restricted_siblings == v_node)))
%            extra_obs = v_node;
%            break;
%         else
%             ascendants = find(full_tree_matrix(v_node,:));
%             intersection = intersect(ascendants, restricted_siblings);
%             if (~isempty(intersection))
%                 assert(length(intersection) == 1);
%                 extra_obs = v_node;
%                 break;
%             end
%         end
    end
end

function U = ComputeU(M, Kh)
  [U, junk1, junk2] = svd(M);
  U = U(:,1:Kh); % take only the first two columns
end