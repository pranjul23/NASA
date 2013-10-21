% Ankur Parikh
% SAILING Lab
% Carnegie Mellon University

% Computes probability of the set of variables test_vars (assumed to be all the leaves), with values test_vals
% Must be given "transform_array" as an argument as this is the spectral parameters

function [spectral_joint_val] = SpectralTest(tree_matrix,  type_indicator, Ko, Kh, transform_array, test_vars, test_vals)
	
	global_root = find_root(tree_matrix);
	for i=1:1:length(test_vars)
		node = test_vars(i);
		assert(type_indicator(node) == 1);
		transform_array{node}.transform = transform_array{node}.pre_transform(test_vals(i),:);
	end

	% compute the joint distribution
	spectral_joint_val = compute_spectral_joint(global_root, transform_array, tree_matrix, type_indicator, global_root);
	if (spectral_joint_val < 0)
		spectral_joint_val = 0;
	elseif (spectral_joint_val > 1)
		spectral_joint_val = 1;
	end
end

% computes the joint recursively
function curr_joint = compute_spectral_joint(node_id, transform_array, tree_matrix, type_indicator, global_root)

    curr_joint = transform_array{node_id}.transform;
    if (type_indicator(node_id) == 1)
        return; % its a leaf
    else
    
        children = find(tree_matrix(:,node_id));

        mult_index = 1;
        for i=1:1:length(children)
            child_node = children(i);
            curr_vector = compute_spectral_joint(child_node, transform_array, tree_matrix, type_indicator, global_root);
            curr_joint = ttv(curr_joint, curr_vector', 1); % tensor vector multiplication
            mult_index = mult_index + 1;
        end
        if (node_id ~= global_root)
            curr_joint = curr_joint.data';
        end
    end
end