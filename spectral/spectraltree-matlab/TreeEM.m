% Ankur Parikh
% SAILING Lab
% Carnegie Mellon University

% Implements basic EM for a tree
% Uses Matricized Belief Propagation in E step to go fast

function best_estimated_CPT = TreeEM(tree_matrix, training_samples, type_indicator, Ko, Kh, PRECISION, NUM_RESTARTS)
    
    best_obj = -Inf;
    best_estimated_CPT = {};
    for n=1:1:NUM_RESTARTS
        total_vars = size(tree_matrix, 1);
	
        estimated_CPT = cell(total_vars, 1);

        for v=1:1:total_vars
            if (sum(tree_matrix(v,:)) == 0)
               if (type_indicator(v) == 0)
                    estimated_CPT{v} = GenerateRandomCPT(Kh);
               else
                   estimated_CPT{v} = GenerateRandomCPT(Ko);
               end
            else
                parent_id = find(tree_matrix(v,:));
                if (type_indicator(v) == 0 && type_indicator(parent_id) == 0)
                    estimated_CPT{v} = GenerateRandomCPT([Kh, Kh]);
                elseif (type_indicator(v) == 0 && type_indicator(parent_id) == 1)
                    estimated_CPT{v} = GenerateRandomCPT([Kh, Ko]);
                elseif (type_indicator(v) == 1 && type_indicator(parent_id) == 0)
                    estimated_CPT{v} = GenerateRandomCPT([Ko, Kh]);
                elseif (type_indicator(v) == 1 && type_indicator(parent_id) == 1)
                    estimated_CPT{v} = GenerateRandomCPT([Ko, Ko]);
                end          
            end
        end

        obj_val = -Inf;
        new_obj_val = -Inf;
       % epsilon = 1;
        prev_likelihood = Inf;
        isConverged = false;
        while(~isConverged)
            [evidence_marginals_array evidence_joint_prob_array] = E_step(tree_matrix, training_samples, type_indicator, estimated_CPT, Ko);
            new_obj_val = sum(log(evidence_joint_prob_array));
            new_estimated_CPT = M_step(tree_matrix, training_samples, type_indicator, evidence_marginals_array, evidence_joint_prob_array, Kh, Ko);
            
            curr_likelihood = sum(log(evidence_joint_prob_array))
          %  epsilon = ComputeEpsilon(estimated_CPT, new_estimated_CPT);
            estimated_CPT = new_estimated_CPT;
          %  epsilon
      %      assert(epsilon > 0);   
            obj_val = new_obj_val;
            if (prev_likelihood ~= Inf)
                isConverged = TestConvergence(prev_likelihood,  curr_likelihood, PRECISION);
            end
            prev_likelihood = curr_likelihood;
        end
        
        if (obj_val > best_obj)
            best_obj = obj_val;
            best_estimated_CPT = estimated_CPT;
        end
    end
 end

function [evidence_marginals_array, evidence_joint_prob_array] = E_step(tree_matrix, training_samples, type_indicator, estimated_CPT, Ko)

    [N, total_vars] = size(training_samples); 
    [all_evidence_marginals, evidence_joint_prob_array] = RunMatrixBeliefPropagation(tree_matrix, estimated_CPT, find(type_indicator), training_samples(:,logical(type_indicator)), Ko);
    evidence_marginals_array = zeros(total_vars, N, Ko, Ko);
    for v=1:1:total_vars
        marginals_matrix = all_evidence_marginals{v};
        parent_id = find(tree_matrix(v,:));
        if (isempty(parent_id))
            z1 = length(marginals_matrix) / N;
            new_marginals_matrix = zeros(N*z1,z1);
            for k=1:1:z1
                z1_indices = (k:z1:(z1*N))';
                new_marginals_matrix((k:z1:end)',k) = marginals_matrix(z1_indices);
            end
            marginals_matrix = new_marginals_matrix;
        end
        [z1 z2] = size(marginals_matrix);
        z1 = z1 / N;
        for k=1:1:z1
            z1_indices = (k:z1:(z1*(N)))';
            evidence_marginals_array(v,:,k,1:z2) = marginals_matrix(z1_indices,:) ./ (evidence_joint_prob_array * ones(1,z2)); 
        end    
    end
    
%     temp_evidence_marginals_array = zeros(total_vars, N, Ko, Ko);
%     temp_evidence_joint_prob_array = zeros(N,1);
%     for n=1:1:N
%         [evidence_marginal, evidence_joint_prob] = RunBeliefPropagation(tree_matrix, estimated_CPT, find(type_indicator), training_samples(n,logical(type_indicator)));
%         total_vars = length(evidence_marginal);
%         for v=1:1:total_vars
%             marginals_matrix = evidence_marginal{v} / evidence_joint_prob;
%             parent_id = find(tree_matrix(v,:));
%             if (isempty(parent_id))
%                 marginals_matrix = diag(marginals_matrix);
%             end
%             [z1 z2] = size(marginals_matrix);
%             temp_evidence_marginals_array(v, n,1:z1,1:z2) = marginals_matrix; 
%         end
%         temp_evidence_joint_prob_array(n) = evidence_joint_prob;
%     end
end

function estimated_CPT = M_step(tree_matrix, training_samples, type_indicator, evidence_marginals_array, evidence_joint_prob_array, Kh, Ko)

    total_vars = size(tree_matrix, 1);
	estimated_CPT = cell(total_vars, 1);
     N = size(training_samples, 1); 
     
	% compute empirical quantities
	for v=1:1:total_vars
        parent_id = find(tree_matrix(v,:));
        if (isempty(parent_id))
            estimated_CPT{v} = zeros(Kh,1);
            if (type_indicator(v) == 0)
                marginal = squeeze(mean(evidence_marginals_array(v,:,:,:),2));
                estimated_CPT{v} = diag(marginal(1:Kh,1:Kh));
            else
             for k=1:1:Ko
                estimated_CPT{v}(k) = sum(training_samples(:,v) == k) / N;
             end
            end
            continue
        end
   
        % we have to do pairwise marginals, compute the joint first,
        % normalize at the end
        
        % we have to pick out the evidence that is correct
        if (type_indicator(v) == 1 && type_indicator(parent_id) == 1)
            % just compute counts
            estimated_CPT{v} = zeros(Ko,Ko);
            for i=1:1:Ko % parent
                for j=1:1:Ko %v
                    subsample_indicator = (training_samples(:,parent_id) == i);
                    estimated_CPT{v}(j,i) = sum(training_samples(subsample_indicator,v) == j);
                end
            end
        elseif (type_indicator(v) == 0 && type_indicator(parent_id) == 1)
           % select samples such that parent has the correct value then use
           % evidence marginals
           estimated_CPT{v} = zeros(Kh,Ko);
           for i=1:1:Ko
               for j=1:1:Kh
                    subsample_indicator = (training_samples(:,parent_id) == i);
                    estimated_CPT{v}(j,i)  = sum(evidence_marginals_array(v,subsample_indicator, j,i));
               end
           end
        elseif (type_indicator(v) == 1 && type_indicator(parent_id) == 0)
            estimated_CPT{v} = zeros(Ko,Kh);
            for i=1:1:Kh
                for j=1:1:Ko
                    subsample_indicator = (training_samples(:,v) == j);
                    estimated_CPT{v}(j,i)  = sum(evidence_marginals_array(v,subsample_indicator, j,i));
                end
            end
        else
            estimated_CPT{v} =  squeeze(sum(evidence_marginals_array(v,:,1:Kh,1:Kh),2));
        end
        
        % normalize
        assert(size(estimated_CPT{v},2) == Kh);
        estimated_CPT{v} = estimated_CPT{v} + .00000001; % just so things we don't get NaN's
        for k=1:1:Kh   % this is a bad hack
            estimated_CPT{v}(:,k) = estimated_CPT{v}(:,k) ./ sum(estimated_CPT{v}(:,k));
        end
    end
end   

function epsilon = ComputeEpsilon(estimated_CPT, new_estimated_CPT)
    total_vars = length(estimated_CPT);
    epsilon = 0;
    for v=1:1:total_vars
        epsilon = epsilon + norm(estimated_CPT{v} - new_estimated_CPT{v}, 'fro');
    end
    epsilon = epsilon / total_vars;
end


function cpt_matrix = GenerateRandomCPT(Karray)
    if (length(Karray) == 1)
        cpt_matrix = rand(Karray(1),1);
        cpt_matrix = cpt_matrix / sum(cpt_matrix);
    elseif (length(Karray) == 2)
       cpt_matrix = rand(Karray(1),Karray(2));
       for k=1:1:Karray(2)
        cpt_matrix(:,k) = cpt_matrix(:,k) / sum(cpt_matrix(:,k));
       end
    end
end

function isConverged = TestConvergence(prev_likelihood,  curr_likelihood, thresh)

	avg = abs(prev_likelihood + curr_likelihood)  / 2.0;
	diff = abs(curr_likelihood - prev_likelihood);

	if ((diff / avg) < thresh)
		isConverged = true;
	else
		isConverged = false;
    end
end



