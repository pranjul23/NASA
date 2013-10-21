% Ankur Parikh
% SAILING Lab
% Carnegie Mellon University

function [best_condition_matrix, EmpiricalProbMap] = ComputeConditionNumbers(training_samples, type_indicator, Ko, Kh, EmpiricalProbMap)

    total_vars = length(type_indicator);
    eigenvalue_matrix = zeros(total_vars, total_vars);
  %  eigenvalue_matrix(logical(type_indicator), logical(type_indicator)) = 1;
    for i=1:1:total_vars
	    for j=(i+1):1:total_vars
           if (type_indicator(i) == 1 && type_indicator(j) == 1) 
	     [prob_matrix, EmpiricalProbMap] = ComputeFastEmpiricalProb([i j], Ko, training_samples, type_indicator, eigenvalue_matrix, cell(total_vars, 1), EmpiricalProbMap);
               mineig = min(abs(eigs(prob_matrix, Kh)));
               eigenvalue_matrix(i,j) = mineig;
           end
        end
    end

    
    eigenvalue_matrix = max(eigenvalue_matrix, eigenvalue_matrix');

    best_condition_matrix = zeros(total_vars, total_vars);
    
    for v=1:1:total_vars
        if (type_indicator(v) == 0)
           continue
        end
        vec = eigenvalue_matrix(v,:);
        [C IX] = sort(vec, 'descend');
        best_condition_matrix(v,:) = IX;
    end

end
