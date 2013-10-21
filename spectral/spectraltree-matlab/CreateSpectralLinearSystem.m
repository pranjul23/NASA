function [Amat bvec] = CreateSpectralLinearSystem(Atensor, Btensor, fakeX, mult_modeA, mult_modeX)

	numMultElements = size(Atensor, mult_modeA);	

    
%     if (mult_modeX == 2)
%         mult_modeX = 1;
%     elseif(mult_modeX == 1)
%         mult_modeX = 2;
%     end
%     if (mult_modeA == 2)
%         mult_modeA = 1;
%     elseif (mult_modeA == 1)
%         mult_modeA = 2;
%     end
    

    num_B_elements = prod(size(Btensor));
    num_X_elements = prod(size(fakeX));
    
    Amat = zeros(num_B_elements, num_X_elements);
    bvec = zeros(num_B_elements, 1);
 
    B_val_matrix = createAllCombos(size(Btensor));
    dims = size(fakeX);
 %   row_dim = dims(1); 
 %   col_dim = dims(2);
 %   dims(1) = col_dim;
 %  dims(2) = row_dim;
    Xoffsets = ComputeOffsets(dims);

    for n=1:1:num_B_elements
   %    n
        B_cell_index = num2cell(B_val_matrix(n,:));
        bvec(n) = Btensor(B_cell_index{:});
    %    bvec(n) = Btensor(n);
    %    val_array = B_val_matrix(n,:);
        
        % hack because matlab switches row/column order
      %  col = val_array(2);
   %     row = val_array(1);
    %    val_array(1) = col;
    %    val_array(2) = row;
        for k=1:1:numMultElements
    %        k
            x_cell_index = B_cell_index;
          %  x_cell_index = val_array;
            x_cell_index{mult_modeX} = k;
            
            A_cell_index =[k B_cell_index{end}];            
            
            xindex = ComputeIndex(x_cell_index, Xoffsets);
            
            Amat(n, xindex) = Atensor(k, B_cell_index{end});
        end
    end
end


function offsets = ComputeOffsets(dims)

    offsets = zeros(length(dims), 1);
   
	for i=length(dims):-1:1
		if (i == length(dims))
			offsets(i) = 1;
        else
			offsets(i) = offsets(i+1) * dims(i + 1);
        end
    end
end

function index = ComputeIndex(val_array, offsets)

    index = 1;
    for i=1:1:length(val_array)
        index = index + (val_array{i} - 1)*offsets(i);
    end
    
end



