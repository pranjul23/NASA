%compute P(O1,O3)

function O1O3 = computeO1O3_ver2(A1, A, D1, D, O)

[o x] = size(O);
d = size(D1,1);

O3X2D2 = zeros(o,x,d);

for ind_o3 = 1:o
    for ind_d2 = 1:d
        for ind_x2 = 1:x
                        
            for ind_x3 = 1:x
                
                O3X2D2(ind_o3, ind_x2, ind_d2) = O3X2D2(ind_o3, ind_x2, ind_d2) + ...
                    O(ind_o3, ind_x3) * A(ind_x3, ind_x2, ind_d2);
            end
            
        end        
    end
end

%normalize
for ind_x2 = 1:x
    for ind_d2 = 1:d        
        O3X2D2(:, ind_x2, ind_d2) = O3X2D2(:, ind_x2, ind_d2)/sum(O3X2D2(:, ind_x2, ind_d2));        
    end
end

F = zeros(o, x*d);
for ind_o3 = 1:o
for ind_x2 = 1:x
for ind_d2 = 1:d 
    
    F(ind_o3, ind_d2 + d*(ind_x2-1)) = O3X2D2(ind_o3, ind_x2, ind_d2);
    
end
end
end

O1X2D2 = zeros(o,x,d);

for ind_o1 = 1:o
    for ind_x2 = 1:x
        for ind_d2 = 1:d
          
            for ind_x1 = 1:x
            for ind_d1 = 1:d
                
                O1X2D2(ind_o1, ind_x2, ind_d2) = O1X2D2(ind_o1, ind_x2, ind_d2) + ...
                    O(ind_o1, ind_x1) * A1(ind_x1) * D1(ind_d1, ind_x1) * D(ind_d2, ind_x2, ind_d1) * A(ind_x2, ind_x1, ind_d1);
                
            end
            end
            
        end
    end
end


%normalize
O1X2D2 = O1X2D2/sum(sum(sum(O1X2D2)));



O1O3 = zeros(o,o);

for ind_o1 = 1:o
    for ind_o3 = 1:o
        
        for ind_x2 = 1:x
        for ind_d2 = 1:d    
            
            O1O3(ind_o1, ind_o3) = O1O3(ind_o1, ind_o3) + ...
                O1X2D2(ind_o1, ind_x2, ind_d2) * O3X2D2(ind_o3, ind_x2, ind_d2);
            
        end
        end
    end
end
    