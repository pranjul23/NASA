%compute P(O1, O6)

function O1O6 = computeO1O6(A1, A, D1, D, O)

[o x] = size(O);
d = size(D1,1);


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


S = zeros(o, x*d);
for ind_o1 = 1:o
for ind_d2 = 1:d 
for ind_x2 = 1:x
    
    S(ind_o1, ind_x2 + x*(ind_d2-1)) = O1X2D2(ind_o1, ind_x2, ind_d2);
    
end
end
end



O6X2D2 = zeros(o, x, d);

for ind_o6 = 1:o
    for ind_d2 = 1:d
       for ind_x2 = 1:x 
        
        for ind_x3 = 1:x
        for ind_x4 = 1:x
        for ind_x5 = 1:x
        for ind_x6 = 1:x
        for ind_d3 = 1:d
        for ind_d4 = 1:d
        for ind_d5 = 1:d    
           
            O6X2D2(ind_o6, ind_x2, ind_d2) = O6X2D2(ind_o6, ind_x2, ind_d2) + ...
                O(ind_o6, ind_x6) * A(ind_x6, ind_x5, ind_d5) * D(ind_d5, ind_x5, ind_d4) * A(ind_x5, ind_x4, ind_d4) * D(ind_d4, ind_x4, ind_d3) * A(ind_x4, ind_x3, ind_d3) * D(ind_d3, ind_x3, ind_d2) * A(ind_x3, ind_x2, ind_d2);
        end
        end
        end
        end
        end
        end
        end
       end
    end
end


%normalize
for ind_x2 = 1:x
    for ind_d2 = 1:d
        O6X2D2(:, ind_x2, ind_d2) = O6X2D2(:, ind_x2, ind_d2)/sum(O6X2D2(:, ind_x2, ind_d2));
    end
end




F = zeros(o, x*d);
for ind_o6 = 1:o
for ind_d2 = 1:d 
for ind_x2 = 1:x
    
    F(ind_o6, ind_x2 + x*(ind_d2-1)) = O6X2D2(ind_o6, ind_x2, ind_d2);
    
end
end
end





O1O6 = zeros(o,o);

for ind_o1 = 1:o
    for ind_o6 = 1:o
        
        for ind_x2 = 1:x
        for ind_d2 = 1:d    
            
            O1O6(ind_o1, ind_o6) = O1O6(ind_o1, ind_o6) + ...
                O1X2D2(ind_o1, ind_x2, ind_d2) * O6X2D2(ind_o6, ind_x2, ind_d2);
            
        end
        end
    end
end
    